// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/ephemeris.hpp"
#include "core/dynamics_rotational.hpp"
#include "core/state.hpp"
#include "core/status.hpp"
#include "core/time.hpp"
#include "util/math.hpp"
#include "util/units.hpp"

#include <algorithm>
#include <cmath>

static StatusCode valid_ephemeris_epoch_metadata(const EphemerisEpochMetadata& epoch) {
    if (!std::isfinite(epoch.ref_epoch.day) || !std::isfinite(epoch.ref_epoch.frac))
        return StatusCode::invalid_input;

    if (epoch.ref_epoch.frac >= 1 || epoch.ref_epoch.frac < 0)
        return StatusCode::invalid_input; // normalize the jd?

    if (epoch.offset_unit != UTime::second) return StatusCode::invalid_input; // TEMP:

    return StatusCode::ok;
}

static StatusCode valid_ephemeris_times(const svec<f64>& dt) {
    if (dt.empty()) return StatusCode::empty_ephemeris;

    for (i32 i = 0; i < dt.size(); ++i) {
        if (!std::isfinite(dt[i])) return StatusCode::invalid_input;
        if (i > 0) {
            if (dt[i] - dt[i - 1] <= 0.0)
                return StatusCode::non_monotonic_time; // non monotonic
        }
    }

    return StatusCode::ok;
}

static StatusCode valid_cartesian_ephemeris_metadata(
    const CartesianEphemerisMetadata& metadata
) {
    StatusCode status;

    status = valid_ephemeris_epoch_metadata(metadata.epoch);
    if (status != StatusCode::ok) return status;

    if (metadata.frame.frame.empty() || metadata.frame.center.empty()
        || metadata.frame.object.empty())
        return StatusCode::invalid_ephemeris_metadata;

    if (metadata.source.source_name.empty() || metadata.source.source_type.empty())
        return StatusCode::invalid_ephemeris_metadata;

    if (metadata.units.length != ULength::kilometer
        || metadata.units.time != UTime::second)
        return StatusCode::invalid_ephemeris_metadata; // TEMP:

    return StatusCode::ok;
}

static StatusCode valid_orientation_ephemeris_metadata(
    const OrientationEphemerisMetadata& metadata
) {
    StatusCode status = valid_ephemeris_epoch_metadata(metadata.epoch);
    if (status != StatusCode::ok) return status;

    if (metadata.frame.object.empty() || metadata.frame.source_frame.empty()
        || metadata.frame.target_frame.empty())
        return StatusCode::invalid_ephemeris_metadata;

    if (metadata.source.source_name.empty() || metadata.source.source_type.empty())
        return StatusCode::invalid_ephemeris_metadata;

    // TEMP: require canonical units until conversion
    bool canonical_units = metadata.units.angular_velocity_angle == UAngle::radian
                           && metadata.units.angular_velocity_time == UTime::second;
    if (!canonical_units) return StatusCode::invalid_ephemeris_metadata;

    // TEMP: enforce StateAtt convention until conversion
    bool canonical_convention
        = metadata.convention.quaternion_order == QuaternionComponentOrder::vector_scalar
          && metadata.convention.rotation == RotationConvention::passive
          && metadata.convention.angular_velocity_frame
                 == AngularVelocityExpressionFrame::target
          && metadata.convention.angular_velocity_direction
                 == AngularVelocityDirection::target_relative_source;
    if (!canonical_convention) return StatusCode::invalid_ephemeris_metadata;

    return StatusCode::ok;
}

StatusCode validate_cartesian_ephemeris_table(const CartesianEphemerisTable& table) {
    StatusCode status;

    if (table.states.empty()) return StatusCode::empty_ephemeris;
    if (table.states.size() != table.dt.size()) return StatusCode::size_mismatch;

    status = valid_ephemeris_times(table.dt);
    if (status != StatusCode::ok) return status;

    status = valid_cartesian_ephemeris_metadata(table.metadata);
    if (status != StatusCode::ok) return status;

    for (i32 i = 0; i < table.states.size(); ++i) {
        if (!finite_state(table.states[i])) return StatusCode::invalid_state;
    }

    return StatusCode::ok;
}

StatusCode validate_orientation_ephemeris_table(
    const OrientationEphemerisTable& table,
    f64 quaternion_tol
) {
    if (!std::isfinite(quaternion_tol) || quaternion_tol <= 0.0)
        return StatusCode::invalid_input;

    if (table.states.empty()) return StatusCode::empty_ephemeris;
    if (table.states.size() != table.dt.size()) return StatusCode::size_mismatch;

    StatusCode status = valid_ephemeris_times(table.dt);
    if (status != StatusCode::ok) return status;

    status = valid_orientation_ephemeris_metadata(table.metadata);
    if (status != StatusCode::ok) return status;

    for (i32 i = 0; i < table.states.size(); ++i) {
        const StateAtt& state = table.states[i];
        if (!finite_state_att(state)) return StatusCode::invalid_att_state;
        if (std::abs(state.q.norm() - 1.0) > quaternion_tol)
            return StatusCode::invalid_att_state;

        if (!table.has_angular_velocity && state.w.norm() > quaternion_tol)
            return StatusCode::invalid_att_state;

        if (i > 0 && table.states[i - 1].q.dot(state.q) < 0.0)
            return StatusCode::invalid_att_state;
    }

    return StatusCode::ok;
}

StatusCode canonicalize_orientation_ephemeris_samples(
    OrientationEphemerisTable& table,
    f64 tol
) {
    if (!finite_pos(tol)) return StatusCode::invalid_input;

    for (i32 i = 0; i < table.states.size(); ++i) {
        StateAtt& state = table.states[i];

        if (!finite_state(state, tol)) return StatusCode::invalid_att_state;
        normalize_quaternion_inplace(state.q, tol);

        if (i > 0) {
            StateAtt& state_prev = table.states[i - 1];
            // prevent large quaternion hemisphere jumps for interpolation
            if (state_prev.q.dot(state.q) < 0.0) {
                state.q = -state.q;
            }
        }

        if (!table.has_angular_velocity) {
            state.w = vec3d0;
        }
    }

    return validate_orientation_ephemeris_table(table);
}

enum struct EphemerisQueryLocation { exact, interior, before_coverage, after_coverage };

struct EphemerisBracket {
    i32 lower = -1;
    i32 upper = -1;
    EphemerisQueryLocation location = EphemerisQueryLocation::interior;
};

static StatusCode find_ephemeris_bracket(
    const svec<f64>& epochs,
    f64 dt,
    f64 tol,
    EphemerisBracket& out
) {
    if (epochs.empty()) return StatusCode::empty_ephemeris;
    if (!std::isfinite(dt) || !finite_pos(tol)) return StatusCode::invalid_input;

    EphemerisBracket temp{};

    auto it = std::lower_bound(epochs.begin(), epochs.end(), dt); // least upper bound

    if (it == epochs.end()) {
        if (equal_tol(epochs.back(), dt, tol)) {
            temp.location = EphemerisQueryLocation::exact;
            temp.lower = epochs.size() - 1;
            temp.upper = epochs.size() - 1;
        } else {
            temp.location = EphemerisQueryLocation::after_coverage;
            temp.lower = epochs.size() - 1;
            temp.upper = epochs.size() - 1;
        }
    } else {
        i32 idx = static_cast<i32>(std::distance(epochs.begin(), it));

        if (equal_tol(epochs[idx], dt, tol)) {
            temp.location = EphemerisQueryLocation::exact;
            temp.lower = idx;
            temp.upper = idx;
        } else if (idx != 0 && equal_tol(epochs[idx - 1], dt, tol)) {
            temp.location = EphemerisQueryLocation::exact;
            temp.lower = idx - 1;
            temp.upper = idx - 1;
        } else if (idx == 0) {
            temp.location = EphemerisQueryLocation::before_coverage;
            temp.lower = idx;
            temp.upper = idx;
        } else {
            temp.location = EphemerisQueryLocation::interior;
            temp.lower = idx - 1;
            temp.upper = idx;
        }
    }

    out = temp;
    return StatusCode::ok;
}

static StatusCode validate_interior_ephemeris_bracket(
    const EphemerisBracket& bracket,
    const i32 n
) {
    if (n < 2 || bracket.location != EphemerisQueryLocation::interior || bracket.lower < 0
        || bracket.upper < 0 || bracket.lower >= n || bracket.upper >= n
        || bracket.upper != bracket.lower + 1) {
        return StatusCode::invalid_input;
    }

    return StatusCode::ok;
}

static StatusCode interpolate_cartesian_linear(
    const CartesianEphemerisTable& table,
    const EphemerisBracket& bracket,
    f64 dt,
    f64 tol,
    StateTr& out
) {
    StatusCode status;

    if (table.dt.size() != table.states.size()) return StatusCode::invalid_input;
    status = validate_interior_ephemeris_bracket(bracket, table.dt.size());
    if (status != StatusCode::ok) return status;

    StateTr temp{};

    f64 t0 = table.dt[bracket.lower];
    f64 t1 = table.dt[bracket.upper];
    f64 h = t1 - t0;

    if (!finite_pos(h, tol)) return StatusCode::non_monotonic_time;

    f64 u = (dt - t0) / h; // normalized interp time in [0,1] within tol
    // essentially u = (t - t0) / (t1 - t0)
    if (!std::isfinite(u) || u < -tol || u > 1.0 + tol) return StatusCode::invalid_input;
    u = std::clamp(u, 0.0, 1.0);

    const StateTr& x0 = table.states[bracket.lower];
    const StateTr& x1 = table.states[bracket.upper];

    temp.r = x0.r + u * (x1.r - x0.r);

    if (table.has_velocity) {
        temp.v = x0.v + u * (x1.v - x0.v);
    } else {
        temp.v = vec3d0;
    }

    if (!finite_state(temp)) return StatusCode::non_finite_result;

    out = temp;
    return StatusCode::ok;
}

static StatusCode interpolate_cartesian_cubic_hermite(
    const CartesianEphemerisTable& table,
    const EphemerisBracket& bracket,
    f64 dt,
    f64 tol,
    StateTr& out
) {
    // https://en.wikipedia.org/wiki/Cubic_Hermite_spline
    if (!table.has_velocity) return StatusCode::unsupported_method;

    StatusCode status;

    if (table.dt.size() != table.states.size()) return StatusCode::invalid_input;
    status = validate_interior_ephemeris_bracket(bracket, table.dt.size());
    if (status != StatusCode::ok) return status;

    StateTr temp{};

    const StateTr& x0 = table.states[bracket.lower];
    const StateTr& x1 = table.states[bracket.upper];

    f64 t0 = table.dt[bracket.lower];
    f64 t1 = table.dt[bracket.upper];
    f64 h = t1 - t0;

    if (!finite_pos(h, tol)) return StatusCode::non_monotonic_time;

    f64 u = (dt - t0) / h;
    if (!std::isfinite(u) || u < -tol || u > 1.0 + tol) return StatusCode::invalid_input;
    u = std::clamp(u, 0.0, 1.0);
    f64 u2 = u * u;
    f64 u3 = u2 * u;

    // cubic hermite polynomials
    f64 h00 = 2.0 * u3 - 3.0 * u2 + 1.0;
    f64 h10 = u3 - 2.0 * u2 + u;
    f64 h01 = -2.0 * u3 + 3.0 * u2;
    f64 h11 = u3 - u2;
    temp.r = h00 * x0.r + h10 * h * x0.v + h01 * x1.r + h11 * h * x1.v;

    // polynomial derivatives
    f64 dh00 = 6.0 * u2 - 6.0 * u;
    f64 dh10 = 3.0 * u2 - 4.0 * u + 1.0;
    f64 dh01 = -6.0 * u2 + 6.0 * u;
    f64 dh11 = 3.0 * u2 - 2.0 * u;
    temp.v = (dh00 * x0.r + dh10 * h * x0.v + dh01 * x1.r + dh11 * h * x1.v) / h;

    if (!finite_state(temp)) return StatusCode::non_finite_result;

    out = temp;
    return StatusCode::ok;
}

static StatusCode extrapolate_cartesian_ephemeris(
    const CartesianEphemerisTable& table,
    const EphemerisBracket& bracket,
    f64 dt,
    ExtrapolationMethod method,
    StateTr& out
) {
    if (bracket.location != EphemerisQueryLocation::before_coverage
        && bracket.location != EphemerisQueryLocation::after_coverage)
        return StatusCode::invalid_input;
    if (bracket.lower != bracket.upper) return StatusCode::invalid_input;

    StateTr temp{};
    i32 idx = bracket.lower;
    if (idx < 0 || idx >= table.dt.size()) return StatusCode::invalid_input;

    switch (method) {
    case ExtrapolationMethod::reject: return StatusCode::sample_not_found;
    case ExtrapolationMethod::hold: {
        temp = table.states[idx];
    } break;
    case ExtrapolationMethod::constant_velocity: {
        if (!table.has_velocity) return StatusCode::unsupported_method;

        temp = table.states[idx];
        f64 delta_t = dt - table.dt[idx];
        temp.r += delta_t * temp.v;
    } break;
    }

    if (!finite_state(temp)) return StatusCode::non_finite_result;

    out = temp;
    return StatusCode::ok;
}

StatusCode sample_cartesian_ephemeris(
    const CartesianEphemerisTable& table,
    f64 dt,
    StateTr& out,
    const StateSampleOptions& opts
) {

    StatusCode status;

    if (table.dt.empty() || table.states.empty()) return StatusCode::empty_ephemeris;
    if (table.dt.size() != table.states.size()) return StatusCode::size_mismatch;
    if (!std::isfinite(dt) || !finite_pos(opts.tol)) return StatusCode::invalid_input;

    EphemerisBracket bracket;
    status = find_ephemeris_bracket(table.dt, dt, opts.tol, bracket);
    if (status != StatusCode::ok) return status;
    i32 idx_l = bracket.lower;
    i32 idx_u = bracket.upper;

    StateTr temp{};
    switch (bracket.location) {
    case EphemerisQueryLocation::exact: {
        temp.r = table.states[idx_l].r;
        if (table.has_velocity) {
            temp.v = table.states[idx_l].v;
        }
    } break;
    case EphemerisQueryLocation::interior: {
        switch (opts.tr_interp) {
        case InterpolationMethod::nearest: {
            f64 delta_l = std::abs(dt - table.dt[idx_l]);
            f64 delta_u = std::abs(dt - table.dt[idx_u]);
            i32 idx;
            if (delta_l < delta_u)
                idx = idx_l;
            else if (delta_u < delta_l)
                idx = idx_u;
            else
                idx = idx_l;

            temp.r = table.states[idx].r;
            if (table.has_velocity) {
                temp.v = table.states[idx].v;
            }

        } break;
        case InterpolationMethod::linear: {
            status = interpolate_cartesian_linear(table, bracket, dt, opts.tol, temp);
            if (status != StatusCode::ok) return status;
        } break;
        case InterpolationMethod::slerp: {
            return StatusCode::unsupported_method;
        } break;
        case InterpolationMethod::cubic_hermite: {
            status
                = interpolate_cartesian_cubic_hermite(table, bracket, dt, opts.tol, temp);
            if (status != StatusCode::ok) return status;
        } break;
        }
    } break;
    case EphemerisQueryLocation::before_coverage: [[fallthrough]];
    case EphemerisQueryLocation::after_coverage: {
        status
            = extrapolate_cartesian_ephemeris(table, bracket, dt, opts.tr_extrap, temp);
        if (status != StatusCode::ok) return status;
    } break;
    }

    if (!finite_state(temp)) return StatusCode::non_finite_result;
    out = temp;
    return StatusCode::ok;
}

static StatusCode interpolate_orientation_slerp(
    const OrientationEphemerisTable& table,
    const EphemerisBracket& bracket,
    f64 dt,
    f64 tol,
    StateAtt& out
) {
    if (table.dt.size() != table.states.size()) return StatusCode::invalid_input;

    StatusCode status = validate_interior_ephemeris_bracket(bracket, table.dt.size());
    if (status != StatusCode::ok) return status;

    f64 t0 = table.dt[bracket.lower];
    f64 t1 = table.dt[bracket.upper];
    f64 h = t1 - t0;
    if (!finite_pos(h, tol)) return StatusCode::non_monotonic_time;

    f64 u = (dt - t0) / h;
    if (!std::isfinite(u) || u < -tol || u > 1.0 + tol) return StatusCode::invalid_input;
    u = std::clamp(u, 0.0, 1.0);

    const StateAtt& x0 = table.states[bracket.lower];
    const StateAtt& x1 = table.states[bracket.upper];
    if (!finite_state_att(x0, tol) || !finite_state_att(x1, tol))
        return StatusCode::invalid_att_state;

    vec4d q0 = x0.q.normalized();
    vec4d q1 = x1.q.normalized();
    f64 q_dot = q0.dot(q1);
    if (q_dot < 0.0) {
        q1 = -q1;
        q_dot = -q_dot;
    }
    q_dot = std::clamp(q_dot, -1.0, 1.0);

    StateAtt temp{};
    if (q_dot > 1.0 - std::max(tol * 1000.0, 1.0e-12)) {
        temp.q = ((1.0 - u) * q0 + u * q1).normalized();
    } else {
        f64 theta = std::acos(q_dot);
        f64 sin_theta = std::sin(theta);
        if (std::abs(sin_theta) <= tol) return StatusCode::interp_failed;

        f64 scale0 = std::sin((1.0 - u) * theta) / sin_theta;
        f64 scale1 = std::sin(u * theta) / sin_theta;
        temp.q = (scale0 * q0 + scale1 * q1).normalized();
    }

    if (table.has_angular_velocity) {
        temp.w = x0.w + u * (x1.w - x0.w);
    }

    if (!finite_state_att(temp, tol)) return StatusCode::non_finite_result;

    out = temp;
    return StatusCode::ok;
}

static StatusCode extrapolate_orientation_ephemeris(
    const OrientationEphemerisTable& table,
    const EphemerisBracket& bracket,
    f64 dt,
    ExtrapolationMethod method,
    StateAtt& out
) {
    if (bracket.location != EphemerisQueryLocation::before_coverage
        && bracket.location != EphemerisQueryLocation::after_coverage)
        return StatusCode::invalid_input;
    if (bracket.lower != bracket.upper) return StatusCode::invalid_input;

    i32 idx = bracket.lower;
    if (idx < 0 || idx >= table.dt.size() || idx >= table.states.size())
        return StatusCode::invalid_input;

    StateAtt temp{};
    switch (method) {
    case ExtrapolationMethod::reject: return StatusCode::sample_not_found;
    case ExtrapolationMethod::hold: {
        temp = table.states[idx];
    } break;
    case ExtrapolationMethod::constant_velocity: {
        if (!table.has_angular_velocity) return StatusCode::unsupported_method;

        temp = table.states[idx];
        temp.q = step_q_simple_spin(temp, dt - table.dt[idx]);
    } break;
    }

    if (!finite_state_att(temp)) return StatusCode::non_finite_result;
    normalize_quaternion_inplace(temp.q);

    out = temp;
    return StatusCode::ok;
}

StatusCode sample_orientation_ephemeris(
    const OrientationEphemerisTable& table,
    f64 dt,
    StateAtt& out,
    const StateSampleOptions& opts
) {
    if (table.dt.empty() || table.states.empty()) return StatusCode::empty_ephemeris;
    if (table.dt.size() != table.states.size()) return StatusCode::size_mismatch;
    if (!std::isfinite(dt) || !finite_pos(opts.tol)) return StatusCode::invalid_input;

    EphemerisBracket bracket;
    StatusCode status = find_ephemeris_bracket(table.dt, dt, opts.tol, bracket);
    if (status != StatusCode::ok) return status;

    i32 idx_l = bracket.lower;
    i32 idx_u = bracket.upper;
    StateAtt temp{};

    switch (bracket.location) {
    case EphemerisQueryLocation::exact: {
        temp = table.states[idx_l];
        if (!table.has_angular_velocity) temp.w = vec3d0;
    } break;
    case EphemerisQueryLocation::interior: {
        switch (opts.att_interp) {
        case InterpolationMethod::nearest: {
            f64 delta_l = std::abs(dt - table.dt[idx_l]);
            f64 delta_u = std::abs(dt - table.dt[idx_u]);
            i32 idx = delta_l <= delta_u ? idx_l : idx_u;
            temp = table.states[idx];
            if (!table.has_angular_velocity) temp.w = vec3d0;
        } break;
        case InterpolationMethod::slerp: {
            status = interpolate_orientation_slerp(table, bracket, dt, opts.tol, temp);
            if (status != StatusCode::ok) return status;
        } break;
        case InterpolationMethod::linear: [[fallthrough]];
        case InterpolationMethod::cubic_hermite: return StatusCode::unsupported_method;
        }
    } break;
    case EphemerisQueryLocation::before_coverage: [[fallthrough]];
    case EphemerisQueryLocation::after_coverage: {
        status
            = extrapolate_orientation_ephemeris(table, bracket, dt, opts.att_extrap, temp);
        if (status != StatusCode::ok) return status;
    } break;
    }

    if (!finite_state_att(temp, opts.tol)) return StatusCode::non_finite_result;
    normalize_quaternion_inplace(temp.q, opts.tol);

    out = temp;
    return StatusCode::ok;
}
