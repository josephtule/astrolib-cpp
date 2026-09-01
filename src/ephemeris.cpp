// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/ephemeris.hpp"
#include "core/state.hpp"
#include "core/status.hpp"
#include "core/time.hpp"
#include "util/math.hpp"
#include "util/units.hpp"

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