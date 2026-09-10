// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/measurement_world.hpp"
#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/observations.hpp"
#include "core/state.hpp"
#include "core/status.hpp"
#include <cmath>

StatusCode resolve_observer_measurement_context(
    const World& world,
    EntityId observer_id,
    f64 t,
    ObserverMeasurementContext& out,
    f64 tol_time
) {
    if (observer_id == kInvalidEntityId) return StatusCode::observer_not_found;
    if (!std::isfinite(t) || !std::isfinite(tol_time) || tol_time < 0.0
        || !std::isfinite(world.t_sim())) {
        return StatusCode::invalid_input;
    }

    // TEMP: only allow measurements aligned with world time steps
    if (std::abs(t - world.t_sim()) > tol_time) {
        return StatusCode::time_mismatch;
    }

    const Body* observer = world.body(observer_id);
    if (observer == nullptr) return StatusCode::observer_not_found;

    if (!world.is_active(observer_id)) return StatusCode::inactive_entity;

    ObserverMeasurementContext temp;
    temp.platform_id = observer_id;
    temp.platform_type = observer->body_type;
    switch (observer->body_type) {
    default: return StatusCode::unsupported_type;
    case BodyType::unknown: [[fallthrough]];
    case BodyType::celestial: return StatusCode::unsupported_type;
    case BodyType::satellite: {
        const Satellite* sat = world.satellite(observer_id);
        if (sat == nullptr) return StatusCode::observer_not_found;

        temp.x_tr_observer_I = sat->x_tr;
        temp.x_att_observer = sat->x_att;

        if (!finite_state(temp.x_tr_observer_I)) return StatusCode::invalid_state;
        if (!finite_state(temp.x_att_observer)) return StatusCode::invalid_att_state;
    } break;
    case BodyType::station: {
        const Station* stat = world.station(observer_id);
        if (stat == nullptr) return StatusCode::observer_not_found;

        if (stat->anchored) {
            const Celestial* anchor = world.celestial(stat->anchor_id);
            if (anchor == nullptr) return StatusCode::anchor_not_found;
            temp.x_tr_observer_I = world.stat_x_tr_inertial(observer_id);
            temp.x_att_observer = world.stat_x_att_inertial(observer_id);
        } else {
            temp.x_tr_observer_I = stat->x_tr;
            temp.x_att_observer = stat->x_att;
        }

        if (!finite_state(temp.x_tr_observer_I)) return StatusCode::invalid_state;
        if (!finite_state(temp.x_att_observer)) return StatusCode::invalid_att_state;
    } break;
    }

    // TEMP: see note in struct definition
    temp.r_sensor_I = temp.x_tr_observer_I.r;
    temp.q_S_I = temp.x_att_observer.q;

    out = temp;
    return StatusCode::ok;
}

StatusCode resolve_instrument_measurement_context(
    const World& world,
    EntityId observer_id,
    InstrumentId instrument_id,
    f64 t,
    ObserverMeasurementContext& out,
    f64 tol_time
) {
    ObserverMeasurementContext temp;
    StatusCode status = resolve_observer_measurement_context(
        world, observer_id, t, temp, tol_time
    );
    if (status != StatusCode::ok) return status;

    status = resolve_platform_instrument(world, observer_id, instrument_id, temp.instrument);
    if (status != StatusCode::ok) return status;
    temp.instrument_id = instrument_id;
    out = temp;
    return StatusCode::ok;
}

StatusCode resolve_platform_instrument(
    const World& world,
    EntityId observer_id,
    InstrumentId instrument_id,
    const PlatformInstrument*& out
) {
    const Body* body = world.body(observer_id);
    if (body == nullptr) return StatusCode::observer_not_found;
    if (!world.is_active(observer_id)) return StatusCode::inactive_entity;
    const InstrumentSuite* suite = nullptr;
    StatusCode status = instrument_suite_from_body(*body, suite);
    if (status != StatusCode::ok) return status;
    auto it = suite->instruments.find(instrument_id);
    if (it == suite->instruments.end()) return StatusCode::instrument_not_found;
    if (!it->second.enabled) return StatusCode::instrument_disabled;
    out = &it->second;
    return StatusCode::ok;
}

StatusCode predict_relative_angular_observation(
    const ObserverMeasurementContext& observer,
    const StateTr& reference_target,
    const StateTr& target,
    vec2d& out,
    f64 tol_range,
    f64 tol_pole
) {
    if (!finite_state(observer.x_tr_observer_I) || !finite_state(reference_target)
        || !finite_state(target)) return StatusCode::invalid_state;
    vec3d r_target_observer_I = target.r - observer.x_tr_observer_I.r;
    vec3d r_reference_observer_I = reference_target.r - observer.x_tr_observer_I.r;
    return delta_radec_from_rel(
        r_target_observer_I, r_reference_observer_I, out,
        UAngle::radian, tol_range, tol_pole
    );
}

StatusCode world_predict_relative_angular_observation(
    const World& world,
    EntityId observer_id,
    EntityId reference_target_id,
    EntityId target_id,
    f64 t,
    RelativeAngularObservation& out,
    f64 tol_time,
    f64 tol_range,
    f64 tol_pole
) {
    if (observer_id == target_id || observer_id == reference_target_id) {
        return StatusCode::invalid_input;
    }
    ObserverMeasurementContext observer;
    StatusCode status = resolve_observer_measurement_context(world, observer_id, t, observer, tol_time);
    if (status != StatusCode::ok) return status;
    const Body* target = world.body(target_id);
    const Body* reference = world.body(reference_target_id);
    if (target == nullptr) return StatusCode::target_not_found;
    if (reference == nullptr) return StatusCode::missing_reference;
    if (!world.is_active(target_id) || !world.is_active(reference_target_id)) {
        return StatusCode::inactive_entity;
    }

    // anchored station targets also need their derived world positions
    auto target_state = [&](const Body* body) {
        return body->body_type == BodyType::station
            ? world.stat_x_tr_inertial(body->id) : body->x_tr;
    };
    for (const Body* body : {target, reference}) {
        if (body->body_type == BodyType::station) {
            const Station* station = world.station(body->id);
            if (station->anchored && world.celestial(station->anchor_id) == nullptr) {
                return StatusCode::anchor_not_found;
            }
        }
    }
    RelativeAngularObservation temp = out;
    status = predict_relative_angular_observation(
        observer, target_state(reference), target_state(target), temp.delta_radec,
        tol_range, tol_pole
    );
    if (status != StatusCode::ok) return status;
    temp.t = t;
    temp.observer_id = observer_id;
    temp.reference_target_id = reference_target_id;
    temp.target_id = target_id;
    out = temp;
    return StatusCode::ok;
}

StatusCode relative_angular_residual(
    const vec2d& observed,
    const vec2d& predicted,
    vec2d& out
) {
    if (!observed.allFinite() || !predicted.allFinite()) return StatusCode::invalid_input;
    vec2d temp = observed - predicted;
    if (!temp.allFinite()) return StatusCode::non_finite_result;
    temp(0) = wrap_angle(temp(0), -pi, pi);
    out = temp;
    return StatusCode::ok;
}

StatusCode jacobian_relative_angular_observation(
    const ObserverMeasurementContext& observer,
    const StateTr& reference_target,
    const StateTr& target,
    RelativeAngularJacobians& out,
    f64 tol_range,
    f64 tol_pole
) {
    vec2d angles;
    StatusCode status = predict_relative_angular_observation(
        observer, reference_target, target, angles, tol_range, tol_pole
    );
    if (status != StatusCode::ok) return status;
    // raw wrapped RA is discontinuous at +/- pi
    if (pi - std::abs(angles(0)) <= tol_pole) return StatusCode::invalid_state;
    MeasurementContext ctx;
    ctx.x_tr_observer = observer.x_tr_observer_I;
    ctx.x_tr_target = target;
    // geometry is already checked with separate range and angular tolerances
    matXd H_target = jacobian_radec(ctx, UAngle::radian, 0.0);
    ctx.x_tr_target = reference_target;
    matXd H_reference = jacobian_radec(ctx, UAngle::radian, 0.0);
    if (H_target.rows() != 2 || H_target.cols() != 6
        || H_reference.rows() != 2 || H_reference.cols() != 6) {
        return StatusCode::invalid_state;
    }
    RelativeAngularJacobians temp;
    temp.target = H_target;
    temp.reference = -H_reference;
    temp.observer = H_reference - H_target;
    if (!temp.target.allFinite() || !temp.reference.allFinite() || !temp.observer.allFinite()) {
        return StatusCode::non_finite_result;
    }
    out = temp;
    return StatusCode::ok;
}
