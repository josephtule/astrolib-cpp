// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/measurement_world.hpp"
#include "core/body.hpp"
#include "core/entity.hpp"
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
