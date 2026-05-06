#pragma once

#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/measurement.hpp"
#include "core/state.hpp"
#include "core/world.hpp"
#include "util/units.hpp"

// TODO: split between header/implementation if this file grows

inline vecXd world_predict_measurement(
    const World& world,
    ObservationType type,
    EntityId observer_id,
    EntityId target_id,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol_strict
) {
    const Station* observer = world.station(observer_id);
    const Body* target = world.body(target_id);
    if (observer == nullptr || target == nullptr) return vecXd{};
    if (type == ObservationType::azel
        && (observer->anchored == false || observer->anchor_id == kInvalidEntityId)) {
        return vecXd{};
    }

    vecXd z;
    if (type != ObservationType::azel) {
        z = predict_measurement(
            type,
            target->x_tr,
            world.stat_x_tr_inertial(observer_id), // inertial relative position
            angle_out,
            tol
        );
    } else {
        // azel in enu
        vec3d r_rel_ENU = world.stat_rel_enu(observer_id, target_id);
        MeasurementContext ctx;
        ctx.has_station_local = true;
        ctx.station_llh = observer->llh_BCBF;
        ctx.x_observer.r = vec3d0;
        ctx.x_target.r = r_rel_ENU;
        ctx.azel_frame = AzelInputFrame::enu;
        z = predict_measurement(type, ctx, angle_in, angle_out, tol);
    }

    return z;
}