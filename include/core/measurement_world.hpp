#pragma once

#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/estimation_common.hpp"
#include "core/measurement.hpp"
#include "core/state.hpp"
#include "core/world.hpp"
#include "util/units.hpp"

inline vecXd world_predict_measurement(
    const World& world,
    ObservationType type,
    EntityId observer_id,
    EntityId target_id,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
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
        ctx.x_tr_observer.r = vec3d0;
        ctx.x_tr_target.r = r_rel_ENU;
        ctx.azel_frame = AzelInputFrame::enu;
        z = predict_measurement(type, ctx, angle_in, angle_out, tol);
    }

    return z;
}

inline ODStatus make_world_station_measurement_context(
    const World& world,
    MeasurementContext& ctx,
    EntityId observer_id,
    const StateTr& x_tr_target_pred,
    ObservationType type
) {
    ctx = MeasurementContext{};

    const Station* observer = world.station(observer_id);
    if (observer == nullptr) {
        return ODStatus::observer_not_found;
    }

    ctx.x_tr_target = x_tr_target_pred;

    if (type == ObservationType::azel) {

        if (!observer->anchored || observer->anchor_id == kInvalidEntityId) {
            return ODStatus::observer_not_found;
        }

        const Body* anchor = world.body(observer->anchor_id);
        if (anchor == nullptr) {
            return ODStatus::observer_not_found;
        }

        ctx.has_station_local = true;
        ctx.azel_frame = AzelInputFrame::enu;
        ctx.x_tr_observer.r = vec3d0;
        ctx.x_tr_target.r = world.stat_rel_enu(observer_id, x_tr_target_pred);
        ctx.station_llh = observer->llh_BCBF;
    } else {
        ctx.x_tr_observer = world.stat_x_tr_inertial(observer_id);
    }

    return ODStatus::ok;
}
