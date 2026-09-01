// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/body.hpp"
#include "core/dynamics_rotational.hpp"
#include "core/entity.hpp"
#include "core/estimation_common.hpp"
#include "core/interpolation.hpp"
#include "core/measurement.hpp"
#include "core/observation_type.hpp"
#include "core/od_dynamics.hpp"
#include "core/state.hpp"
#include "core/station_geometry.hpp"
#include "core/transform.hpp"
#include "core/world.hpp"
#include "core/world_history.hpp"
#include "util/units.hpp"

inline vec4d od_q_bcbf_from_inertial(const ODDynamicsConfig& dyn_cfg, f64 dt) {
    vec4d q_BCBF_I = dyn_cfg.q_cb0;
    if (dyn_cfg.update_body_attitude
        && dyn_cfg.att_model == ODAnchorAttModel::simple_spin) {
        StateAtt x_att_cb{.q = dyn_cfg.q_cb0, .w = dyn_cfg.w_cb};
        q_BCBF_I = step_q_simple_spin(x_att_cb, dt);
    }
    return q_BCBF_I;
}

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

inline StatusCode world_predict_measurement(
    const World& world,
    ObservationType type,
    EntityId observer_id,
    EntityId target_id,
    vecXd& z,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
) {
    const Station* observer = world.station(observer_id);
    const Body* target = world.body(target_id);
    if (observer == nullptr) {
        return StatusCode::observer_not_found;
    }
    if (target == nullptr) {
        return StatusCode::target_not_found;
    }

    z = world_predict_measurement(
        world,
        type,
        observer_id,
        target_id,
        angle_in,
        angle_out,
        tol
    );

    if (z.size() == 0) {
        return StatusCode::empty_measurements;
    }

    return StatusCode::ok;
}

inline StatusCode make_world_station_measurement_context(
    const World& world,
    MeasurementContext& ctx,
    EntityId observer_id,
    const StateTr& x_tr_target_pred,
    ObservationType type
) {
    ctx = MeasurementContext{};

    const Station* observer = world.station(observer_id);
    if (observer == nullptr) {
        return StatusCode::observer_not_found;
    }

    ctx.x_tr_target = x_tr_target_pred;

    if (type == ObservationType::azel) {
        if (!observer->anchored || observer->anchor_id == kInvalidEntityId) {
            return StatusCode::observer_not_found;
        }

        const Body* anchor = world.body(observer->anchor_id);
        if (anchor == nullptr) {
            return StatusCode::observer_not_found;
        }

        ctx.has_station_local = true;
        ctx.azel_frame = AzelInputFrame::enu;
        ctx.x_tr_observer.r = vec3d0;
        ctx.x_tr_target.r = world.stat_rel_enu(observer_id, x_tr_target_pred);
        ctx.station_llh = observer->llh_BCBF;
    } else {
        ctx.x_tr_observer = world.stat_x_tr_inertial(observer_id);
    }

    return StatusCode::ok;
}

inline StatusCode world_predict_measurement_from_state(
    const World& world,
    ObservationType type,
    EntityId observer_id,
    const StateTr& x_tr_target_pred,
    const ODDynamicsConfig& dyn_cfg,
    f64 dt,
    vecXd& z,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
) {
    MeasurementContext ctx;
    StatusCode status = make_world_station_measurement_context(
        world,
        ctx,
        observer_id,
        x_tr_target_pred,
        type
    );
    if (!od_status_success(status)) {
        return status;
    }

    if (type != ObservationType::azel) {
        z = predict_measurement(type, ctx, angle_in, angle_out, tol);
    } else {
        const Station* stat = world.station(observer_id);
        if (stat == nullptr) return StatusCode::observer_not_found;
        if (!stat->anchored || stat->anchor_id == kInvalidEntityId) {
            return StatusCode::observer_not_found;
        }

        const Body* anchor = world.body(stat->anchor_id);
        if (anchor == nullptr) return StatusCode::observer_not_found;

        vec4d q_BCBF_I = od_q_bcbf_from_inertial(dyn_cfg, dt);
        mat3d R_BCBF_I = ep_to_dcm(q_BCBF_I);
        mat3d R_ENU_BCBF = world.stat_rot_enu_from_body(observer_id);

        StateTr x_anchor = od_anchor_tr_at_time(dyn_cfg, dt);
        vec3d r_target_body_I = x_tr_target_pred.r - x_anchor.r;
        vec3d r_target_body_BCBF = R_BCBF_I * r_target_body_I;

        switch (ctx.azel_frame) {
        case AzelInputFrame::enu: {
            vec3d r_rel_BCBF = r_target_body_BCBF - stat->r_body_BCBF;
            vec3d r_rel_ENU = R_ENU_BCBF * r_rel_BCBF;
            ctx.x_tr_observer.r = vec3d0;
            ctx.x_tr_target.r = r_rel_ENU;
            ctx.station_llh = stat->llh_BCBF;
            ctx.has_station_local = true;
        } break;
        case AzelInputFrame::bcbf:
            ctx.x_tr_observer.r = stat->r_body_BCBF;
            ctx.x_tr_target.r = r_target_body_BCBF;
            ctx.station_llh = stat->llh_BCBF;
            ctx.has_station_local = true;
            break;
        }

        z = predict_measurement(type, ctx, angle_in, angle_out, tol);
    }

    if (z.size() == 0) {
        return StatusCode::empty_measurements;
    }

    return StatusCode::ok;
}

inline StatusCode world_predict_measurement_history(
    const World& world,
    const WorldHistory& history,
    ObservationType type,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    vecXd& z_pred,
    const StateSampleOptions& sample_opts = StateSampleOptions{},
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
) {
    StateTr x_tr_observer;
       StatusCode status
        = sample_tr_history(world, history, observer_id, t, x_tr_observer, sample_opts);
    if (!od_status_success(status)) return status;

    StateTr x_tr_target;
    status = sample_tr_history(world, history, target_id, t, x_tr_target, sample_opts);
    if (!od_status_success(status)) return status;

    MeasurementContext ctx;

    ctx.x_tr_target = x_tr_target;

    if (type == ObservationType::azel) {
        const Station* observer = world.station(observer_id);
        if (observer == nullptr) {
            return StatusCode::observer_not_found;
        }
        const Celestial* anchor = world.celestial(observer->anchor_id);
        if (anchor == nullptr) return StatusCode::observer_not_found;

        StateTr x_tr_anchor;
        status = sample_tr_history(world, history, anchor->id, t, x_tr_anchor, sample_opts);
        if (!od_status_success(status)) return status;

        StateAtt x_att_anchor;
        status = sample_att_history(world, history, anchor->id, t, x_att_anchor, sample_opts);
        if (!od_status_success(status)) return status;
        if (!observer->anchored || observer->anchor_id == kInvalidEntityId) {
            return StatusCode::observer_not_found;
        }

        vec3d r_target_body_I = x_tr_target.r - x_tr_anchor.r;
        vec3d r_target_body_BCBF
            = ep_rotate_fast_passive(x_att_anchor.q, r_target_body_I);

        vec3d r_rel_BCBF = r_target_body_BCBF - observer->r_body_BCBF;
        mat3d R_ENU_BCBF = stat_rot_enu_from_detic(observer->llh_BCBF, angle_in);

        ctx.has_station_local = true;
        ctx.azel_frame = AzelInputFrame::enu;
        ctx.x_tr_observer.r = vec3d0;
        ctx.x_tr_target.r = R_ENU_BCBF * r_rel_BCBF;
        ctx.station_llh = observer->llh_BCBF;
    } else {
        ctx.x_tr_observer = x_tr_observer;
    }

    z_pred = predict_measurement(type, ctx, angle_in, angle_out, tol);
    i32 dim = measurement_dim(type);
    if (z_pred.size() != dim) return StatusCode::size_mismatch;

    return StatusCode::ok;
}

inline StatusCode world_jacobian_measurement(
    const World& world,
    ObservationType type,
    EntityId observer_id,
    const StateTr& x_target_pred,
    const ODDynamicsConfig& dyn_cfg,
    f64 dt,
    matXd& H,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 eps_pos = 1e-3,
    f64 eps_vel = 1e-6,
    f64 tol = tol12
) {
    MeasurementContext ctx;
    StatusCode status = make_world_station_measurement_context(
        world,
        ctx,
        observer_id,
        x_target_pred,
        type
    );
    if (!od_status_success(status)) {
        return status;
    }

    i32 dim = measurement_dim(type);
    if (type != ObservationType::azel) {
        H = measurement_jacobian(type, ctx, angle_in, angle_out, eps_pos, eps_vel, tol);
    } else {
        switch (ctx.azel_frame) {
        case AzelInputFrame::enu: {
            vec4d q_BCBF_I = od_q_bcbf_from_inertial(dyn_cfg, dt);
            mat3d R_ENU_BCBF = world.stat_rot_enu_from_body(observer_id);
            mat3d R_BCBF_I = ep_to_dcm(q_BCBF_I);
            mat3d R_ENU_I = R_ENU_BCBF * R_BCBF_I;
            H = jacobian_azel_inertial_from_enu(ctx, R_ENU_I, angle_out, tol);
        } break;
        case AzelInputFrame::bcbf: {
            vec4d q_BCBF_I = od_q_bcbf_from_inertial(dyn_cfg, dt);
            mat3d R_BCBF_I = ep_to_dcm(q_BCBF_I);
            H = jacobian_azel_inertial_from_bcbf(ctx, R_BCBF_I, angle_in, angle_out, tol);
        } break;
        }
    }
    if (H.cols() != 6 || H.rows() != dim) {
        return StatusCode::empty_measurements;
    }

    return StatusCode::ok;
}
