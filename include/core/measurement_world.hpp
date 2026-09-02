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

struct ObserverMeasurementContext {
    EntityId platform_id = kInvalidEntityId;
    InstrumentId instrument_id = kInvalidInstrumentId;

    BodyType platform_type = BodyType::unknown;

    // inertial state of observer
    StateTr x_tr_observer_I;
    StateAtt x_att_observer;

    // inertial sensor state initial same as host, later allow offset/mounting
    vec3d r_sensor_I = vec3d0;
    vec4d q_S_I = q_identity;

    // borrowed pointer; do not retain this context across instrument edits
    const PlatformInstrument* instrument = nullptr;
};

struct RelativeAngularObservation {
    f64 t = 0.0; // observation epoch
    EntityId observer_id = kInvalidEntityId;
    EntityId reference_target_id = kInvalidEntityId;
    EntityId target_id = kInvalidEntityId;
    // target minus reference, radians; not delta_ra * cos(dec)
    vec2d delta_radec = vec2d0;
    // supplied differential covariance, rad^2; prediction leaves this unchanged
    // from absolute angles: R = R_target + R_reference - C_tr - C_tr.transpose()
    mat2d R = mat2d0;
};

struct RelativeAngularJacobians {
    matd<2, 6> target = matd<2, 6>::Zero();
    matd<2, 6> reference = matd<2, 6>::Zero();
    matd<2, 6> observer = matd<2, 6>::Zero();
};

// explicit states share an epoch and inertial axes; output is radians
StatusCode predict_relative_angular_observation(
    const ObserverMeasurementContext& observer,
    const StateTr& reference_target,
    const StateTr& target,
    vec2d& out,
    f64 tol_range = tol12,
    f64 tol_pole = tol12
);

StatusCode world_predict_relative_angular_observation(
    const World& world,
    EntityId observer_id,
    EntityId reference_target_id,
    EntityId target_id,
    f64 t,
    RelativeAngularObservation& out,
    f64 tol_time = tol12,
    f64 tol_range = tol12,
    f64 tol_pole = tol12
);

StatusCode relative_angular_residual(
    const vec2d& observed,
    const vec2d& predicted,
    vec2d& out
);

StatusCode jacobian_relative_angular_observation(
    const ObserverMeasurementContext& observer,
    const StateTr& reference_target,
    const StateTr& target,
    RelativeAngularJacobians& out,
    f64 tol_range = tol12,
    f64 tol_pole = tol12
);

StatusCode resolve_observer_measurement_context(
    const World& world,
    EntityId observer_id,
    f64 t,
    ObserverMeasurementContext& out,
    f64 tol_time = tol12
);

StatusCode resolve_instrument_measurement_context(
    const World& world,
    EntityId observer_id,
    InstrumentId instrument_id,
    f64 t,
    ObserverMeasurementContext& out,
    f64 tol_time = tol12
);

// Ownership lookup only; history callers resolve the state at the requested epoch.
StatusCode resolve_platform_instrument(
    const World& world,
    EntityId observer_id,
    InstrumentId instrument_id,
    const PlatformInstrument*& out
);

inline vec4d od_q_bcbf_from_inertial(const ODDynamicsConfig& dyn_cfg, f64 dt) {
    vec4d q_BCBF_I = dyn_cfg.q_cb0;
    if (dyn_cfg.update_body_attitude
        && dyn_cfg.att_model == ODAnchorAttModel::simple_spin) {
        StateAtt x_att_cb{.q = dyn_cfg.q_cb0, .w = dyn_cfg.w_cb};
        q_BCBF_I = step_q_simple_spin(x_att_cb, dt);
    }
    return q_BCBF_I;
}

inline StatusCode make_world_measurement_context(
    const World& world,
    MeasurementContext& ctx,
    EntityId observer_id,
    const StateTr& x_tr_target_pred,
    ObservationType type
) {
    ctx = MeasurementContext{};

    ObserverMeasurementContext observer_ctx;
    StatusCode status = resolve_observer_measurement_context(
        world,
        observer_id,
        world.t_sim(),
        observer_ctx
    );
    if (status != StatusCode::ok) return status;
    if (!finite_state(x_tr_target_pred)) return StatusCode::invalid_state;
    const Station* observer = world.station(observer_id);

    ctx.x_tr_target = x_tr_target_pred;

    if (type == ObservationType::azel) {
        if (observer == nullptr) return StatusCode::unsupported_type;
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
        ctx.x_tr_observer = observer_ctx.x_tr_observer_I;
    }

    return StatusCode::ok;
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
    if (observer_id == target_id) return StatusCode::invalid_input;
    const Body* target = world.body(target_id);
    if (target == nullptr) return StatusCode::target_not_found;
    if (!world.is_active(target_id)) return StatusCode::inactive_entity;
    if (measurement_dim(type) <= 0) return StatusCode::unsupported_type;
    MeasurementContext ctx;
    StatusCode status
        = make_world_measurement_context(world, ctx, observer_id, target->x_tr, type);
    if (status != StatusCode::ok) return status;
    vecXd temp = predict_measurement(type, ctx, angle_in, angle_out, tol);
    if (temp.size() != measurement_dim(type)) return StatusCode::size_mismatch;
    if (!temp.allFinite()) return StatusCode::non_finite_result;
    z = std::move(temp);
    return StatusCode::ok;
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
    vecXd z;
    world_predict_measurement(
        world,
        type,
        observer_id,
        target_id,
        z,
        angle_in,
        angle_out,
        tol
    );
    return z;
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
    StatusCode status
        = make_world_measurement_context(world, ctx, observer_id, x_tr_target_pred, type);
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
    if (observer_id == target_id) return StatusCode::invalid_input;
    const Body* platform = world.body(observer_id);
    if (platform == nullptr) return StatusCode::observer_not_found;
    if (world.body(target_id) == nullptr) return StatusCode::target_not_found;
    if (!world.is_active(observer_id) || !world.is_active(target_id)) {
        return StatusCode::inactive_entity;
    }
    if (platform->body_type != BodyType::station
        && platform->body_type != BodyType::satellite) {
        return StatusCode::unsupported_type;
    }
    if (measurement_dim(type) <= 0) return StatusCode::unsupported_type;
    if (type == ObservationType::azel) {
        const Station* station = world.station(observer_id);
        if (station == nullptr || !station->anchored) return StatusCode::unsupported_type;
        if (world.celestial(station->anchor_id) == nullptr)
            return StatusCode::anchor_not_found;
    }

    StateTr x_tr_observer;
    StatusCode status = sample_tr_history(
        world,
        history,
        observer_id,
        t,
        x_tr_observer,
        sample_opts.translation
    );
    if (!od_status_success(status)) return status;

    StateTr x_tr_target;
    status = sample_tr_history(
        world,
        history,
        target_id,
        t,
        x_tr_target,
        sample_opts.translation
    );
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
        status = sample_tr_history(
            world,
            history,
            anchor->id,
            t,
            x_tr_anchor,
            sample_opts.translation
        );
        if (!od_status_success(status)) return status;

        StateAtt x_att_anchor;
        status = sample_att_history(
            world,
            history,
            anchor->id,
            t,
            x_att_anchor,
            sample_opts.orientation
        );
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
    StatusCode status
        = make_world_measurement_context(world, ctx, observer_id, x_target_pred, type);
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
