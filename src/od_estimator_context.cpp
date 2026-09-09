// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0
#include "core/od_estimator_context.hpp"
#include "od_environment_internal.hpp"
#include "core/station_geometry.hpp"
#include <algorithm>
#include <cmath>

StatusCode validate_od_estimator_context(const ODEstimatorContext& ctx) {
    if ((ctx.dynamics == nullptr) == (ctx.world == nullptr)) return StatusCode::invalid_input;
    if (!std::isfinite(ctx.fixed_step_size) || ctx.fixed_step_size < 0.0)
        return StatusCode::invalid_input;
    if (ctx.world && ctx.fixed_step_size != 0.0) return StatusCode::unsupported_method;
    if (ctx.world) return ctx.world->target_id() == kInvalidEntityId
        ? StatusCode::invalid_input : StatusCode::ok;
    switch (ctx.integrator) {
    case IntegratorTypeFixed::rk1: case IntegratorTypeFixed::rk2:
    case IntegratorTypeFixed::rk2_heun: case IntegratorTypeFixed::rk2_ralston:
    case IntegratorTypeFixed::rk3: case IntegratorTypeFixed::rk3_ralston:
    case IntegratorTypeFixed::rk4: case IntegratorTypeFixed::rk4_38:
    case IntegratorTypeFixed::rk5_nystrom: case IntegratorTypeFixed::rk6_butcher: break;
    default: return StatusCode::unsupported_method;
    }
    return validate_od_dynamics_context(*ctx.dynamics);
}

StatusCode propagate_od_estimator(
    const ODEstimatorContext& ctx, f64 t0, const StateTr& x0, f64 tf,
    i32 steps, bool reset_reference, VarStateTr& out
) {
    StatusCode status = validate_od_estimator_context(ctx);
    if (status != StatusCode::ok) return status;
    if (!std::isfinite(t0) || !std::isfinite(tf) || steps <= 0
        || !statetr_to_vec6d(x0).allFinite()) return StatusCode::invalid_input;
    if (ctx.world) {
        if (reset_reference) {
            if (t0 != ctx.world->reference_time()) return StatusCode::time_mismatch;
            status = ctx.world->reset(x0);
            if (status != StatusCode::ok) return status;
        } else {
            if (ctx.world->world().t_sim() != t0) return StatusCode::time_mismatch;
            ctx.world->world().body(ctx.world->target_id())->x_tr = x0;
        }
        auto options = ctx.options;
        // identical adaptive acceptance for linearization and line-search candidates
        options.with_stm = true;
        options.Phi0 = mat6d1;
        auto result = ctx.world->propagate(tf, options);
        if (result.has_stm) out = VarStateTr{.x = result.x, .Phi = result.Phi};
        return result.status;
    }
    VarStateTr y; y.x = x0; y.Phi = mat6d1;
    const f64 interval = tf - t0;
    if (!std::isfinite(interval)) return StatusCode::invalid_input;
    f64 dt = ctx.fixed_step_size > 0.0
        ? std::copysign(std::min(ctx.fixed_step_size, std::abs(interval)), interval)
        : interval / steps;
    if (!std::isfinite(dt) || (tf != t0 && t0 + dt == t0))
        return StatusCode::step_size_underflow;
    auto derivative = [&](f64 t, const VarStateTr& state) {
        VarDerivTr d{};
        if (status != StatusCode::ok) return d;
        status = od_environment_detail::derivative(*ctx.dynamics, t, state.x, d.dx);
        if (status != StatusCode::ok) return d;
        mat6d A = mat6d0; A.block<3,3>(0,3) = mat3d::Identity();
        for (const auto& source : ctx.dynamics->sources) {
            ODSourceState xs;
            status = od_environment_detail::source_state_at_time(*ctx.dynamics, source.id, t, xs);
            if (status != StatusCode::ok) return d;
            StateTr rel; rel.r = state.x.r - xs.x_tr.r;
            mat6d J = source.gravity_model == GravityModel::pointmass
                ? jacobian_tr_two_body(rel, source.mu)
                : jacobian_tr_zonal(rel, source.mu, source.ref_radius, source.degree, xs.x_att.q, source.J);
            A.block<3,3>(3,0) += J.block<3,3>(3,0);
        }
        d.dPhi = A * state.Phi;
        return d;
    };
    if (tf != t0) {
        f64 t = t0;
        i32 i = 0;
        while (ctx.fixed_step_size > 0.0 ? t != tf : i < steps) {
            // Keep the requested step; shorten only the step reaching tf.
            const f64 remaining = tf - t;
            const f64 h = ctx.fixed_step_size > 0.0
                ? std::copysign(std::min(ctx.fixed_step_size, std::abs(remaining)), remaining)
                : dt;
            if (t + h == t) return StatusCode::step_size_underflow;
            y = step_integrator<VarStateTr, VarDerivTr>(derivative, t, y, h, ctx.integrator).second;
            if (status != StatusCode::ok) return status;
            if (!statetr_to_vec6d(y.x).allFinite() || !y.Phi.allFinite())
                return StatusCode::non_finite_result;
            if (ctx.fixed_step_size > 0.0) t = std::abs(h) == std::abs(remaining) ? tf : t + h;
            else t = t0 + (++i) * dt;
        }
    }
    out = y;
    return StatusCode::ok;
}

StatusCode od_estimator_measurement_context(
    const ODEstimatorContext& ctx, const Measurement& measurement,
    const StateTr& target, const StateTr& observer, MeasurementContext& out,
    mat3d& R_measurement_I
) {
    StatusCode valid = validate_od_estimator_context(ctx);
    if (valid != StatusCode::ok) return valid;
    EntityId target_id = ctx.world ? ctx.world->target_id() : ctx.dynamics->propagation.target_id;
    if (measurement.target_id != kInvalidEntityId && measurement.target_id != target_id)
        return StatusCode::invalid_input;
    MeasurementContext temp = make_measurement_context(target, observer);
    mat3d R = mat3d::Identity();
    const World* geometry = ctx.observer_geometry;
    if (!geometry && ctx.world) geometry = &ctx.world->world();
    const Station* station = geometry ? geometry->station(measurement.observer_id) : nullptr;
    if (station && station->anchored) {
        ODSourceState anchor;
        StatusCode status;
        if (ctx.world) {
            if (ctx.world->world().t_sim() != measurement.t) return StatusCode::time_mismatch;
            const Body* body = ctx.world->world().body(station->anchor_id);
            if (!body) return StatusCode::body_not_found;
            if (body->ephemeris_providers.orientation
                && !body->ephemeris_providers.orientation->table.has_angular_velocity)
                return StatusCode::unsupported_method;
            anchor.x_tr = body->x_tr; anchor.x_att = body->x_att;
            if (body->body_type == BodyType::celestial
                && static_cast<const Celestial*>(body)->attitude_model == CelestialAttitudeModel::fixed
                && !body->ephemeris_providers.orientation) anchor.x_att.w = vec3d0;
            temp.x_tr_observer.r = anchor.x_tr.r + ep_rotate_fast_passive(ep_conj(anchor.x_att.q), station->r_body_BCBF);
            temp.x_tr_observer.v = anchor.x_tr.v + ep_rotate_fast_passive(ep_conj(anchor.x_att.q), anchor.x_att.w.cross(station->r_body_BCBF));
        } else {
            status = od_source_state_at_time(*ctx.dynamics, station->anchor_id, measurement.t, anchor);
            if (status != StatusCode::ok) return status;
            status = od_anchored_observer_state_at_time(*ctx.dynamics, station->anchor_id,
                station->r_body_BCBF, measurement.t, temp.x_tr_observer);
            if (status != StatusCode::ok) return status;
        }
        if (measurement.type == ObservationType::azel) {
            R = geometry->stat_rot_enu_from_body(station->id) * ep_to_dcm(anchor.x_att.q);
            temp.x_tr_target.r = R * (target.r-temp.x_tr_observer.r);
            temp.x_tr_target.v = vec3d0;
            temp.x_tr_observer = StateTr{};
            temp.has_station_local = true;
            temp.azel_frame = AzelInputFrame::enu;
        }
    } else if (measurement.type == ObservationType::azel) return StatusCode::unsupported_type;
    out = temp; R_measurement_I = R;
    return StatusCode::ok;
}
