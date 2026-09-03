// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/od_environment.hpp"
#include "od_environment_internal.hpp"
#include <algorithm>
#include <cmath>

StatusCode validate_od_dynamics_context(const ODDynamicsContext& ctx) {
    const auto& cfg = ctx.propagation;
    switch (cfg.mode) {
    case ODPropagationMode::lightweight:
    case ODPropagationMode::provider: break;
    case ODPropagationMode::full_world: return StatusCode::unsupported_method;
    default: return StatusCode::invalid_input;
    }
    if (cfg.target_id == kInvalidEntityId || !std::isfinite(cfg.t0))
        return StatusCode::invalid_input;
    if (cfg.source_ids.size() != ctx.sources.size()) return StatusCode::size_mismatch;
    for (size_t i = 0; i < cfg.source_ids.size(); ++i) {
        EntityId id = cfg.source_ids[i];
        if (id == kInvalidEntityId || id == cfg.target_id
            || std::find(cfg.source_ids.begin(), cfg.source_ids.begin() + i, id)
                != cfg.source_ids.begin() + i) return StatusCode::invalid_input;
        auto count = std::count_if(ctx.sources.begin(), ctx.sources.end(),
            [id](const ODSourceConfig& source) { return source.id == id; });
        if (count != 1) return StatusCode::invalid_input;
    }
    for (const auto& source : ctx.sources) {
        if (!source.x_tr0.r.allFinite() || !source.x_tr0.v.allFinite())
            return StatusCode::invalid_state;
        if (!source.x_att0.q.allFinite() || !source.x_att0.w.allFinite()
            || std::abs(source.x_att0.q.norm() - 1.0) > 1e-8)
            return StatusCode::invalid_att_state;
        if (source.tr_model != ODAnchorTrModel::fixed
            && source.tr_model != ODAnchorTrModel::constant_velocity)
            return StatusCode::unsupported_method;
        if (source.att_model != ODAnchorAttModel::fixed
            && source.att_model != ODAnchorAttModel::simple_spin)
            return StatusCode::unsupported_method;
        if (cfg.mode != ODPropagationMode::provider) continue;
        if (source.providers.translation || source.providers.orientation) {
            if (!ctx.has_epoch || !std::isfinite(ctx.epoch.day)
                || !std::isfinite(ctx.epoch.frac) || ctx.inertial_frame.empty()
                || source.object.empty()) return StatusCode::invalid_input;
            switch (ctx.time_scale) {
            case TimeScale::utc: case TimeScale::ut1: case TimeScale::tai:
            case TimeScale::tt: case TimeScale::tdb: case TimeScale::gps: break;
            default: return StatusCode::invalid_input;
            }
            if (!std::isfinite(ctx.offsets.tai_utc) || !std::isfinite(ctx.offsets.ut1_utc)
                || !std::isfinite(ctx.offsets.tt_tai) || !std::isfinite(ctx.offsets.tai_gps))
                return StatusCode::invalid_input;
        }
        if (source.providers.translation) {
            const auto& provider = *source.providers.translation;
            StatusCode status = validate_ephemeris_provider(provider);
            if (status != StatusCode::ok) return status;
            const auto& frame = provider.table.metadata.frame;
            if (ctx.center.empty() || frame.center != ctx.center
                || frame.frame != ctx.inertial_frame || frame.object != source.object)
                return StatusCode::invalid_input;
        }
        if (source.providers.orientation) {
            const auto& provider = *source.providers.orientation;
            StatusCode status = validate_orientation_provider(provider);
            if (status != StatusCode::ok) return status;
            const auto& frame = provider.table.metadata.frame;
            if (source.body_frame.empty() || frame.source_frame != ctx.inertial_frame
                || frame.target_frame != source.body_frame || frame.object != source.object)
                return StatusCode::invalid_input;
        }
    }
    return StatusCode::ok;
}

StatusCode od_environment_detail::source_state_at_time(
    const ODDynamicsContext& ctx,
    EntityId source_id,
    f64 t,
    ODSourceState& out
) {
    if (!std::isfinite(t)) return StatusCode::invalid_input;
    StatusCode status;
    auto it = std::find_if(ctx.sources.begin(), ctx.sources.end(),
        [source_id](const ODSourceConfig& source) { return source.id == source_id; });
    if (it == ctx.sources.end()) return StatusCode::body_not_found;
    const auto& source = *it;
    f64 dt = t - ctx.propagation.t0;
    if (!std::isfinite(dt)) return StatusCode::non_finite_result;
    ODSourceState temp{source.x_tr0, source.x_att0};
    bool provider_mode = ctx.propagation.mode == ODPropagationMode::provider;
    JulianDate epoch = ctx.epoch;
    if (provider_mode && (source.providers.translation || source.providers.orientation)) {
        epoch.frac += dt / 86400.0;
        epoch = normalize_jd(epoch);
    }
    if (provider_mode && source.providers.translation) {
        status = query_ephemeris_provider(
            *source.providers.translation, epoch, ctx.time_scale, ctx.offsets, temp.x_tr
        );
        if (status != StatusCode::ok) return status;
    } else if (source.tr_model == ODAnchorTrModel::constant_velocity) {
        temp.x_tr.r += dt * temp.x_tr.v;
    }
    if (provider_mode && source.providers.orientation) {
        status = query_orientation_provider(
            *source.providers.orientation, epoch, ctx.time_scale, ctx.offsets, temp.x_att
        );
        if (status != StatusCode::ok) return status;
    } else if (source.att_model == ODAnchorAttModel::simple_spin) {
        temp.x_att.q = step_q_simple_spin(source.x_att0, dt);
    }
    if (!temp.x_tr.r.allFinite() || !temp.x_tr.v.allFinite()
        || !temp.x_att.q.allFinite() || !temp.x_att.w.allFinite())
        return StatusCode::non_finite_result;
    out = temp;
    return StatusCode::ok;
}

StatusCode od_environment_detail::derivative(
    const ODDynamicsContext& ctx, f64 t, const StateTr& x_target_I, DerivTr& out
) {
    if (!x_target_I.r.allFinite() || !x_target_I.v.allFinite() || !std::isfinite(t))
        return StatusCode::invalid_state;
    StatusCode status;
    DerivTr temp;
    temp.dr = x_target_I.v;
    temp.dv = vec3d0;
    for (EntityId id : ctx.propagation.source_ids) {
        const auto& source = *std::find_if(ctx.sources.begin(), ctx.sources.end(),
            [id](const ODSourceConfig& s) { return s.id == id; });
        if (!std::isfinite(source.mu) || source.mu <= 0.0) return StatusCode::invalid_input;
        if (source.gravity_model != GravityModel::pointmass
            && source.gravity_model != GravityModel::zonal) return StatusCode::unsupported_method;
        if (source.gravity_model == GravityModel::zonal
            && (!std::isfinite(source.ref_radius) || source.ref_radius <= 0.0
                || source.degree < 0 || source.degree > 6 || !source.J.allFinite()))
            return StatusCode::invalid_input;
        ODSourceState state;
        status = source_state_at_time(ctx, id, t, state);
        if (status != StatusCode::ok) return status;
        vec3d r_rel_I = x_target_I.r - state.x_tr.r;
        if (!r_rel_I.allFinite() || r_rel_I.norm() <= tol12) return StatusCode::invalid_state;
        if (source.gravity_model == GravityModel::pointmass) {
            temp.dv += accel_gravity_pointmass(r_rel_I, source.mu);
        } else {
            vec3d r_rel_B = ep_rotate_fast_passive(state.x_att.q, r_rel_I);
            vec3d a_B = accel_gravity_zonal(
                r_rel_B, source.mu, source.ref_radius, source.degree, source.J
            );
            temp.dv += ep_rotate_fast_passive(ep_conj(state.x_att.q), a_B);
        }
    }
    if (!temp.dv.allFinite()) return StatusCode::non_finite_result;
    out = temp;
    return StatusCode::ok;
}

StatusCode od_anchored_observer_state_at_time(
    const ODDynamicsContext& ctx, EntityId anchor_id, const vec3d& r_observer_body_B,
    f64 t, StateTr& out
) {
    if (!r_observer_body_B.allFinite()) return StatusCode::invalid_state;
    ODSourceState anchor;
    StatusCode status = od_source_state_at_time(ctx, anchor_id, t, anchor);
    if (status != StatusCode::ok) return status;
    StateTr temp;
    vec4d q_IB = ep_conj(anchor.x_att.q);
    temp.r = anchor.x_tr.r + ep_rotate_fast_passive(q_IB, r_observer_body_B);
    // fixed orientation contributes no rotation velocity
    const auto& source = *std::find_if(ctx.sources.begin(), ctx.sources.end(),
        [anchor_id](const ODSourceConfig& s) { return s.id == anchor_id; });
    bool rotating = source.att_model == ODAnchorAttModel::simple_spin
        || (ctx.propagation.mode == ODPropagationMode::provider && source.providers.orientation);
    if (ctx.propagation.mode == ODPropagationMode::provider && source.providers.orientation
        && !source.providers.orientation->table.has_angular_velocity)
        return StatusCode::unsupported_method;
    temp.v = anchor.x_tr.v;
    if (rotating) temp.v += ep_rotate_fast_passive(q_IB, anchor.x_att.w.cross(r_observer_body_B));
    if (!temp.r.allFinite() || !temp.v.allFinite()) return StatusCode::non_finite_result;
    out = temp;
    return StatusCode::ok;
}

StatusCode od_source_state_at_time(
    const ODDynamicsContext& ctx, EntityId id, f64 t, ODSourceState& out
) {
    StatusCode status = validate_od_dynamics_context(ctx);
    if (status != StatusCode::ok) return status;
    return od_environment_detail::source_state_at_time(ctx, id, t, out);
}

StatusCode derivtr_od_context(
    const ODDynamicsContext& ctx, f64 t, const StateTr& x, DerivTr& out
) {
    StatusCode status = validate_od_dynamics_context(ctx);
    if (status != StatusCode::ok) return status;
    return od_environment_detail::derivative(ctx, t, x, out);
}
