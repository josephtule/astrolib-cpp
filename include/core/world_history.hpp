// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/estimation_common.hpp"
#include "core/interpolation.hpp"
#include "core/state.hpp"
#include "core/station_geometry.hpp"
#include "core/status.hpp"
#include "core/transform.hpp"
#include "core/world.hpp"
#include "util/typedefs.hpp"
#include <cmath>
#include <cstdlib>

// TODO: change _linear to general, create interpolation type

struct WorldHistorySample {
    // separate from WorldStateSnapshot, only for OD
    f64 t;
    umap<EntityId, StateTr> x_tr;
    umap<EntityId, StateAtt> x_att;
};

struct WorldHistory {
    dque<WorldHistorySample> samples;
    i32 max_samples = 2;
};

inline WorldHistorySample capture_world_history_sample(const World& world) {
    WorldHistorySample sample;

    sample.t = world.t_sim();
    for (EntityId id : world.active_entity_ids()) {
        const Body* body = world.body(id);
        if (body == nullptr) {
            continue;
        }
        sample.x_tr[id] = body->x_tr;
        sample.x_att[id] = body->x_att;
    }

    return sample;
}

inline void push_world_history_sample(
    WorldHistory& history,
    const WorldHistorySample& sample
) {
    history.samples.push_back(sample);
    while (history.samples.size() > history.max_samples) {
        history.samples.pop_front();
    }
}

inline void clear_world_history(WorldHistory& history) { history.samples.clear(); }

inline StatusCode sample_nearest_idx(
    const WorldHistory& history,
    f64 t,
    i32& out,
    f64 tol = tol12
) {
    if (history.samples.empty()) {
        return StatusCode::empty_history;
    }

    f64 closest_dt = INFINITY;
    i32 closest_i = 0;
    for (i32 i = 0; i < history.samples.size(); ++i) {
        f64 dt = std::abs(history.samples[i].t - t);

        if (dt <= tol) {
            closest_dt = dt;
            closest_i = i;
            break;
        }

        if (dt < closest_dt) {
            closest_dt = dt;
            closest_i = i;
        } else {
            // increasing dt, ordered by time so break
            break;
        }
    }

    out = closest_i;

    return StatusCode::ok;
}

inline StatusCode sample_nearest(
    const WorldHistory& history,
    f64 t,
    WorldHistorySample& out,
    f64 tol = tol12
) {
    i32 idx;

    auto status = sample_nearest_idx(history, t, idx, tol);
    if (!od_status_success(status)) {
        return status;
    }

    out = history.samples[idx];

    return StatusCode::ok;
}

inline StatusCode sample_tr_nearest(
    const WorldHistory& history,
    EntityId id,
    f64 t,
    StateTr& out,
    f64 tol = tol12
) {
    WorldHistorySample sample;
    StatusCode status = sample_nearest(history, t, sample, tol);
    if (!od_status_success(status)) {
        return status;
    }

    auto it = sample.x_tr.find(id);
    if (it == sample.x_tr.end()) {
        return StatusCode::sample_not_found;
    }
    out = it->second;

    return StatusCode::ok;
}

inline StatusCode sample_station_tr_nearest(
    const World& world,
    const WorldHistory& history,
    EntityId station_id,
    f64 t,
    StateTr& out,
    f64 tol = tol12
) {
    WorldHistorySample sample;
    StatusCode status = sample_nearest(history, t, sample, tol);
    if (!od_status_success(status)) {
        return status;
    }

    const Station* stat = world.station(station_id);
    if (!stat->anchored) {
        return sample_tr_nearest(history, station_id, t, out, tol);
    }

    EntityId anchor_id = stat->anchor_id;

    const auto it_tr = sample.x_tr.find(anchor_id);
    if (it_tr == sample.x_tr.end()) {
        return StatusCode::sample_not_found;
    }
    const auto it_att = sample.x_att.find(anchor_id);
    if (it_att == sample.x_att.end()) {
        return StatusCode::sample_not_found;
    }

    StateTr x_tr_anchor = it_tr->second;
    StateAtt x_att_anchor = it_att->second;

    vec4d q_NB = ep_conj(x_att_anchor.q);
    vec3d r_offset_inertial = ep_rotate_fast_passive(q_NB, stat->r_body_BCBF);
    out.r = x_tr_anchor.r + r_offset_inertial;
    vec3d w_inertial = ep_rotate_fast_passive(q_NB, x_att_anchor.w);
    out.v = x_tr_anchor.v + w_inertial.cross(r_offset_inertial);

    return StatusCode::ok;
}

inline StatusCode sample_att_nearest(
    const WorldHistory& history,
    EntityId id,
    f64 t,
    StateAtt& out,
    f64 tol = tol12
) {
    WorldHistorySample sample;
    StatusCode status = sample_nearest(history, t, sample, tol);
    if (!od_status_success(status)) {
        return status;
    }

    auto it = sample.x_att.find(id);
    if (it == sample.x_att.end()) {
        return StatusCode::sample_not_found;
    }
    out = it->second;

    return StatusCode::ok;
}

inline StatusCode sample_station_att_nearest(
    const World& world,
    const WorldHistory& history,
    EntityId station_id,
    f64 t,
    StateAtt& out,
    f64 tol = tol12
) {
    WorldHistorySample sample;
    StatusCode status = sample_nearest(history, t, sample, tol);
    if (!od_status_success(status)) {
        return status;
    }

    const Station* stat = world.station(station_id);
    if (!stat->anchored) {
        return sample_att_nearest(history, station_id, t, out, tol);
    }
    EntityId anchor_id = stat->anchor_id;

    const auto it_tr = sample.x_tr.find(anchor_id);
    if (it_tr == sample.x_tr.end()) {
        return StatusCode::sample_not_found;
    }
    const auto it_att = sample.x_att.find(anchor_id);
    if (it_att == sample.x_att.end()) {
        return StatusCode::sample_not_found;
    }

    StateAtt x_att_anchor;
    status = sample_att_nearest(history, anchor_id, t, x_att_anchor, tol);
    if (!od_status_success(status)) return status;

    out.q = stat_att_enu_from_detic(x_att_anchor, stat->llh_BCBF);

    return StatusCode::ok;
}

inline StatusCode sample_tr_interp_linear(
    const WorldHistorySample& sample1,
    const WorldHistorySample& sample2,
    EntityId id,
    f64 t,
    StateTr& out,
    f64 tol = tol12
) {
    // check valid id
    auto it = sample1.x_tr.find(id);
    if (it == sample1.x_tr.end()) {
        return StatusCode::sample_not_found;
    }
    StateTr x_tr_1 = it->second;

    it = sample2.x_tr.find(id);
    if (it == sample2.x_tr.end()) {
        return StatusCode::sample_not_found;
    }
    StateTr x_tr_2 = it->second;

    vec6d x_tr_1_vec = statetr_to_vec6d(x_tr_1);
    vec6d x_tr_2_vec = statetr_to_vec6d(x_tr_2);

    if (std::abs(sample1.t - sample2.t) <= tol) {
        if (std::abs(sample1.t - t) > tol) {
            return StatusCode::time_mismatch;
        }
        out = x_tr_1;
        return StatusCode::ok;
    }

    vecXd x_tr_t = interp_linear(t, x_tr_1_vec, sample1.t, x_tr_2_vec, sample2.t, tol);
    if (x_tr_t.size() != 6) {
        return StatusCode::interp_failed;
    }
    out = vec6d_to_statetr(vec6d(x_tr_t));

    return StatusCode::ok;
}

inline StatusCode sample_att_interp_linear(
    const WorldHistorySample& sample1,
    const WorldHistorySample& sample2,
    EntityId id,
    f64 t,
    StateAtt& out,
    f64 tol = tol12
) {
    // check valid id
    auto it = sample1.x_att.find(id);
    if (it == sample1.x_att.end()) {
        return StatusCode::sample_not_found;
    }
    StateAtt x_att_1 = it->second;

    it = sample2.x_att.find(id);
    if (it == sample2.x_att.end()) {
        return StatusCode::sample_not_found;
    }
    StateAtt x_att_2 = it->second;

    if (std::abs(sample1.t - sample2.t) <= tol) {
        if (std::abs(sample1.t - t) > tol) {
            return StatusCode::time_mismatch;
        }
        out = x_att_1;
        return StatusCode::ok;
    }

    vec4d q_t = interp_quat_linear(t, x_att_1.q, sample1.t, x_att_2.q, sample2.t, tol);
    if (q_t.norm() <= tol || !q_t.allFinite()) {
        return StatusCode::interp_failed;
    }
    out.q = q_t;

    vecXd w_t = interp_linear(t, x_att_1.w, sample1.t, x_att_2.w, sample2.t);
    if (w_t.size() != 3) {
        return StatusCode::interp_failed;
    }
    out.w = vec3d(w_t);

    return StatusCode::ok;
}

inline StatusCode sample_bracket(
    const WorldHistory& history,
    f64 t,
    WorldHistorySample& sample1,
    WorldHistorySample& sample2,
    f64 tol = tol12
) {
    const dque<WorldHistorySample>& samples = history.samples;

    if (samples.empty()) {
        return StatusCode::empty_history;
    }

    if (samples.size() < 2) {
        sample1 = samples[0];
        sample2 = samples[0];
        if (std::abs(sample1.t - t) <= tol) {
            return StatusCode::ok;
        }
        return StatusCode::time_mismatch;
    }

    bool found = false;
    for (i32 i = 0; i < static_cast<i32>(samples.size()) - 1; ++i) {
        // assumes increasing time by index
        const WorldHistorySample& s1 = samples[i];
        const WorldHistorySample& s2 = samples[i + 1];

        if (s1.t <= t + tol && t <= s2.t + tol) {
            sample1 = s1;
            sample2 = s2;
            found = true;
            break;
        }
    }

    if (found) {
        return StatusCode::ok;
    }
    return StatusCode::time_mismatch;
}

inline StatusCode sample_tr_interp_linear(
    const WorldHistory& history,
    EntityId id,
    f64 t,
    StateTr& out,
    f64 tol = tol12
) {
    WorldHistorySample sample1;
    WorldHistorySample sample2;
    auto status = sample_bracket(history, t, sample1, sample2, tol);
    if (!od_status_success(status)) {
        return status;
    }

    status = sample_tr_interp_linear(sample1, sample2, id, t, out, tol);
    if (!od_status_success(status)) {
        return status;
    }

    return StatusCode::ok;
}

inline StatusCode sample_att_interp_linear(
    const WorldHistory& history,
    EntityId id,
    f64 t,
    StateAtt& out,
    f64 tol = tol12
) {
    WorldHistorySample sample1;
    WorldHistorySample sample2;
    auto status = sample_bracket(history, t, sample1, sample2, tol);
    if (!od_status_success(status)) {
        return status;
    }

    status = sample_att_interp_linear(sample1, sample2, id, t, out, tol);
    if (!od_status_success(status)) {
        return status;
    }

    return StatusCode::ok;
}

inline StatusCode sample_station_tr_interp_linear(
    const World& world,
    const WorldHistory& history,
    EntityId station_id,
    f64 t,
    StateTr& out,
    f64 tol = tol12
) {
    const Station* stat = world.station(station_id);
    if (stat == nullptr) {
        return StatusCode::observer_not_found;
    }

    if (!stat->anchored) {
        return sample_tr_interp_linear(history, station_id, t, out, tol);
    }

    EntityId anchor_id = stat->anchor_id;
    StateTr x_tr_anchor;
    StatusCode status = sample_tr_interp_linear(history, anchor_id, t, x_tr_anchor, tol);
    if (!od_status_success(status)) return status;
    StateAtt x_att_anchor;
    status = sample_att_interp_linear(history, anchor_id, t, x_att_anchor, tol);
    if (!od_status_success(status)) return status;

    vec4d q_NB = ep_conj(x_att_anchor.q);
    vec3d r_offset_inertial = ep_rotate_fast_passive(q_NB, stat->r_body_BCBF);
    out.r = x_tr_anchor.r + r_offset_inertial;
    vec3d w_inertial = ep_rotate_fast_passive(q_NB, x_att_anchor.w);
    out.v = x_tr_anchor.v + w_inertial.cross(r_offset_inertial);

    return StatusCode::ok;
}

inline StatusCode sample_station_att_interp_linear(
    const World& world,
    const WorldHistory& history,
    EntityId station_id,
    f64 t,
    StateAtt& out,
    f64 tol = tol12
) {
    const Station* stat = world.station(station_id);
    if (stat == nullptr) {
        return StatusCode::observer_not_found;
    }

    if (!stat->anchored) {
        return sample_att_interp_linear(history, station_id, t, out, tol);
    }

    EntityId anchor_id = stat->anchor_id;
    StateAtt x_att_anchor;
    StatusCode status
        = sample_att_interp_linear(history, anchor_id, t, x_att_anchor, tol);
    if (!od_status_success(status)) return status;

    out.q = stat_att_enu_from_detic(x_att_anchor, stat->llh_BCBF);

    return StatusCode::ok;
}

inline StatusCode sample_tr_history(
    const World& world,
    const WorldHistory& history,
    EntityId id,
    f64 t,
    StateTr& out,
    const CartesianSampleOptions& opts
) {
    StatusCode status = StatusCode::ok;
    const Body* body = world.body(id);
    if (body == nullptr) {
        return StatusCode::body_not_found;
    }
    switch (opts.interpolation) {
    case CartesianInterpolationMethod::nearest: {
        if (body->body_type == BodyType::station) {
            return sample_station_tr_nearest(world, history, id, t, out, opts.tol);
        }
        return sample_tr_nearest(history, id, t, out, opts.tol);
    } break;
    case CartesianInterpolationMethod::linear: {
        if (body->body_type == BodyType::station) {
            return sample_station_tr_interp_linear(world, history, id, t, out, opts.tol);
        }
        return sample_tr_interp_linear(history, id, t, out, opts.tol);
    } break;
    case CartesianInterpolationMethod::cubic_hermite:
        return StatusCode::unsupported_method;
    }

    return status;
}

inline StatusCode sample_att_history(
    const World& world,
    const WorldHistory& history,
    EntityId id,
    f64 t,
    StateAtt& out,
    const OrientationSampleOptions& opts
) {
    StatusCode status = StatusCode::ok;

    const Body* body = world.body(id);
    if (body == nullptr) {
        return StatusCode::body_not_found;
    }

    switch (opts.interpolation) {
    case OrientationInterpolationMethod::nearest: {
        if (body->body_type == BodyType::station) {
            status = sample_station_att_nearest(world, history, id, t, out, opts.tol);
            return status;
        }
        return sample_att_nearest(history, id, t, out, opts.tol);
    } break;
    case OrientationInterpolationMethod::slerp: {
        if (body->body_type == BodyType::station) {
            status
                = sample_station_att_interp_linear(world, history, id, t, out, opts.tol);
            return status;
        }
        status = sample_att_interp_linear(history, id, t, out, opts.tol);
        return status;
    } break;
    }

    return StatusCode::ok;
}

struct WorldHistoryProvider {
    const World* world = nullptr;
    const WorldHistory* history = nullptr;
};

inline StatusCode validate_world_history_provider(const WorldHistoryProvider& provider) {
    if (provider.world == nullptr) {
        return StatusCode::invalid_input;
    }
    if (provider.history == nullptr || provider.history->samples.empty()) {
        return StatusCode::empty_history;
    }

    return StatusCode::ok;
}

inline StatusCode provider_body_tr(
    const WorldHistoryProvider& provider,
    EntityId body_id,
    f64 t,
    StateTr& out,
    f64 tol = tol12
) {
    StatusCode status = validate_world_history_provider(provider);
    if (!od_status_success(status)) {
        return status;
    }

    status = sample_tr_interp_linear(*provider.history, body_id, t, out, tol);
    if (!od_status_success(status)) {
        return status;
    }

    return StatusCode::ok;
}

inline StatusCode provider_body_att(
    const WorldHistoryProvider& provider,
    EntityId body_id,
    f64 t,
    StateAtt& out,
    f64 tol = tol12
) {
    StatusCode status = validate_world_history_provider(provider);
    if (!od_status_success(status)) {
        return status;
    }

    status = sample_att_interp_linear(*provider.history, body_id, t, out, tol);
    if (!od_status_success(status)) {
        return status;
    }

    return StatusCode::ok;
}

inline StatusCode provider_station_tr(
    const WorldHistoryProvider& provider,
    EntityId station_id,
    f64 t,
    StateTr& out,
    f64 tol = tol12
) {
    StatusCode status = validate_world_history_provider(provider);
    if (!od_status_success(status)) {
        return status;
    }

    status = sample_station_tr_interp_linear(
        *provider.world,
        *provider.history,
        station_id,
        t,
        out,
        tol
    );
    if (!od_status_success(status)) {
        return status;
    }

    return StatusCode::ok;
}

inline StatusCode provider_station_att(
    const WorldHistoryProvider& provider,
    EntityId station_id,
    f64 t,
    StateAtt& out,
    f64 tol = tol12
) {
    StatusCode status = validate_world_history_provider(provider);
    if (!od_status_success(status)) {
        return status;
    }

    status = sample_station_att_interp_linear(
        *provider.world,
        *provider.history,
        station_id,
        t,
        out,
        tol
    );
    if (!od_status_success(status)) {
        return status;
    }

    return StatusCode::ok;
}
