#pragma once

#include "core/entity.hpp"
#include "core/estimation_common.hpp"
#include "core/interpolation.hpp"
#include "core/state.hpp"
#include "core/world.hpp"
#include "util/typedefs.hpp"
#include <cstddef>

struct WorldHistorySample {
    // separate from WorldStateSnapshot, only for
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

inline ODStatus sample_nearest(
    const WorldHistory& history,
    f64 t,
    WorldHistorySample& out,
    f64 tol = tol12
) {
    if (history.samples.empty()) {
        return ODStatus::empty_history;
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

    out = history.samples[closest_i];

    return ODStatus::ok;
}

inline ODStatus sample_tr_nearest(
    const WorldHistory& history,
    EntityId id,
    f64 t,
    StateTr& out,
    f64 tol = tol12
) {
    WorldHistorySample sample;
    ODStatus status = sample_nearest(history, t, sample, tol);
    if (!od_status_success(status)) {
        return status;
    }

    auto it = sample.x_tr.find(id);
    if (it == sample.x_tr.end()) {
        return ODStatus::sample_not_found;
    }
    out = it->second;

    return ODStatus::ok;
}

inline ODStatus sample_att_nearest(
    const WorldHistory& history,
    EntityId id,
    f64 t,
    StateAtt& out,
    f64 tol = tol12
) {
    WorldHistorySample sample;
    ODStatus status = sample_nearest(history, t, sample, tol);
    if (!od_status_success(status)) {
        return status;
    }

    auto it = sample.x_att.find(id);
    if (it == sample.x_att.end()) {
        return ODStatus::sample_not_found;
    }
    out = it->second;

    return ODStatus::ok;
}

inline ODStatus sample_tr_interp_linear(
    f64 t,
    EntityId id,
    const WorldHistorySample& sample1,
    const WorldHistorySample& sample2,
    StateTr& out
) {
    // check valid id
    auto it = sample1.x_tr.find(id);
    if (it == sample1.x_tr.end()) {
        return ODStatus::sample_not_found;
    }
    StateTr x_tr_1 = it->second;

    it = sample2.x_tr.find(id);
    if (it == sample2.x_tr.end()) {
        return ODStatus::sample_not_found;
    }
    StateTr x_tr_2 = it->second;

    vec6d x_tr_1_vec = statetr_to_vec6d(x_tr_1);
    vec6d x_tr_2_vec = statetr_to_vec6d(x_tr_2);

    vec6d x_tr_t = interp_linear(t, x_tr_1_vec, sample1.t, x_tr_2_vec, sample2.t);
    out = vec6d_to_statetr(x_tr_t);

    return ODStatus::ok;
}

inline ODStatus sample_att_interp_linear(
    f64 t,
    EntityId id,
    const WorldHistorySample& sample1,
    const WorldHistorySample& sample2,
    StateAtt& out,
    f64 tol = tol12
) {
    // check valid id
    auto it = sample1.x_att.find(id);
    if (it == sample1.x_att.end()) {
        return ODStatus::sample_not_found;
    }
    StateAtt x_att_1 = it->second;

    it = sample2.x_att.find(id);
    if (it == sample2.x_att.end()) {
        return ODStatus::sample_not_found;
    }
    StateAtt x_att_2 = it->second;

    vec4d q_t = interp_quat_linear(t, x_att_1.q, sample1.t, x_att_2.q, sample2.t, tol);
    out.q = q_t;

    vec3d w_t = interp_linear(t, x_att_1.w, sample1.t, x_att_2.w, sample2.t);
    out.w = w_t;

    return ODStatus::ok;
}