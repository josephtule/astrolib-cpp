// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/entity.hpp"
#include "core/integrator_adaptive.hpp"
#include "core/integrator_fixed.hpp"
#include "core/state.hpp"
#include "core/world.hpp"

#include <algorithm>

struct WorldAdaptiveConfig {
    AdaptiveIntegratorConfig opts{};
    bool use_substeps = false;
};

struct WorldStepperConfig {
    IntegratorType integrator_tr{IntegratorTypeFixed::rk4};
    IntegratorType integrator_att{IntegratorTypeFixed::rk4};

    i32 substeps = 1;   // subdivisions per tick
    i32 ticks = 1;      // repeated integration ticks per call
    f64 dt_scale = 1.0; // simulated-time multiplier applied to input dt
    bool step_tr = true;
    bool step_att = true;
    bool paused = false;

    WorldAdaptiveConfig adaptive{};
};

struct WorldStepperWorkspace {
    svec<EntityId> propagated_tr_ids;  // all propagated tr
    svec<EntityId> propagated_att_ids; // satellites, free-stations
    svec<EntityId> celestial_att_ids;  // celestials with simple-spin/provided
    svec<EntityId> gravity_source_ids; // gravity sources
    svec<EntityId> source_att_ids;     // sources requiring attitude
    svec<EntityId> staged_att_ids; // union of propagated_att_ids and celestial_att_ids
    bool dirty = true;
};

void rebuild_world_stepper_workspace(
    const World& world,
    WorldStepperWorkspace& workspace
);

inline void invalidate_stepper_wksp(WorldStepperWorkspace& wksp) {
    wksp.dirty = true;
}

struct WorldStepperStats {
    i32 ticks_completed = 0;
    i32 substeps_completed = 0;
    f64 dt_sim_advanced = 0.0;
    AdaptiveIntegratorStats adaptive{};
};

struct WorldStepResult {
    StatusCode status = StatusCode::invalid_state;
    f64 t = 0.0;
    f64 final_error_norm = 0.0;
    WorldStepperStats stats{};
};

inline WorldStepperStats operator+(
    const WorldStepperStats& stats1,
    const WorldStepperStats& stats2
) {
    AdaptiveIntegratorStats adaptive = stats1.adaptive;
    i64 accepted_before = adaptive.accepted_steps;

    adaptive.attempted_steps += stats2.adaptive.attempted_steps;
    adaptive.accepted_steps += stats2.adaptive.accepted_steps;
    adaptive.rejected_steps += stats2.adaptive.rejected_steps;
    adaptive.deriv_evals += stats2.adaptive.deriv_evals;

    if (stats2.adaptive.accepted_steps > 0) {
        if (accepted_before == 0) {
            adaptive.min_accepted_dt = stats2.adaptive.min_accepted_dt;
            adaptive.max_accepted_dt = stats2.adaptive.max_accepted_dt;
        } else {
            adaptive.min_accepted_dt
                = std::min(adaptive.min_accepted_dt, stats2.adaptive.min_accepted_dt);
            adaptive.max_accepted_dt
                = std::max(adaptive.max_accepted_dt, stats2.adaptive.max_accepted_dt);
        }
        adaptive.final_accepted_dt = stats2.adaptive.final_accepted_dt;
    }

    return WorldStepperStats{
        .ticks_completed = stats1.ticks_completed + stats2.ticks_completed,
        .substeps_completed = stats1.substeps_completed + stats2.substeps_completed,
        .dt_sim_advanced = stats1.dt_sim_advanced + stats2.dt_sim_advanced,
        .adaptive = adaptive
    };
}
inline WorldStepperStats& operator+=(
    WorldStepperStats& stats1,
    const WorldStepperStats& stats2
) {
    stats1 = stats1 + stats2;
    return stats1;
}

DerivTr derivtr_world(const World& world, EntityId id, const StateTr& x);

bool step_tr_world(World& world, EntityId id, f64 dt, const WorldStepperConfig& cfg);

WorldStepResult step_world_legacy(
    World& world,
    f64 dt,
    const WorldStepperConfig& cfg
);

WorldStepResult step_world_legacy(
    World& world,
    f64 dt,
    const WorldStepperConfig& cfg,
    WorldStepperWorkspace& wksp
);


WorldStepResult step_world(
    World& world,
    f64 dt,
    const WorldStepperConfig& cfg,
    WorldStepperWorkspace& wksp
);

WorldStepResult step_world(
    World& world,
    f64 dt,
    const WorldStepperConfig& cfg
);