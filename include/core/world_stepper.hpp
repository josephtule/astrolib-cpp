// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/entity.hpp"
#include "core/integrator_adaptive.hpp"
#include "core/integrator_fixed.hpp"
#include "core/state.hpp"
#include "core/world.hpp"

struct WorldStepperConfig {
    IntegratorTypeFixed integrator_tr = IntegratorTypeFixed::rk2;
    IntegratorTypeFixed integrator_att = IntegratorTypeFixed::rk4;
    i32 substeps = 1;   // subdivisions per tick
    i32 ticks = 1;      // repeated integration ticks per call
    f64 dt_scale = 1.0; // simulated-time multiplier applied to input dt
    bool step_tr = true;
    bool step_att = true;
    bool paused = false;
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

struct WorldStepperStats {
    bool success = true;
    i32 ticks_completed = 0;
    i32 substeps_completed = 0;
    f64 dt_sim_advanced = 0.0;
};
inline WorldStepperStats operator+(
    const WorldStepperStats& stats1,
    const WorldStepperStats& stats2
) {
    return WorldStepperStats{
        .success = stats1.success && stats2.success,
        .ticks_completed = stats1.ticks_completed + stats2.ticks_completed,
        .substeps_completed = stats1.substeps_completed + stats2.substeps_completed,
        .dt_sim_advanced = stats1.dt_sim_advanced + stats2.dt_sim_advanced
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

WorldStepperStats step_world(World& world, f64 dt, const WorldStepperConfig& cfg);

WorldStepperStats step_world(
    World& world,
    f64 dt,
    const WorldStepperConfig& cfg,
    WorldStepperWorkspace& wksp
);

WorldStepperStats step_world_tableau(
    World& world,
    f64 dt,
    const WorldStepperConfig& cfg
);

WorldStepperStats step_world_tableau(
    World& world,
    f64 dt,
    const WorldStepperConfig& cfg,
    WorldStepperWorkspace& wksp
);

struct WorldAdaptiveIntegratorConfig {
    IntegratorTypeAdaptive integrator_tr = IntegratorTypeAdaptive::dopri54;
    // IntegratorTypeAdaptive integrator_tr = IntegratorTypeAdaptive::bosha32;
    IntegratorTypeAdaptive integrator_att = IntegratorTypeAdaptive::dopri54;

    AdaptiveIntegratorConfig opts{};
    bool use_substeps = false;
};

struct WorldAdaptiveStepResult {
    StatusCode status = StatusCode::invalid_input;

    f64 t = 0.0;
    f64 dt_sim_advanced = 0.0;

    f64 final_error_norm = 0.0;

    i32 ticks_completed = 0;
    i32 substeps_completed = 0;
    AdaptiveIntegratorStats stats{};
};

WorldAdaptiveStepResult step_world_adaptive(
    World& world,
    f64 dt,
    const WorldStepperConfig& stepper_cfg,
    const WorldAdaptiveIntegratorConfig& adaptive_cfg,
    WorldStepperWorkspace& wksp
);

WorldAdaptiveStepResult step_world_adaptive(
    World& world,
    f64 dt,
    const WorldStepperConfig& stepper_cfg,
    const WorldAdaptiveIntegratorConfig& adaptive_cfg
);
