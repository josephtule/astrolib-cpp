#pragma once

#include "core/entity.hpp"
#include "core/integrator.hpp"
#include "core/state.hpp"
#include "core/world.hpp"

struct WorldStepperConfig {
    IntegratorType integrator = IntegratorType::rk1;
    i32 substeps = 1;     // subdivisions per tick
    i32 ticks = 1;        // repeated integration ticks per call
    f64 time_scale = 1.0; // simulated-time multiplier applied to input dt
    bool step_translation = true;
    bool step_attitude = false;
};

struct WorldStepperStats {
    bool success = false;
    i32 ticks_completed;
    i32 substeps_completed;
    f64 dt_sim_advanced;
};

DerivTr derivtr_world(const World& world, EntityId id, const StateTr& x);

bool step_tr_world(World& world, EntityId id, f64 dt, const WorldStepperConfig& cfg);

WorldStepperStats step_world(World& world, f64 dt, const WorldStepperConfig& cfg);
