#pragma once

#include "core/estimation_common.hpp"
#include "core/estimation_recursive.hpp"
#include "core/world.hpp"

struct ODWorldMeasurementEvent {
    Measurement measurement;
    EntityId observer_id = kInvalidEntityId;
    EntityId target_id = kInvalidEntityId;
};

ODEKFStepResult od_ekf_step_world(
    const World& world,
    const ODEKFState& filter,
    const ODWorldMeasurementEvent& event,
    const ODDynamicsConfig& dyn_config,
    i32 prop_steps,
    const mat6d& Q,
    f64 tol_time = tol12
);

ODStatus ekf_observer_state_from_world(
    const World& world,
    EntityId observer_id,
    StateTr& x_tr_observer
);

struct ODRealtimeEKFInput {
    const World* world;
     ODEKFState filter;
    const ODWorldMeasurementEvent* event;
    f64 t_target = 0.0;
    ODDynamicsConfig dyn_config;
    i32 prop_steps = 100;
    mat6d Q = mat6d0;
    f64 tol_time = tol12;
};

ODEKFStepResult od_ekf_predict_step(
    const ODEKFState& filter,
    f64 t_target,
    const ODDynamicsConfig& dyn_config,
    i32 prop_steps,
    const mat6d& Q,
    f64 tol_time = tol12
);

ODEKFStepResult od_ekf_update_world(const ODRealtimeEKFInput& input);
