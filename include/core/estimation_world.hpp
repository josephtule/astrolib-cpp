#pragma once

#include "core/estimation_common.hpp"
#include "core/estimation_recursive.hpp"
#include "core/measurement.hpp"
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
    f64 tol_time = tol12,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 eps_pos = 1e-3,
    f64 eps_vel = 1e-6,
    f64 tol = tol12
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

ODStatus make_world_measurement_event(
    const World& world,
    ObservationType type,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    const matXd& R,
    ODWorldMeasurementEvent& event,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
);

struct ODRealtimeEvent {
    f64 t = 0.0;
    bool has_measurement = false;
    ODWorldMeasurementEvent event;
};

struct ODRealtimeEKFResult {
    ODStatus status = ODStatus::invalid_input;
    ODEKFState filter;
    i32 processed_events = 0;
    i32 measurement_updates = 0;
    i32 prediction_updates = 0;
    f64 residual_norm = 0.0;
    f64 raw_residual_norm = 0.0;
};

ODRealtimeEKFResult od_ekf_update_world_events(
    const World& world,
    const ODEKFState& initial_filter,
    const svec<ODRealtimeEvent>& events,
    const ODDynamicsConfig& dyn_config,
    i32 prop_steps,
    const mat6d& Q,
    f64 tol_time = tol12
);

ODStatus validate_realtime_ekf_events(
    const svec<ODRealtimeEvent>& events,
    f64 t_prev,
    f64 tol_time = tol12
);

struct ODRealtimeScheduleItem {
    InstrumentId instrument_id = kInvalidInstrumentId;
    f64 t = 0.0;
    bool has_measurement = false;
    ObservationType type = ObservationType::range;
    EntityId observer_id = kInvalidEntityId;
    EntityId target_id = kInvalidEntityId;
};

ODStatus make_world_measurement_event(
    const World& world,
    ObservationType type,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    ODWorldMeasurementEvent& event,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
);

ODStatus make_world_measurement_event_instrument(
    const World& world,
    InstrumentId instrument_id,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    ODWorldMeasurementEvent& event,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
);

ODStatus make_noisy_world_measurement_event_instrument(
    const World& world,
    InstrumentId instrument_id,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    ODWorldMeasurementEvent& event,
    const MeasurementNoiseOptions& noise_opts,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
);

