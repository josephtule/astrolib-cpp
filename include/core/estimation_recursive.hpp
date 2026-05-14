#pragma once

#include "core/estimation_common.hpp"
#include "core/measurement.hpp"
#include "core/od_dynamics.hpp"
#include "core/state.hpp"
#include "util/constants.hpp"

struct ODEKFState {
    StateTr x; // state
    mat6d P;   // state estimate covariance
    f64 t = 0.0;
};

// step used for real-time kf, also used by offline kf
struct ODEKFStepInput {
    ODEKFState filter;
    Measurement measurement;
    StateTr x_tr_observer;
    ODDynamicsConfig dyn_config;
    i32 prop_steps = 100;
    mat6d Q = mat6d0; // process noise covariance (process uncertainty)
    f64 tol_time = tol12;
};

struct ODEKFStepResult {
    ODStatus status = ODStatus::invalid_input;
    ODEKFState filter;
    vecXd residual;
    f64 residual_norm = 0.0;
    f64 raw_residual_norm = 0.0;
};

// offline ekf
struct ODEKFOfflineInput {
    ODEKFState initial_filter;
    svec<Measurement> measurements;
    svec<StateTr> observer_states;
    ODDynamicsConfig dyn_config;
    i32 prop_steps = 100;
    mat6d Q = mat6d0;
    f64 tol_time = tol12;
};

struct ODEKFResult {
    ODStatus status = ODStatus::invalid_input;
    ODEKFState filter;
    i32 processed_measurements = 0;
    f64 residual_norm = 0.0;
    f64 raw_residual_norm = 0.0;
};

ODStatus od_ekf_step_validate_input(const ODEKFStepInput& input);
ODStatus od_ekf_validate_input(const ODEKFOfflineInput& input);
ODEKFStepResult od_ekf_step(const ODEKFStepInput& input);
ODEKFResult od_ekf_offline(const ODEKFOfflineInput& input);

struct ODEKFPredictResult {
    VarStateTr y;
    mat6d P = mat6d1;
    f64 t = 0.0;
    ODStatus status = ODStatus::invalid_input;
};

ODEKFPredictResult od_ekf_predict(
    const ODEKFState& filter,
    f64 t_target,
    const ODDynamicsConfig& dyn_config,
    i32 prop_steps,
    const mat6d& Q,
    f64 tol = tol12
);

