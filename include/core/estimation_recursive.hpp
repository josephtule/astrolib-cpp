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
    StateTr observer_state;
    ODDynamicsConfig dyn_config;
    i32 prop_steps = 100;
    mat6d Q = mat6d0; // process noise covariance (process uncertainty)
    f64 tol_time = tol_strict;
};

struct ODEKFStepResult {
    bool success = false;
    ODStatus status = ODStatus::invalid_input;
    ODEKFState filter;
    vecXd residual;
    f64 residual_norm = 0.0;
    f64 raw_residual_norm = 0.0;
};

// offline ekf
struct ODEKFInput {
    ODEKFState initial_filter;
    svec<Measurement> measurements;
    svec<StateTr> observer_states;
    ODDynamicsConfig dyn_config;
    i32 prop_steps = 100;
    mat6d Q = mat6d0;
    f64 tol_time = tol_strict;
};

struct ODEKFResult {
    bool success = false;
    ODStatus status = ODStatus::invalid_input;
    ODEKFState filter;
    i32 processed_measurements = 0;
    f64 residual_norm = 0.0;
    f64 raw_residual_norm = 0.0;
};

ODStatus od_ekf_step_validate_input(const ODEKFStepInput& input);
ODStatus od_ekf_validate_input(const ODEKFInput& input);
ODEKFStepResult od_ekf_step(const ODEKFStepInput& input);
ODEKFResult od_ekf_offline(const ODEKFInput& input);
