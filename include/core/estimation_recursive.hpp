#pragma once

#include "core/measurement.hpp"
#include "core/od_dynamics.hpp"
#include "core/state.hpp"

struct ODKFInput {
    StateTr x0_est; // estimated state
    mat6d P0; // state estimation covariance matrix
    svec<Measurement> measurements;
    svec<StateTr> observer_states;
    ODDynamicsConfig dyn_config;
    f64 t0;
    i32 prop_steps;
    mat6d Q; // process noise covariance matrix
};

enum struct ODKFStatus {};

struct ODKFResult {};