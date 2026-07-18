// Copyright 2025-2026 Joseph Tu Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Eigen/Core"
#include "core/estimation_common.hpp"
#include "core/measurement.hpp"
#include "core/od_dynamics.hpp"
#include "core/state.hpp"
#include "util/constants.hpp"
#include "util/vecdefs.hpp"

struct ODBatchInput {
    StateTr x0_guess;
    svec<Measurement> measurements;
    svec<StateTr> observer_states;
    ODDynamicsConfig dyn_config;
    f64 t0 = 0.0;
    i32 max_iters = 5;
    i32 prop_steps = 100; // TODO: this is a temp count
    f64 tol_dx = 1e-8;
    f64 tol_residual = 1e-10;
    // max dx to avoid overcorrection (in state-length units)
    f64 max_dx_r_norm = 5e4;
    f64 max_dx_v_norm = 15.0;
    bool use_line_search = true;
};

struct ODBatchResult {
    StatusCode status = StatusCode::invalid_input;
    StateTr x0_est;
    mat6d covariance = mat6d0;
    mat6d normal_inv = mat6d0;
    f64 residual_norm = 0.0;
    f64 raw_residual_norm = 0.0;
    f64 dx_norm = 0.0;
    i32 iterations = 0;
};

struct ODBatchResidualEval {
    StatusCode status = StatusCode::invalid_input;
    f64 weighted_norm = 0.0;
    f64 raw_norm = 0.0;
};

StatusCode od_batch_validate_input(const ODBatchInput& input);
ODBatchResidualEval od_batch_eval_residual_norm(
    const ODBatchInput& input,
    const StateTr& x0_ref
);
ODBatchResult od_batch_lumve(const ODBatchInput& input);
