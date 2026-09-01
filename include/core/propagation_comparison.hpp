// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/integrator_adaptive.hpp"
#include "core/integrator_common.hpp"
#include "core/status.hpp"
#include "core/world_stepper.hpp"
#include "util/vecdefs.hpp"

struct PropagationRunConfig {
    string label;
    IntegratorType integrator{IntegratorTypeFixed::rk4};

    f64 nominal_dt = 1.0;
    svec<f64> output_times;
    svec<string> body_config_ids;

    bool step_tr = true;
    bool step_att = true;
    bool include_attitude_output = true;

    AdaptiveIntegratorConfig adaptive{};
    string frame_label = "simulation_inertial";
};

struct PropagationStateSample {
    f64 t = 0.0;
    string body_config_id;
    StateTr x_tr{};
    StateAtt x_att{};
    bool has_attitude = false;
};

struct PropagationRunResult {
    StatusCode status = StatusCode::invalid_state;
    string label;
    IntegratorType integrator{IntegratorTypeFixed::rk4};

    svec<PropagationStateSample> samples;
    WorldStepperStats stats{};
    f64 runtime_ms = 0.0;
};