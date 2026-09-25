// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/integrator_adaptive.hpp"
#include "core/integrator_common.hpp"
#include "core/status.hpp"
#include "core/world_stepper.hpp"
#include "util/constants.hpp"
#include "util/typedefs.hpp"
#include "util/vecdefs.hpp"

enum struct PropagationBackend {
    world_stepper_legacy,
    world_stepper,
    compiled_problem,
};

inline string propagation_backend_string(PropagationBackend backend) {
    switch (backend) {
    case PropagationBackend::world_stepper_legacy: return "world_stepper_legacy";
    case PropagationBackend::world_stepper: return "world_stepper";
    case PropagationBackend::compiled_problem: return "compiled_problem";
    }
    return "unknown";
}

struct PropagationComparisonTolerance {
    f64 time = tol12;
    f64 position = tol9;
    f64 velocity = tol12;
    f64 quaternion = tol12;
    f64 angular_velocity = tol12;
    f64 invariant = tol9;
};

struct PropagationInvariantConfig {
    string central_body_config_id;
    bool include_orbital_invariants = false;
};

struct PropagationRunConfig {
    string label;
    IntegratorType integrator{IntegratorTypeFixed::rk4};

    PropagationBackend backend = PropagationBackend::world_stepper;

    f64 t0 = 0.0;
    f64 tf = 0.0;
    svec<f64> sample_epochs; // ordered
    f64 nominal_dt = 1.0;
    svec<string> body_config_ids;

    bool step_tr = true;
    bool step_att = true;
    bool include_attitude_output = true;

    AdaptiveIntegratorConfig adaptive{};
    string frame_label = "simulation_inertial";
    PropagationInvariantConfig invariants{};
    PropagationComparisonTolerance tolerance{};
};

struct PropagationInvariantSample {
    bool has_orbital_invariants = false;
    f64 specific_energy = 0.0;
    vec3d specific_angular_momentum = vec3d0;
    f64 specific_angular_momentum_norm = 0.0;
};

struct PropagationStateSample {
    f64 requested_t = 0.0;
    f64 t = 0.0;
    string body_config_id;
    StateTr x_tr{};
    StateAtt x_att{};
    bool has_attitude = false;
    f64 quaternion_norm = 0.0;
    PropagationInvariantSample invariants{};
};

struct PropagationRunMetrics {
    i64 derivative_evaluations = 0;
    i64 provider_translation_queries = 0;
    i64 provider_orientation_queries = 0;
    i64 workspace_rebuilds = 0;
    i64 stage_builds = 0;
    i64 allocation_count = 0;
    bool allocation_count_available = false;
    f64 runtime_ms = 0.0;
};

struct PropagationRunResult {
    StatusCode status = StatusCode::invalid_state;
    string label;
    string frame_label;
    PropagationBackend backend = PropagationBackend::world_stepper;
    IntegratorType integrator{IntegratorTypeFixed::rk4};

    f64 requested_t0 = 0.0;
    f64 requested_tf = 0.0;
    f64 actual_t0 = 0.0;
    f64 actual_tf = 0.0;
    f64 failure_t = 0.0;

    svec<PropagationStateSample> samples;
    WorldStepperStats stats{};
    PropagationRunMetrics metrics{};
};

struct PropagationSampleDifference {
    f64 requested_t = 0.0;
    string body_config_id;

    f64 time = 0.0;
    f64 position = 0.0;
    f64 velocity = 0.0;
    f64 quaternion = 0.0;
    f64 quaternion_norm = 0.0;
    f64 angular_velocity = 0.0;
    f64 specific_energy = 0.0;
    f64 specific_angular_momentum = 0.0;

    bool passed = false;
};

struct PropagationComparisonResult {
    StatusCode status = StatusCode::invalid_state;
    bool passed = false;

    string reference_label;
    string candidate_label;
    f64 first_failed_t = 0.0;
    string first_failed_body_config_id;

    f64 max_time = 0.0;
    f64 max_position = 0.0;
    f64 max_velocity = 0.0;
    f64 max_quaternion = 0.0;
    f64 max_quaternion_norm = 0.0;
    f64 max_angular_velocity = 0.0;
    f64 max_specific_energy = 0.0;
    f64 max_specific_angular_momentum = 0.0;

    svec<PropagationSampleDifference> samples;
};

StatusCode run_propagation_case(
    World& world,
    const umap<string, EntityId>& body_ids,
    const WorldStepperConfig& stepper_cfg,
    const PropagationRunConfig& run_cfg,
    PropagationRunResult& out
);

StatusCode compare_propagation_runs(
    const PropagationRunResult& reference,
    const PropagationRunResult& candidate,
    const PropagationComparisonTolerance& tolerance,
    PropagationComparisonResult& out
);

StatusCode save_propagation_run_csv(
    const string& filepath,
    const PropagationRunResult& result
);
