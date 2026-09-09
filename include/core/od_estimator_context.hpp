// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "core/od_environment.hpp"
#include "core/od_workspace.hpp"
#include "core/measurement.hpp"

struct ODEstimatorContext {
    const ODDynamicsContext* dynamics = nullptr;
    ODWorldWorkspace* world = nullptr; // exactly one of dynamics/world; borrowed per run
    const World* observer_geometry = nullptr; // station IDs and fixed attachment geometry
    IntegratorTypeFixed integrator = IntegratorTypeFixed::rk4;
    f64 fixed_step_size = 0.0; // dynamics path: 0 uses steps; positive uses fixed steps + final remainder
    ODWorldPropagationOptions options;
};
StatusCode validate_od_estimator_context(const ODEstimatorContext& ctx);
StatusCode propagate_od_estimator(
    const ODEstimatorContext& ctx, f64 t0, const StateTr& x0, f64 tf,
    i32 steps, bool reset_reference, VarStateTr& out
);
StatusCode od_estimator_measurement_context(
    const ODEstimatorContext& ctx, const Measurement& measurement,
    const StateTr& target, const StateTr& observer, MeasurementContext& out,
    mat3d& R_measurement_I
);
