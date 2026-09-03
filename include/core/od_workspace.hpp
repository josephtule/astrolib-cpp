// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/scenario_io.hpp"
#include "core/world_stepper.hpp"
#include <variant>

struct ODWorldPropagationOptions {
    f64 max_step = 60.0; // maximum interval passed to the world stepper
    i32 max_steps = 100000;
    bool with_stm = true;
    mat6d Phi0 = mat6d1; // identity gives the transition from this call's start epoch
    mat6d stm_abs_tol = mat6d::Constant(1e-9);
    f64 stm_rel_tol = 1e-9;
};

struct ODWorldPropagationResult {
    StatusCode status = StatusCode::invalid_input;
    StateTr x;
    mat6d Phi = mat6d1;
    bool has_stm = false;
    f64 t = 0.0;
    WorldStepperStats stats;
};

// one private world per independent full-world estimator session
class ODWorldWorkspace {
    using ReferenceBody = std::variant<Celestial, Satellite, Station>;
    World world_;
    WorldStepperWorkspace workspace_;
    WorldStepperConfig stepper_;
    WorldStepperConfig reference_stepper_;
    ScenarioBuildResult bindings_;
    svec<ReferenceBody> reference_bodies_;
    svec<EntityId> reference_active_ids_;
    EntityId target_id_ = kInvalidEntityId;
    f64 t0_ = 0.0;
    JulianDate epoch_;
    TimeScale time_scale_ = TimeScale::tt;
    TimeOffsets offsets_;
    bool initialized_ = false;

  public:
    StatusCode initialize(
        const ScenarioConfig& scenario, const string& target,
        const svec<string>& sources
    );
    StatusCode reset(const StateTr& x_target0_I);
    ODWorldPropagationResult propagate(
        f64 t_target, const ODWorldPropagationOptions& options = {}
    );
    World& world() { return world_; }
    const World& world() const { return world_; }
    WorldStepperWorkspace& workspace() { return workspace_; }
    WorldStepperConfig& stepper() { return stepper_; }
    const ScenarioBuildResult& bindings() const { return bindings_; }
    EntityId target_id() const { return target_id_; }
    f64 reference_time() const { return t0_; }
};
