// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/observation_type.hpp"
#include "core/status.hpp"
#include "util/constants.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

struct MeasurementContext;

struct ObserverUncertainty {
    bool enabled = false;
    // [r, v] covariance at measurement epoch, in observer-state inertial axes
    mat6d P = mat6d0;
};

// first-order independent observer error
// does not modify the sensor covariance
StatusCode effective_measurement_covariance(
    const matXd& R_sensor,
    const matXd& H_observer,
    const mat6d& P_observer,
    matXd& out
);

StatusCode measurement_observer_jacobian(
    ObservationType type,
    const MeasurementContext& ctx,
    matXd& out,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 eps_pos = 1e-3,
    f64 eps_vel = 1e-6,
    f64 tol = tol12
);
