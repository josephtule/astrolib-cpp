// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/measurement_types.hpp"
#include "util/constants.hpp"

vecXd predict_measurement(
    ObservationType type,
    const MeasurementContext& ctx,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
);

vecXd predict_measurement(
    ObservationType type,
    const StateTr& x_target,
    const StateTr& x_observer,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
);

vecXd measurement_residual(
    ObservationType type,
    ecref<vecXd> z_obs,
    ecref<vecXd> z_pred,
    UAngle angle_in = UAngle::radian
);
