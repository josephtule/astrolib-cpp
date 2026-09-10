// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/measurement_types.hpp"
#include "util/constants.hpp"

matXd jacobian_fd_measurement(
    ObservationType type,
    const MeasurementContext& ctx,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 eps_pos = 1e-3,
    f64 eps_vel = 1e-6,
    f64 tol = tol12
);

matXd jacobian_radec(
    const MeasurementContext& ctx,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
);

matXd jacobian_azel_enu(
    const MeasurementContext& ctx,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
);

matXd jacobian_azel_inertial_from_enu(
    const MeasurementContext& ctx,
    const mat3d& R_ENU_I,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
);

matXd jacobian_azel_bcbf(
    const MeasurementContext& ctx,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
);

matXd jacobian_azel_inertial_from_bcbf(
    const MeasurementContext& ctx,
    const mat3d& R_BCBF_I,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
);

matXd measurement_jacobian(
    ObservationType type,
    const MeasurementContext& ctx,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 eps_pos = 1e-3,
    f64 eps_vel = 1e-6,
    f64 tol = tol12
);
