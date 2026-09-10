// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/measurement_jacobian.hpp"

#include "core/measurement_model.hpp"
#include "core/station_geometry.hpp"
#include "util/units.hpp"

#include <cmath>

matXd jacobian_fd_measurement(
    ObservationType type,
    const MeasurementContext& ctx,
    UAngle angle_in,
    UAngle angle_out,
    f64 eps_pos,
    f64 eps_vel,
    f64 tol
) {
    i32 rows = measurement_dim(type);
    i32 cols = 6;
    matXd H(rows, cols);

    vecXd z0 = predict_measurement(type, ctx, angle_in, angle_out, tol);
    if (z0.size() != rows) return matXd{};

    MeasurementContext ctx_pert = ctx;
    for (i32 i = 0; i < cols; ++i) {
        // forward differencing
        f64 eps_i = i < 3 ? eps_pos : eps_vel;

        vec6d x_target_pert_vec = statetr_to_vec6d(ctx.x_tr_target);
        x_target_pert_vec(i) += eps_i;
        StateTr x_target_pert = vec6d_to_statetr(x_target_pert_vec);
        ctx_pert.x_tr_target = x_target_pert;
        vecXd z_pert = predict_measurement(type, ctx_pert, angle_in, angle_out, tol);
        if (z_pert.size() != rows) return matXd{};
        vecXd residual = measurement_residual(type, z_pert, z0, angle_out);
        H.col(i) = residual / eps_i;
    }

    return H;
}

matXd jacobian_radec(const MeasurementContext& ctx, UAngle angle_out, f64 tol) {
    vec3d r_rel = ctx.x_tr_target.r - ctx.x_tr_observer.r;
    f64 rho2 = r_rel.squaredNorm();
    f64 rho = std::sqrt(rho2);
    f64 rhoxy2 = r_rel.segment<2>(0).squaredNorm();
    f64 rhoxy = std::sqrt(rhoxy2);
    f64 denom = rhoxy * rho2;

    if (rho2 <= tol || rhoxy2 <= tol || denom <= tol) return matXd{};

    matd<2, 6> H;
    H.row(0) = vec6d{-r_rel(1) / rhoxy2, r_rel(0) / rhoxy2, 0.0, 0.0, 0.0, 0.0};
    H.row(1) = vec6d{
        -r_rel(2) * r_rel(0) / denom,
        -r_rel(2) * r_rel(1) / denom,
        rhoxy / rho2,
        0.0,
        0.0,
        0.0
    };
    H *= convert_angle(1.0, UAngle::radian, angle_out);
    return H;
}

matXd jacobian_azel_enu(const MeasurementContext& ctx, UAngle angle_out, f64 tol) {
    vec3d r_rel_enu = ctx.x_tr_target.r - ctx.x_tr_observer.r;
    f64 e = r_rel_enu(0);
    f64 n = r_rel_enu(1);
    f64 u = r_rel_enu(2);

    f64 rho2 = r_rel_enu.squaredNorm();
    f64 rhoen2 = r_rel_enu.segment<2>(0).squaredNorm();
    f64 rhoen = std::sqrt(rhoen2);
    f64 denom = rhoen * rho2;

    if (rho2 <= tol || rhoen2 <= tol || denom <= tol) return matXd{};

    matd<2, 6> H;
    H.row(0) = vec6d{n / rhoen2, -e / rhoen2, 0.0, 0.0, 0.0, 0.0};
    H.row(1) = vec6d{-u * e / denom, -u * n / denom, rhoen / rho2, 0.0, 0.0, 0.0};
    H *= convert_angle(1.0, UAngle::radian, angle_out);
    return H;
}

matXd jacobian_azel_inertial_from_enu(
    const MeasurementContext& ctx,
    const mat3d& R_ENU_I,
    UAngle angle_out,
    f64 tol
) {
    matXd H = jacobian_azel_enu(ctx, angle_out, tol);
    if (H.rows() != 2 || H.cols() != 6) return H;
    H.block<2, 3>(0, 0) = H.block<2, 3>(0, 0) * R_ENU_I;
    return H;
}

matXd jacobian_azel_bcbf(
    const MeasurementContext& ctx,
    UAngle angle_in,
    UAngle angle_out,
    f64 tol
) {
    f64 lat = ctx.station_llh(0);
    f64 lon = ctx.station_llh(1);
    mat3d R_ENU_BCBF = rot_enu_from_bcbf(lat, lon, angle_in);
    MeasurementContext enu_ctx = ctx;
    enu_ctx.x_tr_observer = R_ENU_BCBF * enu_ctx.x_tr_observer;
    enu_ctx.x_tr_target = R_ENU_BCBF * enu_ctx.x_tr_target;

    matXd H = jacobian_azel_enu(enu_ctx, angle_out, tol);
    if (H.rows() != 2 || H.cols() != 6) return H;
    H.block<2, 3>(0, 0) = H.block<2, 3>(0, 0) * R_ENU_BCBF;
    return H;
}

matXd jacobian_azel_inertial_from_bcbf(
    const MeasurementContext& ctx,
    const mat3d& R_BCBF_I,
    UAngle angle_in,
    UAngle angle_out,
    f64 tol
) {
    matXd H = jacobian_azel_bcbf(ctx, angle_in, angle_out, tol);
    if (H.rows() != 2 || H.cols() != 6) return H;
    H.block<2, 3>(0, 0) = H.block<2, 3>(0, 0) * R_BCBF_I;
    return H;
}

matXd measurement_jacobian(
    ObservationType type,
    const MeasurementContext& ctx,
    UAngle angle_in,
    UAngle angle_out,
    f64 eps_pos,
    f64 eps_vel,
    f64 tol
) {
    if (measurement_dim(type) == 0) return matXd{};

    switch (type) {
    case ObservationType::radec: return jacobian_radec(ctx, angle_out, tol);
    // TODO: add jacobians for the below if possible
    case ObservationType::azel: [[fallthrough]]; // only works in world context
    case ObservationType::range: [[fallthrough]];
    case ObservationType::range_rate: [[fallthrough]];
    case ObservationType::pos: [[fallthrough]];
    case ObservationType::pos_vel: [[fallthrough]];
    case ObservationType::rel_pos: [[fallthrough]];
    case ObservationType::rel_pos_vel: [[fallthrough]];
    default:
        return jacobian_fd_measurement(
            type,
            ctx,
            angle_in,
            angle_out,
            eps_pos,
            eps_vel,
            tol
        );
    }
}
