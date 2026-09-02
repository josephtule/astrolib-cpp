// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/measurement_uncertainty.hpp"
#include "core/measurement.hpp"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <cmath>

namespace {

// normalize by standard deviations so pos/vel units do not set the PSD tolerance
StatusCode checked_covariance(const matXd& input, bool positive_definite, matXd& out) {
    if (input.rows() == 0 || input.rows() != input.cols() || !input.allFinite()) {
        return StatusCode::invalid_covariance;
    }
    const auto n = input.rows();
    vecXd scale(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        if (input(i, i) < 0.0 || (positive_definite && input(i, i) == 0.0)) {
            return StatusCode::invalid_covariance;
        }
        scale(i) = std::sqrt(input(i, i));
    }
    matXd normalized = matXd::Zero(n, n);
    for (Eigen::Index i = 0; i < n; ++i) {
        for (Eigen::Index j = 0; j < n; ++j) {
            if (scale(i) == 0.0 || scale(j) == 0.0) {
                // zero variance requires zero covariance with every component
                if (input(i, j) != 0.0) return StatusCode::invalid_covariance;
            } else {
                normalized(i, j) = input(i, j) / scale(i) / scale(j);
            }
        }
    }
    if (!normalized.allFinite()) return StatusCode::invalid_covariance;
    constexpr f64 tol = 1e-12;
    if ((normalized - normalized.transpose()).cwiseAbs().maxCoeff() > tol) {
        return StatusCode::invalid_covariance;
    }
    normalized = (0.5 * normalized + 0.5 * normalized.transpose()).eval();
    if (positive_definite) {
        Eigen::LLT<matXd> factor(normalized);
        if (factor.info() != Eigen::Success) return StatusCode::invalid_covariance;
        out = (0.5 * input + 0.5 * input.transpose()).eval();
    } else {
        Eigen::SelfAdjointEigenSolver<matXd> factor(normalized);
        if (factor.info() != Eigen::Success
            || !factor.eigenvalues().allFinite()
            || factor.eigenvalues().minCoeff() < -tol * n) {
            return StatusCode::invalid_covariance;
        }

        if (factor.eigenvalues().minCoeff() < 0.0) {
            // remove only roundoff-sized negative modes before projecting uncertainty
            normalized = factor.eigenvectors()
                * factor.eigenvalues().cwiseMax(0.0).asDiagonal()
                * factor.eigenvectors().transpose();
            out = scale.asDiagonal() * normalized * scale.asDiagonal();
        } else {
            out = (0.5 * input + 0.5 * input.transpose()).eval();
        }
    }
    return out.allFinite() ? StatusCode::ok : StatusCode::invalid_covariance;
}

} // namespace

StatusCode measurement_observer_jacobian(
    ObservationType type,
    const MeasurementContext& ctx,
    matXd& out,
    UAngle angle_in,
    UAngle angle_out,
    f64 eps_pos,
    f64 eps_vel,
    f64 tol
) {
    matXd temp;
    switch (type) {
        case ObservationType::pos:
        case ObservationType::pos_vel:
            temp = matXd::Zero(measurement_dim(type), 6);
            break;
        case ObservationType::range:
        case ObservationType::range_rate:
        case ObservationType::radec:
        case ObservationType::rel_pos:
        case ObservationType::rel_pos_vel:
            temp = -measurement_jacobian(
                type, ctx, angle_in, angle_out, eps_pos, eps_vel, tol
            );
            break;
        default:
            // azel also depends on the observer's local frame
            return StatusCode::unsupported_type;
    }
    if (temp.rows() != measurement_dim(type) || temp.cols() != 6) {
        return StatusCode::size_mismatch;
    }
    if (!temp.allFinite()) return StatusCode::non_finite_result;
    out = std::move(temp);
    return StatusCode::ok;
}

StatusCode effective_measurement_covariance(
    const matXd& R_sensor,
    const matXd& H_observer,
    const mat6d& P_observer,
    matXd& out
) {
    if (R_sensor.rows() == 0 || R_sensor.rows() != R_sensor.cols()) {
        return StatusCode::invalid_covariance;
    }
    if (H_observer.rows() != R_sensor.rows() || H_observer.cols() != 6) {
        return StatusCode::size_mismatch;
    }
    if (!H_observer.allFinite()) return StatusCode::invalid_input;

    matXd R, P;
    StatusCode status = checked_covariance(R_sensor, true, R);
    if (status != StatusCode::ok) return status;

    status = checked_covariance(P_observer, false, P);
    if (status != StatusCode::ok) return status;

    matXd temp = R + H_observer * P * H_observer.transpose();
    if (!temp.allFinite()) return StatusCode::non_finite_result;

    temp = (0.5 * temp + 0.5 * temp.transpose()).eval();
    matXd validated;
    status = checked_covariance(temp, true, validated);
    if (status != StatusCode::ok) return status;

    out = std::move(validated);
    return StatusCode::ok;
}
