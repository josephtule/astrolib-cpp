// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/measurement_noise.hpp"

#include "core/measurement_model.hpp"
#include "util/constants.hpp"
#include "util/math.hpp"
#include "util/units.hpp"

#include "Eigen/Cholesky"

#include <cmath>
#include <random>

StatusCode apply_measurement_noise_diagonal(
    Measurement& meas,
    const MeasurementNoiseOptions& opts,
    ObservationType type
) {
    i32 dim = measurement_dim(type);
    if (dim <= 0 || meas.z.size() != dim || !meas.z.allFinite()) {
        return dim <= 0 ? StatusCode::unsupported_type : StatusCode::invalid_state;
    }
    if (meas.R.rows() != dim || meas.R.cols() != dim || !meas.R.allFinite()) {
        return StatusCode::invalid_covariance;
    }

    std::normal_distribution<f64> noise_unit(0.0, 1.0);
    for (i32 i = 0; i < dim; ++i) {
        f64 var_i = meas.R(i, i);
        if (var_i < 0.0) return StatusCode::invalid_covariance;
        meas.z(i) += noise_unit(opts.rng) * std::sqrt(var_i);
    }

    if (type == ObservationType::radec || type == ObservationType::azel) {
        f64 wrap_min = convert_angle(0.0, UAngle::radian, opts.u_angle);
        f64 wrap_max = convert_angle(twopi, UAngle::radian, opts.u_angle);
        meas.z(0) = wrap_angle(meas.z(0), wrap_min, wrap_max, opts.u_angle, opts.u_angle);
    }
    return StatusCode::ok;
}

StatusCode apply_measurement_noise_cholesky(
    Measurement& meas,
    const MeasurementNoiseOptions& opts,
    ObservationType type
) {
    i32 dim = measurement_dim(type);
    if (dim <= 0 || meas.z.size() != dim || !meas.z.allFinite()) {
        return dim <= 0 ? StatusCode::unsupported_type : StatusCode::invalid_state;
    }
    if (meas.R.rows() != dim || meas.R.cols() != dim || !meas.R.allFinite()) {
        return StatusCode::invalid_covariance;
    }

    Eigen::LLT<matXd> R_llt(meas.R);
    if (R_llt.info() != Eigen::Success) return StatusCode::invalid_covariance;

    vecXd dz(dim);
    std::normal_distribution<f64> noise_unit(0.0, 1.0);
    for (i32 i = 0; i < dim; ++i) dz(i) = noise_unit(opts.rng);

    meas.z += R_llt.matrixL() * dz;
    if (!meas.z.allFinite()) return StatusCode::invalid_state;

    if (type == ObservationType::radec || type == ObservationType::azel) {
        f64 wrap_min = convert_angle(0.0, UAngle::radian, opts.u_angle);
        f64 wrap_max = convert_angle(twopi, UAngle::radian, opts.u_angle);
        meas.z(0) = wrap_angle(meas.z(0), wrap_min, wrap_max, opts.u_angle, opts.u_angle);
    }
    return StatusCode::ok;
}
