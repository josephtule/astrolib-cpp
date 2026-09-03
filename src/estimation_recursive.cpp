// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/estimation_recursive.hpp"
#include "core/estimation_common.hpp"
#include "core/measurement.hpp"
#include "core/od_dynamics.hpp"
#include "core/state.hpp"
#include "core/od_estimator_context.hpp"

#include <cmath>

StatusCode od_ekf_step_validate_input(const ODEKFStepInput& input) {
    if ((!input.propagation && input.dyn_config.mu <= 0.0) || input.prop_steps <= 0) {
        return StatusCode::validation_failed;
    }
    if (!std::isfinite(input.tol_time) || input.tol_time < 0.0) {
        return StatusCode::validation_failed;
    }
    if (input.observer_uncertainty.enabled && input.measurement.R.size() == 0) {
        return StatusCode::invalid_covariance;
    }
    i32 dim = measurement_dim(input.measurement.type);
    if (dim <= 0) {
        return StatusCode::unsupported_type;
    }
    if (input.measurement.R.size() != 0
        && (input.measurement.R.cols() != dim || input.measurement.R.rows() != dim
            || !input.measurement.R.allFinite())) {
        return StatusCode::invalid_covariance;
    }
    if (!input.Q.allFinite() || !input.filter.P.allFinite()) {
        return StatusCode::invalid_covariance;
    }
    if (input.measurement.z.size() != dim || !input.measurement.z.allFinite()) {
        return StatusCode::invalid_state;
    }

    if (!input.propagation && (input.dyn_config.zonal_degree < 0 || input.dyn_config.zonal_degree > 6)) {
        return StatusCode::validation_failed;
    }

    return input.propagation ? validate_od_estimator_context(*input.propagation) : StatusCode::ok;
}

StatusCode od_ekf_validate_input(const ODEKFOfflineInput& input) {
    if (input.measurements.size() == 0) return StatusCode::empty_measurements;
    if (!input.observer_uncertainties.empty()
        && input.observer_uncertainties.size() != input.measurements.size()) {
        return StatusCode::size_mismatch;
    }
    if (input.measurements.size() != input.observer_states.size()) {
        return StatusCode::size_mismatch;
    }
    if ((!input.propagation && input.dyn_config.mu <= 0.0) || input.prop_steps <= 0) {
        return StatusCode::validation_failed;
    }
    if (!std::isfinite(input.tol_time) || input.tol_time < 0.0) {
        return StatusCode::validation_failed;
    }
    if (!input.Q.allFinite() || !input.initial_filter.P.allFinite()) {
        return StatusCode::invalid_covariance;
    }
    for (i32 i = 0; i < input.measurements.size(); ++i) {
        const Measurement& meas = input.measurements[i];
        if (!input.observer_uncertainties.empty()
            && input.observer_uncertainties[i].enabled && meas.R.size() == 0) {
            return StatusCode::invalid_covariance;
        }
        i32 dim = measurement_dim(meas.type);
        if (dim <= 0) {
            return StatusCode::unsupported_type;
        }
        if (meas.z.size() != dim || !meas.z.allFinite()) {
            return StatusCode::invalid_state;
        }
        if (meas.R.size() != 0
            && (meas.R.cols() != dim || meas.R.rows() != dim || !meas.R.allFinite())) {
            return StatusCode::invalid_covariance;
        }
    }
    return input.propagation ? validate_od_estimator_context(*input.propagation) : StatusCode::ok;
}

ODEKFPredictResult od_ekf_predict(
    const ODEKFState& filter,
    f64 t_target,
    const ODDynamicsConfig& dyn_config,
    i32 prop_steps,
    const mat6d& Q,
    f64 tol,
    const ODEstimatorContext* propagation
) {
    ODEKFPredictResult result;
    result.y.x = filter.x;
    result.y.Phi = mat6d1;
    result.P = filter.P;
    result.t = filter.t;
    f64 dt = t_target - filter.t;
    if (dt < -tol) {
        result.status = StatusCode::propagation_failed;
        return result;
    } else if (std::abs(dt) <= tol) {
        // same epoch
        dt = 0.0;
    }
    bool propagate = dt != 0.0;
    mat6d Q_eff;
    if (std::abs(dt) <= tol) { // sensor fusion mode
        Q_eff = mat6d0;
    } else {
        Q_eff = Q;
    }
    VarStateTr y0;
    y0.x = filter.x;
    y0.Phi = mat6d1;

    // propagate prediction state and STM
    VarStateTr yf = y0;
    if (propagation) {
        result.status = propagate_od_estimator(*propagation, filter.t, filter.x,
            filter.t + dt, prop_steps, false, yf);
        if (result.status != StatusCode::ok) {
            if (propagation->world && propagation->world->world().t_sim() != filter.t) {
                result.y = yf;
                result.t = propagation->world->world().t_sim();
                // partial interval: no full-interval Q injection on failed propagation
                result.P = yf.Phi * filter.P * yf.Phi.transpose();
            }
            return result;
        }
    } else if (propagate) {
        yf = propagate_var_tr_od(filter.t, y0, dt, prop_steps, dyn_config);
    } else {
        yf = y0;
    }

    // state and covariance prediction
    // process noise belongs to elapsed dynamics, not each sensor update
    result.y.x = yf.x;
    result.y.Phi = yf.Phi;
    result.P = result.y.Phi * filter.P * result.y.Phi.transpose() + Q_eff;
    if (!statetr_to_vec6d(result.y.x).allFinite() || !result.y.Phi.allFinite()
        || !result.P.allFinite()) {
        result.status = StatusCode::propagation_failed;
        return result;
    }

    result.t = filter.t + dt;

    result.status = StatusCode::ok;
    return result;
}

ODEKFStepResult od_ekf_step(const ODEKFStepInput& input) {
    // NOTE: for sensor fusion, add multiple measurements at the same timestamp
    // can be from same or different sources
    ODEKFStepResult result;
    const ODEKFState& filter = input.filter;
    const Measurement& meas = input.measurement;
    const StateTr& x_tr_obsv = input.x_tr_observer;
    i32 dim = measurement_dim(meas.type);

    result.status = od_ekf_step_validate_input(input);
    if (result.status != StatusCode::ok) {
        return result;
    } else {
        result.filter = filter;
        result.status = StatusCode::ok;
        result.residual_norm = 0.0;
        result.raw_residual_norm = 0.0;
    }

    ODEKFPredictResult prediction = od_ekf_predict(
        filter,
        meas.t,
        input.dyn_config,
        input.prop_steps,
        input.Q,
        input.tol_time,
        input.propagation
    );

    if (prediction.status != StatusCode::ok) {
        result.filter.x = prediction.y.x;
        result.filter.P = prediction.P;
        result.filter.t = prediction.t;
        result.status = prediction.status;
        return result;
    }

    StateTr& x_pred = prediction.y.x;
    mat6d& P_pred = prediction.P;
    f64& t_pred = prediction.t;

    result.filter.x = x_pred;
    result.filter.P = P_pred;
    result.filter.t = t_pred;


    if (input.observer_uncertainty.enabled && meas.type == ObservationType::azel) {
        result.status = StatusCode::unsupported_type;
        return result;
    }

    // measurement prediction
    MeasurementContext ctx = make_measurement_context(x_pred, x_tr_obsv);
    mat3d R_measurement_I = mat3d::Identity();
    if (input.propagation) {
        Measurement geometry_measurement = meas;
        geometry_measurement.t = t_pred;
        result.status = od_estimator_measurement_context(*input.propagation, geometry_measurement,
            x_pred, x_tr_obsv, ctx, R_measurement_I);
        if (result.status != StatusCode::ok) return result;
    }
    vecXd z_pred = predict_measurement(meas.type, ctx, input.angle_in, input.angle_out, input.tol_measurement);
    if (z_pred.size() != dim) {
        result.status = StatusCode::size_mismatch;
        return result;
    }

    // residuals
    vecXd res = measurement_residual(meas.type, meas.z, z_pred, input.angle_out);
    if (!res.allFinite()) {
        result.status = StatusCode::invalid_state;
        return result;
    }
    result.raw_residual_norm = res.norm();

    // measurement jacobian
    matXd H = measurement_jacobian(meas.type, ctx, input.angle_in, input.angle_out,
        input.eps_pos, input.eps_vel, input.tol_measurement);
    if (H.cols() != 6 || H.rows() != dim || !H.allFinite()) {
        result.status = StatusCode::size_mismatch;
        return result;
    }

    H.leftCols(3) = (H.leftCols(3) * R_measurement_I).eval();

    // measurement covariance
    matXd R;
    if (meas.R.size() == 0) {
        R = matXd::Identity(dim, dim);
    } else {
        R = meas.R;
    }

    if (input.observer_uncertainty.enabled) {
        matXd H_observer;
        result.status = measurement_observer_jacobian(meas.type, ctx, H_observer,
            input.angle_in, input.angle_out, input.eps_pos, input.eps_vel, input.tol_measurement);
        if (result.status != StatusCode::ok) return result;
        result.status = effective_measurement_covariance(
            meas.R, H_observer, input.observer_uncertainty.P, R
        );
        if (result.status != StatusCode::ok) return result;
    }

    // innovation covariance
    matXd S = H * P_pred * H.transpose() + R;
    if (!S.allFinite()) {
        result.status = StatusCode::invalid_covariance;
        return result;
    }

    // solve for Kalman gain
    Eigen::LDLT<matXd> S_ldlt(S);
    if (S_ldlt.info() != Eigen::Success) {
        result.status = StatusCode::singular_innovation;
        return result;
    }
    matXd X = S_ldlt.solve(H * P_pred); // solve S * X = H * P_pred for X
    // matXd K = P_pred * H.transpose() * S.inverse();
    matXd K = X.transpose();
    if (!X.allFinite() || !K.allFinite()) {
        result.status = StatusCode::singular_innovation;
        return result;
    }

    // a posteriori state and STM, updated estimate
    vec6d dx_vec = K * res;
    DerivTr dx = vec6d_to_derivtr(dx_vec);
    StateTr x_post = x_pred + dx;
    // mat6d P_post = (mat6d1 - K * H) * P_pred;
    mat6d P_post = (mat6d1 - K * H) * P_pred * (mat6d1 - K * H).transpose()
                   + K * R * K.transpose();       // Joseph form, more stable
    P_post = 0.5 * (P_post + P_post.transpose()); // force symmetry
    if (!dx_vec.allFinite() || !statetr_to_vec6d(x_post).allFinite()
        || !P_post.allFinite()) {
        result.status = StatusCode::correction_rejected;
        return result;
    }

    vecXd weighted_res = S_ldlt.solve(res);
    if (!weighted_res.allFinite()) {
        result.status = StatusCode::singular_innovation;
        return result;
    }
    f64 res_norm2 = res.dot(weighted_res);
    if (!std::isfinite(res_norm2) || res_norm2 < 0.0) {
        result.status = StatusCode::singular_innovation;
        return result;
    }

    // store results
    result.filter.x = x_post;
    result.filter.P = P_post;
    result.filter.t = t_pred;
    result.residual = res;
    result.residual_norm = std::sqrt(res_norm2);
    // result.residual_norm = std::sqrt(res.transpose() * S.inverse() * res);
    result.raw_residual_norm = res.norm();
    if (input.propagation && input.propagation->world)
        input.propagation->world->world().body(input.propagation->world->target_id())->x_tr = x_post;
    result.status = StatusCode::ok;
    return result;
}

ODEKFResult od_ekf_offline(const ODEKFOfflineInput& input) {
    ODEKFResult result;
    ODEKFState filter = input.initial_filter;

    result.status = od_ekf_validate_input(input);
    if (result.status != StatusCode::ok) {
        return result;
    }

    if (input.propagation && input.propagation->world) {
        if (filter.t != input.propagation->world->reference_time()) {
            result.status = StatusCode::time_mismatch;
            return result;
        }
        result.status = input.propagation->world->reset(filter.x);
        if (result.status != StatusCode::ok) return result;
    }
    // EKF loop
    for (i32 i = 0; i < input.measurements.size(); ++i) {
        const Measurement& meas = input.measurements[i];
        ODEKFStepInput step_input{
            .filter = filter,
            .measurement = meas,
            .x_tr_observer = input.observer_states[i],
            .dyn_config = input.dyn_config,
            .prop_steps = input.prop_steps,
            .Q = input.Q,
            .tol_time = input.tol_time
        };

        if (!input.observer_uncertainties.empty()) {
            step_input.observer_uncertainty = input.observer_uncertainties[i];
        }
        step_input.propagation = input.propagation;
        ODEKFStepResult step_result = od_ekf_step(step_input);
        if (!od_status_success(step_result.status)) {
            result.filter = step_result.filter;
            result.status = step_result.status;
            result.processed_measurements = i;
            return result;
        }
        filter = step_result.filter;
        result.processed_measurements = i + 1;
        result.residual_norm = step_result.residual_norm;
        result.raw_residual_norm = step_result.raw_residual_norm;
    }

    result.filter = filter;
    result.status = StatusCode::ok;
    result.filter = filter;

    return result;
}
