#pragma once

#include "core/estimation_common.hpp"
#include "core/measurement.hpp"
#include "core/od_dynamics.hpp"
#include "core/state.hpp"
#include "util/constants.hpp"
#include <cmath>

struct ODEKFState {
    StateTr x; // state
    mat6d P;   // state estimate covariance
    f64 t = 0.0;
};

// step used for real-time kf, also used by offline kf

struct ODEKFStepInput {
    ODEKFState filter;
    Measurement measurement;
    StateTr observer_state;
    ODDynamicsConfig dyn_config;
    i32 prop_steps = 100;
    mat6d Q = mat6d0; // process noise covariance (process uncertainty)
    f64 tol_time = tol_strict;
};

struct ODEKFStepResult {
    bool success = false;
    ODStatus status = ODStatus::invalid_input;
    ODEKFState filter;
    vecXd residual;
    f64 residual_norm = 0.0;
    f64 raw_residual_norm = 0.0;
};

// offline ekf
struct ODEKFInput {
    ODEKFState initial_filter;
    svec<Measurement> measurements;
    svec<StateTr> observer_states;
    ODDynamicsConfig dyn_config;
    i32 prop_steps = 100;
    mat6d Q = mat6d0;
    f64 tol_time = tol_strict;
};
struct ODEKFResult {
    bool success = false;
    ODStatus status = ODStatus::invalid_input;
    ODEKFState filter;
    i32 processed_measurements = 0;
    f64 residual_norm = 0.0;
    f64 raw_residual_norm = 0.0;
};

inline ODStatus od_ekf_step_validate_input(const ODEKFStepInput& input) {
    if (input.dyn_config.mu <= 0.0 || input.prop_steps <= 0) {
        return ODStatus::invalid_input;
    }
    if (!std::isfinite(input.tol_time) || input.tol_time < 0.0) {
        return ODStatus::invalid_input;
    }
    i32 dim = measurement_dim(input.measurement.type);
    if (dim <= 0) {
        return ODStatus::invalid_input;
    }
    if (input.measurement.R.size() != 0
        && (input.measurement.R.cols() != dim || input.measurement.R.rows() != dim
            || !input.measurement.R.allFinite())) {
        return ODStatus::invalid_covariance;
    }
    if (!input.Q.allFinite() || !input.filter.P.allFinite()) {
        return ODStatus::invalid_covariance;
    }
    if (input.measurement.z.size() != dim || !input.measurement.z.allFinite()) {
        return ODStatus::invalid_input;
    }
    return ODStatus::ok;
}

inline ODStatus od_ekf_validate_input(const ODEKFInput& input) {
    if (input.measurements.size() == 0) return ODStatus::empty_measurements;
    if (input.measurements.size() != input.observer_states.size()) {
        return ODStatus::size_mismatch;
    }
    if (input.dyn_config.mu <= 0.0 || input.prop_steps <= 0) {
        return ODStatus::invalid_input;
    }
    if (!std::isfinite(input.tol_time) || input.tol_time < 0.0) {
        return ODStatus::invalid_input;
    }
    if (!input.Q.allFinite() || !input.initial_filter.P.allFinite()) {
        return ODStatus::invalid_covariance;
    }
    for (i32 i = 0; i < input.measurements.size(); ++i) {
        const Measurement& meas = input.measurements[i];
        i32 dim = measurement_dim(meas.type);
        if (dim <= 0) {
            return ODStatus::invalid_input;
        }
        if (meas.z.size() != dim || !meas.z.allFinite()) {
            return ODStatus::invalid_input;
        }
        if (meas.R.size() != 0
            && (meas.R.cols() != dim || meas.R.rows() != dim || !meas.R.allFinite())) {
            return ODStatus::invalid_covariance;
        }
    }
    return ODStatus::ok;
}

inline ODEKFStepResult od_ekf_step(const ODEKFStepInput& input) {
    // NOTE: for sensor fusion, add three measurements at the same timestamp
    ODEKFStepResult result;
    const ODEKFState& filter = input.filter;
    const Measurement& meas = input.measurement;
    const StateTr& x_obsv = input.observer_state;
    i32 dim = measurement_dim(meas.type);

    result.status = od_ekf_step_validate_input(input);
    if (result.status != ODStatus::ok) {
        return result;
    } else {
        result.filter = filter;
        result.status = ODStatus::ok;
        result.residual_norm = 0.0;
        result.raw_residual_norm = 0.0;
    }

    f64 dt = meas.t - filter.t;
    if (dt < -input.tol_time) {
        // reject negative for now, TODO: create separate for smoothing/backward filtering
        result.status = ODStatus::invalid_input;
        return result;
    } else if (std::abs(dt) <= input.tol_time) {
        // same epoch
        dt = 0.0;
    }
    bool propagate = dt != 0.0;
    mat6d Q_eff;
    if (std::abs(dt) <= input.tol_time) { // sensor fusion mode
        Q_eff = mat6d0;
    } else {
        Q_eff = input.Q;
    }
    VarStateTr y0;
    y0.x = input.filter.x;
    y0.Phi = mat6d1;

    // propagate prediction state and STM
    VarStateTr yf;
    if (propagate) {
        yf = propagate_var_tr_od_rk4(
            filter.t,
            y0,
            dt,
            input.prop_steps,
            input.dyn_config
        );
    } else {
        yf = y0;
    }

    // state and covariance prediction
    // process noise belongs to elapsed dynamics, not each sensor update
    StateTr x_pred = yf.x;
    mat6d Phi = yf.Phi;
    mat6d P_pred = Phi * filter.P * Phi.transpose() + Q_eff;
    f64 t_pred = filter.t + dt;
    if (!statetr_to_vec6(x_pred).allFinite() || !Phi.allFinite() || !P_pred.allFinite()) {
        result.status = ODStatus::propagation_failed;
        return result;
    }
    // predicted state fallback
    result.filter.x = x_pred;
    result.filter.P = P_pred;
    result.filter.t = t_pred;

    // measurement prediction
    MeasurementContext ctx;
    ctx.x_target = x_pred;
    ctx.x_observer = x_obsv;
    vecXd z_pred = predict_measurement(meas.type, ctx);
    if (z_pred.size() != dim) {
        result.status = ODStatus::invalid_input;
        return result;
    }

    // residuals
    vecXd res = measurement_residual(meas.type, meas.z, z_pred);
    if (!res.allFinite()) {
        result.status = ODStatus::invalid_input;
        return result;
    }
    result.raw_residual_norm = res.norm();

    // measurement jacobian
    matXd H = measurement_jacobian(meas.type, ctx);
    if (H.cols() != 6 || H.rows() != dim || !H.allFinite()) {
        result.status = ODStatus::invalid_input;
        return result;
    }

    // measurement covariance
    matXd R;
    if (meas.R.size() == 0) {
        R = matXd::Identity(dim, dim);
    } else {
        R = meas.R;
    }

    // innovation covariance
    matXd S = H * P_pred * H.transpose() + R;
    if (!S.allFinite()) {
        result.status = ODStatus::invalid_input;
        return result;
    }

    // solve for Kalman gain
    Eigen::LDLT<matXd> S_ldlt(S);
    if (S_ldlt.info() != Eigen::Success) {
        result.status = ODStatus::singular_innovation;
        return result;
    }
    matXd X = S_ldlt.solve(H * P_pred); // solve S * X = H * P_pred for X
    // matXd K = P_pred * H.transpose() * S.inverse();
    matXd K = X.transpose();
    if (!X.allFinite() || !K.allFinite()) {
        result.status = ODStatus::singular_innovation;
        return result;
    }

    // a posteriori state and STM, updated estimate
    vec6d dx_vec = K * res;
    DerivTr dx = vec6_to_derivtr(dx_vec);
    StateTr x_post = x_pred + dx;
    // mat6d P_post = (mat6d1 - K * H) * P_pred;
    mat6d P_post = (mat6d1 - K * H) * P_pred * (mat6d1 - K * H).transpose()
                   + K * R * K.transpose();       // Joseph form, more stable
    P_post = 0.5 * (P_post + P_post.transpose()); // force symmetry
    if (!dx_vec.allFinite() || !statetr_to_vec6(x_post).allFinite()
        || !P_post.allFinite()) {
        result.status = ODStatus::correction_rejected;
        return result;
    }

    vecXd weighted_res = S_ldlt.solve(res);
    if (!weighted_res.allFinite()) {
        result.status = ODStatus::singular_innovation;
        return result;
    }
    f64 res_norm2 = res.dot(weighted_res);
    if (!std::isfinite(res_norm2) || res_norm2 < 0.0) {
        result.status = ODStatus::singular_innovation;
        return result;
    }

    // store results
    result.filter.x = x_post;
    result.filter.P = P_post;
    result.filter.t = meas.t;
    result.residual = res;
    result.residual_norm = std::sqrt(res_norm2);
    // result.residual_norm = std::sqrt(res.transpose() * S.inverse() * res);
    result.raw_residual_norm = res.norm();
    result.success = true;
    result.status = ODStatus::ok;
    return result;
}

inline ODEKFResult od_ekf_offline(const ODEKFInput& input) {
    ODEKFResult result;
    ODEKFState filter = input.initial_filter;

    result.status = od_ekf_validate_input(input);
    if (result.status != ODStatus::ok) {
        return result;
    }

    // EKF loop
    for (i32 i = 0; i < input.measurements.size(); ++i) {
        const Measurement& meas = input.measurements[i];
        ODEKFStepInput step_input{
            .filter = filter,
            .measurement = meas,
            .observer_state = input.observer_states[i],
            .dyn_config = input.dyn_config,
            .prop_steps = input.prop_steps,
            .Q = input.Q,
            .tol_time = input.tol_time
        };

        ODEKFStepResult step_result = od_ekf_step(step_input);
        if (!step_result.success) {
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
    result.status = ODStatus::ok;
    result.success = true;
    result.filter = filter;

    return result;
}
