#pragma once

#include "core/measurement.hpp"
#include "core/od_dynamics.hpp"
#include "core/state.hpp"
#include "util/constants.hpp"
#include "util/vecdefs.hpp"
#include <functional>

enum struct ODBatchStatus {
    ok,
    invalid_input,
    empty_measurements,
    size_mismatch,
    propagation_failed,
    singular_normal_matrix,
    max_iters_reached,
};

struct ODBatchInput {
    StateTr x0_guess;
    svec<Measurement> measurements;
    svec<StateTr> observer_states;
    ODDynamicsConfig dyn_config;
    f64 t0 = 0.0;
    i32 max_iters = 5;
    i32 prop_steps = 100; // TODO: this is a temp count
    f64 tol_dx = 1e-8;
    f64 tol_residual = 1e-10;
};

struct ODBatchResult {
    bool success = false;
    ODBatchStatus status = ODBatchStatus::invalid_input;
    StateTr x0_est;
    mat6d covariance = mat6d0;
    mat6d normal_inv = mat6d0;
    f64 residual_norm = 0.0;
    f64 dx_norm = 0.0;
    i32 iterations = 0;
};

inline ODBatchStatus od_batch_validate_input(const ODBatchInput& input) {
    if (input.measurements.size() == 0) return ODBatchStatus::empty_measurements;
    if (input.measurements.size() != input.observer_states.size())
        return ODBatchStatus::size_mismatch;
    if (input.dyn_config.mu <= 0.0 || input.prop_steps <= 0 || input.max_iters <= 0)
        return ODBatchStatus::invalid_input;
    return ODBatchStatus::ok;
}

inline ODBatchResult od_batch_lumve(const ODBatchInput& input) {
    ODBatchResult result;
    result.x0_est = input.x0_guess;

    result.status = od_batch_validate_input(input);
    if (result.status != ODBatchStatus::ok) {
        result.success = false;
        return result;
    } else {
        result.success = true;
        result.status = ODBatchStatus::ok;
        result.iterations = 0;
        result.residual_norm = 0;
        result.dx_norm = 0;
    }

    // iteration variables
    StateTr x0_ref = input.x0_guess;
    result.x0_est = x0_ref;
    result.status = ODBatchStatus::max_iters_reached;
    result.success = false;

    vecXd residual_stack;
    matXd H_stack;
    mat6d normal;
    vec6d rhs;
    vec6d dx0;
    for (i32 iter = 0; iter < input.max_iters; ++iter) {
        result.iterations = iter + 1;

        // build batch system
        i32 total_dim = 0;

        for (i32 i = 0; i < input.measurements.size(); ++i) {
            Measurement meas = input.measurements[i];
            i32 dim = measurement_dim(meas.type);
            if (dim <= 0) {
                result.status = ODBatchStatus::invalid_input;
                return result;
            }
            if (meas.z.size() != dim) {
                result.status = ODBatchStatus::invalid_input;
                return result;
            }
            // measurement covariance guards
            if (meas.R.size() != 0 && (meas.R.cols() != dim || meas.R.rows() != dim)) {
                result.status = ODBatchStatus::invalid_input;
                return result;
            }
            if (!meas.R.allFinite()) {
                result.status = ODBatchStatus::invalid_input;
                return result;
            }
            total_dim += dim;
        }
        if (total_dim < 6) { // less than state size
            result.status = ODBatchStatus::singular_normal_matrix;
            return result;
        } // end measurement dimension loop

        // stacked residuals and design matrix (observation matrix)
        // vecXd y = vecXd::Zero(total_dim);
        // matXd H = matXd::Zero(total_dim, 6);
        vec6d lambda = vec6d0; // normal matrices
        mat6d Lambda = mat6d0;
        f64 residual_norm2 = 0.0;

        i32 row0 = 0;
        for (i32 i = 0; i < input.measurements.size(); ++i) {
            const Measurement& meas = input.measurements[i];
            const StateTr& x_obsv = input.observer_states[i];
            i32 dim = measurement_dim(meas.type);

            // TODO: use solver instead of invert
            matXd W_i; // weight matrix
            if (meas.R.size() == 0) {
                W_i = matXd::Identity(dim, dim);
            } else {
                W_i = meas.R.inverse(); // weight matrix
            }

            // initiate variational state and propagate
            f64 dt = meas.t - input.t0;
            VarStateTr y0;
            y0.x = x0_ref;
            y0.Phi = mat6d::Identity();
            VarStateTr yf = propagate_var_tr_od_rk4(
                input.t0,
                y0,
                dt,
                input.prop_steps,
                input.dyn_config
            );

            vec6d xf_vec = statetr_to_vec6d(yf.x);
            if (!xf_vec.allFinite() || !yf.Phi.allFinite()) {
                result.status = ODBatchStatus::propagation_failed;
                return result;
            }

            // measurement prediction
            MeasurementContext ctx;
            ctx.x_target = yf.x;
            ctx.x_observer = x_obsv;
            vecXd z_pred = predict_measurement(meas.type, ctx);
            if (z_pred.size() != dim) {
                result.status = ODBatchStatus::invalid_input;
                return result;
            }

            // residuals
            vecXd res_i = measurement_residual(meas.type, meas.z, z_pred);
            if (!res_i.allFinite()) {
                result.status = ODBatchStatus::invalid_input;
                return result;
            }

            // jacobian and initial state jacobian
            matXd G_i = measurement_jacobian(meas.type, ctx);
            // Htilde in matlab code (see scratch/homework 06)
            if (G_i.rows() != dim || G_i.cols() != 6 || !G_i.allFinite()) {
                result.status = ODBatchStatus::invalid_input;
                return result;
            }
            matXd H_i = G_i * yf.Phi;

            lambda += H_i.transpose() * W_i * res_i;
            Lambda += H_i.transpose() * W_i * H_i;
            residual_norm2 += res_i.transpose() * W_i * res_i;

            // store residual and jacobian
            // y.segment(row0, dim) = res_i;
            // H.block(row0, 0, dim, 6) = H_i;
            row0 += dim;
        } // end measurement loop

        if (row0 != total_dim) {
            result.status = ODBatchStatus::invalid_input;
            return result;
        }

        // normal matrices
        if (!Lambda.allFinite() || Lambda.diagonal().cwiseAbs().minCoeff() <= tol_tight) {
            result.status = ODBatchStatus::singular_normal_matrix;
            return result;
        }

        // TODO: use CompleteOrthogonalDecomposition later for robustness
        // use ldlh decomp for solving (see sOPT)
        Eigen::LDLT<mat6d> ldlt(Lambda);
        if (ldlt.info() != Eigen::Success) {
            result.status = ODBatchStatus::singular_normal_matrix;
            return result;
        }
        vec6d dx = ldlt.solve(lambda);
        if (!dx.allFinite()) {
            result.status = ODBatchStatus::singular_normal_matrix;
            return result;
        }

        // TODO: add line search
        // f64 alpha = 0.25;
        // StateTr x0_cand = x0_ref + alpha * vec6_to_derivtr(dx);
        // x0_ref = x0_cand;
        x0_ref += vec6_to_statetr(dx);

        // store result
        result.iterations = iter + 1;
        result.dx_norm = dx.norm();
        result.x0_est = x0_ref;
        result.residual_norm = sqrt(residual_norm2);
        result.normal_inv = Lambda.inverse();
        result.covariance = result.normal_inv;

        if (result.dx_norm < input.tol_dx || result.residual_norm < input.tol_residual) {
            // tolerances should be looser for angle-only measurements
            result.success = true;
            result.status = ODBatchStatus::ok;
            return result;
        }

    } // end lumve iterations

    // max iterations reached
    result.status = ODBatchStatus::max_iters_reached;
    result.x0_est = x0_ref;

    return result;
}