// Copyright 2025-2026 Joseph Tu Le
// SPDX-License-Identifier: Apache-2.0

#include "core/estimation_batch.hpp"
#include "core/estimation_common.hpp"
#include "core/measurement.hpp"
#include "core/od_dynamics.hpp"

#include <Eigen/Core>
#include <cmath>

StatusCode od_batch_validate_input(const ODBatchInput& input) {
    if (input.measurements.size() == 0) return StatusCode::empty_measurements;
    if (input.measurements.size() != input.observer_states.size()) {
        return StatusCode::size_mismatch;
    }
    if (input.dyn_config.mu <= 0.0 || input.prop_steps <= 0 || input.max_iters <= 0)
        return StatusCode::validation_failed;

    if (input.dyn_config.zonal_degree < 0 || input.dyn_config.zonal_degree > 6) {
        return StatusCode::validation_failed;
    }

    return StatusCode::ok;
}

ODBatchResidualEval od_batch_eval_residual_norm(
    const ODBatchInput& input,
    const StateTr& x0_ref
) {
    // helper to compute residual norms (weighted & raw)
    ODBatchResidualEval eval;
    eval.status = od_batch_validate_input(input);
    if (eval.status != StatusCode::ok) {
        return eval;
    }
    f64 weighted_norm2 = 0.0;
    f64 raw_norm2 = 0.0;

    for (i32 i = 0; i < input.measurements.size(); ++i) {
        const Measurement& meas = input.measurements[i];
        const StateTr& x_obsv = input.observer_states[i];
        i32 dim = measurement_dim(meas.type);
        if (dim <= 0 || meas.z.size() != dim) {
            eval.status = dim <= 0 ? StatusCode::unsupported_type
                                    : StatusCode::size_mismatch;
            return eval;
        }

        if (meas.R.size() != 0 && (meas.R.cols() != dim || meas.R.rows() != dim)) {
            eval.status = StatusCode::invalid_covariance;
            return eval;
        }
        if (!meas.R.allFinite()) {
            eval.status = StatusCode::invalid_covariance;
            return eval;
        }

        // propagate candidate state to measurement time
        f64 dt = meas.t - input.t0;
        StateTr x_pred_i
            = propagate_tr_od(input.t0, x0_ref, dt, input.prop_steps, input.dyn_config);
        if (!statetr_to_vec6d(x_pred_i).allFinite()) {
            eval.status = StatusCode::propagation_failed;
            return eval;
        }

        // measurement prediction
        MeasurementContext ctx = make_measurement_context(x_pred_i, x_obsv);
        vecXd z_pred = predict_measurement(meas.type, ctx);
        if (z_pred.size() != dim) {
            eval.status = StatusCode::size_mismatch;
            return eval;
        }

        // residuals
        vecXd res_i = measurement_residual(meas.type, meas.z, z_pred);
        if (!res_i.allFinite()) {
            eval.status = StatusCode::invalid_state;
            return eval;
        }

        raw_norm2 += res_i.squaredNorm();

        // weight residuals using measurement covariance
        // empty measurement covariance uses identity
        if (meas.R.size() == 0) {
            weighted_norm2 += res_i.squaredNorm();
        } else {
            Eigen::LDLT<matXd> R_ldlt(meas.R);
            if (R_ldlt.info() != Eigen::Success) {
                eval.status = StatusCode::invalid_covariance;
                return eval;
            }

            vecXd weighted_res = R_ldlt.solve(res_i);
            if (!weighted_res.allFinite()) {
                eval.status = StatusCode::invalid_covariance;
                return eval;
            }

            f64 weighted_norm2_i = res_i.dot(weighted_res);
            if (weighted_norm2_i < 0.0) {
                eval.status = StatusCode::invalid_covariance;
                return eval;
            }
            weighted_norm2 += weighted_norm2_i;
        }
    }

    eval.weighted_norm = std::sqrt(weighted_norm2);
    eval.raw_norm = std::sqrt(raw_norm2);
    eval.status = StatusCode::ok;
    return eval;
}

ODBatchResult od_batch_lumve(const ODBatchInput& input) {
    // TODO: switch for different measurement types
    ODBatchResult result;
    result.x0_est = input.x0_guess;

    result.status = od_batch_validate_input(input);
    if (result.status != StatusCode::ok) {
        return result;
    } else {
        result.status = StatusCode::ok;
        result.iterations = 0.0;
        result.residual_norm = 0.0;
        result.dx_norm = 0.0;
    }

    // iteration variables
    StateTr x0_ref = input.x0_guess;
    result.x0_est = x0_ref;
    result.status = StatusCode::max_iters_reached;

    for (i32 iter = 0; iter < input.max_iters; ++iter) {
        result.iterations = iter + 1;

        // build batch system
        i32 total_dim = 0;

        for (i32 i = 0; i < input.measurements.size(); ++i) {
            const Measurement& meas = input.measurements[i];
            i32 dim = measurement_dim(meas.type);
            if (dim <= 0) {
                result.status = StatusCode::invalid_input;
                return result;
            }
            if (meas.z.size() != dim) {
                result.status = StatusCode::invalid_input;
                return result;
            }
            // measurement covariance guards
            if (meas.R.size() != 0 && (meas.R.cols() != dim || meas.R.rows() != dim)) {
                result.status = StatusCode::invalid_input;
                return result;
            }
            if (!meas.R.allFinite()) {
                result.status = StatusCode::invalid_input;
                return result;
            }
            total_dim += dim;
        }
        if (total_dim < 6) { // less than state size
            result.status = StatusCode::singular_normal_matrix;
            return result;
        } // end measurement dimension loop

        vec6d lambda = vec6d0; // normal matrices
        mat6d Lambda = mat6d0;
        f64 residual_norm2 = 0.0;
        f64 raw_residual_norm2 = 0.0;

        i32 row0 = 0;
        for (i32 i = 0; i < input.measurements.size(); ++i) {
            const Measurement& meas = input.measurements[i];
            const StateTr& x_obsv = input.observer_states[i];
            i32 dim = measurement_dim(meas.type);

            // initiate variational state and propagate
            f64 dt = meas.t - input.t0;
            VarStateTr y0;
            y0.x = x0_ref;
            y0.Phi = mat6d1;
            VarStateTr yf = propagate_var_tr_od(
                input.t0,
                y0,
                dt,
                input.prop_steps,
                input.dyn_config
            );

            vec6d xf_vec = statetr_to_vec6d(yf.x);
            if (!xf_vec.allFinite() || !yf.Phi.allFinite()) {
                result.status = StatusCode::propagation_failed;
                return result;
            }

            // measurement prediction
            MeasurementContext ctx = make_measurement_context(yf.x, x_obsv);

            vecXd z_pred = predict_measurement(meas.type, ctx);
            if (z_pred.size() != dim) {
                result.status = StatusCode::invalid_input;
                return result;
            }

            // residuals
            vecXd res_i = measurement_residual(meas.type, meas.z, z_pred);
            if (!res_i.allFinite()) {
                result.status = StatusCode::invalid_input;
                return result;
            }

            // jacobian and initial state jacobian
            matXd G_i = measurement_jacobian(meas.type, ctx);
            // Htilde in matlab code (see scratch/homework 06)
            if (G_i.rows() != dim || G_i.cols() != 6 || !G_i.allFinite()) {
                result.status = StatusCode::invalid_input;
                return result;
            }
            matXd H_i = G_i * yf.Phi;

            vecXd weighted_res;
            matXd weighted_H;
            if (meas.R.size() == 0) {
                weighted_res = res_i;
                weighted_H = H_i;
            } else {
                Eigen::LDLT<matXd> R_ldlt(meas.R);
                if (R_ldlt.info() != Eigen::Success) {
                    result.status = StatusCode::invalid_covariance;
                    return result;
                }
                weighted_res = R_ldlt.solve(res_i);
                weighted_H = R_ldlt.solve(H_i);
            }
            if (!weighted_res.allFinite() || !weighted_H.allFinite()) {
                result.status = StatusCode::invalid_covariance;
                return result;
            }
            f64 weighted_norm2_i = res_i.dot(weighted_res);
            if (weighted_norm2_i < 0.0) {
                result.status = StatusCode::invalid_covariance;
                return result;
            }

            lambda += H_i.transpose() * weighted_res;
            Lambda += H_i.transpose() * weighted_H;
            residual_norm2 += weighted_norm2_i;
            raw_residual_norm2 += res_i.squaredNorm();
            // store residual and jacobian
            // y.segment(row0, dim) = res_i;
            // H.block(row0, 0, dim, 6) = H_i;
            row0 += dim;
        } // end measurement loop

        if (row0 != total_dim) {
            result.status = StatusCode::invalid_input;
            return result;
        }

        f64 current_weighted_norm = std::sqrt(residual_norm2);
        f64 current_raw_norm = std::sqrt(raw_residual_norm2);

        // normal matrices
        if (!Lambda.allFinite() || Lambda.diagonal().cwiseAbs().minCoeff() <= tol9) {
            result.status = StatusCode::singular_normal_matrix;
            return result;
        }

        // TODO: use CompleteOrthogonalDecomposition later for robustness
        // use ldlh decomp for solving (see sOPT)
        Eigen::LDLT<mat6d> ldlt(Lambda);
        if (ldlt.info() != Eigen::Success) {
            result.status = StatusCode::singular_normal_matrix;
            return result;
        }
        vec6d dx_vec = ldlt.solve(lambda);
        if (!dx_vec.allFinite()) {
            result.status = StatusCode::singular_normal_matrix;
            return result;
        }
        StateTr dx = vec6d_to_statetr(dx_vec);
        if (dx.r.norm() > input.max_dx_r_norm || dx.v.norm() > input.max_dx_v_norm) {
            result.status = StatusCode::correction_rejected;
            return result;
        }

        StateTr x0_next = x0_ref;
        ODBatchResidualEval next_eval{
            .status = StatusCode::ok,
            .weighted_norm = current_weighted_norm,
            .raw_norm = current_raw_norm
        };
        vec6d dx_accepted_vec = dx_vec;
        if (input.use_line_search) {
            bool accepted = false;
            for (i32 attempt = 0; attempt < 8; ++attempt) {
                // TODO: add attempts as option
                f64 alpha = 1.0 / std::pow(2.0, attempt);
                vec6d dx_cand_vec = alpha * dx_vec;
                StateTr x0_cand = x0_ref + vec6d_to_statetr(dx_cand_vec);
                ODBatchResidualEval cand_eval
                    = od_batch_eval_residual_norm(input, x0_cand);
                if (!od_status_success(cand_eval.status)) {
                    continue;
                }
                if (cand_eval.weighted_norm < current_weighted_norm) {
                    accepted = true;
                    x0_next = x0_cand;
                    next_eval = cand_eval;
                    dx_accepted_vec = dx_cand_vec;
                    break;
                }
            }
            if (!accepted) {
                result.residual_norm = current_weighted_norm;
                result.raw_residual_norm = current_raw_norm;
                result.x0_est = x0_ref;
                result.normal_inv = Lambda.inverse();
                result.covariance = result.normal_inv;
                if (iter > 0 && statetr_to_vec6d(result.x0_est).allFinite()
                    && result.covariance.allFinite()) {
                    // residual probably already in noise floor so accept
                    // TODO: NOTE: this isn't totally valid, fix later, add stalled or
                    // local min status
                    result.status = StatusCode::ok;
                    return result;
                }
                result.status = StatusCode::correction_rejected;
                return result;
            }
        } else {
            x0_next = x0_ref + dx;
        }
        x0_ref = x0_next;

        // store result
        result.iterations = iter + 1;
        result.dx_norm = dx_accepted_vec.norm();
        result.x0_est = x0_ref;
        result.residual_norm = next_eval.weighted_norm;
        result.raw_residual_norm = next_eval.raw_norm;
        result.normal_inv = Lambda.inverse();
        result.covariance = result.normal_inv;

        if (result.dx_norm < input.tol_dx || result.residual_norm < input.tol_residual) {
            // tolerances should be looser for angle-only measurements
            result.status = StatusCode::ok;
            return result;
        }

    } // end lumve iterations

    // max iterations reached
    result.status = StatusCode::max_iters_reached;
    result.x0_est = x0_ref;

    return result;
}
