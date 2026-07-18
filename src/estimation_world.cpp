// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/estimation_world.hpp"
#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/estimation_common.hpp"
#include "core/estimation_recursive.hpp"
#include "core/interpolation.hpp"
#include "core/measurement.hpp"
#include "core/measurement_world.hpp"
#include "core/observation_type.hpp"
#include "core/state.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

StatusCode ekf_observer_state_from_world(
    const World& world,
    EntityId observer_id,
    StateTr& x_tr_observer
) {

    const Body* body = world.body(observer_id);
    if (body == nullptr) {
        return StatusCode::observer_not_found;
    }

    if (body->body_type != BodyType::station) {
        x_tr_observer = body->x_tr;
        return StatusCode::ok;
    }

    const Station* stat = world.station(observer_id);
    if (stat == nullptr) {
        return StatusCode::observer_not_found;
    }
    x_tr_observer = world.stat_x_tr_inertial(observer_id);

    return StatusCode::ok;
}

ODEKFStepResult od_ekf_step_world(
    const World& world,
    const ODEKFState& filter,
    const ODWorldMeasurementEvent& event,
    const ODDynamicsConfig& dyn_config,
    i32 prop_steps,
    const mat6d& Q,
    f64 tol_time,
    UAngle angle_in,
    UAngle angle_out,
    f64 eps_pos,
    f64 eps_vel,
    f64 tol
) {
    ODEKFStepResult result;

    const Body* target = world.body(event.target_id);
    if (target == nullptr) {
        result.status = StatusCode::target_not_found;
        return result;
    }

    StateTr x_tr_observer;
    result.status
        = ekf_observer_state_from_world(world, event.observer_id, x_tr_observer);
    if (result.status != StatusCode::ok) {
        return result;
    }

    ODEKFStepInput input{
        .filter = filter,
        .measurement = event.measurement,
        .x_tr_observer = x_tr_observer,
        .dyn_config = dyn_config,
        .prop_steps = prop_steps,
        .Q = Q,
        .tol_time = tol_time
    };
    result.status = od_ekf_step_validate_input(input);

    // result = od_ekf_step(input);
    if (!od_status_success(result.status)) {
        return result;
    } else {
        result.filter = filter;
        result.status = StatusCode::ok;
        result.residual_norm = 0.0;
        result.raw_residual_norm = 0.0;
    }

    ODEKFPredictResult prediction = od_ekf_predict(
        filter,
        event.measurement.t,
        input.dyn_config,
        input.prop_steps,
        input.Q,
        input.tol_time
    );
    if (!od_status_success(prediction.status)) {
        result.status = prediction.status;
        return result;
    }

    StateTr& x_pred = prediction.y.x;
    mat6d& P_pred = prediction.P;

    // predicted measurement
    vecXd z_pred;
    StatusCode pred_status = world_predict_measurement_from_state(
        world,
        event.measurement.type,
        event.observer_id,
        x_pred,
        dyn_config,
        event.measurement.t - dyn_config.t0,
        z_pred,
        angle_in,
        angle_out,
        tol
    );
    if (!od_status_success(pred_status)) {
        result.status = pred_status;
        return result;
    }

    // residuals
    vecXd res = measurement_residual(
        event.measurement.type,
        event.measurement.z,
        z_pred,
        angle_out
    );
    if (!res.allFinite()) {
        result.status = StatusCode::invalid_input;
        return result;
    }

    // measurement jacobian
    i32 dim = measurement_dim(event.measurement.type);
    matXd H;
    StatusCode status = world_jacobian_measurement(
        world,
        event.measurement.type,
        event.observer_id,
        x_pred,
        dyn_config,
        event.measurement.t - dyn_config.t0,
        H,
        angle_in,
        angle_out,
        eps_pos,
        eps_vel,
        tol
    );
    if (!od_status_success(status)) {
        result.status = status;
        return result;
    }
    if (H.cols() != 6 || H.rows() != dim || !H.allFinite()) {
        result.status = StatusCode::invalid_input;
        return result;
    }

    // measurement covariance
    matXd R;
    if (event.measurement.R.size() == 0) {
        R = matXd::Identity(dim, dim);
    } else {
        R = event.measurement.R;
    }

    // innovation covariance
    matXd S = H * P_pred * H.transpose() + R;
    if (!S.allFinite()) {
        result.status = StatusCode::invalid_input;
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
    result.filter.t = event.measurement.t;
    result.residual = res;
    result.residual_norm = std::sqrt(res_norm2);
    // result.residual_norm = std::sqrt(res.transpose() * S.inverse() * res);
    result.raw_residual_norm = res.norm();
    result.status = StatusCode::ok;

    return result;
}

ODEKFStepResult od_ekf_predict_step(
    const ODEKFState& filter,
    f64 t_target,
    const ODDynamicsConfig& dyn_config,
    i32 prop_steps,
    const mat6d& Q,
    f64 tol_time
) {

    ODEKFStepResult step_result;
    step_result.filter = filter;

    ODEKFPredictResult prediction
        = od_ekf_predict(filter, t_target, dyn_config, prop_steps, Q, tol_time);
    if (prediction.status != StatusCode::ok) {
        step_result.status = prediction.status;
        return step_result;
    }

    step_result.filter.x = prediction.y.x;
    step_result.filter.P = prediction.P;
    step_result.filter.t = prediction.t;
    step_result.status = StatusCode::prediction_only;

    return step_result;
}

ODEKFStepResult od_ekf_update_world(const ODRealtimeEKFInput& input) {
    ODEKFStepResult result;

    if (input.event == nullptr) {
        // prediction only
        result = od_ekf_predict_step(
            input.filter,
            input.t_target,
            input.dyn_config,
            input.prop_steps,
            input.Q,
            input.tol_time
        );
    } else {
        // estimate
        if (input.world == nullptr) {
            result.status = StatusCode::invalid_input;
            return result;
        }
        if (std::abs(input.t_target - input.event->measurement.t) > input.tol_time) {
            result.status = StatusCode::time_mismatch;
            return result;
        }
        result = od_ekf_step_world(
            *input.world,
            input.filter,
            *input.event,
            input.dyn_config,
            input.prop_steps,
            input.Q,
            input.tol_time
        );
    }

    return result;
}

StatusCode make_world_measurement_event(
    const World& world,
    ObservationType type,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    const matXd& R,
    ODWorldMeasurementEvent& event,
    UAngle angle_in,
    UAngle angle_out,
    f64 tol
) {
    vecXd z;
    StatusCode meas_status = world_predict_measurement(
        world,
        type,
        observer_id,
        target_id,
        z,
        angle_in,
        angle_out,
        tol
    );
    if (!od_status_success(meas_status)) {
        return meas_status;
    }

    Measurement meas;
    meas.z = z;
    meas.t = t;
    meas.R = R;
    meas.type = type;
    meas.observer_id = observer_id;
    meas.target_id = target_id;

    event.measurement = meas;
    event.observer_id = observer_id;
    event.target_id = target_id;

    return StatusCode::ok;
}

StatusCode validate_realtime_ekf_events(
    const svec<ODRealtimeEvent>& events,
    f64 t_prev,
    f64 tol_time
) {
    if (events.empty()) {
        return StatusCode::empty_events;
    }
    for (const auto event : events) {
        if (event.t < t_prev - tol_time) {
            return StatusCode::time_mismatch;
        }
        if (event.has_measurement) {
            if (event.event.observer_id == kInvalidEntityId) {
                return StatusCode::observer_not_found;
            }
            if (event.event.target_id == kInvalidEntityId) {
                return StatusCode::target_not_found;
            }
            if (event.event.measurement.z.size()
                != measurement_dim(event.event.measurement.type)) {
                return StatusCode::invalid_input;
            }
            if (std::abs(event.event.measurement.t - event.t) > tol_time) {
                return StatusCode::time_mismatch;
            }
        }
    }

    return StatusCode::ok;
}

StatusCode make_world_measurement_event(
    const World& world,
    ObservationType type,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    ODWorldMeasurementEvent& event,
    UAngle angle_in,
    UAngle angle_out,
    f64 tol
) {
    const Station* observer = world.station(observer_id);
    if (observer == nullptr) {
        return StatusCode::observer_not_found;
    }

    matXd R;
    StatusCode status = station_measurement_covariance(*observer, type, R);
    if (status != StatusCode::ok) {
        return status;
    }

    return make_world_measurement_event(
        world,
        type,
        observer_id,
        target_id,
        t,
        R,
        event,
        angle_in,
        angle_out,
        tol
    );
}

StatusCode make_world_measurement_event_instrument(
    const World& world,
    InstrumentId instrument_id,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    ODWorldMeasurementEvent& event,
    UAngle angle_in,
    UAngle angle_out,
    f64 tol
) {
    const Station* observer = world.station(observer_id);
    if (observer == nullptr) {
        return StatusCode::observer_not_found;
    }

    StationInstrument instrument;
    StatusCode status = get_station_instrument(*observer, instrument, instrument_id);
    if (status != StatusCode::ok) {
        return status;
    }

    return make_world_measurement_event(
        world,
        instrument.type,
        observer_id,
        target_id,
        t,
        instrument.R,
        event,
        angle_in,
        angle_out,
        tol
    );
}

StatusCode make_noisy_world_measurement_event_instrument(
    const World& world,
    InstrumentId instrument_id,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    ODWorldMeasurementEvent& event,
    const MeasurementNoiseOptions& noise_opts,
    UAngle angle_in,
    UAngle angle_out,
    f64 tol
) {
    const Station* observer = world.station(observer_id);
    if (observer == nullptr) {
        return StatusCode::observer_not_found;
    }

    StationInstrument instrument;
    StatusCode status = get_station_instrument(*observer, instrument, instrument_id);
    if (status != StatusCode::ok) {
        return status;
    }

    status = make_world_measurement_event(
        world,
        instrument.type,
        observer_id,
        target_id,
        t,
        instrument.R,
        event,
        angle_in,
        angle_out,
        tol
    );
    if (status != StatusCode::ok) {
        return status;
    }

    if (noise_opts.enabled) {
        if (noise_opts.diagonal) {
            return apply_measurement_noise_diagonal(
                event.measurement,
                noise_opts,
                instrument.type
            );
        } else {
            return apply_measurement_noise_cholesky(
                event.measurement,
                noise_opts,
                instrument.type
            );
        }
    }

    return StatusCode::ok;
}

StatusCode make_world_measurement_event_history(
    const World& world,
    const WorldHistory& history,
    ObservationType type,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    const matXd& R,
    ODWorldMeasurementEvent& event,
    const HistorySampleOptions& sample_opts,
    UAngle angle_in,
    UAngle angle_out,
    f64 tol
) {
    vecXd z_pred;
    StatusCode status = world_predict_measurement_history(
        world,
        history,
        type,
        observer_id,
        target_id,
        t,
        z_pred,
        sample_opts,
        angle_in,
        angle_out,
        tol
    );
    if (!od_status_success(status)) return status;
    i32 dim = measurement_dim(type);
    if (z_pred.size() != dim) return StatusCode::size_mismatch;

    const Station* observer = world.station(observer_id);
    if (observer == nullptr) {
        return StatusCode::observer_not_found;
    }

    Measurement meas;
    meas.t = t;
    meas.z = z_pred;
    meas.R = R;
    meas.observer_id = observer_id;
    meas.target_id = target_id;
    meas.type = type;

    event.measurement = meas;
    event.observer_id = observer_id;
    event.target_id = target_id;

    return StatusCode::ok;
}

StatusCode make_world_measurement_event_history_instrument(
    const World& world,
    const WorldHistory& history,
    InstrumentId instrument_id,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    ODWorldMeasurementEvent& event,
    const HistorySampleOptions& sample_opts,
    UAngle angle_in,
    UAngle angle_out,
    f64 tol
) {
    const Station* observer = world.station(observer_id);
    if (observer == nullptr) {
        return StatusCode::observer_not_found;
    }

    StationInstrument instrument;
    StatusCode status = get_station_instrument(*observer, instrument, instrument_id);
    if (status != StatusCode::ok) {
        return status;
    }

    status = make_world_measurement_event_history(
        world,
        history,
        instrument.type,
        observer_id,
        target_id,
        t,
        instrument.R,
        event,
        sample_opts,
        angle_in,
        angle_out,
        tol
    );
    if (status != StatusCode::ok) {
        return status;
    }

    return StatusCode::ok;
}

StatusCode make_noisy_world_measurement_event_history_instrument(
    const World& world,
    const WorldHistory& history,
    InstrumentId instrument_id,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    ODWorldMeasurementEvent& event,
    const MeasurementNoiseOptions& noise_opts,
    const HistorySampleOptions& sample_opts,
    UAngle angle_in,
    UAngle angle_out,
    f64 tol
) {
    StatusCode status = make_world_measurement_event_history_instrument(
        world,
        history,
        instrument_id,
        observer_id,
        target_id,
        t,
        event,
        sample_opts,
        angle_in,
        angle_out,
        tol
    );
    if (status != StatusCode::ok) {
        return status;
    }

    const Station* observer = world.station(observer_id);
    if (observer == nullptr) {
        return StatusCode::observer_not_found;
    }

    StationInstrument instrument;
    status = get_station_instrument(*observer, instrument, instrument_id);
    if (status != StatusCode::ok) {
        return status;
    }

    if (noise_opts.enabled) {
        if (noise_opts.diagonal) {
            return apply_measurement_noise_diagonal(
                event.measurement,
                noise_opts,
                instrument.type
            );
        } else {
            return apply_measurement_noise_cholesky(
                event.measurement,
                noise_opts,
                instrument.type
            );
        }
    }

    return StatusCode::ok;
}

ODRealtimeEKFResult od_ekf_update_world_events(
    const World& world,
    const ODEKFState& initial_filter,
    const svec<ODRealtimeEvent>& events,
    const ODDynamicsConfig& dyn_config,
    i32 prop_steps,
    const mat6d& Q,
    f64 tol_time
) {
    ODRealtimeEKFResult result;
    result.filter = initial_filter;

    if (events.empty()) {
        result.status = StatusCode::empty_events;
        return result;
    }

    for (const ODRealtimeEvent& event : events) {
        ODRealtimeEKFInput input{
            .world = &world,
            .event = event.has_measurement ? &event.event : nullptr
        };
        input.dyn_config = dyn_config;
        input.Q = Q;
        input.prop_steps = prop_steps;
        input.filter = result.filter;
        input.t_target = event.t;
        input.tol_time = tol_time;

        ODEKFStepResult step_result = od_ekf_update_world(input);
        if (!od_status_success(step_result.status)) {
            result.status = step_result.status;
            result.raw_residual_norm = step_result.raw_residual_norm;
            result.residual_norm = step_result.residual_norm;
            result.filter = step_result.filter;
            return result;
        }

        if (step_result.status == StatusCode::ok) {
            ++result.measurement_updates;
        } else if (step_result.status == StatusCode::prediction_only) {
            ++result.prediction_updates;
        }

        result.filter = step_result.filter;
        result.status = step_result.status;
        result.raw_residual_norm = step_result.raw_residual_norm;
        result.residual_norm = step_result.residual_norm;
        ++result.processed_events;
    }

    return result;
}
