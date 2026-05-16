#include "core/estimation_world.hpp"
#include "core/body.hpp"
#include "core/estimation_common.hpp"
#include "core/estimation_recursive.hpp"
#include "core/measurement.hpp"
#include "core/measurement_world.hpp"
#include "core/state.hpp"
#include "util/units.hpp"

ODStatus ekf_observer_state_from_world(
    const World& world,
    EntityId observer_id,
    StateTr& x_tr_observer
) {

    const Body* body = world.body(observer_id);
    if (body == nullptr) {
        return ODStatus::observer_not_found;
    }

    if (body->body_type != BodyType::station) {
        x_tr_observer = body->x_tr;
        return ODStatus::ok;
    }

    const Station* stat = world.station(observer_id);
    if (stat == nullptr) {
        return ODStatus::observer_not_found;
    }
    x_tr_observer = world.stat_x_tr_inertial(observer_id);

    return ODStatus::ok;
}

ODEKFStepResult od_ekf_step_world(
    const World& world,
    const ODEKFState& filter,
    const ODWorldMeasurementEvent& event,
    const ODDynamicsConfig& dyn_config,
    i32 prop_steps,
    const mat6d& Q,
    f64 tol_time
) {
    ODEKFStepResult result;

    const Body* target = world.body(event.target_id);
    if (target == nullptr) {
        result.status = ODStatus::target_not_found;
        return result;
    }

    StateTr x_tr_observer;
    result.status
        = ekf_observer_state_from_world(world, event.observer_id, x_tr_observer);
    if (result.status != ODStatus::ok) {
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

    result = od_ekf_step(input);

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
    if (prediction.status != ODStatus::ok) {
        step_result.status = prediction.status;
        return step_result;
    }

    step_result.filter.x = prediction.y.x;
    step_result.filter.P = prediction.P;
    step_result.filter.t = prediction.t;
    step_result.status = ODStatus::prediction_only;

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
            result.status = ODStatus::invalid_input;
            return result;
        }
        if (std::abs(input.t_target - input.event->measurement.t) > input.tol_time) {
            result.status = ODStatus::time_mismatch;
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

ODStatus make_world_measurement_event(
    const World& world,
    ObservationType type,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    const matXd& R,
    ODWorldMeasurementEvent& event,
    UAngle angle_out,
    f64 tol
) {
    vecXd z;
    ODStatus meas_status = world_predict_measurement(
        world,
        type,
        observer_id,
        target_id,
        z,
        UAngle::radian,
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

    return ODStatus::ok;
}

ODStatus validate_realtime_ekf_events(
    const svec<ODRealtimeEvent>& events,
    f64 t_prev,
    f64 tol_time
) {
    if (events.empty()) {
        return ODStatus::empty_events;
    }
    for (const auto event : events) {
        if (event.t < t_prev - tol_time) {
            return ODStatus::time_mismatch;
        }
        if (event.has_measurement) {
            if (event.event.observer_id == kInvalidEntityId) {
                return ODStatus::observer_not_found;
            }
            if (event.event.target_id == kInvalidEntityId) {
                return ODStatus::target_not_found;
            }
            if (event.event.measurement.z.size()
                != measurement_dim(event.event.measurement.type)) {
                return ODStatus::invalid_input;
            }
            if (std::abs(event.event.measurement.t - event.t) > tol_time) {
                return ODStatus::time_mismatch;
            }
        }
    }

    return ODStatus::ok;
}

ODStatus make_world_measurement_event(
    const World& world,
    ObservationType type,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    ODWorldMeasurementEvent& event,
    UAngle angle_out,
    f64 tol
) {
    const Station* observer = world.station(observer_id);
    if (observer == nullptr) {
        return ODStatus::observer_not_found;
    }

    matXd R;
    ODStatus status = stat_meas_cov(*observer, type, R);
    if (status != ODStatus::ok) {
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
        angle_out,
        tol
    );
}

ODStatus make_world_measurement_event_instrument(
    const World& world,
    InstrumentId instrument_id,
    EntityId observer_id,
    EntityId target_id,
    f64 t,
    ODWorldMeasurementEvent& event,
    UAngle angle_out,
    f64 tol
) {
    const Station* observer = world.station(observer_id);
    if (observer == nullptr) {
        return ODStatus::observer_not_found;
    }

    StationInstrument instrument;
    ODStatus status = get_station_instrument(*observer, instrument, instrument_id);
    if (status != ODStatus::ok) {
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
        angle_out,
        tol
    );
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
        result.status = ODStatus::empty_events;
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

        if (step_result.status == ODStatus::ok) {
            ++result.measurement_updates;
        } else if (step_result.status == ODStatus::prediction_only) {
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