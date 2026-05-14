#include "core/estimation_world.hpp"
#include "core/body.hpp"
#include "core/estimation_common.hpp"
#include "core/estimation_recursive.hpp"
#include "core/measurement.hpp"
#include "core/state.hpp"

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

    if (input.world == nullptr) {
        result.status = ODStatus::invalid_input;
        return result;
    }

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