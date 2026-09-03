// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/od_workspace.hpp"
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <typeinfo>

StatusCode ODWorldWorkspace::initialize(
    const ScenarioConfig& scenario, const string& target, const svec<string>& sources
) {
    StatusCode status = validate_scenario_config(scenario);
    if (status != StatusCode::ok) return status;
    ODWorldWorkspace temp;
    status = build_world_from_scenario_config(scenario, temp.world_, temp.bindings_, temp.stepper_);
    if (status != StatusCode::ok) return status;
    auto target_it = temp.bindings_.body_ids.find(target);
    if (target_it == temp.bindings_.body_ids.end()) return StatusCode::target_not_found;
    temp.target_id_ = target_it->second;
    auto* sat = temp.world_.satellite(temp.target_id_);
    if (!sat || !temp.world_.is_active(temp.target_id_) || sat->emits_gravity)
        return StatusCode::invalid_input;
    svec<EntityId> source_ids;
    for (const auto& name : sources) {
        auto it = temp.bindings_.celestial_ids.find(name);
        if (it == temp.bindings_.celestial_ids.end()) return StatusCode::body_not_found;
        if (!temp.world_.is_active(it->second)
            || std::find(source_ids.begin(), source_ids.end(), it->second) != source_ids.end())
            return StatusCode::invalid_input;
        const auto* cel = temp.world_.celestial(it->second);
        if (cel->gravity_model != GravityModel::pointmass && cel->gravity_model != GravityModel::zonal)
            return StatusCode::unsupported_method;
        source_ids.push_back(it->second);
    }
    // estimated translation must not be overwritten by a prescribed target ephemeris
    sat->ephemeris_providers.translation.reset();
    sat->propagate_tr = true;
    temp.stepper_.paused = false;
    temp.stepper_.ticks = 1;
    temp.stepper_.dt_scale = 1.0;
    temp.stepper_.step_tr = true;
    for (EntityId id : temp.world_.all_entity_ids()) {
        Body* body = temp.world_.body(id);
        body->emits_gravity = std::find(source_ids.begin(), source_ids.end(), id) != source_ids.end();
        switch (body->body_type) {
        case BodyType::celestial: temp.reference_bodies_.emplace_back(*temp.world_.celestial(id)); break;
        case BodyType::satellite: temp.reference_bodies_.emplace_back(*temp.world_.satellite(id)); break;
        case BodyType::station: temp.reference_bodies_.emplace_back(*temp.world_.station(id)); break;
        default: return StatusCode::unsupported_type;
        }
    }
    temp.reference_active_ids_ = temp.world_.active_entity_ids();
    temp.reference_stepper_ = temp.stepper_;
    temp.t0_ = temp.world_.t_sim();
    temp.epoch_ = temp.world_.get_date_jd();
    temp.time_scale_ = temp.world_.get_time_scale();
    temp.offsets_ = temp.world_.get_time_offsets();
    temp.initialized_ = true;
    rebuild_world_stepper_workspace(temp.world_, temp.workspace_);
    *this = std::move(temp);
    return StatusCode::ok;
}

StatusCode ODWorldWorkspace::reset(const StateTr& x_target0_I) {
    if (!initialized_ || !x_target0_I.r.allFinite() || !x_target0_I.v.allFinite())
        return StatusCode::invalid_input;
    if (world_.all_entity_ids().size() != reference_bodies_.size())
        return StatusCode::invalid_input;
    // reject topology/type edits before touching any state
    for (const auto& entry : reference_bodies_) {
        bool valid = std::visit([&](const auto& ref) {
            const Body* body = world_.body(ref.id);
            return body && typeid(*body) == typeid(ref);
        }, entry);
        if (!valid) return StatusCode::invalid_input;
    }
    for (const auto& entry : reference_bodies_) {
        std::visit([&](const auto& ref) {
            using T = std::decay_t<decltype(ref)>;
            *static_cast<T*>(world_.body(ref.id)) = ref;
        }, entry);
    }
    // recreate the original active ordering, independent of intervening toggles
    auto active = world_.active_entity_ids();
    for (EntityId id : active) world_.make_inactive(id);
    for (EntityId id : reference_active_ids_) world_.make_active(id);
    world_.reset_time(t0_);
    world_.set_date(epoch_, time_scale_);
    world_.set_time_offsets(offsets_);
    world_.body(target_id_)->x_tr = x_target0_I;
    stepper_ = reference_stepper_;
    rebuild_world_stepper_workspace(world_, workspace_);
    return StatusCode::ok;
}

ODWorldPropagationResult ODWorldWorkspace::propagate(
    f64 t_target, const ODWorldPropagationOptions& options
) {
    ODWorldPropagationResult result;
    result.t = world_.t_sim();
    if (!initialized_) return result;
    const Body* target = world_.body(target_id_);
    if (!target) { result.status = StatusCode::target_not_found; return result; }
    result.x = target->x_tr;
    if (!std::isfinite(t_target) || !std::isfinite(options.max_step)
        || options.max_step <= 0.0 || options.max_steps < 1
        || stepper_.paused || !stepper_.step_tr || stepper_.ticks != 1
        || stepper_.dt_scale != 1.0 || workspace_.target_stm != nullptr)
        return result;
    if (!world_.is_active(target_id_) || !target->propagate_tr
        || target->ephemeris_providers.translation || target->emits_gravity) {
        result.status = StatusCode::unsupported_method;
        return result;
    }
    WorldTargetSTM stm;
    stm.target_id = target_id_;
    stm.Phi = options.Phi0;
    stm.abs_tol = options.stm_abs_tol;
    stm.rel_tol = options.stm_rel_tol;
    if (options.with_stm && (!stm.Phi.allFinite() || !stm.abs_tol.allFinite()
        || stm.abs_tol.minCoeff() <= 0.0 || !std::isfinite(stm.rel_tol) || stm.rel_tol < 0.0))
        return result;
    // pointer is borrowed only for this call, including all error exits
    struct ClearSTM {
        WorldStepperWorkspace& workspace;
        ~ClearSTM() { workspace.target_stm = nullptr; }
    } clear{workspace_};
    workspace_.target_stm = options.with_stm ? &stm : nullptr;
    result.Phi = stm.Phi;
    result.has_stm = options.with_stm;
    for (i32 step = 0; result.t != t_target; ++step) {
        if (step >= options.max_steps) {
            result.status = StatusCode::max_steps_reached;
            return result;
        }
        f64 remaining = t_target - result.t;
        f64 dt = std::copysign(std::min(std::abs(remaining), options.max_step), remaining);
        if (!std::isfinite(dt) || result.t + dt == result.t) {
            result.status = StatusCode::step_size_underflow;
            return result;
        }
        f64 before = result.t;
        WorldStepResult stepped = step_world(world_, dt, stepper_, workspace_);
        result.stats += stepped.stats;
        result.t = world_.t_sim();
        result.x = world_.body(target_id_)->x_tr;
        result.Phi = stm.Phi;
        if (stepped.status != StatusCode::ok) {
            result.status = stepped.status;
            return result;
        }
        if (result.t == before) {
            result.status = StatusCode::step_size_underflow;
            return result;
        }
    }
    result.status = StatusCode::ok;
    return result;
}
