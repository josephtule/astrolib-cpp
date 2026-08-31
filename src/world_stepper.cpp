// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/world_stepper.hpp"
#include "core/body.hpp"
#include "core/dynamics_rotational.hpp"
#include "core/entity.hpp"
#include "core/integrator_adaptive.hpp"
#include "core/integrator_fixed.hpp"
#include "core/state.hpp"
#include "core/status.hpp"
#include "core/world.hpp"
#include "util/math.hpp"
#include "util/typedefs.hpp"

#include <cmath>

// NOTE: this is a preliminary implementation

static svec<EntityId> propagated_tr_ids(const World& world) {
    svec<EntityId> ids;

    for (EntityId id : world.active_entity_ids()) {
        const Body* body = world.body(id);
        if (body == nullptr) continue;
        if (!body->propagate_tr) continue;
        if (body->body_type == BodyType::station) {
            const Station* stat = world.station(id);
            if (stat->anchored || !stat->mass_properties.active)
                continue; // only propagate free stations
        }
        ids.push_back(id);
    }

    return ids;
}

static svec<EntityId> propagated_att_ids(const World& world) {
    svec<EntityId> ids;

    for (EntityId id : world.active_entity_ids()) {
        const Body* body = world.body(id);
        if (body == nullptr) continue;
        if (!body->propagate_att) continue;
        switch (body->body_type) {
        case BodyType::celestial: continue;
        case BodyType::station: {
            const Station* stat = world.station(id);
            if (stat->anchored || !stat->mass_properties.active) continue;
        } break;
        case BodyType::satellite: {
            {
                const Satellite* sat = world.satellite(id);
                if (!sat->mass_properties.active) continue;
            }
            break;
        }
        default: continue;
        }
        ids.push_back(id);
    }

    return ids;
}

static svec<EntityId> celestial_att_ids(const World& world) {
    svec<EntityId> ids;

    for (EntityId id : world.active_entity_ids()) {
        if (world.is_celestial(id)) {
            const Celestial* cel = world.celestial(id);
            if (cel->propagate_att
                && cel->attitude_model != CelestialAttitudeModel::fixed) {
                ids.push_back(id);
            }
        }
    }

    return ids;
}

static svec<EntityId> gravity_source_ids(const World& world) {
    svec<EntityId> ids;

    for (EntityId id : world.active_entity_ids()) {
        const Body* body = world.body(id);
        if (body == nullptr) continue;
        if (!body->emits_gravity) continue;
        ids.push_back(id);
    }
    return ids;
}

static svec<EntityId> source_att_ids(const World& world) {
    svec<EntityId> ids;

    for (EntityId id : world.active_entity_ids()) {
        const Body* body = world.body(id);
        if (body == nullptr) continue;
        // NOTE: atmosphere currently always false
        if (!body->emits_gravity && !body->has_atmosphere) continue;
        const Celestial* cel = world.celestial(id);
        if (cel == nullptr) continue;
        if (!cel->propagate_att) continue;
        if (cel->attitude_model == CelestialAttitudeModel::fixed) continue;
        ids.push_back(id);
    }

    return ids;
}

static svec<EntityId> staged_att_ids(
    const svec<EntityId>& propagated_ids,
    const svec<EntityId>& celestial_ids
) {
    svec<EntityId> ids = propagated_ids;

    for (EntityId id : celestial_ids) {
        if (std::find(ids.begin(), ids.end(), id) == ids.end()) {
            ids.push_back(id);
        }
    }

    return ids;
}

DerivTr derivtr_world(const World& world, EntityId id, const StateTr& x) {
    // TODO: source states still read from world state
    // not fully staged yet for moving sources
    DerivTr dx;
    dx.dr = x.v;
    dx.dv = world.gravity_accel_on(id, x);

    return dx;
}

DerivAtt derivatt_world(const World& world, EntityId id, const StateAtt& x) {
    DerivAtt dx;

    const Body* body = world.body(id);
    switch (body->body_type) {
    case BodyType::celestial: break;
    case BodyType::satellite: {
        const Satellite* sat = world.satellite(id);
        if (sat == nullptr || !sat->mass_properties.active) break;
        if (sat->mass_properties.principal_axes) {
            dx = d_rigidbody(world.t_sim(), x, sat->mass_properties.I);
        } else {
            dx = d_rigidbody(
                world.t_sim(),
                x,
                sat->mass_properties.I,
                sat->mass_properties.I_inv
            );
        }
    } break;
    case BodyType::station: {
        const Station* stat = world.station(id);
        if (stat == nullptr || stat->anchored || !stat->mass_properties.active) break;
        if (stat->mass_properties.principal_axes) {
            dx = d_rigidbody(world.t_sim(), x, stat->mass_properties.I);
        } else {
            dx = d_rigidbody(
                world.t_sim(),
                x,
                stat->mass_properties.I,
                stat->mass_properties.I_inv
            );
        }
    } break;

    case BodyType::unknown: break;
    }

    return dx;
}

bool step_tr_world(World& world, EntityId id, f64 dt, const WorldStepperConfig& cfg) {
    Body* body = world.body(id);
    if (body == nullptr) return false;

    auto f = [&](f64 t, StateTr x) -> DerivTr { return derivtr_world(world, id, x); };
    auto tx = step_integrator<StateTr, DerivTr>(
        f,
        world.t_sim(),
        body->x_tr,
        dt,
        cfg.integrator_tr
    );
    body->x_tr = tx.second;

    return true;
}

bool step_att_world(World& world, EntityId id, f64 dt, const WorldStepperConfig& cfg) {
    Body* body = world.body(id);
    if (body == nullptr) return false;

    switch (body->body_type) {
    case BodyType::celestial: return false;
    case BodyType::station: {
        Station* stat = world.station(id);
        if (stat->anchored || !stat->mass_properties.active) return false;
    } break;
    case BodyType::satellite: {
        Satellite* sat = world.satellite(id);
        if (!sat->mass_properties.active) return false;
    } break;
    default: return false;
    }

    auto f = [&](f64 t, StateAtt x) -> DerivAtt { return derivatt_world(world, id, x); };
    auto tx = step_integrator<StateAtt, DerivAtt>(
        f,
        world.t_sim(),
        body->x_att,
        dt,
        cfg.integrator_att
    );
    body->x_att = tx.second;
    normalize_quaternion_inplace<f64>(body->x_att.q);

    return true;
}

bool step_cel_att_world(World& world, EntityId id, f64 dt) {
    Celestial* cel = world.celestial(id);
    if (cel == nullptr) return false;

    switch (cel->attitude_model) {
    case CelestialAttitudeModel::fixed: return false;
    case CelestialAttitudeModel::simple_spin: {
        cel->x_att.q = step_q_simple_spin(cel->x_att, dt);
    } break;
    case CelestialAttitudeModel::provider: {
        return false; // TODO: add provider later
    } break;
    }

    return true;
}

using WorldTrDerivMap = umap<EntityId, DerivTr>;
using WorldAttDerivMap = umap<EntityId, DerivAtt>;

struct WorldTrStage {
    svec<EntityId> ids;
    umap<EntityId, StateTr> x;
};
struct DerivTrWeight {
    // K_i in RK integrators (translational)
    const umap<EntityId, DerivTr>* k = nullptr;
    f64 scale = 0.0;
};

struct WorldAttStage {
    svec<EntityId> ids;
    umap<EntityId, StateAtt> x;
};
struct DerivAttWeight {
    const umap<EntityId, DerivAtt>* k = nullptr;
    f64 scale = 0.0;
};

struct WorldStage {
    WorldTrStage tr;
    WorldAttStage att;
};

struct WorldStageDeriv {
    WorldTrDerivMap tr;
    WorldAttDerivMap att;
};

static StatusCode build_tr_stage(
    const World& world,
    const WorldStepperWorkspace& wksp,
    WorldTrStage& stage
) {
    WorldTrStage stage_new;

    for (EntityId id : wksp.propagated_tr_ids) {
        const Body* body = world.body(id);
        if (body == nullptr) return StatusCode::body_not_found;
        if (!finite_state_tr(body->x_tr)) return StatusCode::invalid_state;

        stage_new.ids.push_back(id);
        stage_new.x.emplace(id, body->x_tr);
    }

    stage = stage_new;
    return StatusCode::ok;
}

static WorldTrStage build_tr_stage(
    const World& world,
    const WorldStepperWorkspace& wksp
) {
    WorldTrStage stage;
    // TODO: update fixed-step RK paths to propagate stage-build statuses.
    if (build_tr_stage(world, wksp, stage) != StatusCode::ok) return {};
    return stage;
}

static StatusCode build_att_stage(
    const World& world,
    const WorldStepperWorkspace& wksp,
    WorldAttStage& stage
) {
    WorldAttStage stage_new;

    for (EntityId id : wksp.staged_att_ids) {
        const Body* body = world.body(id);
        if (body == nullptr) return StatusCode::body_not_found;
        if (!finite_state_att(body->x_att)) {
            return StatusCode::invalid_att_state;
        }

        stage_new.ids.push_back(id);
        stage_new.x.emplace(id, body->x_att);
    }

    stage = stage_new;
    return StatusCode::ok;
}

static WorldAttStage build_att_stage(
    const World& world,
    const WorldStepperWorkspace& wksp
) {
    WorldAttStage stage;
    // TODO: update fixed-step RK paths to propagate stage-build statuses.
    if (build_att_stage(world, wksp, stage) != StatusCode::ok) return {};
    return stage;
}

struct WorldStageBuildResult {
    StatusCode status = StatusCode::invalid_state;
    WorldStage stage;
};

struct WorldAdaptiveTrialResult {
    StatusCode status = StatusCode::invalid_state;

    // candidates, separate for quaternion error
    WorldStage stage_high;
    WorldStage stage_low;

    i32 deriv_evals = 0;
};

struct WorldFixedTrialResult {
    StatusCode status = StatusCode::invalid_state;
    WorldStage stage;
    i32 deriv_evals = 0;
};

static WorldStageBuildResult build_world_stage(
    const World& world,
    const WorldStepperWorkspace& wksp
) {
    WorldStageBuildResult result;

    result.status = build_tr_stage(world, wksp, result.stage.tr);
    if (result.status != StatusCode::ok) return result;

    result.status = build_att_stage(world, wksp, result.stage.att);
    if (result.status != StatusCode::ok) return result;

    result.status = StatusCode::ok;
    return result;
}

static StatusCode celestial_att_at_stage(
    const World& world,
    EntityId id,
    const StateAtt& x_att_base,
    f64 t_base,
    f64 t_stage,
    StateAtt& x_att_stage
) {
    const Celestial* cel = world.celestial(id);
    if (cel == nullptr) return StatusCode::body_not_found;
    x_att_stage = x_att_base;

    switch (cel->attitude_model) {
    case CelestialAttitudeModel::fixed: break;
    case CelestialAttitudeModel::simple_spin: {
        x_att_stage.q = step_q_simple_spin(x_att_base, t_stage - t_base);
        normalize_quaternion_inplace<f64>(x_att_stage.q);
    } break;
    case CelestialAttitudeModel::provider: {
        // TODO: do later
        x_att_stage = x_att_base;
        return StatusCode::unsupported_method;
    } break;
    }

    if (!finite_state_att(x_att_stage)) return StatusCode::invalid_att_state;
    return StatusCode::ok;
}

template <size_t Stages>
static StatusCode build_att_tableau_stage(
    const World& world,
    const WorldAttStage& base_stage,
    const array<WorldStageDeriv, Stages>& k,
    size_t stage_index,
    f64 t,
    f64 dt,
    const RKTableau<Stages>& tableau,
    const WorldStepperWorkspace& wksp,
    WorldAttStage& stage
) {
    WorldAttStage trial = base_stage;

    trial.ids = base_stage.ids;

    if (stage_index >= Stages) return StatusCode::invalid_input;

    for (EntityId id : wksp.propagated_att_ids) {
        StateAtt x_stage = base_stage.x.at(id);

        for (i32 j = 0; j < stage_index; ++j) {
            x_stage += dt * tableau.a[stage_index][j] * k[j].att.at(id);
            if (!finite_state_att(x_stage)) return StatusCode::invalid_state;
        }
        trial.x.at(id) = x_stage;
    }

    for (EntityId id : wksp.celestial_att_ids) {
        StateAtt x_stage = base_stage.x.at(id);

        f64 t_stage = t + tableau.c[stage_index] * dt;
        StatusCode status
            = celestial_att_at_stage(world, id, x_stage, t, t_stage, x_stage);
        if (!finite_state_att(x_stage)) return StatusCode::invalid_state;
        if (status != StatusCode::ok) return status;
        trial.x.at(id) = x_stage;
    }

    stage = trial;
    return StatusCode::ok;
}

static WorldAttStage build_source_att_stage(
    const World& world,
    const WorldStepperWorkspace& wksp,
    bool stage_source_att
) {
    WorldAttStage stage;
    if (!stage_source_att) return stage;

    svec<EntityId> ids = wksp.source_att_ids;
    for (EntityId id : ids) {
        const Celestial* cel = world.celestial(id);
        if (cel == nullptr || !cel->emits_gravity) {
            // || cel->gravity_model == GravityModel::pointmass){
            continue;
        }
        // TODO: add has atmosphere and radiation
        // only bodies that affect force via their attitude used

        stage.ids.push_back(id);
        stage.x.emplace(id, cel->x_att);
    }

    return stage;
}

// Gravity sources use their staged translational state when propagated
// fixed sources fall back to their stored world state.
static StateTr source_tr_from_stage_or_world(
    const World& world,
    const WorldTrStage& stage,
    EntityId source_id
) {
    // NOTE: currently only celestials are valid gravity emitters.
    const Celestial* cel = world.celestial(source_id);
    if (cel == nullptr) return StateTr{};

    auto it = stage.x.find(source_id);
    if (it == stage.x.end()) return cel->x_tr; // world fallback

    return it->second;
}

static StatusCode source_tr_from_stage_or_world(
    const World& world,
    const WorldTrStage& stage,
    EntityId source_id,
    StateTr& x_tr_source
) {
    const Celestial* cel = world.celestial(source_id);
    if (cel == nullptr) return StatusCode::body_not_found;

    auto it = stage.x.find(source_id);
    if (it == stage.x.end())
        x_tr_source = cel->x_tr; // return body not found instead?
    else
        x_tr_source = it->second;

    if (!finite_state_tr(x_tr_source)) return StatusCode::invalid_state;

    return StatusCode::ok;
}

static StateAtt source_att_from_stage_or_world(
    const World& world,
    const WorldAttStage& stage,
    EntityId source_id
) {
    const Body* body = world.body(source_id);
    if (body == nullptr) return StateAtt{};

    auto it = stage.x.find(source_id);
    if (it == stage.x.end()) return body->x_att;

    return it->second;
}

static StatusCode source_att_from_stage_or_world(
    const World& world,
    const WorldAttStage& stage,
    EntityId source_id,
    StateAtt& x_att_source
) {
    const Celestial* cel = world.celestial(source_id);
    if (cel == nullptr) return StatusCode::body_not_found;

    auto it = stage.x.find(source_id);
    if (it == stage.x.end())
        x_att_source = cel->x_att;
    else
        x_att_source = it->second;

    if (!finite_state_att(x_att_source)) return StatusCode::invalid_att_state;

    return StatusCode::ok;
}

static WorldTrStage build_tr_trial_stage(
    const WorldTrStage& base_stage,
    std::initializer_list<DerivTrWeight> weights
) {
    WorldTrStage trial;
    trial.ids = base_stage.ids;

    for (EntityId id : base_stage.ids) {
        StateTr x = base_stage.x.at(id);

        for (const DerivTrWeight& weight : weights) {
            if (weight.k == nullptr) continue;
            x += weight.scale * weight.k->at(id);
        }
        trial.x.emplace(id, x);
    }

    return trial;
}

template <size_t Stages>
static StatusCode build_tr_tableau_stage(
    const WorldTrStage& base_stage,
    const array<WorldStageDeriv, Stages>& k,
    size_t stage_index,
    f64 dt,
    const RKTableau<Stages>& tableau,
    WorldTrStage& stage
) {
    WorldTrStage trial = base_stage;
    trial.ids = base_stage.ids;

    if (stage_index >= Stages) return StatusCode::invalid_input;

    for (EntityId id : base_stage.ids) {
        StateTr x_stage = base_stage.x.at(id);

        for (i32 j = 0; j < stage_index; ++j) {
            x_stage += dt * tableau.a[stage_index][j] * k[j].tr.at(id);
            if (!finite_state_tr(x_stage)) return StatusCode::invalid_state;
        }
        trial.x.at(id) = x_stage;
    }

    stage = trial;

    return StatusCode::ok;
}

static WorldAttStage build_source_att_trial_stage(
    const World& world,
    const WorldAttStage& base_stage,
    const f64 dt_stage
) {
    WorldAttStage trial;
    trial.ids = base_stage.ids;

    for (EntityId id : base_stage.ids) {
        const Celestial* cel = world.celestial(id);
        if (cel == nullptr) continue;

        StateAtt x = base_stage.x.at(id);

        switch (cel->attitude_model) {
        case CelestialAttitudeModel::fixed: break; // no change
        case CelestialAttitudeModel::simple_spin: {
            x.q = step_q_simple_spin(x, dt_stage);
        } break;
        case CelestialAttitudeModel::provider: {
            // TODO: add provider here attitude here
        } break;
        }

        trial.x.emplace(id, x);
    }

    return trial;
}

static vec3d staged_gravity_accel_on(
    const World& world,
    const WorldTrStage& stage_tr,
    const WorldAttStage& stage_att,
    EntityId target_id,
    const StateTr& x_tr_target,
    const WorldStepperWorkspace& wksp
) {
    vec3d a = vec3d0;
    const Body* target = world.body(target_id);
    if (target == nullptr) return a;

    const svec<EntityId>& source_ids = wksp.gravity_source_ids;
    for (EntityId source_id : source_ids) {
        if (target_id == source_id) continue;
        const Body* body = world.body(source_id);
        if (body == nullptr || !body->emits_gravity) continue;
        const Celestial* source = world.celestial(source_id);
        if (source == nullptr) continue;
        StateTr x_tr_source = source_tr_from_stage_or_world(world, stage_tr, source_id);
        StateAtt x_att_source
            = source_att_from_stage_or_world(world, stage_att, source_id);
        a += world.gravity_accel_from(
            target_id,
            x_tr_target,
            source_id,
            x_tr_source,
            x_att_source
        );
    }

    return a;
}

static StatusCode staged_gravity_accel_on(
    const World& world,
    const WorldTrStage& stage_tr,
    const WorldAttStage& stage_att,
    EntityId target_id,
    const StateTr& x_tr_target,
    const WorldStepperWorkspace& wksp,
    vec3d& a_grav
) {
    StatusCode status = StatusCode::invalid_input;
    vec3d a_temp = vec3d0;
    const Body* target = world.body(target_id);
    if (target == nullptr) return StatusCode::body_not_found;

    const svec<EntityId>& source_ids = wksp.gravity_source_ids;
    for (EntityId source_id : source_ids) {
        if (target_id == source_id) continue;
        const Celestial* source = world.celestial(source_id);
        if (source == nullptr) return StatusCode::body_not_found;

        StateTr x_tr_source;
        status = source_tr_from_stage_or_world(world, stage_tr, source_id, x_tr_source);
        if (status != StatusCode::ok) return status;

        StateAtt x_att_source;
        status
            = source_att_from_stage_or_world(world, stage_att, source_id, x_att_source);
        if (status != StatusCode::ok) return status;

        a_temp += world.gravity_accel_from(
            target_id,
            x_tr_target,
            source_id,
            x_tr_source,
            x_att_source
        );
    }
    if (!finite_vec(a_temp)) return StatusCode::non_finite_result;

    a_grav = a_temp;
    return StatusCode::ok;
}

static DerivTr derivtr_world_staged(
    const World& world,
    const WorldTrStage& stage_tr,
    const WorldAttStage& stage_att,
    EntityId target_id,
    const StateTr& x_tr_target,
    const WorldStepperWorkspace& wksp
) {
    // TODO: add other forces later
    DerivTr dx;
    dx.dr = x_tr_target.v;
    dx.dv = staged_gravity_accel_on(
        world,
        stage_tr,
        stage_att,
        target_id,
        x_tr_target,
        wksp
    );
    return dx;
}

static StatusCode derivtr_world_staged(
    const World& world,
    f64 t_stage,
    const WorldTrStage& stage_tr,
    const WorldAttStage& stage_att,
    EntityId target_id,
    const StateTr& x_tr_target,
    const WorldStepperWorkspace& wksp,
    DerivTr& dx
) {
    DerivTr dx_temp;
    dx_temp.dr = x_tr_target.v;
    StatusCode status = staged_gravity_accel_on(
        world,
        stage_tr,
        stage_att,
        target_id,
        x_tr_target,
        wksp,
        dx_temp.dv
    );
    if (status != StatusCode::ok) return status;
    if (!finite_deriv_tr(dx_temp)) return StatusCode::invalid_state;

    dx = dx_temp;
    return StatusCode::ok;
}

static DerivTr derivtr_world_staged(
    const World& world,
    f64 t_stage,
    const WorldTrStage& stage_tr,
    const WorldAttStage& stage_att,
    EntityId target_id,
    const StateTr& x_tr_target,
    const WorldStepperWorkspace& wksp
) {
    // TODO: providers will use t_stage, but not yet implemented
    return derivtr_world_staged(world, stage_tr, stage_att, target_id, x_tr_target, wksp);
}

static StatusCode derivatt_world_staged(
    const World& world,
    f64 t_stage,
    const WorldTrStage& stage_tr,
    const WorldAttStage& stage_att,
    EntityId target_id,
    const StateAtt& x_att_target,
    const WorldStepperWorkspace& wksp,
    DerivAtt& dx
) {
    DerivAtt dx_temp;
    // TODO: add other torques later
    const Body* body = world.body(target_id);
    if (body == nullptr) return StatusCode::body_not_found;

    if (body->body_type == BodyType::satellite) {
        const Satellite* sat = world.satellite(target_id);
        if (sat == nullptr) return StatusCode::body_not_found;
        if (sat->mass_properties.principal_axes) {
            dx_temp = d_rigidbody(t_stage, x_att_target, sat->mass_properties.I);
        } else {
            dx_temp = d_rigidbody(
                t_stage,
                x_att_target,
                sat->mass_properties.I,
                sat->mass_properties.I_inv
            );
        }
    }

    if (body->body_type == BodyType::station) {
        const Station* stat = world.station(target_id);
        if (stat == nullptr) return StatusCode::body_not_found;
        if (!stat->anchored) {
            if (stat->mass_properties.principal_axes) {
                dx_temp = d_rigidbody(t_stage, x_att_target, stat->mass_properties.I);
            } else {
                dx_temp = d_rigidbody(
                    t_stage,
                    x_att_target,
                    stat->mass_properties.I,
                    stat->mass_properties.I_inv
                );
            }
        }
    }

    if (!finite_deriv_att(dx_temp)) return StatusCode::invalid_att_state;

    dx = dx_temp;
    return StatusCode::ok;
}

static StatusCode evaluate_world_stage_derivatives(
    const World& world,
    f64 t_stage,
    const WorldStage& stage,
    const WorldStepperWorkspace& wksp,
    WorldStageDeriv& dx
) {
    StatusCode status;
    WorldStageDeriv dx_temp{};

    for (EntityId id : wksp.propagated_att_ids) {
        auto x_att_it = stage.att.x.find(id);
        if (x_att_it == stage.att.x.end()) return StatusCode::invalid_att_state;
        StateAtt x_att_target = x_att_it->second;

        DerivAtt dx_att;
        status = derivatt_world_staged(
            world,
            t_stage,
            stage.tr,
            stage.att,
            id,
            x_att_target,
            wksp,
            dx_att
        );
        if (status != StatusCode::ok) return status;
        dx_temp.att.emplace(id, dx_att);
    }

    for (EntityId id : wksp.propagated_tr_ids) {
        auto x_tr_it = stage.tr.x.find(id);
        if (x_tr_it == stage.tr.x.end()) return StatusCode::invalid_state;
        StateTr x_tr_target = x_tr_it->second;

        DerivTr dx_tr;
        status = derivtr_world_staged(
            world,
            t_stage,
            stage.tr,
            stage.att,
            id,
            x_tr_target,
            wksp,
            dx_tr
        );
        if (status != StatusCode::ok) return status;
        dx_temp.tr.emplace(id, dx_tr);
    }

    dx = std::move(dx_temp);
    return StatusCode::ok;
}

template <size_t Stages>
static StatusCode build_world_tableau_stage(
    const World& world,
    const WorldStage& base_stage,
    const array<WorldStageDeriv, Stages>& k,
    size_t stage_index,
    f64 t,
    f64 dt,
    const RKTableau<Stages>& tableau,
    const WorldStepperWorkspace& wksp,
    WorldStage& stage
) {
    StatusCode status = StatusCode::invalid_input;

    status = build_att_tableau_stage(
        world,
        base_stage.att,
        k,
        stage_index,
        t,
        dt,
        tableau,
        wksp,
        stage.att
    );
    if (status != StatusCode::ok) return status;

    status = build_tr_tableau_stage(base_stage.tr, k, stage_index, dt, tableau, stage.tr);
    if (status != StatusCode::ok) return status;

    return StatusCode::ok;
}

template <size_t Stages>
static WorldFixedTrialResult step_world_fixed_rk_trial(
    const World& world,
    f64 t,
    f64 dt,
    const RKTableau<Stages>& tableau,
    const WorldStepperWorkspace& wksp
) {
    WorldFixedTrialResult result;

    if (!isfinite(t) || !isfinite(dt) || dt == 0.0 || tableau.embedded) {
        result.status = StatusCode::invalid_input;
        return result;
    }

    WorldStageBuildResult base_result = build_world_stage(world, wksp);
    if (base_result.status != StatusCode::ok) {
        result.status = base_result.status;
        return result;
    }
    const WorldStage& base_stage = base_result.stage;

    array<WorldStageDeriv, Stages> k{};
    for (size_t i = 0; i < Stages; ++i) {
        WorldStage stage_i;

        StatusCode status = build_world_tableau_stage(
            world,
            base_stage,
            k,
            i,
            t,
            dt,
            tableau,
            wksp,
            stage_i
        );
        if (status != StatusCode::ok) {
            result.status = status;
            return result;
        }

        f64 t_stage = t + tableau.c[i] * dt;
        status = evaluate_world_stage_derivatives(world, t_stage, stage_i, wksp, k[i]);
        ++result.deriv_evals;

        if (status != StatusCode::ok) {
            result.status = status;
            return result;
        }
    }

    WorldStage stage_next = base_stage;

    for (EntityId id : wksp.propagated_tr_ids) {
        auto x_base_it = base_stage.tr.x.find(id);
        if (x_base_it == base_stage.tr.x.end()) {
            result.status = StatusCode::invalid_state;
            return result;
        }
        StateTr x_next = x_base_it->second;

        for (size_t i = 0; i < Stages; ++i) {
            auto dx_it = k[i].tr.find(id);
            if (dx_it == k[i].tr.end()) {
                result.status = StatusCode::invalid_state;
                return result;
            }
            x_next += dt * tableau.b_high[i] * dx_it->second;
        }
        if (!finite_state(x_next)) {
            result.status = StatusCode::invalid_state;
            return result;
        }

        stage_next.tr.x.at(id) = x_next;
    }

    for (EntityId id : wksp.propagated_att_ids) {
        auto x_base_it = base_stage.att.x.find(id);
        if (x_base_it == base_stage.att.x.end()) {
            result.status = StatusCode::invalid_att_state;
            return result;
        }
        StateAtt x_next = x_base_it->second;

        for (size_t i = 0; i < Stages; ++i) {
            auto dx_it = k[i].att.find(id);
            if (dx_it == k[i].att.end()) {
                result.status = StatusCode::invalid_att_state;
                return result;
            }
            x_next += dt * tableau.b_high[i] * dx_it->second;
        }
        if (!finite_state(x_next)) {
            result.status = StatusCode::invalid_att_state;
            return result;
        }
        normalize_quaternion_inplace<f64>(x_next.q);

        stage_next.att.x.at(id) = x_next;
    }

    for (EntityId id : wksp.celestial_att_ids) {
        auto x_base_it = base_stage.att.x.find(id);
        if (x_base_it == base_stage.att.x.end()) {
            result.status = StatusCode::invalid_att_state;
            return result;
        }

        StateAtt x_end;
        result.status
            = celestial_att_at_stage(world, id, x_base_it->second, t, t + dt, x_end);
        if (result.status != StatusCode::ok) return result;
        if (!finite_state(x_end)) {
            result.status = StatusCode::invalid_att_state;
            return result;
        }
        normalize_quaternion_inplace<f64>(x_end.q);

        stage_next.att.x.at(id) = x_end;
    }

    result.stage = std::move(stage_next);
    result.status = StatusCode::ok;
    return result;
}

template <size_t Stages>
static WorldAdaptiveTrialResult step_world_embedded_rk_trial(
    const World& world,
    f64 t,
    f64 dt,
    const RKTableau<Stages>& tableau,
    const WorldStepperWorkspace& wksp
) {
    WorldAdaptiveTrialResult result;

    if (!isfinite(t) || !isfinite(dt) || dt == 0.0) {
        result.status = StatusCode::invalid_input;
        return result;
    }

    // base stage
    WorldStageBuildResult base_result = build_world_stage(world, wksp);
    if (base_result.status != StatusCode::ok) {
        result.status = base_result.status;
        return result;
    }
    const WorldStage& base_stage = base_result.stage;

    array<WorldStageDeriv, Stages> k{};
    for (size_t i = 0; i < Stages; ++i) {
        WorldStage stage_i;

        StatusCode status = build_world_tableau_stage(
            world,
            base_stage,
            k,
            i,
            t,
            dt,
            tableau,
            wksp,
            stage_i
        );
        if (status != StatusCode::ok) {
            result.status = status;
            return result;
        }

        f64 t_stage = t + tableau.c[i] * dt;

        status = evaluate_world_stage_derivatives(world, t_stage, stage_i, wksp, k[i]);
        ++result.deriv_evals;

        if (status != StatusCode::ok) {
            result.status = status;
            return result;
        }
    }

    WorldStage stage_high = base_stage;
    WorldStage stage_low = base_stage;

    for (EntityId id : wksp.propagated_tr_ids) {
        auto x_base_it = base_stage.tr.x.find(id);
        if (x_base_it == base_stage.tr.x.end()) {
            result.status = StatusCode::invalid_state;
            return result;
        }
        StateTr x_base = x_base_it->second;
        StateTr x_high = x_base;
        StateTr x_low = x_base;

        for (i32 i = 0; i < Stages; ++i) {
            auto dx_it = k[i].tr.find(id);
            if (dx_it == k[i].tr.end()) {
                result.status = StatusCode::invalid_state;
                return result;
            }
            DerivTr dx = dx_it->second;

            x_high += dt * tableau.b_high[i] * dx;
            x_low += dt * tableau.b_low[i] * dx;
        }
        if (!finite_state(x_high) || !finite_state(x_low)) {
            result.status = StatusCode::invalid_state;
            return result;
        }

        stage_high.tr.x.at(id) = x_high;
        stage_low.tr.x.at(id) = x_low;
    }

    for (EntityId id : wksp.propagated_att_ids) {
        auto x_base_it = base_stage.att.x.find(id);
        if (x_base_it == base_stage.att.x.end()) {
            result.status = StatusCode::invalid_state;
            return result;
        }
        StateAtt x_base = x_base_it->second;
        StateAtt x_high = x_base;
        StateAtt x_low = x_base;

        for (i32 i = 0; i < Stages; ++i) {
            auto dx_it = k[i].att.find(id);
            if (dx_it == k[i].att.end()) {
                result.status = StatusCode::invalid_state;
                return result;
            }
            DerivAtt dx = dx_it->second;

            x_high += dt * tableau.b_high[i] * dx;
            x_low += dt * tableau.b_low[i] * dx;
        }

        if (!finite_state(x_high) || !finite_state(x_low)) {
            result.status = StatusCode::invalid_state;
            return result;
        }
        normalize_quaternion_inplace<f64>(x_high.q);
        normalize_quaternion_inplace<f64>(x_low.q);

        stage_high.att.x.at(id) = x_high;
        stage_low.att.x.at(id) = x_low;
    }

    for (EntityId id : wksp.celestial_att_ids) {
        auto x_base_it = base_stage.att.x.find(id);
        if (x_base_it == base_stage.att.x.end()) {
            result.status = StatusCode::invalid_state;
            return result;
        }
        StateAtt x_base = x_base_it->second;

        StateAtt x_end;
        result.status = celestial_att_at_stage(world, id, x_base, t, t + dt, x_end);
        if (result.status != StatusCode::ok) return result;

        if (!finite_state(x_end)) {
            result.status = StatusCode::invalid_state;
            return result;
        }
        normalize_quaternion_inplace<f64>(x_end.q);

        stage_high.att.x.at(id) = x_end;
        stage_low.att.x.at(id) = x_end;
    }

    result.stage_high = std::move(stage_high);
    result.stage_low = std::move(stage_low);
    result.status = StatusCode::ok;
    return result;
}

static StatusCode validate_adaptive_world_integrator_config(
    WorldAdaptiveIntegratorConfig cfg,
    WorldStepperConfig stepper
) {
    StatusCode status = validate_adaptive_integrator_config(cfg.opts);
    if (status != StatusCode::ok) return status;

    // TEMP: temporarily only allow dopri54
    if (cfg.integrator_tr != IntegratorTypeAdaptive::dopri54
        || cfg.integrator_att != IntegratorTypeAdaptive::dopri54) {
        return StatusCode::unsupported_method;
    }
    // TEMP:
    if (cfg.integrator_tr != cfg.integrator_att) return StatusCode::unsupported_method;

    return StatusCode::ok;
}

static f64 adaptive_error_norm_world_tr(
    const StateTr& x_base,
    const StateTr& x_high,
    const StateTr& x_low,
    const WorldAdaptiveIntegratorConfig& cfg
) {
    if (!finite_state(x_base) || !finite_state(x_high) || !finite_state(x_low)) {
        return inf<f64>;
    }

    vec6d base_vec = statetr_to_vec6d(x_base);
    vec6d high_vec = statetr_to_vec6d(x_high);
    vec6d low_vec = statetr_to_vec6d(x_low);

    vec6d error = high_vec - low_vec;

    vec6d abs_tol;
    abs_tol << cfg.opts.abs_tol_r, cfg.opts.abs_tol_r, cfg.opts.abs_tol_r,
        cfg.opts.abs_tol_v, cfg.opts.abs_tol_v, cfg.opts.abs_tol_v;

    vec6d mag = base_vec.cwiseAbs().cwiseMax(high_vec.cwiseAbs());
    vec6d scale = abs_tol + cfg.opts.rel_tol * mag;
    if (!finite_pos(scale)) return inf<f64>;
    vec6d weighted_error = error.cwiseQuotient(scale);

    f64 error_norm = std::sqrt(weighted_error.squaredNorm() / 6.0);
    if (!finite_nonneg(error_norm)) return inf<f64>;

    return error_norm;
}

// TODO: move this later
template <typename T>
static T quaternion_angle_between(const vec4<T>& q1, const vec4<T>& q2) {
    if (!finite_vec(q1) || !finite_vec(q2)) return inf<f64>;
    T q1_norm = q1.norm();
    T q2_norm = q2.norm();
    if (!finite_pos<T>(q1_norm) || !finite_pos<T>(q2_norm)) return inf<f64>;

    T cos_half_angle = std::abs(q1.dot(q2)) / (q1_norm * q2_norm);
    cos_half_angle = std::clamp(cos_half_angle, 0.0, 1.0);
    T angle = 2.0 * std::acos(cos_half_angle);

    return angle;
}

static f64 adaptive_error_norm_world_att(
    const StateAtt& x_base,
    const StateAtt& x_high,
    const StateAtt& x_low,
    const WorldAdaptiveIntegratorConfig& cfg
) {
    if (!finite_state(x_base) || !finite_state(x_high) || !finite_state(x_low)) {
        return inf<f64>;
    }

    f64 angle_error = quaternion_angle_between(x_high.q, x_low.q);
    f64 angle_high = quaternion_angle_between(x_base.q, x_high.q);
    f64 angle_low = quaternion_angle_between(x_base.q, x_low.q);

    if (!isfinite(angle_error) || !isfinite(angle_high) || !isfinite(angle_low)) {
        return inf<f64>;
    }

    f64 angle_scale
        = cfg.opts.abs_tol_angle + cfg.opts.rel_tol * std::max(angle_high, angle_low);
    f64 weighted_angle = angle_error / angle_scale;

    vec3d w_error = x_high.w - x_low.w;
    vec3d w_mag = x_base.w.cwiseAbs().cwiseMax(x_high.w.cwiseAbs());
    vec3d w_scale = cfg.opts.abs_tol_w * vec3d1 + cfg.opts.rel_tol * w_mag;

    if (!finite_pos(w_scale)) return inf<f64>;
    vec3d weighted_w_error = w_error.cwiseQuotient(w_scale);

    f64 error_norm = std::sqrt(
        (weighted_angle * weighted_angle + weighted_w_error.squaredNorm()) / 4.0
    );
    if (!isfinite(error_norm)) return inf<f64>;

    return error_norm;
}

static f64 adaptive_error_norm_world(
    const World& world,
    const WorldAdaptiveTrialResult& trial,
    const WorldStepperConfig& stepper_cfg,
    const WorldAdaptiveIntegratorConfig& adaptive_cfg,
    const WorldStepperWorkspace& wksp
) {
    if (trial.status != StatusCode::ok) return inf<f64>;

    f64 world_norm = 0.0;
    bool evaluated_block = false;

    if (stepper_cfg.step_tr) {
        for (EntityId id : wksp.propagated_tr_ids) {
            const Body* body = world.body(id);
            if (body == nullptr) return inf<f64>;

            auto x_high_it = trial.stage_high.tr.x.find(id);
            if (x_high_it == trial.stage_high.tr.x.end()) return inf<f64>;
            auto x_low_it = trial.stage_low.tr.x.find(id);
            if (x_low_it == trial.stage_low.tr.x.end()) return inf<f64>;

            f64 block_norm = adaptive_error_norm_world_tr(
                body->x_tr,
                x_high_it->second,
                x_low_it->second,
                adaptive_cfg
            );

            if (!std::isfinite(block_norm)) return inf<f64>;
            world_norm = std::max(world_norm, block_norm);
            evaluated_block = true;
        }
    }

    if (stepper_cfg.step_att) {
        for (EntityId id : wksp.propagated_att_ids) {
            const Body* body = world.body(id);
            if (body == nullptr) return inf<f64>;

            auto x_high_it = trial.stage_high.att.x.find(id);
            if (x_high_it == trial.stage_high.att.x.end()) return inf<f64>;
            auto x_low_it = trial.stage_low.att.x.find(id);
            if (x_low_it == trial.stage_low.att.x.end()) return inf<f64>;

            f64 block_norm = adaptive_error_norm_world_att(
                body->x_att,
                x_high_it->second,
                x_low_it->second,
                adaptive_cfg
            );

            if (!std::isfinite(block_norm)) return inf<f64>;
            world_norm = std::max(world_norm, block_norm);
            evaluated_block = true;
        }
    }

    if (!evaluated_block) return 0.0;
    if (!std::isfinite(world_norm)) return inf<f64>;

    return world_norm;
}

static StatusCode commit_world_stage(
    World& world,
    const WorldStage& accepted_stage,
    const WorldStepperConfig& stepper_cfg,
    const WorldStepperWorkspace& wksp
) {
    // validate
    if (stepper_cfg.step_tr) {
        for (EntityId id : wksp.propagated_tr_ids) {
            const Body* body = world.body(id);
            if (body == nullptr) return StatusCode::body_not_found;

            auto x_it = accepted_stage.tr.x.find(id);
            if (x_it == accepted_stage.tr.x.end()) return StatusCode::invalid_state;

            if (!finite_state(x_it->second)) return StatusCode::invalid_state;
        }
    }
    if (stepper_cfg.step_att) {
        for (EntityId id : wksp.staged_att_ids) {
            const Body* body = world.body(id);
            if (body == nullptr) return StatusCode::body_not_found;

            auto x_it = accepted_stage.att.x.find(id);
            if (x_it == accepted_stage.att.x.end()) return StatusCode::invalid_att_state;

            if (!finite_state(x_it->second)) return StatusCode::invalid_att_state;
        }
    }

    // commit
    if (stepper_cfg.step_tr) {
        for (EntityId id : wksp.propagated_tr_ids) {
            world.body(id)->x_tr = accepted_stage.tr.x.at(id);
        }
    }
    if (stepper_cfg.step_att) {
        for (EntityId id : wksp.staged_att_ids) {
            world.body(id)->x_att = accepted_stage.att.x.at(id);
        }
    }

    return StatusCode::ok;
}

template <size_t Stages>
static StatusCode step_world_fixed_tableau(
    World& world,
    f64 t,
    f64 dt,
    const RKTableau<Stages>& tableau,
    const WorldStepperConfig& cfg,
    const WorldStepperWorkspace& wksp
) {
    WorldFixedTrialResult trial
        = step_world_fixed_rk_trial(world, t, dt, tableau, wksp);
    if (trial.status != StatusCode::ok) return trial.status;
    return commit_world_stage(world, trial.stage, cfg, wksp);
}

static StatusCode dispatch_world_fixed_tableau(
    World& world,
    f64 t,
    f64 dt,
    IntegratorTypeFixed integrator,
    const WorldStepperConfig& cfg,
    const WorldStepperWorkspace& wksp
) {
    switch (integrator) {
    case IntegratorTypeFixed::rk1:
        return step_world_fixed_tableau(world, t, dt, rk1_tableau, cfg, wksp);
    case IntegratorTypeFixed::rk2:
        return step_world_fixed_tableau(world, t, dt, rk2_tableau, cfg, wksp);
    case IntegratorTypeFixed::rk2_heun:
        return step_world_fixed_tableau(
            world,
            t,
            dt,
            rk2_heun_tableau,
            cfg,
            wksp
        );
    case IntegratorTypeFixed::rk2_ralston:
        return step_world_fixed_tableau(
            world,
            t,
            dt,
            rk2_ralston_tableau,
            cfg,
            wksp
        );
    case IntegratorTypeFixed::rk3:
        return step_world_fixed_tableau(world, t, dt, rk3_tableau, cfg, wksp);
    case IntegratorTypeFixed::rk4:
        return step_world_fixed_tableau(world, t, dt, rk4_tableau, cfg, wksp);
    }

    return StatusCode::unsupported_method;
}

template <size_t Stages>
static WorldAdaptiveStepResult propagate_world_embedded_rk(
    World& world,
    f64 tf,
    const RKTableau<Stages>& tableau,
    const WorldStepperConfig& stepper_cfg,
    const WorldAdaptiveIntegratorConfig& adaptive_cfg,
    const WorldStepperWorkspace& wksp
) {
    WorldAdaptiveStepResult result;
    result.t = world.t_sim();
    if (!isfinite(result.t) || !isfinite(tf)) {
        result.status = StatusCode::invalid_input;
        return result;
    }
    if (result.t == tf) {
        result.status = StatusCode::ok;
        return result;
    }

    result.status = validate_adaptive_world_integrator_config(adaptive_cfg, stepper_cfg);
    if (result.status != StatusCode::ok) return result;

    if (!tableau.embedded) {
        result.status = StatusCode::invalid_input;
        return result;
    }

    f64 dir = tf > world.t_sim() ? 1.0 : -1.0;
    f64 dt_trial = dir
                   * std::clamp(
                       adaptive_cfg.opts.dt_initial,
                       adaptive_cfg.opts.dt_min,
                       adaptive_cfg.opts.dt_max
                   );

    while (true) {
        if (result.stats.attempted_steps >= adaptive_cfg.opts.max_attempts) {
            result.status = StatusCode::max_steps_reached;
            return result;
        }

        f64 t = result.t;
        f64 dt_remaining = tf - t;
        bool final_attempt = std::abs(dt_trial) >= std::abs(dt_remaining);

        if (final_attempt) {
            dt_trial = dt_remaining;
        }

        if (!isfinite(dt_trial)) {
            result.status = StatusCode::invalid_input;
            return result;
        }
        WorldAdaptiveTrialResult trial
            = step_world_embedded_rk_trial(world, t, dt_trial, tableau, wksp);
        ++result.stats.attempted_steps;
        result.stats.deriv_evals += trial.deriv_evals;
        if (trial.status != StatusCode::ok) {
            result.status = trial.status;
            return result;
        }

        f64 error_norm
            = adaptive_error_norm_world(world, trial, stepper_cfg, adaptive_cfg, wksp);

        if (!isfinite(error_norm)) {
            result.status = StatusCode::non_finite_result;
            return result;
        }

        bool accepted = error_norm <= 1.0;
        f64 scale = adaptive_step_scale(error_norm, adaptive_cfg.opts);
        if (!isfinite(scale)) {
            result.status = StatusCode::non_finite_result;
            return result;
        }

        if (!accepted) {
            ++result.stats.rejected_steps;
            if (result.stats.rejected_steps >= adaptive_cfg.opts.max_rejections) {
                result.status = StatusCode::max_rejections_reached;
                return result;
            }
            scale = std::min(scale, 1.0);
            dt_trial *= scale;

            if (std::abs(dt_trial) < adaptive_cfg.opts.dt_min) {
                result.status = StatusCode::step_size_underflow;
                return result;
            }
            continue;
        } else {
            StatusCode status
                = commit_world_stage(world, trial.stage_high, stepper_cfg, wksp);
            if (status != StatusCode::ok) {
                result.status = status;
                return result;
            }

            world.advance_time(dt_trial);
            result.dt_sim_advanced += dt_trial;
            result.t = world.t_sim();
            result.final_error_norm = error_norm;

            f64 dt_abs = std::abs(dt_trial);
            if (result.stats.accepted_steps == 0) {
                result.stats.min_accepted_dt = dt_abs;
                result.stats.max_accepted_dt = dt_abs;
            } else {
                result.stats.min_accepted_dt
                    = std::min(dt_abs, result.stats.min_accepted_dt);
                result.stats.max_accepted_dt
                    = std::max(dt_abs, result.stats.max_accepted_dt);
            }

            ++result.stats.accepted_steps;
            result.stats.final_accepted_dt = dt_abs;

            if (final_attempt) {
                result.t = tf;
                result.status = StatusCode::ok;
                return result;
            }

            f64 dt_next = std::clamp(
                std::abs(dt_trial) * scale,
                adaptive_cfg.opts.dt_min,
                adaptive_cfg.opts.dt_max
            );
            dt_trial = dir * dt_next;
        }
    }
}

static bool step_tr_world_staged_rk1(
    World& world,
    f64 t,
    f64 dt,
    const WorldStepperWorkspace& wksp,
    bool stage_source_att
) {
    // rk1 stage 1
    WorldTrStage stage0_tr = build_tr_stage(world, wksp);
    WorldAttStage stage0_att = build_source_att_stage(world, wksp, stage_source_att);
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0_tr.ids) {
        const StateTr& x0 = stage0_tr.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0_tr, stage0_att, id, x0, wksp);
    }

    for (EntityId id : stage0_tr.ids) {
        StateTr x0 = stage0_tr.x.at(id);
        StateTr x_next = x0 + dt * k1.at(id);

        Body* body = world.body(id);
        if (body == nullptr) return false;
        body->x_tr = x_next;
    }

    return true;
}

static bool step_tr_world_staged_rk2(
    World& world,
    f64 t,
    f64 dt,
    const WorldStepperWorkspace& wksp,
    bool stage_source_att
) {
    // rk2 stage 1
    WorldTrStage stage0_tr = build_tr_stage(world, wksp);
    WorldAttStage stage0_att = build_source_att_stage(world, wksp, stage_source_att);
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0_tr.ids) {
        const StateTr& x0 = stage0_tr.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0_tr, stage0_att, id, x0, wksp);
    }

    // rk2 stage 2
    f64 k1_scale = 0.5 * dt;
    f64 dt_stage = 0.5 * dt;
    WorldTrStage stage1_tr
        = build_tr_trial_stage(stage0_tr, {{.k = &k1, .scale = k1_scale}});
    WorldAttStage stage1_att = build_source_att_trial_stage(world, stage0_att, dt_stage);
    umap<EntityId, DerivTr> k2;
    for (EntityId id : stage1_tr.ids) {
        const StateTr& x1 = stage1_tr.x.at(id);
        k2[id] = derivtr_world_staged(world, stage1_tr, stage1_att, id, x1, wksp);
    }

    for (EntityId id : stage0_tr.ids) {
        StateTr x0 = stage0_tr.x.at(id);
        StateTr x_next = x0 + dt * k2.at(id);

        Body* body = world.body(id);
        if (body == nullptr) return false;
        body->x_tr = x_next;
    }

    return true;
}

static bool step_tr_world_staged_rk2heun(
    World& world,
    f64 t,
    f64 dt,
    const WorldStepperWorkspace& wksp,
    bool stage_source_att
) {
    WorldTrStage stage0_tr = build_tr_stage(world, wksp);
    WorldAttStage stage0_att = build_source_att_stage(world, wksp, stage_source_att);

    // rk2 stage 1
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0_tr.ids) {
        const StateTr& x0 = stage0_tr.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0_tr, stage0_att, id, x0, wksp);
    }

    // rk2 stage 2
    f64 k1_scale = dt;
    f64 dt_stage = dt;
    WorldTrStage stage1_tr
        = build_tr_trial_stage(stage0_tr, {{.k = &k1, .scale = k1_scale}});
    WorldAttStage stage1_att = build_source_att_trial_stage(world, stage0_att, dt_stage);
    umap<EntityId, DerivTr> k2;
    for (EntityId id : stage1_tr.ids) {
        const StateTr& x1 = stage1_tr.x.at(id);
        k2[id] = derivtr_world_staged(world, stage1_tr, stage1_att, id, x1, wksp);
    }

    for (EntityId id : stage0_tr.ids) {
        StateTr x0 = stage0_tr.x.at(id);
        StateTr x_next = x0 + dt * (k1.at(id) + k2.at(id)) / 2.0;

        Body* body = world.body(id);
        if (body == nullptr) return false;
        body->x_tr = x_next;
    }

    return true;
}

static bool step_tr_world_staged_rk2ralston(
    World& world,
    f64 t,
    f64 dt,
    const WorldStepperWorkspace& wksp,
    bool stage_source_att
) {
    // rk2 stage 1
    WorldTrStage stage0_tr = build_tr_stage(world, wksp);
    WorldAttStage stage0_att = build_source_att_stage(world, wksp, stage_source_att);
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0_tr.ids) {
        const StateTr& x0 = stage0_tr.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0_tr, stage0_att, id, x0, wksp);
    }

    // rk2 stage 2
    f64 k1_scale = (2.0 / 3.0) * dt;
    f64 dt_stage = (2.0 / 3.0) * dt;
    WorldTrStage stage1_tr
        = build_tr_trial_stage(stage0_tr, {{.k = &k1, .scale = k1_scale}});
    WorldAttStage stage1_att = build_source_att_trial_stage(world, stage0_att, dt_stage);
    umap<EntityId, DerivTr> k2;
    for (EntityId id : stage1_tr.ids) {
        const StateTr& x1 = stage1_tr.x.at(id);
        k2[id] = derivtr_world_staged(world, stage1_tr, stage1_att, id, x1, wksp);
    }

    for (EntityId id : stage0_tr.ids) {
        StateTr x0 = stage0_tr.x.at(id);
        StateTr x_next = x0 + dt * (k1.at(id) / 4.0 + 3.0 / 4.0 * k2.at(id));

        Body* body = world.body(id);
        if (body == nullptr) return false;
        body->x_tr = x_next;
    }

    return true;
}

static bool step_tr_world_staged_rk3(
    World& world,
    f64 t,
    f64 dt,
    const WorldStepperWorkspace& wksp,
    bool stage_source_att
) {
    // rk3 stage 1
    WorldTrStage stage0_tr = build_tr_stage(world, wksp);
    WorldAttStage stage0_att = build_source_att_stage(world, wksp, stage_source_att);
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0_tr.ids) {
        const StateTr& x0 = stage0_tr.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0_tr, stage0_att, id, x0, wksp);
    }

    // rk3 stage 2
    f64 k1_scale = 0.5 * dt;
    f64 dt_stage = 0.5 * dt;
    WorldTrStage stage1_tr
        = build_tr_trial_stage(stage0_tr, {{.k = &k1, .scale = k1_scale}});
    WorldAttStage stage1_att = build_source_att_trial_stage(world, stage0_att, dt_stage);
    umap<EntityId, DerivTr> k2;
    for (EntityId id : stage1_tr.ids) {
        const StateTr& x1 = stage1_tr.x.at(id);
        k2[id] = derivtr_world_staged(world, stage1_tr, stage1_att, id, x1, wksp);
    }

    // rk3 stage 3
    f64 k1_stage3_scale = -dt;
    f64 k2_stage3_scale = 2.0 * dt;
    dt_stage = dt;
    WorldTrStage stage2_tr = build_tr_trial_stage(
        stage0_tr,
        {{.k = &k1, .scale = k1_stage3_scale}, {.k = &k2, .scale = k2_stage3_scale}}
    );
    WorldAttStage stage2_att = build_source_att_trial_stage(world, stage0_att, dt_stage);
    umap<EntityId, DerivTr> k3;
    for (EntityId id : stage2_tr.ids) {
        const StateTr& x2 = stage2_tr.x.at(id);
        k3[id] = derivtr_world_staged(world, stage2_tr, stage2_att, id, x2, wksp);
    }

    for (EntityId id : stage0_tr.ids) {
        StateTr x0 = stage0_tr.x.at(id);
        StateTr x_next = x0 + dt / 6.0 * (k1.at(id) + 4.0 * k2.at(id) + k3.at(id));

        Body* body = world.body(id);
        if (body == nullptr) return false;
        body->x_tr = x_next;
    }

    return true;
}

static bool step_tr_world_staged_rk4(
    World& world,
    f64 t,
    f64 dt,
    const WorldStepperWorkspace& wksp,
    bool stage_source_att
) {
    // rk4 stage 1
    WorldTrStage stage0_tr = build_tr_stage(world, wksp);
    WorldAttStage stage0_att = build_source_att_stage(world, wksp, stage_source_att);
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0_tr.ids) {
        const StateTr& x0 = stage0_tr.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0_tr, stage0_att, id, x0, wksp);
    }

    // rk4 stage 2
    f64 k1_scale = 0.5 * dt;
    f64 dt_stage = 0.5 * dt;
    WorldTrStage stage1_tr
        = build_tr_trial_stage(stage0_tr, {{.k = &k1, .scale = k1_scale}});
    WorldAttStage stage1_att = build_source_att_trial_stage(world, stage0_att, dt_stage);
    umap<EntityId, DerivTr> k2;
    for (EntityId id : stage1_tr.ids) {
        const StateTr& x1 = stage1_tr.x.at(id);
        k2[id] = derivtr_world_staged(world, stage1_tr, stage1_att, id, x1, wksp);
    }

    // rk4 stage 3
    f64 k2_scale = 0.5 * dt;
    dt_stage = 0.5 * dt;
    WorldTrStage stage2_tr
        = build_tr_trial_stage(stage0_tr, {{.k = &k2, .scale = k2_scale}});
    WorldAttStage stage2_att = build_source_att_trial_stage(world, stage0_att, dt_stage);
    umap<EntityId, DerivTr> k3;
    for (EntityId id : stage2_tr.ids) {
        const StateTr& x2 = stage2_tr.x.at(id);
        k3[id] = derivtr_world_staged(world, stage2_tr, stage2_att, id, x2, wksp);
    }

    // rk4 stage 4
    f64 k3_scale = dt;
    dt_stage = dt;
    WorldTrStage stage3_tr
        = build_tr_trial_stage(stage0_tr, {{.k = &k3, .scale = k3_scale}});
    WorldAttStage stage3_att = build_source_att_trial_stage(world, stage0_att, dt_stage);
    umap<EntityId, DerivTr> k4;
    for (EntityId id : stage3_tr.ids) {
        const StateTr& x3 = stage3_tr.x.at(id);
        k4[id] = derivtr_world_staged(world, stage3_tr, stage3_att, id, x3, wksp);
    }

    for (EntityId id : stage0_tr.ids) {
        StateTr x0 = stage0_tr.x.at(id);
        StateTr x_next
            = x0 + dt / 6.0 * (k1.at(id) + 2.0 * k2.at(id) + 2.0 * k3.at(id) + k4.at(id));

        Body* body = world.body(id);
        if (body == nullptr) return false;
        body->x_tr = x_next;
    }

    return true;
}

static bool step_tr_world_staged(
    World& world,
    f64 t,
    f64 dt,
    const WorldStepperConfig& cfg,
    const WorldStepperWorkspace& wksp
) {
    bool stage_source_att = cfg.step_att;

    switch (cfg.integrator_tr) {
    case IntegratorTypeFixed::rk1:
        return step_tr_world_staged_rk1(world, t, dt, wksp, stage_source_att);
    case IntegratorTypeFixed::rk2:
        return step_tr_world_staged_rk2(world, t, dt, wksp, stage_source_att);
    case IntegratorTypeFixed::rk2_heun:
        return step_tr_world_staged_rk2heun(world, t, dt, wksp, stage_source_att);
    case IntegratorTypeFixed::rk2_ralston:
        return step_tr_world_staged_rk2ralston(world, t, dt, wksp, stage_source_att);
    case IntegratorTypeFixed::rk3:
        return step_tr_world_staged_rk3(world, t, dt, wksp, stage_source_att);
    case IntegratorTypeFixed::rk4:
        return step_tr_world_staged_rk4(world, t, dt, wksp, stage_source_att);
    }

    return false;
}

WorldStepperStats step_world(World& world, f64 dt, const WorldStepperConfig& cfg) {
    WorldStepperWorkspace wksp;
    return step_world(world, dt, cfg, wksp);
}

WorldStepperStats step_world(
    World& world,
    f64 dt,
    const WorldStepperConfig& cfg,
    WorldStepperWorkspace& wksp
) {
    WorldStepperStats stats{.success = false};
    if (wksp.dirty) {
        rebuild_world_stepper_workspace(world, wksp);
    }

    if (cfg.substeps < 1) return stats;
    if (cfg.ticks < 1) return stats;
    if (!std::isfinite(cfg.dt_scale) || cfg.dt_scale <= 0.0) return stats;

    f64 dt_tick = dt * cfg.dt_scale;
    f64 dt_sub = dt_tick / cfg.substeps;

    svec<EntityId> att_ids = wksp.propagated_att_ids;
    svec<EntityId> cel_att_ids = wksp.celestial_att_ids;

    for (i32 tick = 0; tick < cfg.ticks; ++tick) {
        for (i32 substep = 0; substep < cfg.substeps; ++substep) {
            // translational
            if (cfg.step_tr) {
                bool step_ok
                    = step_tr_world_staged(world, world.t_sim(), dt_sub, cfg, wksp);
                if (!step_ok) {
                    stats.success = step_ok;
                    return stats;
                }
            }

            // attitude
            for (EntityId id : att_ids) {
                if (cfg.step_att) {
                    bool step_ok = step_att_world(world, id, dt_sub, cfg);
                    if (!step_ok) {
                        stats.success = step_ok;
                        return stats;
                    }
                }
            }

            // celestial attitude
            for (EntityId id : cel_att_ids) {
                if (cfg.step_att) {
                    bool step_ok = step_cel_att_world(world, id, dt_sub);
                    if (!step_ok) {
                        stats.success = step_ok;
                        return stats;
                    }
                }
            }

            world.advance_time(dt_sub);
            stats.dt_sim_advanced += dt_sub;
            stats.substeps_completed += 1;
        }
        stats.ticks_completed += 1;
    }

    stats.success = true;
    return stats;
}

WorldStepperStats step_world_tableau(
    World& world,
    f64 dt,
    const WorldStepperConfig& cfg
) {
    WorldStepperWorkspace wksp;
    return step_world_tableau(world, dt, cfg, wksp);
}

WorldStepperStats step_world_tableau(
    World& world,
    f64 dt,
    const WorldStepperConfig& cfg,
    WorldStepperWorkspace& wksp
) {
    WorldStepperStats stats{.success = false};

    if (!isfinite(dt) || cfg.substeps < 1 || cfg.ticks < 1) return stats;
    if (!finite_pos(cfg.dt_scale)) return stats;
    if (cfg.step_tr && cfg.step_att && cfg.integrator_tr != cfg.integrator_att) {
        return stats;
    }

    if (wksp.dirty) rebuild_world_stepper_workspace(world, wksp);

    f64 dt_tick = dt * cfg.dt_scale;
    f64 dt_sub = dt_tick / cfg.substeps;

    IntegratorTypeFixed integrator = cfg.step_att && !cfg.step_tr
                                         ? cfg.integrator_att
                                         : cfg.integrator_tr;

    for (i32 tick = 0; tick < cfg.ticks; ++tick) {
        for (i32 substep = 0; substep < cfg.substeps; ++substep) {
            if (dt_sub != 0.0 && (cfg.step_tr || cfg.step_att)) {
                StatusCode status = dispatch_world_fixed_tableau(
                    world,
                    world.t_sim(),
                    dt_sub,
                    integrator,
                    cfg,
                    wksp
                );
                if (status != StatusCode::ok) return stats;
            }

            world.advance_time(dt_sub);
            stats.dt_sim_advanced += dt_sub;
            ++stats.substeps_completed;
        }
        ++stats.ticks_completed;
    }

    stats.success = true;
    return stats;
}

static void accumulate_adaptive_stats(
    AdaptiveIntegratorStats& total,
    const AdaptiveIntegratorStats& current
) {
    i64 accepted_before = total.accepted_steps;

    total.attempted_steps += current.attempted_steps;
    total.accepted_steps += current.accepted_steps;
    total.rejected_steps += current.rejected_steps;
    total.deriv_evals += current.deriv_evals;

    if (current.accepted_steps == 0) return;

    if (accepted_before == 0) {
        total.min_accepted_dt = current.min_accepted_dt;
        total.max_accepted_dt = current.max_accepted_dt;
    } else {
        total.min_accepted_dt = std::min(total.min_accepted_dt, current.min_accepted_dt);
        total.max_accepted_dt = std::max(total.max_accepted_dt, current.max_accepted_dt);
    }
    total.final_accepted_dt = current.final_accepted_dt;
}

WorldAdaptiveStepResult step_world_adaptive(
    World& world,
    f64 dt,
    const WorldStepperConfig& stepper_cfg,
    const WorldAdaptiveIntegratorConfig& adaptive_cfg,
    WorldStepperWorkspace& wksp
) {
    WorldAdaptiveStepResult result;

    result.status = validate_adaptive_world_integrator_config(adaptive_cfg, stepper_cfg);
    if (result.status != StatusCode::ok) return result;

    if (!isfinite(dt) || !finite_pos(stepper_cfg.dt_scale)) {
        result.status = StatusCode::invalid_input;
        return result;
    }
    if (stepper_cfg.ticks < 1 || stepper_cfg.substeps < 1) {
        result.status = StatusCode::invalid_input;
        return result;
    }

    result.t = world.t_sim();
    if (wksp.dirty) {
        rebuild_world_stepper_workspace(world, wksp);
    }

    f64 dt_tick = dt * stepper_cfg.dt_scale;
    f64 dt_sub = 0.0;
    i32 subs;
    if (adaptive_cfg.use_substeps) {
        subs = stepper_cfg.substeps;
    } else {
        subs = 1;
    }
    dt_sub = dt_tick / subs;

    for (i32 tick = 0; tick < stepper_cfg.ticks; ++tick) {
        for (i32 substep = 0; substep < subs; ++substep) {
            f64 tf_sub = world.t_sim() + dt_sub;

            WorldAdaptiveStepResult sub_result;

            switch (adaptive_cfg.integrator_tr) {
            // TEMP: only dopri for now
            case IntegratorTypeAdaptive::dopri54: {
                sub_result = propagate_world_embedded_rk(
                    world,
                    tf_sub,
                    dopri54_tableau,
                    stepper_cfg,
                    adaptive_cfg,
                    wksp
                );
            } break;
            default: {
                sub_result.status = StatusCode::unsupported_method;
                sub_result.t = world.t_sim();
            } break;
            }

            result.t = sub_result.t;
            result.dt_sim_advanced += sub_result.dt_sim_advanced;
            accumulate_adaptive_stats(result.stats, sub_result.stats);

            if (sub_result.stats.accepted_steps > 0) {
                result.final_error_norm = sub_result.final_error_norm;
            }

            if (sub_result.status != StatusCode::ok) {
                result.status = sub_result.status;
                result.t = world.t_sim();
                return result;
            }

            ++result.substeps_completed;
        }

        ++result.ticks_completed;
    }

    result.status = StatusCode::ok;
    result.t = world.t_sim();
    return result;
}

WorldAdaptiveStepResult step_world_adaptive(
    World& world,
    f64 dt,
    const WorldStepperConfig& stepper_cfg,
    const WorldAdaptiveIntegratorConfig& adaptive_cfg
) {
    WorldStepperWorkspace wksp;
    return step_world_adaptive(world, dt, stepper_cfg, adaptive_cfg, wksp);
}

void rebuild_world_stepper_workspace(const World& world, WorldStepperWorkspace& wksp) {
    // TODO: add conditional/caching
    wksp.propagated_tr_ids = propagated_tr_ids(world);
    wksp.propagated_att_ids = propagated_att_ids(world);
    wksp.celestial_att_ids = celestial_att_ids(world);
    wksp.gravity_source_ids = gravity_source_ids(world);
    wksp.source_att_ids = source_att_ids(world);
    wksp.staged_att_ids = staged_att_ids(wksp.propagated_att_ids, wksp.celestial_att_ids);
    wksp.dirty = false;
}
