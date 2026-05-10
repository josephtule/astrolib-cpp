#include "core/world_stepper.hpp"
#include "core/body.hpp"
#include "core/dynamics_rotational.hpp"
#include "core/entity.hpp"
#include "core/integrator.hpp"
#include "core/state.hpp"
#include "core/transform.hpp"
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

static WorldTrStage build_tr_stage(const World& world) {
    WorldTrStage stage;

    svec<EntityId> ids = propagated_tr_ids(world);
    for (EntityId id : ids) {
        const Body* body = world.body(id);
        if (body == nullptr) continue;

        stage.ids.push_back(id);
        stage.x.emplace(id, body->x_tr);
    }
    return stage;
}

static WorldAttStage build_source_att_stage(const World& world) {
    WorldAttStage stage;

    svec<EntityId> ids = world.active_entity_ids();
    for (EntityId id : ids) {
        const Body* body = world.body(id);
        if (body == nullptr || !body->emits_gravity) continue; // TODO: add has atmosphere
        // only bodies that affect force via their attitude used

        stage.ids.push_back(id);
        stage.x.emplace(id, body->x_att);
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
    const StateTr& x_tr_target
) {
    vec3d a = vec3d0;
    const Body* target = world.body(target_id);
    if (target == nullptr) return a;

    svec<EntityId> ids = world.active_entity_ids();
    for (EntityId id : ids) {
        if (target_id == id) continue;
        const Body* body = world.body(id);
        if (body == nullptr || !body->emits_gravity) continue;
        const Celestial* source = world.celestial(id);
        if (source == nullptr) continue;
        StateTr x_tr_source = source_tr_from_stage_or_world(world, stage_tr, id);
        StateAtt x_att_source = source_att_from_stage_or_world(world, stage_att, id);
        a += world.gravity_accel_from(
            target_id,
            x_tr_target,
            id,
            x_tr_source,
            x_att_source
        );
    }

    return a;
}

static DerivTr derivtr_world_staged(
    const World& world,
    const WorldTrStage& stage_tr,
    const WorldAttStage& stage_att,
    EntityId target_id,
    const StateTr& x_tr_target
) {
    DerivTr dx;
    dx.dr = x_tr_target.v;
    dx.dv = staged_gravity_accel_on(world, stage_tr, stage_att, target_id, x_tr_target);
    return dx;
}

static bool step_tr_world_staged_rk1(World& world, f64 t, f64 dt) {
    // rk1 stage 1
    WorldTrStage stage0_tr = build_tr_stage(world);
    WorldAttStage stage0_att = build_source_att_stage(world);
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0_tr.ids) {
        const StateTr& x0 = stage0_tr.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0_tr, stage0_att, id, x0);
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

static bool step_tr_world_staged_rk2(World& world, f64 t, f64 dt) {
    // rk2 stage 1
    WorldTrStage stage0_tr = build_tr_stage(world);
    WorldAttStage stage0_att = build_source_att_stage(world);
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0_tr.ids) {
        const StateTr& x0 = stage0_tr.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0_tr, stage0_att, id, x0);
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
        k2[id] = derivtr_world_staged(world, stage1_tr, stage1_att, id, x1);
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

static bool step_tr_world_staged_rk2heun(World& world, f64 t, f64 dt) {
    WorldTrStage stage0_tr = build_tr_stage(world);
    WorldAttStage stage0_att = build_source_att_stage(world);

    // rk2 stage 1
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0_tr.ids) {
        const StateTr& x0 = stage0_tr.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0_tr, stage0_att, id, x0);
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
        k2[id] = derivtr_world_staged(world, stage1_tr, stage1_att, id, x1);
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

static bool step_tr_world_staged_rk2ralston(World& world, f64 t, f64 dt) {
    // rk2 stage 1
    WorldTrStage stage0_tr = build_tr_stage(world);
    WorldAttStage stage0_att = build_source_att_stage(world);
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0_tr.ids) {
        const StateTr& x0 = stage0_tr.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0_tr, stage0_att, id, x0);
    }

    // rk2 stage 2
    f64 k1_scale = 0.75 * dt;
    f64 dt_stage = 0.75 * dt;
    WorldTrStage stage1_tr
        = build_tr_trial_stage(stage0_tr, {{.k = &k1, .scale = k1_scale}});
    WorldAttStage stage1_att = build_source_att_trial_stage(world, stage0_att, dt_stage);
    umap<EntityId, DerivTr> k2;
    for (EntityId id : stage1_tr.ids) {
        const StateTr& x1 = stage1_tr.x.at(id);
        k2[id] = derivtr_world_staged(world, stage1_tr, stage1_att, id, x1);
    }

    for (EntityId id : stage0_tr.ids) {
        StateTr x0 = stage0_tr.x.at(id);
        StateTr x_next = x0 + dt * (k1.at(id) / 3.0 + 2.0 / 3.0 * k2.at(id));

        Body* body = world.body(id);
        if (body == nullptr) return false;
        body->x_tr = x_next;
    }

    return true;
}

static bool step_tr_world_staged_rk3(World& world, f64 t, f64 dt) {
    // rk3 stage 1
    WorldTrStage stage0_tr = build_tr_stage(world);
    WorldAttStage stage0_att = build_source_att_stage(world);
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0_tr.ids) {
        const StateTr& x0 = stage0_tr.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0_tr, stage0_att, id, x0);
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
        k2[id] = derivtr_world_staged(world, stage1_tr, stage1_att, id, x1);
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
        k3[id] = derivtr_world_staged(world, stage2_tr, stage2_att, id, x2);
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

static bool step_tr_world_staged_rk4(World& world, f64 t, f64 dt) {
    // rk4 stage 1
    WorldTrStage stage0_tr = build_tr_stage(world);
    WorldAttStage stage0_att = build_source_att_stage(world);
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0_tr.ids) {
        const StateTr& x0 = stage0_tr.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0_tr, stage0_att, id, x0);
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
        k2[id] = derivtr_world_staged(world, stage1_tr, stage1_att, id, x1);
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
        k3[id] = derivtr_world_staged(world, stage2_tr, stage2_att, id, x2);
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
        k4[id] = derivtr_world_staged(world, stage3_tr, stage3_att, id, x3);
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

static bool step_tr_world_staged(World& world, f64 t, f64 dt, IntegratorType integrator) {
    switch (integrator) {
    case IntegratorType::rk1: return step_tr_world_staged_rk1(world, t, dt);
    case IntegratorType::rk2: return step_tr_world_staged_rk2(world, t, dt);
    case IntegratorType::rk2_heun: return step_tr_world_staged_rk2heun(world, t, dt);
    case IntegratorType::rk2_ralston:
        return step_tr_world_staged_rk2ralston(world, t, dt);
    case IntegratorType::rk3: return step_tr_world_staged_rk3(world, t, dt);
    case IntegratorType::rk4: return step_tr_world_staged_rk4(world, t, dt);
    }

    return false;
}

WorldStepperStats step_world(World& world, f64 dt, const WorldStepperConfig& cfg) {
    WorldStepperStats stats{.success = false};
    if (cfg.substeps < 1) return stats;
    if (cfg.ticks < 1) return stats;
    if (!std::isfinite(cfg.time_scale) || cfg.time_scale <= 0.0) return stats;

    f64 dt_tick = dt * cfg.time_scale;
    f64 dt_sub = dt_tick / cfg.substeps;

    svec<EntityId> att_ids = propagated_att_ids(world);
    svec<EntityId> cel_att_ids = celestial_att_ids(world);

    for (i32 tick = 0; tick < cfg.ticks; ++tick) {
        for (i32 substep = 0; substep < cfg.substeps; ++substep) {
            // translational
            if (cfg.step_translation) {
                bool step_ok = step_tr_world_staged(
                    world,
                    world.t_sim(),
                    dt_sub,
                    cfg.integrator_tr
                );
                if (!step_ok) {
                    stats.success = step_ok;
                    return stats;
                }
            }

            // attitude
            for (EntityId id : att_ids) {
                if (cfg.step_attitude) {
                    bool step_ok = step_att_world(world, id, dt_sub, cfg);
                    if (!step_ok) {
                        stats.success = step_ok;
                        return stats;
                    }
                }
            }

            // celestial attitude
            for (EntityId id : cel_att_ids) {
                if (cfg.step_attitude) {
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
