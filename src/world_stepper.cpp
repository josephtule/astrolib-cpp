#include "core/world_stepper.hpp"
#include "Eigen/Core"
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
        f64 w_mag = cel->x_att.w.norm();
        if (w_mag <= tol12) break;

        vec3d axis = cel->x_att.w / w_mag;
        f64 theta = w_mag * dt;
        vec4d dq;
        dq << axis * std::sin(theta / 2.0), std::cos(theta / 2.0);

        cel->x_att.q = ep_mult(dq, cel->x_att.q);
        normalize_quaternion_inplace<f64>(cel->x_att.q);
    } break;
    case CelestialAttitudeModel::provider: {
        return false; // TODO: add provider later
    } break;
    }

    return true;
}

template <typename State>
struct WorldTrStage {
    svec<EntityId> ids;
    umap<EntityId, State> x;
    // TODO: add attitude later?
};
template <typename Deriv>
struct DerivWeight {
    // K_i in RK integrators (translational)
    const umap<EntityId, Deriv>* k = nullptr;
    f64 scale = 0.0;
};

template <typename State>
static WorldTrStage<State> build_tr_stage(const World& world) {
    WorldTrStage<State> stage;

    svec<EntityId> ids = propagated_tr_ids(world);
    for (EntityId id : ids) {
        const Body* body = world.body(id);
        if (body == nullptr) continue;

        stage.ids.push_back(id);
        stage.x.emplace(id, body->x_tr);
    }
    return stage;
}

// Gravity sources use their staged translational state when propagated
// fixed sources fall back to their stored world state.
template <typename State>
static StateTr source_tr_from_stage_or_world(
    const World& world,
    const WorldTrStage<State>& stage,
    EntityId source_id
) {
    // NOTE: currently only celestials are valid gravity emitters.
    const Celestial* cel = world.celestial(source_id);
    if (cel == nullptr) return StateTr{};

    auto it = stage.x.find(source_id);
    if (it == stage.x.end()) return cel->x_tr; // world fallback

    return it->second;
}

template <typename State>
static WorldTrStage<State> build_tr_trial_stage(
    const WorldTrStage<State>& base_stage,
    std::initializer_list<DerivWeight<DerivTr>> weights
) {
    WorldTrStage<State> trial;
    trial.ids = base_stage.ids;

    for (EntityId id : base_stage.ids) {
        StateTr x = base_stage.x.at(id);

        for (const DerivWeight<DerivTr>& weight : weights) {
            if (weight.k == nullptr) continue;
            x += weight.scale * weight.k->at(id);
        }

        trial.x.emplace(id, x);
    }

    return trial;
}

template <typename State>
static vec3d staged_gravity_accel_on(
    const World& world,
    const WorldTrStage<State>& stage,
    EntityId target_id,
    const StateTr& x_target
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
        StateTr x_source = source_tr_from_stage_or_world(world, stage, id);
        a += world.gravity_accel_from(target_id, x_target, id, x_source);
    }

    return a;
}

template <typename State>
static DerivTr derivtr_world_staged(
    const World& world,
    const WorldTrStage<State>& stage,
    EntityId target_id,
    const StateTr& x_target
) {
    DerivTr dx;
    dx.dr = x_target.v;
    dx.dv = staged_gravity_accel_on(world, stage, target_id, x_target);
    return dx;
}

template <typename State>
static bool step_tr_world_staged_rk1(World& world, f64 t, f64 dt) {
    WorldTrStage stage0 = build_tr_stage<State>(world);

    // rk1 stage 1
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0.ids) {
        const StateTr& x0 = stage0.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0, id, x0);
    }

    for (EntityId id : stage0.ids) {
        StateTr x0 = stage0.x.at(id);
        StateTr x_next = x0 + dt * k1.at(id);

        Body* body = world.body(id);
        if (body == nullptr) return false;
        body->x_tr = x_next;
    }

    return true;
}

template <typename State>
static bool step_tr_world_staged_rk2(World& world, f64 t, f64 dt) {
    WorldTrStage stage0 = build_tr_stage<State>(world);
    umap<EntityId, StateTr> x_next;

    // rk2 stage 1
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0.ids) {
        const StateTr& x0 = stage0.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0, id, x0);
    }

    // rk2 stage 2
    WorldTrStage stage1 = build_tr_trial_stage(stage0, {{.k = &k1, .scale = 0.5 * dt}});
    umap<EntityId, DerivTr> k2;
    for (EntityId id : stage1.ids) {
        const StateTr& x1 = stage1.x.at(id);
        k2[id] = derivtr_world_staged(world, stage1, id, x1);
    }

    for (EntityId id : stage0.ids) {
        StateTr x0 = stage0.x.at(id);
        StateTr x_next = x0 + dt * k2.at(id);

        Body* body = world.body(id);
        if (body == nullptr) return false;
        body->x_tr = x_next;
    }

    return true;
}

template <typename State>
static bool step_tr_world_staged_rk2heun(World& world, f64 t, f64 dt) {
    WorldTrStage stage0 = build_tr_stage<State>(world);
    umap<EntityId, StateTr> x_next;

    // rk2 stage 1
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0.ids) {
        const StateTr& x0 = stage0.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0, id, x0);
    }

    // rk2 stage 2
    WorldTrStage stage1 = build_tr_trial_stage(stage0, {{.k = &k1, .scale = dt}});
    umap<EntityId, DerivTr> k2;
    for (EntityId id : stage1.ids) {
        const StateTr& x1 = stage1.x.at(id);
        k2[id] = derivtr_world_staged(world, stage1, id, x1);
    }

    for (EntityId id : stage0.ids) {
        StateTr x0 = stage0.x.at(id);
        StateTr x_next = x0 + dt * (k1.at(id) + k2.at(id)) / 2.0;

        Body* body = world.body(id);
        if (body == nullptr) return false;
        body->x_tr = x_next;
    }

    return true;
}

template <typename State>
static bool step_tr_world_staged_rk2ralston(World& world, f64 t, f64 dt) {
    WorldTrStage stage0 = build_tr_stage<State>(world);
    umap<EntityId, StateTr> x_next;

    // rk2 stage 1
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0.ids) {
        const StateTr& x0 = stage0.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0, id, x0);
    }

    // rk2 stage 2
    WorldTrStage stage1
        = build_tr_trial_stage(stage0, {{.k = &k1, .scale = 3.0 / 4.0 * dt}});
    umap<EntityId, DerivTr> k2;
    for (EntityId id : stage1.ids) {
        const StateTr& x1 = stage1.x.at(id);
        k2[id] = derivtr_world_staged(world, stage1, id, x1);
    }

    for (EntityId id : stage0.ids) {
        StateTr x0 = stage0.x.at(id);
        StateTr x_next = x0 + dt * (k1.at(id) / 3.0 + 2.0 / 3.0 * k2.at(id));

        Body* body = world.body(id);
        if (body == nullptr) return false;
        body->x_tr = x_next;
    }

    return true;
}

template <typename State>
static bool step_tr_world_staged_rk3(World& world, f64 t, f64 dt) {
    umap<EntityId, StateTr> x_next;

    // rk3 stage 1
    WorldTrStage stage0 = build_tr_stage<State>(world);
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0.ids) {
        const StateTr& x0 = stage0.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0, id, x0);
    }

    // rk3 stage 2
    WorldTrStage stage1 = build_tr_trial_stage(stage0, {{.k = &k1, .scale = dt / 2.0}});
    umap<EntityId, DerivTr> k2;
    for (EntityId id : stage1.ids) {
        const StateTr& x1 = stage1.x.at(id);
        k2[id] = derivtr_world_staged(world, stage1, id, x1);
    }

    // rk3 stage 3
    WorldTrStage stage2 = build_tr_trial_stage(
        stage0,
        {{.k = &k1, .scale = -dt}, {.k = &k2, .scale = 2.0 * dt}}
    );
    umap<EntityId, DerivTr> k3;
    for (EntityId id : stage2.ids) {
        const StateTr& x2 = stage2.x.at(id);
        k3[id] = derivtr_world_staged(world, stage2, id, x2);
    }

    for (EntityId id : stage0.ids) {
        StateTr x0 = stage0.x.at(id);
        StateTr x_next = x0 + dt / 6.0 * (k1.at(id) + 4.0 * k2.at(id) + k3.at(id));

        Body* body = world.body(id);
        if (body == nullptr) return false;
        body->x_tr = x_next;
    }

    return true;
}

template <typename State>
static bool step_tr_world_staged_rk4(World& world, f64 t, f64 dt) {
    umap<EntityId, StateTr> x_next;

    // rk4 stage 1
    WorldTrStage stage0 = build_tr_stage<State>(world);
    umap<EntityId, DerivTr> k1;
    for (EntityId id : stage0.ids) {
        const StateTr& x0 = stage0.x.at(id);
        k1[id] = derivtr_world_staged(world, stage0, id, x0);
    }

    // rk4 stage 2
    WorldTrStage stage1 = build_tr_trial_stage(stage0, {{.k = &k1, .scale = dt / 2.0}});
    umap<EntityId, DerivTr> k2;
    for (EntityId id : stage1.ids) {
        const StateTr& x1 = stage1.x.at(id);
        k2[id] = derivtr_world_staged(world, stage1, id, x1);
    }

    // rk4 stage 3
    WorldTrStage stage2 = build_tr_trial_stage(stage0, {{.k = &k2, .scale = dt / 2.0}});
    umap<EntityId, DerivTr> k3;
    for (EntityId id : stage2.ids) {
        const StateTr& x2 = stage2.x.at(id);
        k3[id] = derivtr_world_staged(world, stage2, id, x2);
    }

    // rk4 stage 4
    WorldTrStage stage3 = build_tr_trial_stage(stage0, {{.k = &k3, .scale = dt}});
    umap<EntityId, DerivTr> k4;
    for (EntityId id : stage3.ids) {
        const StateTr& x3 = stage3.x.at(id);
        k4[id] = derivtr_world_staged(world, stage3, id, x3);
    }

    for (EntityId id : stage0.ids) {
        StateTr x0 = stage0.x.at(id);
        StateTr x_next
            = x0 + dt / 6.0 * (k1.at(id) + 2.0 * k2.at(id) + 2.0 * k3.at(id) + k4.at(id));

        Body* body = world.body(id);
        if (body == nullptr) return false;
        body->x_tr = x_next;
    }

    return true;
}

template <typename State>
static bool step_tr_world_staged(World& world, f64 t, f64 dt, IntegratorType integrator) {
    switch (integrator) {
    case IntegratorType::rk1: return step_tr_world_staged_rk1<State>(world, t, dt);
    case IntegratorType::rk2: return step_tr_world_staged_rk2<State>(world, t, dt);
    case IntegratorType::rk2_heun: return step_tr_world_staged_rk2heun<State>(world, t, dt);
    case IntegratorType::rk2_ralston:
        return step_tr_world_staged_rk2ralston<State>(world, t, dt);
    case IntegratorType::rk3: return step_tr_world_staged_rk3<State>(world, t, dt);
    case IntegratorType::rk4: return step_tr_world_staged_rk4<State>(world, t, dt);
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
                bool step_ok = step_tr_world_staged<StateTr>(
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
