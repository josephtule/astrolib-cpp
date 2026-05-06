#include "core/world_stepper.hpp"
#include "Eigen/Core"
#include "core/body.hpp"
#include "core/dynamics_rotational.hpp"
#include "core/entity.hpp"
#include "core/integrator.hpp"
#include "core/state.hpp"
#include "core/world.hpp"
#include "util/typedefs.hpp"

#include <cmath>

// NOTE: this is a preliminary implementation

static svec<EntityId> propagated_tr_ids(const World& world) {
    svec<EntityId> ids;

    for (EntityId id : world.active_entity_ids()) {
        const Body* body = world.body(id);
        if (body == nullptr) continue;
        if (!body->propagate_tr) continue;
        if (body->body_type == BodyType::celestial)
            continue; // TODO: remove this, keep stationary gravity sources for now
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
        if (body->body_type == BodyType::celestial) continue; // TODO: remove this later
        ids.push_back(id);
    }

    return ids;
}

DerivTr derivtr_world(const World& world, EntityId id, const StateTr& x) {
    // TODO: x currently unused, higher order world stepping requires staging
    DerivTr dx;
    dx.dr = x.v;
    dx.dv = world.gravity_accel_on(id);

    return dx;
}

DerivAtt derivatt_world(const World& world, EntityId id, const StateAtt& x) {
    DerivAtt dx;

    const Body* body = world.body(id);
    switch (body->body_type) {
    case BodyType::celestial: break;
    case BodyType::satellite: {
        const Satellite* sat = world.satellite(id);
        if (sat->principal_axes) {
            dx = d_rigidbody(world.t_sim(), x, sat->I);
        } else {
            dx = d_rigidbody(world.t_sim(), x, sat->I, sat->I_inv);
        }
        // TODO: add non principal axes later
    } break;
    case BodyType::station:
        break;
        // TODO: not sure how to deal with stations, only allow in if free station
        // otherwise update with celestial body (always align with relative position
        // vector?)
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
        cfg.integrator
    );
    body->x_tr = tx.second;

    return true;
}

bool step_att_world(World& world, EntityId id, f64 dt, const WorldStepperConfig& cfg) {
    Body* body = world.body(id);
    if (body == nullptr) return false;
    if (body->body_type == BodyType::celestial || body->body_type == BodyType::station)
        return false; // TODO: this only allows satellites for now

    auto f = [&](f64 t, StateAtt x) -> DerivAtt { return derivatt_world(world, id, x); };
    auto tx = step_integrator<StateAtt, DerivAtt>(f, world.t_sim(), body->x_att, dt);
    body->x_att = tx.second;

    return true;
}

WorldStepperStats step_world(World& world, f64 dt, const WorldStepperConfig& cfg) {
    WorldStepperStats stats{.success = false};
    if (cfg.integrator != IntegratorType::rk1) return stats;
    if (cfg.substeps < 1) return stats;
    if (cfg.ticks < 1) return stats;
    if (!std::isfinite(cfg.time_scale) || cfg.time_scale <= 0.0) return stats;

    f64 dt_tick = dt * cfg.time_scale;
    f64 dt_sub = dt_tick / cfg.substeps;
    svec<EntityId> ids = propagated_tr_ids(world);

    for (i32 tick = 0; tick < cfg.ticks; ++tick) {
        for (i32 substep = 0; substep < cfg.substeps; ++substep) {
            for (EntityId id : ids) {
                if (cfg.step_translation) {
                    bool step_ok = step_tr_world(world, id, dt_sub, cfg);
                    if (!step_ok) {
                        stats.success = step_ok;
                        return stats;
                    }
                }

                if (cfg.step_attitude) {
                    bool step_ok = step_att_world(world, id, dt_sub, cfg);
                    if (!step_ok) {
                        stats.success = step_ok;
                        return stats;
                    }
                }
            }
            world.advance_time(dt_sub);
            stats.dt_sim_advanced += dt;
            stats.substeps_completed += 1;
        }
        stats.ticks_completed += 1;
    }

    return stats;
}