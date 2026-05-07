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
        if (body->body_type == BodyType::celestial)
            continue; // TODO: remove this, keep stationary gravity sources for now
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
        if (w_mag <= tol_strict) break;

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

WorldStepperStats step_world(World& world, f64 dt, const WorldStepperConfig& cfg) {
    WorldStepperStats stats{.success = false};
    if (cfg.integrator_tr != IntegratorType::rk1) return stats;
    if (cfg.substeps < 1) return stats;
    if (cfg.ticks < 1) return stats;
    if (!std::isfinite(cfg.time_scale) || cfg.time_scale <= 0.0) return stats;

    f64 dt_tick = dt * cfg.time_scale;
    f64 dt_sub = dt_tick / cfg.substeps;

    svec<EntityId> tr_ids = propagated_tr_ids(world);
    svec<EntityId> att_ids = propagated_att_ids(world);
    svec<EntityId> cel_att_ids = celestial_att_ids(world);

    for (i32 tick = 0; tick < cfg.ticks; ++tick) {
        for (i32 substep = 0; substep < cfg.substeps; ++substep) {
            // translational
            for (EntityId id : tr_ids) {
                if (cfg.step_translation) {
                    bool step_ok = step_tr_world(world, id, dt_sub, cfg);
                    if (!step_ok) {
                        stats.success = step_ok;
                        return stats;
                    }
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
