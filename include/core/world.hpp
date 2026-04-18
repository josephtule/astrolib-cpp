#pragma once

#include <memory>

#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/state.hpp"
#include "util/typedefs.hpp"

class SystemStepper;

class World {
  private:
    // Entity storage
    EntityId next_id = 1;
    svec<EntityId> active_ids;
    umap<EntityId, std::unique_ptr<Body>> bodies;

    // Simulation state
    f64 t_sim_ = 0.0;
    bool paused = false;

    EntityId allocate_id();
    EntityId insert_body(std::unique_ptr<Body> body);

  public:
    f64 t_sim() const;
    void reset_time(f64 t0 = 0.0);

    // Entity lifecycle and roles
    bool is_active(EntityId id) const;

    // Entity getters
    Body* body(EntityId id);
    const Body* body(EntityId id) const;
    Celestial* celestial(EntityId id);
    const Celestial* celestial(EntityId id) const;
    Satellite* satellite(EntityId id);
    const Satellite* satellite(EntityId id) const;
    Station* station(EntityId id);
    const Station* station(EntityId id) const;

    // Entity spawning and adding
    EntityId spawn_celestial();
    EntityId spawn_satellite();
    EntityId spawn_station();
    EntityId insert_celestial(std::unique_ptr<Celestial> cel);
    EntityId insert_satellite(std::unique_ptr<Satellite> sat);
    EntityId insert_station(std::unique_ptr<Station> stat);

    // Type and emission queries
    bool is_celestial(EntityId id) const;
    bool is_satellite(EntityId id) const;
    bool is_station(EntityId id) const;
    bool emits_gravity(EntityId id) const;
    bool emits_radiation(EntityId id) const;

    // Force queries
    vec3d gravity_accel_on(EntityId target_id) const;
    vec3d gravity_accel_from(EntityId target_id, EntityId source_id) const;

    // Station
    vec3d stat_r_inertial(EntityId station_id) const;
    vec3d stat_v_inertial(EntityId station_id) const;
    StateTr stat_x_tr_inertial(EntityId station_id) const;

    // Inertial frame helpers
    vec3d body_z_inertial(EntityId body_id) const;
    vec3d body_w_inertial(EntityId body_id) const;

  private:
    friend class SystemStepper; // Allow private access
};