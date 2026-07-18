// Copyright 2025-2026 Joseph Tu Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/state.hpp"
#include "core/time.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"

class SystemStepper;

struct BodySnapshot {
    EntityId id = kInvalidEntityId;
    std::string name;
    BodyType body_type = BodyType::unknown;

    StateTr x_tr;
    StateAtt x_att;

    bool propagate_tr = true;
    bool propagate_att = false;

    bool emits_gravity = false;
    bool emits_radiation = false;
};

struct WorldStateSnapshot {
    // State-only checkpoint: typed body config/geometry is intentionally not captured
    // yet. Full hard restore later needs derived body snapshots and exact body storage
    // rebuilds.

    // Entity storage
    EntityId next_id = 1;
    svec<EntityId> active_ids;
    svec<BodySnapshot> bodies;

    // Simulation state
    f64 t_sim_ = 0.0;
};

enum struct BodyFilterMode { all = 0, active = 1, inactive = 2 };

class World {
  private:
    // Entity storage
    EntityId next_id = 1;
    svec<EntityId> all_ids;
    svec<EntityId> active_ids;
    svec<EntityId> inactive_ids;
    umap<EntityId, uptr<Body>> bodies;

    // Simulation state
    f64 t_sim_ = 0.0;
    bool paused = false;

    JulianDate jd{};
    TimeScale time_scale = TimeScale::utc;
    bool date_active = false;

    EntityId allocate_id();
    EntityId insert_body(uptr<Body> body);

  public:
    f64 t_sim() const;

    bool is_date_active() const { return date_active; }
    void set_date(const JulianDate jd);
    void set_date(const ModifiedJulianDate mjd);
    void set_date(const CalendarTime cal);
    JulianDate get_date_jd() const;
    ModifiedJulianDate get_date_mjd() const;
    CalendarTime get_date_cal() const;

    void reset_time(f64 t0 = 0.0);
    void advance_time(f64 dt);

    // World Snapshot saving and loading
    WorldStateSnapshot capture_checkpoint() const;
    bool restore_checkpoint_state(const WorldStateSnapshot& snapshot); // soft restore
    // bool restore_checkpoint(const WorldStateSnapshot& snapshot); // hard restore

    // Entity lifecycle and roles
    bool is_active(EntityId id) const;
    bool make_inactive(EntityId id);
    bool make_active(EntityId id);
    const svec<EntityId>& active_entity_ids() const;
    const svec<EntityId>& all_entity_ids() const;

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
    EntityId insert_celestial(uptr<Celestial> cel);
    EntityId insert_satellite(uptr<Satellite> sat);
    EntityId insert_station(uptr<Station> stat);
    // void insert_celestials(svec<uptr<Celestial>> cels);
    void insert_satellites(svec<uptr<Satellite>> sats);
    // void insert_stations(svec<uptr<Station>>stats);
    i32 num_active_celestials() const;
    i32 num_active_satellites() const;
    i32 num_active_stations() const;
    i32 num_celestials() const;
    i32 num_satellites() const;
    i32 num_stations() const;
    i32 num_inactive_celestials() const;
    i32 num_inactive_satellites() const;
    i32 num_inactive_stations() const;
    // TODO: add duplicate helpers

    // Type and emission queries
    bool is_celestial(EntityId id) const;
    bool is_satellite(EntityId id) const;
    bool is_station(EntityId id) const;
    svec<EntityId> celestial_ids(BodyFilterMode mode = BodyFilterMode::all) const;
    svec<EntityId> satellite_ids(BodyFilterMode mode = BodyFilterMode::all) const;
    svec<EntityId> station_ids(BodyFilterMode mode = BodyFilterMode::all) const;
    svec<EntityId> active_celestial_ids() const;
    svec<EntityId> active_satellite_ids() const;
    svec<EntityId> active_station_ids() const;
    svec<EntityId> all_celestial_ids() const;
    svec<EntityId> all_satellite_ids() const;
    svec<EntityId> all_station_ids() const;
    svec<EntityId> inactive_celestial_ids() const;
    svec<EntityId> inactive_satellite_ids() const;
    svec<EntityId> inactive_station_ids() const;
    bool emits_gravity(EntityId id) const;
    bool emits_radiation(EntityId id) const;

    // Force queries
    vec3d gravity_accel_on(EntityId target_id) const;
    vec3d gravity_accel_on(EntityId target_id, const StateTr& x_target) const;
    vec3d gravity_accel_from(EntityId target_id, EntityId source_id) const;
    vec3d gravity_accel_from(
        EntityId target_id,
        const StateTr& target_x,
        EntityId source_id
    ) const;
    vec3d gravity_accel_from(
        EntityId target_id,
        const StateTr& x_target,
        EntityId source_id,
        const StateTr& x_source
    ) const;
    vec3d gravity_accel_from(
        EntityId target_id,
        const StateTr& x_target,
        EntityId source_id,
        const StateTr& x_source,
        const StateAtt& x_source_att
    ) const;

    // Station
    vec3d stat_r_inertial(EntityId station_id) const;
    vec3d stat_v_inertial(EntityId station_id) const;
    StateTr stat_x_tr_inertial(EntityId station_id) const;
    vec4d stat_q_inertial(EntityId station_id) const;
    StateAtt stat_x_att_inertial(EntityId station_id) const;
    bool set_stat_anchor_detic(
        EntityId station_id,
        EntityId anchor_id,
        const vec3d& llh,
        UAngle angle_in = UAngle::degree
    );
    mat3d stat_rot_enu_from_body(EntityId station_id) const;
    vec3d stat_rel_enu(EntityId station_id, EntityId target_id) const;
    vec3d stat_rel_enu(EntityId station_id, const StateTr& x_tr_target) const;

    // Inertial frame helpers
    vec3d body_z_inertial(EntityId body_id) const;
    vec3d body_w_inertial(EntityId body_id) const;

  private:
    friend class SystemStepper; // Allow private access
};
