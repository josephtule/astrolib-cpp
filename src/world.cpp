// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/world.hpp"
#include "core/body.hpp"
#include "core/dynamics_translational.hpp"
#include "core/entity.hpp"
#include "core/state.hpp"
#include "core/station_geometry.hpp"
#include "core/time.hpp"
#include "core/transform.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include <cstddef>
#include <memory>

f64 World::t_sim() const { return t_sim_; }

void World::reset_time(f64 t0) {
    if (date_active) {
        jd.frac = jd.frac - (t_sim_ - t0) / 86400.0;
        jd = normalize_jd(jd);
    }
    t_sim_ = t0;
}

void World::advance_time(f64 dt) {
    if (date_active) {
        jd.frac = jd.frac + dt / 86400.0;
        jd = normalize_jd(jd);
    }
    t_sim_ += dt;
}

void World::set_date(JulianDate jd_in, TimeScale scale) {
    date_active = true;
    time_scale = scale;
    jd = normalize_jd(jd_in);
}
void World::set_date(ModifiedJulianDate mjd, TimeScale scale) {
    date_active = true;
    time_scale = scale;
    jd = normalize_jd(mjd_to_jd(mjd));
}
void World::set_date(CalendarTime cal, TimeScale scale) {
    date_active = true;
    time_scale = scale;
    jd = normalize_jd(cal_to_jd(cal));
}

JulianDate World::get_date_jd() const { return jd; }
ModifiedJulianDate World::get_date_mjd() const { return jd_to_mjd(jd); }
CalendarTime World::get_date_cal() const { return jd_to_cal(jd); }

WorldStateSnapshot World::capture_checkpoint() const {
    WorldStateSnapshot snapshot;
    snapshot.next_id = this->next_id;
    snapshot.active_ids = active_ids;
    snapshot.t_sim_ = t_sim_;
    snapshot.bodies.reserve(bodies.size());
    for (const i32 id : active_ids) {
        const Body* body = this->body(id);
        if (body != nullptr) {
            BodySnapshot body_snapshot;
            body_snapshot.id = body->id;
            body_snapshot.name = body->name;
            body_snapshot.body_type = body->body_type;
            body_snapshot.x_tr = body->x_tr;
            body_snapshot.x_att = body->x_att;
            body_snapshot.propagate_tr = body->propagate_tr;
            body_snapshot.propagate_att = body->propagate_att;
            body_snapshot.emits_gravity = body->emits_gravity;
            body_snapshot.emits_radiation = body->emits_radiation;
            snapshot.bodies.push_back(body_snapshot);
        }
    }
    return snapshot;
}

bool World::restore_checkpoint_state(const WorldStateSnapshot& snapshot) {
    // this is a soft restore, only restores the states of the bodies
    // does not support adding/removing bodies in between capture point and restore point

    if (active_ids.size() != snapshot.active_ids.size() || next_id != snapshot.next_id) {
        return false;
    }
    for (i32 i = 0; i < active_ids.size(); ++i) {
        if (active_ids[i] != snapshot.active_ids[i]) return false;
    }

    // validate bodies
    for (const BodySnapshot& body_snapshot : snapshot.bodies) {
        Body* temp = body(body_snapshot.id);
        if (temp == nullptr) return false;
        if (body_snapshot.body_type != temp->body_type) return false;
    }

    this->next_id = snapshot.next_id;
    this->active_ids = snapshot.active_ids;
    this->t_sim_ = snapshot.t_sim_;

    for (const BodySnapshot& body_snapshot : snapshot.bodies) {
        Body* temp = body(body_snapshot.id);
        temp->id = body_snapshot.id;
        temp->name = body_snapshot.name;
        temp->body_type = body_snapshot.body_type;
        temp->x_tr = body_snapshot.x_tr;
        temp->x_att = body_snapshot.x_att;
        temp->propagate_tr = body_snapshot.propagate_tr;
        temp->propagate_att = body_snapshot.propagate_att;
        temp->emits_gravity = body_snapshot.emits_gravity;
        temp->emits_radiation = body_snapshot.emits_radiation;
    }

    return true;
}

bool World::is_active(EntityId id) const {
    for (auto temp_id : active_ids) {
        if (id == temp_id) return true;
    }
    return false;
}

const svec<EntityId>& World::active_entity_ids() const { return this->active_ids; }
const svec<EntityId>& World::all_entity_ids() const { return this->all_ids; }

Body* World::body(EntityId id) {
    auto it = bodies.find(id);
    if (it == bodies.end()) return nullptr;
    return it->second.get();
}
const Body* World::body(EntityId id) const {
    auto it = bodies.find(id);
    if (it == bodies.end()) return nullptr;
    return it->second.get();
}

Celestial* World::celestial(EntityId id) {
    Body* ptr = body(id);
    if (ptr == nullptr) return nullptr;
    return dynamic_cast<Celestial*>(ptr);
}
const Celestial* World::celestial(EntityId id) const {
    const Body* ptr = body(id);
    if (ptr == nullptr) return nullptr;
    return dynamic_cast<const Celestial*>(ptr);
}

Satellite* World::satellite(EntityId id) {
    Body* ptr = body(id);
    if (ptr == nullptr) return nullptr;
    return dynamic_cast<Satellite*>(ptr);
}
const Satellite* World::satellite(EntityId id) const {
    const Body* ptr = body(id);
    if (ptr == nullptr) return nullptr;
    return dynamic_cast<const Satellite*>(ptr);
}

Station* World::station(EntityId id) {
    Body* ptr = body(id);
    if (ptr == nullptr) return nullptr;
    return dynamic_cast<Station*>(ptr);
}
const Station* World::station(EntityId id) const {
    const Body* ptr = body(id);
    if (ptr == nullptr) return nullptr;
    return dynamic_cast<const Station*>(ptr);
}

EntityId World::allocate_id() {
    active_ids.push_back(next_id);
    all_ids.push_back(next_id);
    return next_id++;
}
EntityId World::insert_body(uptr<Body> body) {
    if (body == nullptr) return kInvalidEntityId;
    EntityId id = allocate_id();
    body->id = id;
    bodies.emplace(id, std::move(body));
    return id;
}

// TODO: add overloads with data
EntityId World::spawn_celestial() { return insert_body(std::make_unique<Celestial>()); }
EntityId World::spawn_satellite() { return insert_body(std::make_unique<Satellite>()); }
EntityId World::spawn_station() { return insert_body(std::make_unique<Station>()); }

EntityId World::insert_celestial(uptr<Celestial> cel) {
    return insert_body(std::move(cel));
}
EntityId World::insert_satellite(uptr<Satellite> sat) {
    return insert_body(std::move(sat));
}
EntityId World::insert_station(uptr<Station> stat) {
    return insert_body(std::move(stat));
}

void World::insert_satellites(svec<uptr<Satellite>> sats) {
    for (auto& sat : sats) {
        insert_satellite(std::move(sat));
    }
}

i32 World::num_celestials() const {
    i32 count = 0;
    for (EntityId id : all_ids) {
        if (is_celestial(id)) count++;
    }
    return count;
}
i32 World::num_satellites() const {
    i32 count = 0;
    for (EntityId id : all_ids) {
        if (is_satellite(id)) count++;
    }
    return count;
}
i32 World::num_stations() const {
    i32 count = 0;
    for (EntityId id : all_ids) {
        if (is_station(id)) count++;
    }
    return count;
}

i32 World::num_active_celestials() const {
    i32 count = 0;
    for (EntityId id : active_ids) {
        if (is_celestial(id)) count++;
    }
    return count;
}

i32 World::num_active_satellites() const {
    i32 count = 0;
    for (EntityId id : active_ids) {
        if (is_satellite(id)) count++;
    }
    return count;
}

i32 World::num_active_stations() const {
    i32 count = 0;
    for (EntityId id : active_ids) {
        if (is_station(id)) count++;
    }
    return count;
}

i32 World::num_inactive_celestials() const {
    i32 count = 0;
    for (EntityId id : inactive_ids) {
        if (is_celestial(id)) count++;
    }
    return count;
}
i32 World::num_inactive_satellites() const {
    i32 count = 0;
    for (EntityId id : inactive_ids) {
        if (is_satellite(id)) count++;
    }
    return count;
}
i32 World::num_inactive_stations() const {
    i32 count = 0;
    for (EntityId id : inactive_ids) {
        if (is_station(id)) count++;
    }
    return count;
}

bool World::is_celestial(EntityId id) const {
    const auto ptr = celestial(id);
    if (ptr == nullptr) return false;
    if (ptr->body_type != BodyType::celestial) return false;
    return true;
}
bool World::is_satellite(EntityId id) const {
    const auto ptr = satellite(id);
    if (ptr == nullptr) return false;
    if (ptr->body_type != BodyType::satellite) return false;
    return true;
}
bool World::is_station(EntityId id) const {
    const auto ptr = station(id);
    if (ptr == nullptr) return false;
    if (ptr->body_type != BodyType::station) return false;
    return true;
}

svec<EntityId> World::celestial_ids(BodyFilterMode mode) const {
    svec<EntityId> ids;
    switch (mode) {
    case BodyFilterMode::all: {
        ids.reserve(num_celestials());
        for (EntityId id : all_ids) {
            if (is_celestial(id)) ids.push_back(id);
        }
    } break;
    case BodyFilterMode::active: {
        ids.reserve(num_active_celestials());
        for (EntityId id : active_ids) {
            if (is_celestial(id)) ids.push_back(id);
        }
    } break;
    case BodyFilterMode::inactive: {
        ids.reserve(num_inactive_celestials());
        for (EntityId id : inactive_ids) {
            if (is_celestial(id)) ids.push_back(id);
        }
    } break;
    }
    return ids;
}

svec<EntityId> World::satellite_ids(BodyFilterMode mode) const {
    svec<EntityId> ids;
    switch (mode) {
    case BodyFilterMode::all: {
        ids.reserve(num_satellites());
        for (EntityId id : all_ids) {
            if (is_satellite(id)) ids.push_back(id);
        }
    } break;
    case BodyFilterMode::active: {
        ids.reserve(num_active_satellites());
        for (EntityId id : active_ids) {
            if (is_satellite(id)) ids.push_back(id);
        }
    } break;
    case BodyFilterMode::inactive: {
        ids.reserve(num_inactive_satellites());
        for (EntityId id : inactive_ids) {
            if (is_satellite(id)) ids.push_back(id);
        }
    } break;
    }
    return ids;
}

svec<EntityId> World::station_ids(BodyFilterMode mode) const {
    svec<EntityId> ids;
    switch (mode) {
    case BodyFilterMode::all: {
        ids.reserve(num_stations());
        for (EntityId id : all_ids) {
            if (is_station(id)) ids.push_back(id);
        }
    } break;
    case BodyFilterMode::active: {
        ids.reserve(num_active_stations());
        for (EntityId id : active_ids) {
            if (is_station(id)) ids.push_back(id);
        }
    } break;
    case BodyFilterMode::inactive: {
        ids.reserve(num_inactive_stations());
        for (EntityId id : inactive_ids) {
            if (is_station(id)) ids.push_back(id);
        }
    } break;
    }
    return ids;
}

svec<EntityId> World::active_celestial_ids() const {
    svec<EntityId> ids;
    ids.reserve(num_active_celestials());
    for (EntityId id : active_ids) {
        if (is_celestial(id)) ids.push_back(id);
    }
    return ids;
}

svec<EntityId> World::active_satellite_ids() const {
    svec<EntityId> ids;
    ids.reserve(num_active_satellites());
    for (EntityId id : active_ids) {
        if (is_satellite(id)) ids.push_back(id);
    }
    return ids;
}

svec<EntityId> World::active_station_ids() const {
    svec<EntityId> ids;
    ids.reserve(num_active_stations());
    for (EntityId id : active_ids) {
        if (is_station(id)) ids.push_back(id);
    }
    return ids;
}

svec<EntityId> World::all_celestial_ids() const {
    svec<EntityId> ids;
    ids.reserve(num_celestials());
    for (EntityId id : all_ids) {
        if (is_celestial(id)) ids.push_back(id);
    }
    return ids;
}
svec<EntityId> World::all_satellite_ids() const {
    svec<EntityId> ids;
    ids.reserve(num_satellites());
    for (EntityId id : all_ids) {
        if (is_satellite(id)) ids.push_back(id);
    }
    return ids;
}
svec<EntityId> World::all_station_ids() const {
    svec<EntityId> ids;
    ids.reserve(num_stations());
    for (EntityId id : all_ids) {
        if (is_station(id)) ids.push_back(id);
    }
    return ids;
}
svec<EntityId> World::inactive_celestial_ids() const {
    svec<EntityId> ids;
    ids.reserve(num_inactive_celestials());
    for (EntityId id : inactive_ids) {
        if (is_celestial(id)) ids.push_back(id);
    }
    return ids;
}
svec<EntityId> World::inactive_satellite_ids() const {
    svec<EntityId> ids;
    ids.reserve(num_inactive_satellites());
    for (EntityId id : inactive_ids) {
        if (is_satellite(id)) ids.push_back(id);
    }
    return ids;
}
svec<EntityId> World::inactive_station_ids() const {
    svec<EntityId> ids;
    ids.reserve(num_inactive_stations());
    for (EntityId id : inactive_ids) {
        if (is_station(id)) ids.push_back(id);
    }
    return ids;
}

bool World::emits_gravity(EntityId id) const {
    const auto ptr = body(id);
    if (ptr == nullptr) return false;
    return ptr->emits_gravity;
}
bool World::emits_radiation(EntityId id) const {
    const auto ptr = body(id);
    if (ptr == nullptr) return false;
    return ptr->emits_radiation;
}

vec3d World::gravity_accel_from(EntityId target_id, EntityId source_id) const {
    const auto target = body(target_id);
    const auto source = celestial(source_id);
    if (target == nullptr || source == nullptr) return vec3d0;

    return gravity_accel_from(
        target_id,
        target->x_tr,
        source_id,
        source->x_tr,
        source->x_att
    );
}

vec3d World::gravity_accel_from(
    EntityId target_id,
    const StateTr& x_target,
    EntityId source_id
) const {
    const auto source = celestial(source_id);
    if (source == nullptr) return vec3d0;

    return gravity_accel_from(
        target_id,
        x_target,
        source_id,
        source->x_tr,
        source->x_att
    );
}

vec3d World::gravity_accel_from(
    EntityId target_id,
    const StateTr& x_target,
    EntityId source_id,
    const StateTr& x_source
) const {
    const auto source = celestial(source_id);
    if (source == nullptr) return vec3d0;

    return gravity_accel_from(target_id, x_target, source_id, x_source, source->x_att);
}

vec3d World::gravity_accel_from(
    EntityId target_id,
    const StateTr& x_tr_target,
    EntityId source_id,
    const StateTr& x_tr_source,
    const StateAtt& x_att_source
) const {
    vec3d a = vec3d0;

    if (!this->emits_gravity(source_id)) return vec3d0;
    if (target_id == source_id) return vec3d0;
    const auto target = body(target_id);
    const auto source = celestial(source_id);
    if (target == nullptr || source == nullptr) return vec3d0;

    vec3d r_rel_inertial = x_tr_target.r - x_tr_source.r;
    switch (source->gravity_model) {
    case GravityModel::spherical_harmonics: {
        vec4d q_BN = x_att_source.q; // [BN]: N -> B
        vec4d q_NB = ep_conj(q_BN);  // [NB]: B -> N
        // Body fixed relative position for spherical harmonic gravity perturbations
        vec3d r_rel_body = ep_rotate_fast_passive(q_BN, r_rel_inertial);
        vec3d a_body = accel_gravity_spherical_harmonics(
            r_rel_body,
            source->mu,
            source->ref_radius,
            source->degree,
            source->order,
            source->C,
            source->S
        ); // Acceleration in body fixed frame
        a = ep_rotate_fast_passive(q_NB, a_body);
        break;
    }
    case GravityModel::zonal: {
        vec4d q_BN = x_att_source.q; // [BN]: N -> B
        vec4d q_NB = ep_conj(q_BN);  // [NB]: B -> N
        // Body fixed relative position for zonal gravity perturbations
        vec3d r_rel_body = ep_rotate_fast_passive(q_BN, r_rel_inertial);
        vec3d a_body = accel_gravity_zonal(
            r_rel_body,
            source->mu,
            source->ref_radius,
            source->degree,
            source->J
        );                                        // Acceleration in body fixed frame
        a = ep_rotate_fast_passive(q_NB, a_body); // rotate back to inertial
        break;
    }
    case GravityModel::pointmass: {
        a = accel_gravity_pointmass(r_rel_inertial, source->mu);
        break;
    }
    default:
    }

    return a;
}

// TODO: optimize two way accelerations somehow
vec3d World::gravity_accel_on(EntityId target_id) const {
    vec3d a = vec3d0;
    const Body* target = body(target_id);
    if (target == nullptr) return a;

    for (const auto& [source_id, source] : bodies) {
        if (target_id == source_id) continue;
        if (!source->emits_gravity) continue;
        a += gravity_accel_from(target_id, source_id);
    }

    return a;
}

vec3d World::gravity_accel_on(EntityId target_id, const StateTr& x_target) const {
    vec3d a = vec3d0;
    const Body* target = body(target_id);
    if (target == nullptr) return a;

    for (const auto& [source_id, source] : bodies) {
        if (target_id == source_id) continue;
        if (!source->emits_gravity) continue;
        a += gravity_accel_from(target_id, x_target, source_id);
    }

    return a;
}

vec3d World::stat_r_inertial(EntityId station_id) const {
    // gets absolute inertial position
    // (anchored stations not propagated, free stations return internal state)
    vec3d r_inertial = vec3d0;
    const Station* stat = station(station_id);
    if (stat == nullptr) return r_inertial;
    if (!stat->anchored) return stat->x_tr.r;
    if (stat->anchor_id == kInvalidEntityId) return r_inertial;

    const Body* anchor = body(stat->anchor_id);
    if (anchor == nullptr) return r_inertial;

    // r_inertial = anchor->x_tr.r + ep_to_dcm<f64>(anchor->x_att.q).transpose() *
    // stat->r_body;
    vec4d q_NB = ep_conj(anchor->x_att.q);
    r_inertial = anchor->x_tr.r + ep_rotate_fast_passive(q_NB, stat->r_body_BCBF);
    return r_inertial;
}

vec3d World::stat_v_inertial(EntityId station_id) const {
    // inertial velocity of station
    vec3d v_inertial = vec3d0;
    const Station* stat = station(station_id);
    if (stat == nullptr) return v_inertial;
    if (!stat->anchored) return stat->x_tr.v;
    if (stat->anchor_id == kInvalidEntityId) return v_inertial;

    const Body* anchor = body(stat->anchor_id);
    if (anchor == nullptr) return v_inertial;

    vec4d q_NB = ep_conj(anchor->x_att.q);
    vec3d r_offset_inertial = ep_rotate_fast_passive(q_NB, stat->r_body_BCBF);
    vec3d w_inertial = ep_rotate_fast_passive(q_NB, anchor->x_att.w);
    v_inertial = anchor->x_tr.v + w_inertial.cross(r_offset_inertial);
    return v_inertial;
}

StateTr World::stat_x_tr_inertial(EntityId station_id) const {
    StateTr x_tr;
    const Station* stat = station(station_id);
    if (stat == nullptr) return x_tr;
    if (!stat->anchored) return stat->x_tr;
    if (stat->anchor_id == kInvalidEntityId) return x_tr;

    const Body* anchor = body(stat->anchor_id);
    if (anchor == nullptr) return x_tr;

    x_tr.r = stat_r_inertial(station_id);
    x_tr.v = stat_v_inertial(station_id);
    return x_tr;
}

vec4d World::stat_q_inertial(EntityId station_id) const {
    // passive attitude from sim inertial to station-local ENU or inertial orientation
    vec4d q_inertial = q_identity;
    const Station* stat = station(station_id);
    if (stat == nullptr) return q_inertial;
    if (!stat->anchored) return stat->x_att.q;
    if (stat->anchor_id == kInvalidEntityId) return q_inertial;

    const Body* anchor = body(stat->anchor_id);
    if (anchor == nullptr) return q_inertial;

    return stat_att_enu_from_detic(anchor->x_att, stat->llh_BCBF);
}

StateAtt World::stat_x_att_inertial(EntityId station_id) const {
    StateAtt x_att;
    const Station* stat = station(station_id);
    if (stat == nullptr) return x_att;
    if (!stat->anchored) return stat->x_att;
    if (stat->anchor_id == kInvalidEntityId) return x_att;

    const Body* anchor = body(stat->anchor_id);
    if (anchor == nullptr) return x_att;

    // anchored stations don't need x_att.w, TODO: maybe add later for sensor
    x_att.q = stat_q_inertial(station_id);
    return x_att;
}

bool World::set_stat_anchor_detic(
    EntityId station_id,
    EntityId anchor_id,
    const vec3d& llh,
    UAngle angle_in
) {
    Station* stat = station(station_id);
    const Celestial* cel = celestial(anchor_id);
    if (stat == nullptr || cel == nullptr) return false;
    if (cel->semimajor_axis <= 0.0 || cel->semiminor_axis <= 0.0) return false;

    vec3d r_body = stat_r_bcbf_from_detic(llh, *cel, angle_in);
    if (!r_body.allFinite()) return false;

    stat->anchored = true;
    stat->anchor_id = anchor_id;
    stat->r_body_BCBF = r_body;
    f64 lat = llh(0), lon = llh(1);
    if (angle_in != UAngle::radian) {
        lat = convert_angle(lat, angle_in, UAngle::radian);
        lon = convert_angle(lon, angle_in, UAngle::radian);
    }
    stat->llh_BCBF = vec3d{lat, lon, llh(2)};
    // don't propagate anchored stations, just assign
    stat->propagate_tr = false;
    stat->propagate_att = false;

    return true;
}

mat3d World::stat_rot_enu_from_body(EntityId station_id) const {
    const Station* stat = station(station_id);
    if (stat == nullptr) return mat3d1; // safe fallback
    return stat_rot_enu_from_detic(stat->llh_BCBF, UAngle::radian);
}

vec3d World::stat_rel_enu(EntityId station_id, EntityId target_id) const {
    const Station* stat = station(station_id);
    if (stat == nullptr) return vec3d0;
    if (!stat->anchored) return vec3d0;
    if (stat->anchor_id == kInvalidEntityId) return vec3d0;

    const Body* target = body(target_id);
    if (target == nullptr) return vec3d0;

    const Body* anchor = body(stat->anchor_id);
    if (anchor == nullptr) return vec3d0;

    // target relative to body, inertial then bcbf
    vec3d r_target_body_BCI = target->x_tr.r - anchor->x_tr.r;
    vec3d r_target_body_BCBF = ep_rotate_fast_passive(anchor->x_att.q, r_target_body_BCI);

    // target relative to station
    vec3d r_rel_bcbf = r_target_body_BCBF - stat->r_body_BCBF;
    mat3d R_ENU_BCBF = stat_rot_enu_from_body(station_id);

    return R_ENU_BCBF * r_rel_bcbf;
}

vec3d World::stat_rel_enu(EntityId station_id, const StateTr& x_tr_target) const {
    const Station* stat = station(station_id);
    if (stat == nullptr) return vec3d0;
    if (!stat->anchored) return vec3d0;
    if (stat->anchor_id == kInvalidEntityId) return vec3d0;

    const Body* anchor = body(stat->anchor_id);
    if (anchor == nullptr) return vec3d0;

    // target relative to body, inertial then bcbf
    vec3d r_target_body_I = x_tr_target.r - anchor->x_tr.r;
    vec3d r_target_body_BCBF = ep_rotate_fast_passive(anchor->x_att.q, r_target_body_I);

    // target relative to station
    vec3d r_rel_bcbf = r_target_body_BCBF - stat->r_body_BCBF;
    mat3d R_ENU_BCBF = stat_rot_enu_from_body(station_id);

    return R_ENU_BCBF * r_rel_bcbf;
}

vec3d World::body_z_inertial(EntityId body_id) const {
    vec3d z_inertial = vec3d0;
    const Body* body = this->body(body_id);
    if (body == nullptr) return z_inertial;

    vec4d q_NB = ep_conj(body->x_att.q);
    z_inertial = ep_rotate_fast_passive(q_NB, axis_z);
    return z_inertial;
}

vec3d World::body_w_inertial(EntityId body_id) const {
    vec3d w_inertial = vec3d0;
    const Body* body = this->body(body_id);
    if (body == nullptr) return w_inertial;

    vec4d q_NB = ep_conj(body->x_att.q);
    w_inertial = ep_rotate_fast_passive(q_NB, body->x_att.w);
    return w_inertial;
}

bool World::make_inactive(EntityId body_id) {
    for (i32 i = 0; i < active_ids.size(); ++i) {
        if (body_id == active_ids[i]) {
            active_ids.erase(active_ids.begin() + i);
            inactive_ids.push_back(body_id);
            return true;
        }
    }

    return false;
}

bool World::make_active(EntityId body_id) {
    for (i32 i = 0; i < inactive_ids.size(); ++i) {
        if (body_id == inactive_ids[i]) {
            inactive_ids.erase(inactive_ids.begin() + i);
            active_ids.push_back(body_id);
            return true;
        }
    }

    return false;
}
