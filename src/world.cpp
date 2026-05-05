#include "core/world.hpp"
#include "core/body.hpp"
#include "core/dynamics_translational.hpp"
#include "core/entity.hpp"
#include "core/state.hpp"
#include "core/transform.hpp"
#include "util/vecdefs.hpp"

#include <cstddef>
#include <memory>

f64 World::t_sim() const { return t_sim_; }

void World::reset_time(f64 t0) { t_sim_ = t0; }

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
    return next_id++;
}
EntityId World::insert_body(std::unique_ptr<Body> body) {
    if (body == nullptr) return kInvalidEntityId;
    EntityId id = allocate_id();
    body->id = id;
    bodies.emplace(id, std::move(body));
    return id;
}

EntityId World::spawn_celestial() { return insert_body(std::make_unique<Celestial>()); }
EntityId World::spawn_satellite() { return insert_body(std::make_unique<Satellite>()); }
EntityId World::spawn_station() { return insert_body(std::make_unique<Station>()); }

EntityId World::insert_celestial(std::unique_ptr<Celestial> cel) {
    return insert_body(std::move(cel));
}
EntityId World::insert_satellite(std::unique_ptr<Satellite> sat) {
    return insert_body(std::move(sat));
}
EntityId World::insert_station(std::unique_ptr<Station> stat) {
    return insert_body(std::move(stat));
}

bool World::is_celestial(EntityId id) const {
    const auto ptr = body(id);
    if (ptr == nullptr) return false;
    if (ptr->body_type != BodyType::celestial) return false;
    return true;
}
bool World::is_satellite(EntityId id) const {
    const auto ptr = body(id);
    if (ptr == nullptr) return false;
    if (ptr->body_type != BodyType::satellite) return false;
    return true;
}
bool World::is_station(EntityId id) const {
    const auto ptr = body(id);
    if (ptr == nullptr) return false;
    if (ptr->body_type != BodyType::station) return false;
    return true;
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
    vec3d a = vec3d::Zero();

    if (!emits_gravity(source_id)) return a;
    if (target_id == source_id) return a;
    const auto target = body(target_id);
    // const auto source = dynamic_cast<const Celestial*>(body(source_id));
    const auto source = celestial(source_id);
    if (target == nullptr || source == nullptr) return a;

    vec3d r_rel_inertial = target->x_tr.r - source->x_tr.r;
    switch (source->gravity_model) {
    case GravityModel::spherical_harmonics: {
        vec4d q_BN = source->x_att.q; // [BN]: N -> B
        vec4d q_NB = ep_conj(q_BN);   // [NB]: B -> N
        // Body fixed relative position for spherical harmonic gravity perturbations
        vec3d r_rel_body = ep_rotate_fast_passive(q_BN, r_rel_inertial);
        vec3d a_body = accel_gravity_spherical_harmonics(
            r_rel_body,
            source->mu,
            source->mean_radius,
            source->degree,
            source->order,
            source->C,
            source->S
        ); // Acceleration in body fixed frame
        a = ep_rotate_fast_passive(q_NB, a_body);
        break;
    }
    case GravityModel::zonal: {
        vec4d q_BN = source->x_att.q; // [BN]: N -> B
        vec4d q_NB = ep_conj(q_BN);   // [NB]: B -> N
        // Body fixed relative position for zonal gravity perturbations
        vec3d r_rel_body = ep_rotate_fast_passive(q_BN, r_rel_inertial);
        vec3d a_body = accel_gravity_zonal(
            r_rel_body,
            source->mu,
            source->mean_radius,
            source->degree,
            source->J
        ); // Acceleration in body fixed frame
        a = ep_rotate_fast_passive(q_NB, a_body);
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

vec3d World::gravity_accel_on(EntityId target_id) const {
    vec3d a = vec3d::Zero();
    const Body* target = body(target_id);
    if (target == nullptr) return a;

    for (const auto& [source_id, source] : bodies) {
        if (target_id == source_id) continue;
        if (!source->emits_gravity) continue;
        a += gravity_accel_from(target_id, source_id);
    }

    return a;
}

vec3d World::stat_r_inertial(EntityId station_id) const {
    vec3d r_inertial = vec3d::Zero();
    const Station* stat = station(station_id);
    if (stat == nullptr) return r_inertial;
    if (!stat->anchored) return r_inertial;
    if (stat->anchor_id == kInvalidEntityId) return r_inertial;

    const Body* anchor = body(stat->anchor_id);
    if (anchor == nullptr) return r_inertial;

    // r_inertial = anchor->x_tr.r + ep_to_dcm<f64>(anchor->x_att.q).transpose() *
    // stat->r_body;
    vec4d q_NB = ep_conj(anchor->x_att.q);
    r_inertial = anchor->x_tr.r + ep_rotate_fast_passive(q_NB, stat->r_body);
    return r_inertial;
}

vec3d World::stat_v_inertial(EntityId station_id) const {
    vec3d v_inertial = vec3d::Zero();
    const Station* stat = station(station_id);
    if (stat == nullptr) return v_inertial;
    if (!stat->anchored) return v_inertial;
    if (stat->anchor_id == kInvalidEntityId) return v_inertial;

    const Body* anchor = body(stat->anchor_id);
    if (anchor == nullptr) return v_inertial;

    vec4d q_NB = ep_conj(anchor->x_att.q);
    vec3d r_offset_inertial = ep_rotate_fast_passive(q_NB, stat->r_body);
    vec3d w_inertial = ep_rotate_fast_passive(q_NB, anchor->x_att.w);
    v_inertial = anchor->x_tr.v + w_inertial.cross(r_offset_inertial);
    return v_inertial;
}

StateTr World::stat_x_tr_inertial(EntityId station_id) const {
    StateTr x_tr;
    const Station* stat = station(station_id);
    if (stat == nullptr) return x_tr;
    if (!stat->anchored) return x_tr;
    if (stat->anchor_id == kInvalidEntityId) return x_tr;

    const Body* anchor = body(stat->anchor_id);
    if (anchor == nullptr) return x_tr;

    x_tr.r = stat_r_inertial(station_id);
    x_tr.v = stat_v_inertial(station_id);
    return x_tr;
}

vec3d World::body_z_inertial(EntityId body_id) const {
    vec3d z_inertial = vec3d::Zero();
    const Body* body = this->body(body_id);
    if (body == nullptr) return z_inertial;

    vec4d q_NB = ep_conj(body->x_att.q);
    z_inertial = ep_rotate_fast_passive(q_NB, axis_z);
    return z_inertial;
}

vec3d World::body_w_inertial(EntityId body_id) const {
    vec3d w_inertial = vec3d::Zero();
    const Body* body = this->body(body_id);
    if (body == nullptr) return w_inertial;

    vec4d q_NB = ep_conj(body->x_att.q);
    w_inertial = ep_rotate_fast_passive(q_NB, body->x_att.w);
    return w_inertial;
}
