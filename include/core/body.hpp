#pragma once

#include "core/entity.hpp"
#include "core/state.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"
#include <string>

enum struct BodyType { unknown, celestial, satellite, station };

struct Body {
    EntityId id = kInvalidEntityId;
    std::string name;
    BodyType body_type = BodyType::unknown;

    StateTr x_tr;
    StateAtt x_att;

    virtual ~Body() = default;

    bool propagate_tr = true;
    bool propagate_att = false;

    bool emits_gravity = false;
    bool emits_radiation = false;
};

enum struct GravityModel { pointmass, zonal, spherical_harmonics };
enum struct RadiationModel { none, isotropic };

struct Celestial : public Body {
    // Gravity
    GravityModel gravity_model = GravityModel::pointmass;
    f64 mu = 0.0;
    i32 degree = 0, order = 0; // spherical harmonics degree (n) and order (m)
    vec7d J = vec7d::Zero();   // zonal coefs
    matXd C, S;                // sph harmonic coefs

    // Radiation
    RadiationModel radiation_model = RadiationModel::none;

    // Shape
    f64 mean_radius = 0.0;
    f64 semimajor_axis = 0.0;
    f64 semiminor_axis = 0.0;
    f64 eccentricity = 0.0, flattening = 0.0;

    // Attitude
    bool use_simple_spin = false;
    f64 spin_rate = 0.0;

    // Constructor(s)
    Celestial() {
        body_type = BodyType::celestial;
        emits_gravity = true;
    }

    void set_spin_rate(f64 spin_rate_) {
        // in rad/s
        if (use_simple_spin) {
            spin_rate = spin_rate_;
            x_att.w(2) = spin_rate;
        }
    }
};

struct Satellite : public Body {
    // Mass Properties
    f64 mass = 0.0;

    bool principal_axes = false;
    mat3d I = mat3d::Identity();
    mat3d I_inv = mat3d::Identity();

    // Constructor(s)
    Satellite() { body_type = BodyType::satellite; }
};

struct Station : public Body {
    bool anchored = true;
    EntityId anchor_id = kInvalidEntityId;
    vec3d r_body_BCBF = vec3d::Zero(); // Position of station relative to anchor in bcbf
    vec3d llh_BCBF = vec3d::Zero();    // Planetodetic coordinates
    // [lat, lon, h] - [rad, rad, sim units]

    // Constructor(s)
    Station() {
        body_type = BodyType::station;
        propagate_att = false;
        propagate_tr = false;
    }
};