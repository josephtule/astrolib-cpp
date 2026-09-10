// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/entity.hpp"
#include "core/ephemeris_provider.hpp"
#include "core/instrument.hpp"
#include "core/state.hpp"
#include "core/status.hpp"
#include "util/vecdefs.hpp"

#include <string>

enum struct BodyType { unknown = 0, celestial = 1, satellite = 2, station = 3 };

inline std::string body_type_str(BodyType type) {
    switch (type) {
    case BodyType::unknown: return "Unknown";
    case BodyType::celestial: return "Celestial";
    case BodyType::satellite: return "Satellite";
    case BodyType::station: return "Station";
    }
    return "Unknown";
}

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
    bool has_atmosphere = false; // currently does nothing

    BodyEphemerisProviders ephemeris_providers{};
};

enum struct GravityModel { pointmass = 0, zonal = 1, spherical_harmonics = 2 };

inline std::string gravity_model_str(GravityModel model) {
    switch (model) {
    case GravityModel::pointmass: return "pointmass";
    case GravityModel::zonal: return "zonal";
    case GravityModel::spherical_harmonics: return "spherical harmonics";
    }
    return "unknown";
}

enum struct RadiationModel { none = 0, isotropic = 1 };

enum struct CelestialAttitudeModel {
    fixed = 0,
    simple_spin = 1,
    provider = 2
}; // TODO: add normal attitude propagation

struct Celestial : public Body {
    // Gravity
    GravityModel gravity_model = GravityModel::pointmass;
    f64 mu = 0.0;
    i32 degree = 0, order = 0; // spherical harmonics degree (n) and order (m)
    vec7d J = vec7d0;          // zonal coefs
    matXd C, S;                // sph harmonic coefs
    std::string gravity_provider;
    std::string gravity_provider_format;
    std::string gravity_provider_filepath;
    i32 gravity_provider_lineskips = 0;
    bool gravity_provider_normalized = true;

    // Attitude/Orientation
    CelestialAttitudeModel attitude_model = CelestialAttitudeModel::fixed;

    // Radiation
    RadiationModel radiation_model = RadiationModel::none;

    // Shape
    f64 ref_radius = 0.0; // used for gravity computations
    f64 semimajor_axis = 0.0;
    f64 semiminor_axis = 0.0;
    f64 mean_radius = 0.0;
    f64 eccentricity = 0.0;
    f64 flattening = 0.0;

    // Spin
    f64 spin_rate() const { return x_att.w.norm(); }

    // Constructor(s)
    Celestial() {
        body_type = BodyType::celestial;
        emits_gravity = true;
    }

    void set_spin_rate(f64 spin_rate_) {
        // in rad/s
        if (attitude_model == CelestialAttitudeModel::simple_spin) {
            x_att.w(2) = spin_rate_;
        }
    }
};

StatusCode set_celestial_ephemeris_providers(
    Celestial& celestial,
    BodyEphemerisProviders providers
);

inline string celestial_attitude_model_str(const CelestialAttitudeModel model) {
    switch (model) {
    case CelestialAttitudeModel::fixed: return "fixed";
    case CelestialAttitudeModel::simple_spin: return "simple spin";
    case CelestialAttitudeModel::provider: return "provider";
    }
    return "unknown";
}

inline string radiation_model_str(const RadiationModel model) {
    switch (model) {
    case RadiationModel::none: return "none";
    case RadiationModel::isotropic: return "isotropic";
    }
    return "unknown";
}

struct MassProperties {
    f64 mass = 0.0;
    mat3d I = mat3d1;
    mat3d I_inv = mat3d1;
    bool principal_axes = true;

    vec3d offset_body = vec3d0; // center for I (input)

    bool active = false;
};

struct Satellite : public Body {
    MassProperties mass_properties;
    InstrumentSuite instrument_suite;

    // Constructor(s)
    Satellite() {
        body_type = BodyType::satellite;
        mass_properties.active = true;
    }

    void set_I(const mat3d& I) {
        this->mass_properties.I = I;
        this->mass_properties.I_inv = I.inverse();
    }
};

struct Station : public Body {
    bool anchored = true;
    EntityId anchor_id = kInvalidEntityId;
    vec3d r_body_BCBF = vec3d0; // Position of station relative to anchor in bcbf
    vec3d llh_BCBF = vec3d0;    // Planetodetic coordinates
    // [lat, lon, h] - [rad, rad, sim units]

    MassProperties mass_properties;
    InstrumentSuite instrument_suite;

    // Constructor(s)
    Station() {
        body_type = BodyType::station;
        propagate_att = false;
        propagate_tr = false;
    }
};

StatusCode instrument_suite_from_body(Body& body, InstrumentSuite*& out);

StatusCode instrument_suite_from_body(const Body& body, const InstrumentSuite*& out);
