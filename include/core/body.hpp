// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/entity.hpp"
#include "core/ephemeris_provider.hpp"
#include "core/observation_type.hpp"
#include "core/state.hpp"
#include "core/status.hpp"
#include "util/printing.hpp"
#include "util/vecdefs.hpp"

#include <print>
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
};
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

// Printing --------------------------------------------------------------------

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

using InstrumentId = u32;
constexpr InstrumentId kInvalidInstrumentId = 0;
struct PlatformInstrument {
    InstrumentId id = kInvalidInstrumentId;
    std::string name;
    ObservationType type = ObservationType::radec;
    matXd R;
    bool enabled = true;
};

struct InstrumentSuite {
    InstrumentId next_id = 1;
    umap<InstrumentId, PlatformInstrument> instruments;
    svec<InstrumentId> enabled_ids;
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

StatusCode measurement_covariance(
    const InstrumentSuite& suite,
    InstrumentId instrument_id,
    matXd& R
);

StatusCode measurement_covariance(
    const InstrumentSuite& suite,
    ObservationType type,
    matXd& R
);

StatusCode set_instrument(
    InstrumentSuite& suite,
    const PlatformInstrument& instrument
);

StatusCode add_instrument(
    InstrumentSuite& suite,
    const PlatformInstrument& instrument,
    InstrumentId& out_id
);

StatusCode get_instrument(
    const InstrumentSuite& suite,
    InstrumentId id,
    PlatformInstrument& out
);

svec<InstrumentId> enabled_instrument_ids(
    const InstrumentSuite& suite
);

StatusCode enable_instrument(
    InstrumentSuite& suite,
    InstrumentId id
);

StatusCode disable_instrument(
    InstrumentSuite& suite,
    InstrumentId id
);

// TODO: remove station specific instrument functions 
StatusCode station_measurement_covariance(
    const Station& station,
    InstrumentId instrument_id,
    matXd& R
);
StatusCode station_measurement_covariance(
    const Station& station,
    ObservationType type,
    matXd& R
);
StatusCode set_station_instrument(Station& station, const PlatformInstrument& instrument);
StatusCode add_station_instrument(
    Station& station,
    const PlatformInstrument& instrument,
    InstrumentId& out_id
);
StatusCode add_station_instrument(Station& station, const PlatformInstrument& instrument);
StatusCode get_station_instrument(
    const Station& station,
    PlatformInstrument& instrument,
    InstrumentId id
);

StatusCode add_radec_instrument(
    Station& station,
    const mat2d& R,
    std::string name = "Ra/Dec Instrument"
);
StatusCode add_radec_instrument(
    Station& station,
    const mat2d& R,
    InstrumentId& out_id,
    std::string name = "Ra/Dec Instrument"
);

StatusCode add_azel_instrument(
    Station& station,
    const mat2d& R,
    std::string name = "Az/El Instrument"
);
StatusCode add_azel_instrument(
    Station& station,
    const mat2d& R,
    InstrumentId& out_id,
    std::string name = "Az/El Instrument"
);

StatusCode add_range_instrument(
    Station& station,
    const matXd& R,
    std::string name = "Range Instrument"
);
StatusCode add_range_instrument(
    Station& station,
    const matXd& R,
    InstrumentId& out_id,
    std::string name = "Range Instrument"
);

StatusCode add_range_rate_instrument(
    Station& station,
    const matXd& R,
    std::string name = "Range-Rate Instrument"
);
StatusCode add_range_rate_instrument(
    Station& station,
    const matXd& R,
    InstrumentId& out_id,
    std::string name = "Range-Rate Instrument"
);

StatusCode add_pos_instrument(
    Station& station,
    const mat3d& R,
    std::string name = "Simulation Position Instrument"
);
StatusCode add_pos_instrument(
    Station& station,
    const mat3d& R,
    InstrumentId& out_id,
    std::string name = "Inertial Position Instrument"
);

StatusCode add_posvel_instrument(
    Station& station,
    const mat6d& R,
    std::string name = "Simulation State Instrument"
);
StatusCode add_posvel_instrument(
    Station& station,
    const mat6d& R,
    InstrumentId& out_id,
    std::string name = "Inertial State Instrument"
);

StatusCode add_rel_pos_instrument(
    Station& station,
    const mat3d& R,
    std::string name = "Relative Position Instrument"
);
StatusCode add_rel_pos_instrument(
    Station& station,
    const mat3d& R,
    InstrumentId& out_id,
    std::string name = "Relative Position Instrument"
);

StatusCode add_rel_posvel_instrument(
    Station& station,
    const mat6d& R,
    std::string name = "Relative State Instrument"
);
StatusCode add_rel_posvel_instrument(
    Station& station,
    const mat6d& R,
    InstrumentId& out_id,
    std::string name = "Relative State Instrument"
);

svec<InstrumentId> enabled_station_instrument_ids(const Station& station);

StatusCode enable_station_instrument(Station& station, InstrumentId instrument_id);
StatusCode disable_station_instrument(Station& station, InstrumentId instrument_id);

void print_station_instruments(const Station& station);

inline void print_mass_properties(const MassProperties& mp, const string& indent = "") {
    std::println("{}Mass: {}", indent, mp.mass);
    std::println("{}Inertia Tensor: I = {}", indent, vec_string(mp.I));
    std::println("{}Inverse Inertia Tensor: I_inv = {}", indent, vec_string(mp.I_inv));
    std::println("{}Principle Axes: {}", indent, mp.principal_axes);
    std::println("{}Active: {}", indent, mp.active);
}

inline void print_body_common(const Body& body, const string& indent = "") {
    std::println("{}ID: {}", indent, body.id);
    std::println("{}Name: {}", indent, body.name);
    std::println("{}Body Type: {}", indent, body_type_str(body.body_type));
    std::println("{}Propagate Translation: {}", indent, body.propagate_tr);
    std::println("{}Propagate Attitude: {}", indent, body.propagate_att);
    std::println("{}Emits Gravity: {}", indent, body.emits_gravity);
    std::println("{}Emits Radiation: {}", indent, body.emits_radiation);
    std::println("{}Has Atmosphere: {}", indent, body.has_atmosphere);
    print_state_tr(body.x_tr, indent);
    print_state_att(body.x_att, indent);
}

inline void print_celestial(const Celestial& cel, const string& indent = "") {
    const string indent2 = indent + "    ";

    std::println("{}--- Celestial ID: {}", indent, cel.id);
    print_body_common(cel, indent2);
    std::println("{}Gravity Model: {}", indent2, gravity_model_str(cel.gravity_model));
    std::println("{}mu: {}", indent2, cel.mu);
    std::println("{}Degree: {}", indent2, cel.degree);
    std::println("{}Order: {}", indent2, cel.order);
    std::println("{}J: {}", indent2, vec_string(cel.J));
    std::println("{}C Dim: ({}, {})", indent2, cel.C.rows(), cel.C.cols());
    std::println("{}S Dim: ({}, {})", indent2, cel.S.rows(), cel.S.cols());
    std::println(
        "{}Attitude Model: {}",
        indent2,
        celestial_attitude_model_str(cel.attitude_model)
    );
    std::println(
        "{}Radiation Model: {}",
        indent2,
        radiation_model_str(cel.radiation_model)
    );
    std::println("{}Reference Radius: {}", indent2, cel.ref_radius);
    std::println("{}Mean Radius: {}", indent2, cel.mean_radius);
    std::println("{}Semimajor Axis: {}", indent2, cel.semimajor_axis);
    std::println("{}Semiminor Axis: {}", indent2, cel.semiminor_axis);
    std::println("{}Eccentricity: {}", indent2, cel.eccentricity);
    std::println("{}Flattening: {}", indent2, cel.flattening);
    std::println("{}Spin Rate: {}", indent2, cel.spin_rate());
}

inline void print_satellite(const Satellite& sat, const string& indent = "") {
    const string indent2 = indent + "    ";

    std::println("{}--- Satellite ID: {}", indent, sat.id);
    print_body_common(sat, indent2);
    print_mass_properties(sat.mass_properties, indent2);
}

inline void print_station_instrument(
    const PlatformInstrument& instrument,
    const string& indent = ""
) {
    std::println("{}--- Instrument ID: {}", indent, instrument.id);
    std::println("{}Name: {}", indent, instrument.name);
    std::println("{}Instrument Type: {}", indent, observation_type_str(instrument.type));
    std::println("{}Enabled: {}", indent, instrument.enabled);
    std::println(
        "{}Measurement Covariance Dim: ({}, {})",
        indent,
        instrument.R.rows(),
        instrument.R.cols()
    );
    std::println("{}Measurement Covariance: R = {}", indent, vec_string(instrument.R));
}

inline void print_station(const Station& stat, const string& indent = "") {
    const string indent2 = indent + "    ";
    const string indent3 = indent2 + "    ";

    std::println("{}--- Station ID: {}", indent, stat.id);
    print_body_common(stat, indent2);
    std::println("{}Anchored: {}", indent2, stat.anchored);
    std::println("{}Anchor ID: {}", indent2, stat.anchor_id);
    std::println("{}r_body_BCBF: {}", indent2, vec_string(stat.r_body_BCBF));
    std::println("{}llh_BCBF: {}", indent2, vec_string(stat.llh_BCBF));

    if (!stat.anchored) {
        print_mass_properties(stat.mass_properties, indent2);
    }

    std::println(
        "{}Instruments Count: {}",
        indent2,
        stat.instrument_suite.instruments.size()
    );
    std::println(
        "{}Enabled Instrument IDs: {}",
        indent2,
        vec_string(stat.instrument_suite.enabled_ids)
    );
    for (const auto& [id, instrument] : stat.instrument_suite.instruments) {
        print_station_instrument(instrument, indent3);
    }
}


