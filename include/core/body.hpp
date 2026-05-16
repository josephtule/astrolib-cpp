#pragma once

#include "core/entity.hpp"
#include "core/estimation_common.hpp"
#include "core/observation_type.hpp"
#include "core/state.hpp"
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
    bool has_atmosphere = false; // currently does nothing
};

enum struct GravityModel { pointmass, zonal, spherical_harmonics };
inline std::string gravity_model_name(GravityModel model) {
    switch (model) {
    case GravityModel::pointmass: return "pointmass";
    case GravityModel::zonal: return "zonal";
    case GravityModel::spherical_harmonics: return "spherical harmonics";
    }
    return "unknown";
};
enum struct RadiationModel { none, isotropic };
enum struct CelestialAttitudeModel { fixed, simple_spin, provider };
struct Celestial : public Body {
    // Gravity
    GravityModel gravity_model = GravityModel::pointmass;
    f64 mu = 0.0;
    i32 degree = 0, order = 0; // spherical harmonics degree (n) and order (m)
    vec7d J = vec7d0;          // zonal coefs
    matXd C, S;                // sph harmonic coefs

    // Attitude/Orientation
    CelestialAttitudeModel attitude_model = CelestialAttitudeModel::fixed;

    // Radiation
    RadiationModel radiation_model = RadiationModel::none;

    // Shape
    f64 mean_radius = 0.0;
    f64 semimajor_axis = 0.0;
    f64 semiminor_axis = 0.0;
    f64 eccentricity = 0.0, flattening = 0.0;

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

struct MassProperties {
    f64 mass = 0.0;
    mat3d I = mat3d1;
    mat3d I_inv = mat3d1;
    bool principal_axes = true;
    bool active = false;
};
struct Satellite : public Body {
    MassProperties mass_properties;

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

using InstrumentId = u32;
constexpr InstrumentId kInvalidInstrumentId = 0;
struct StationInstrument {
    InstrumentId id = kInvalidInstrumentId;
    std::string name;
    ObservationType type = ObservationType::radec;
    matXd R;
    bool enabled = true;
};

struct Station : public Body {
    bool anchored = true;
    EntityId anchor_id = kInvalidEntityId;
    vec3d r_body_BCBF = vec3d0; // Position of station relative to anchor in bcbf
    vec3d llh_BCBF = vec3d0;    // Planetodetic coordinates
    // [lat, lon, h] - [rad, rad, sim units]
    InstrumentId next_instrument_id = 1;
    umap<u32, StationInstrument> instruments;

    MassProperties mass_properties;

    // Constructor(s)
    Station() {
        body_type = BodyType::station;
        propagate_att = false;
        propagate_tr = false;
    }
};

ODStatus station_instrument_covariance(const Station& station, InstrumentId instrument_id, matXd& R);
ODStatus stat_meas_cov(const Station& station, ObservationType type, matXd& R);
ODStatus station_instrument_covariance(Station& station, const StationInstrument& instrument);
ODStatus add_station_instrument(
    Station& station,
    const StationInstrument& instrument,
    InstrumentId& out_id
);
ODStatus add_station_instrument(Station& station, const StationInstrument& instrument);
ODStatus get_station_instrument(
    const Station& station,
    const StationInstrument* instrument,
    InstrumentId id
);

ODStatus add_radec_instrument(
    Station& station,
    const mat2d& R,
    std::string name = "Ra/Dec Instrument"
);
ODStatus add_radec_instrument(
    Station& station,
    const mat2d& R,
    InstrumentId& out_id,
    std::string name = "Ra/Dec Instrument"
);

ODStatus add_azel_instrument(
    Station& station,
    const mat2d& R,
    std::string name = "Az/El Instrument"
);
ODStatus add_azel_instrument(
    Station& station,
    const mat2d& R,
    InstrumentId& out_id,
    std::string name = "Az/El Instrument"
);

ODStatus add_range_instrument(
    Station& station,
    const matXd& R,
    std::string name = "Range Instrument"
);
ODStatus add_range_instrument(
    Station& station,
    const matXd& R,
    InstrumentId& out_id,
    std::string name = "Range Instrument"
);

ODStatus add_range_rate_instrument(
    Station& station,
    const matXd& R,
    std::string name = "Range-Rate Instrument"
);
ODStatus add_range_rate_instrument(
    Station& station,
    const matXd& R,
    InstrumentId& out_id,
    std::string name = "Range-Rate Instrument"
);

ODStatus add_pos_instrument(
    Station& station,
    const mat3d& R,
    std::string name = "Simulation Position Instrument"
);
ODStatus add_pos_instrument(
    Station& station,
    const mat3d& R,
    InstrumentId& out_id,
    std::string name = "Inertial Position Instrument"
);

ODStatus add_posvel_instrument(
    Station& station,
    const mat6d& R,
    std::string name = "Simulation State Instrument"
);
ODStatus add_posvel_instrument(
    Station& station,
    const mat6d& R,
    InstrumentId& out_id,
    std::string name = "Inertial State Instrument"
);

ODStatus add_rel_pos_instrument(
    Station& station,
    const mat3d& R,
    std::string name = "Relative Position Instrument"
);
ODStatus add_rel_pos_instrument(
    Station& station,
    const mat3d& R,
    InstrumentId& out_id,
    std::string name = "Relative Position Instrument"
);

ODStatus add_rel_posvel_instrument(
    Station& station,
    const mat6d& R,
    std::string name = "Relative State Instrument"
);
ODStatus add_rel_posvel_instrument(
    Station& station,
    const mat6d& R,
    InstrumentId& out_id,
    std::string name = "Relative State Instrument"
);
