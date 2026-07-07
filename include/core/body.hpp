#pragma once

#include "core/entity.hpp"
#include "core/estimation_common.hpp"
#include "core/observation_type.hpp"
#include "core/state.hpp"
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
    svec<InstrumentId> enabled_instrument_ids;

    MassProperties mass_properties;

    // Constructor(s)
    Station() {
        body_type = BodyType::station;
        propagate_att = false;
        propagate_tr = false;
    }
};

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
StatusCode set_station_instrument(Station& station, const StationInstrument& instrument);
StatusCode add_station_instrument(
    Station& station,
    const StationInstrument& instrument,
    InstrumentId& out_id
);
StatusCode add_station_instrument(Station& station, const StationInstrument& instrument);
StatusCode get_station_instrument(
    const Station& station,
    StationInstrument& instrument,
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