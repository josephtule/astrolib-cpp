#pragma once

#include "core/body.hpp"
#include "core/estimation_common.hpp"
#include "core/ingest.hpp"
#include "core/integrator.hpp"
#include "core/observation_type.hpp"
#include "core/od_dynamics.hpp"
#include "core/orbital_elements.hpp"
#include "core/state.hpp"
#include "core/time.hpp"
#include "core/transform.hpp"
#include "core/world.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

struct ScenarioSchemaConfig {
    string name = "astrolib.scenario";
    i32 version = 1;
};
struct ScenarioMetadataConfig {
    string name = "earth_moon_sat_demo";
    u32 rng_seed = 12345;
};

struct ScenarioUnitConfig {
    string length = "km";
    string time = "s";
    string angle = "rad";
    string mass = "kg";
};

enum struct DateType { cal, jd, mjd };
inline string date_type_str(const DateType& type) {
    switch (type) {
    case DateType::cal: return "calendar";
    case DateType::jd: return "julian date";
    case DateType::mjd: return "modified julian date";
    }

    return "unknown";
}
struct ScenarioTimeConfig {
    f64 t0 = 0.0;

    CalendarTime cal{};
    JulianDate jd{};
    ModifiedJulianDate mjd{};

    TimeScale time_scale = TimeScale::utc;
    DateType date_type = DateType::jd;
    CalendarPrintStyle cal_style = CalendarPrintStyle::separate;
};

struct ScenarioFramesConfig {
    string simulation_inertial = "I";
    string attitude_quaternion = "passive_N_to_B";
    string station_local = "ENU";
};

struct ScenarioGravityProviderConfig {
    string id;
    string format;
    // optional for zonal and spherical harmonics
    string filepath;
    i32 lineskips = 0;
    bool normalized = true;
};

enum struct GravityCoefficientSource {
    none,
    provider,
    direct_zonal,
    direct_spherical_harmonics
};
inline string gravity_coefficient_source_str(const GravityCoefficientSource& source) {
    switch (source) {
    case GravityCoefficientSource::none: return "none";
    case GravityCoefficientSource::provider: return "provider";
    case GravityCoefficientSource::direct_zonal: return "direct zonal";
    case GravityCoefficientSource::direct_spherical_harmonics:
        return "direct spherical harmonics";
    }

    return "unknown";
}

struct ScenarioGravityConfig {
    GravityModel model = GravityModel::pointmass;
    f64 mu = 0.0;

    // optional for zonal and spherical harmonics
    f64 radius = 0.0; // NOTE: match to provider if used
    i32 degree = 0;
    i32 order = 0;
    string coefficients;
    vecXd J;
    GravityCoefficientSource coefficient_source = GravityCoefficientSource::none;
};

struct ScenarioPropagationConfig {
    bool translation = false;
    bool attitude = false;
};

enum struct StateTrInputType { pos_vel, classical };
inline string state_tr_type_str(const StateTrInputType& type) {
    switch (type) {
    case StateTrInputType::pos_vel: return "position + velocity";
    case StateTrInputType::classical: return "classical orbital elements";
    }

    return "Invalid StateTr input type";
}
struct ScenarioStateTrConfig {
    vec3d r = vec3d0;
    vec3d v = vec3d0;

    // TODO: add state input type (orbital elements)
    StateTrInputType input_type = StateTrInputType::pos_vel;
    OEClassical coes;
    string central; // central body for orbital elements

    ULength units_length = ULength::kilometer;
    UAngle units_angle = UAngle::degree;
};

inline string state_att_type_str(const AttitudeType& type) {
    switch (type) {
    case AttitudeType::quaternion: return "quaternion";
    case AttitudeType::dcm: return "dcm";
    case AttitudeType::axis_angle: return "axis-angle";
    case AttitudeType::euler_angles: return "euler angles";
    case AttitudeType::crp: return "classical rodriguez parameters";
    case AttitudeType::mrp: return "modified rodriguez parameters";
    }

    return "unknown";
}
struct ScenarioStateAttConfig {
    AttitudeType input_type = AttitudeType::quaternion;
    // quaternion
    vec4d q = q_identity;
    // axis angle
    vec3d axis = vec3d0; // use for crp and mrp as well
    f64 angle = 0.0;
    // euler angles
    vec3d angles{0.0, 0.0, 0.0};
    std::array<i32, 3> sequence{3, 2, 1};
    vec3d w = vec3d0;
    // dcm
    mat3d dcm = mat3d1;

    UAngle units_angle = UAngle::radian;
};

struct ScenarioCelestialModelConfig {
    string id;

    ScenarioGravityConfig gravity_model;

    f64 mean_radius = 0.0;
    f64 semimajor_axis = 0.0;
    f64 semiminor_axis = 0.0;
    f64 eccentricity = 0.0;
    f64 flattening = 0.0;
};
struct ScenarioCelestialConfig {
    string id;
    string name;
    ScenarioCelestialModelConfig model;
    ScenarioStateTrConfig x_tr;
    ScenarioStateAttConfig x_att;
    bool has_attitude_model = false;
    CelestialAttitudeModel attitude_model = CelestialAttitudeModel::fixed;
    // RadiationModel radition_model = RadiationModel::none;
    ScenarioPropagationConfig propagation;
};

struct ScenarioMassPropertiesConfig {
    f64 mass = 0.0;

    string type = "diag";
    bool principle_axes = true;
    mat3d inertia = mat3d1;
};

struct ScenarioSatelliteConfig {
    string id;
    string name;
    ScenarioStateTrConfig x_tr;
    ScenarioStateAttConfig x_att;
    ScenarioPropagationConfig propagation;
    ScenarioMassPropertiesConfig mass_properties;
};

struct ScenarioCovarianceConfig {
    string type = "diagonal";
    matXd covariance;
};

struct ScenarioInstrumentConfig {
    string id;
    ObservationType type = ObservationType::radec;
    bool enabled = true;
    ScenarioCovarianceConfig covariance_cfg;
};

struct ScenarioStationConfig {
    string id;
    string name;

    ScenarioPropagationConfig propagation;
    bool anchored = true;
    string anchor;

    // optional dependent in anchoring
    // anchored
    string coordinate_type = "detic_llh";
    vec3d llh_BCBF = vec3d0;
    vec3d r_body = vec3d0;
    string local_frame;
    // unanchored
    ScenarioStateTrConfig x_tr;
    ScenarioStateAttConfig x_att;
    ScenarioMassPropertiesConfig mass_properties;

    UAngle units_angle = UAngle::radian;
    ULength units_length = ULength::kilometer;

    svec<ScenarioInstrumentConfig> instruments;
};

struct ScenarioWorldStepperConfig {
    IntegratorType integrator_tr = IntegratorType::rk4;
    IntegratorType integrator_att = IntegratorType::rk4;

    u32 substeps = 1;
    u32 ticks = 1;
    f64 time_scale = 1.0;
};

struct ScenarioGraphicsConfig {
    i32 target_fps = 60;
    MyColor background_color = {30, 30, 30, 255};
    bool draw_inertial_axes = true;
};

struct ScenarioConfig {
    ScenarioSchemaConfig schema;
    ScenarioMetadataConfig metadata;

    ScenarioTimeConfig time;

    svec<ScenarioGravityProviderConfig> gravity_providers;
    // svec<ScenarioEphemerisProviderConfig> ephemeris_providers;

    svec<ScenarioCelestialConfig> celestials;
    svec<ScenarioSatelliteConfig> satellites;
    svec<ScenarioStationConfig> stations;
    // svec<ScenarioInstrumentConfig> instrument_templates;

    ScenarioWorldStepperConfig world_stepper;
    ScenarioGraphicsConfig graphics_settings;
};

struct ScenarioBuildResult {
    umap<string, EntityId> body_ids;
    umap<string, EntityId> celestial_ids;
    umap<string, EntityId> satellite_ids;
    umap<string, EntityId> station_ids;
};

StatusCode load_scenario_json(const std::string& filepath, ScenarioConfig& out);
StatusCode validate_scenario_config(const ScenarioConfig& cfg);
StatusCode build_world_from_scenario_config(
    const ScenarioConfig& cfg,
    World& world,
    ScenarioBuildResult& result
);
