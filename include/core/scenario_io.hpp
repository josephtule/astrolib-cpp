#pragma once

#include "core/body.hpp"
#include "core/estimation_common.hpp"
#include "core/integrator.hpp"
#include "core/observation_type.hpp"
#include "core/od_dynamics.hpp"
#include "core/state.hpp"
#include "core/transform.hpp"
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

struct ScenarioTimeConfig {
    f64 t0 = 0.0;
    string time_scale = "sim_seconds";
    string calendar_utc;
    bool has_calendar_utc = false;
};

struct ScenarioFramesConfig {
    string simulation_inertial = "I";
    string attitude_quaternion = "passive_N_to_B";
    string station_local = "ENU";
};

struct ScenarioGravityProviderConfig {
    string id;
    string type;
    // optional for zonal and spherical harmonics
    string filepath;
    bool normalized = true;
};

struct ScenarioGravityConfig {
    GravityModel model = GravityModel::pointmass;
    f64 mu = 0.0;

    // optional for zonal and spherical harmonics
    f64 radius = 0.0;
    i32 degree = 0;
    i32 order = 0;
    string coefficients;
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

struct ScenarioPropagationConfig {
    bool translation = false;
    bool attitude = false;
};

struct ScenarioStateTrConfig {
    vec3d r, v;

    ULength units_length = ULength::kilometer;
};

struct ScenarioStateAttConfig {
    AttitudeType type = AttitudeType::quaternion;
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

struct ScenarioCelestialConfig {
    string id;
    string name;
    ScenarioCelestialModelConfig model;
    ScenarioStateTrConfig x_tr;
    ScenarioStateAttConfig x_att;
    CelestialAttitudeModel attitude_model = CelestialAttitudeModel::fixed;
    // RadiationModel radition_model = RadiationModel::none;
    ScenarioPropagationConfig propagation;
    ScenarioGravityConfig gravity;
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

    svec<ScenarioGravityProviderConfig> gravity_providers;
    // svec<ScenarioEphemerisProviderConfig> ephemeris_providers;

    svec<ScenarioCelestialConfig> celestials;
    svec<ScenarioSatelliteConfig> satellites;
    svec<ScenarioStationConfig> stations;
    // svec<ScenarioInstrumentConfig> instrument_templates;

    ScenarioWorldStepperConfig world_stepper;
    ScenarioGraphicsConfig graphics_settings;
};

StatusCode load_scenario_json(const std::string& filepath, ScenarioConfig& out);
StatusCode validate_scenario_config(const ScenarioConfig& cfg);