#pragma once

#include "core/estimation_common.hpp"
#include "core/state.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

using string = std::string;

struct ScenarioSchemaConfig {
    string name = "astrolib.scenario";
    i32 version = 1;
};
struct ScenarioMetadataConfig {
    string name = "earth_moon_sat_demo";
    u64 rng_seed = 12345;
    bool has_rng_seed = true;
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
    string path;
    bool normalized = true;
};

struct ScenarioGravityConfig {
    string model = "pointmass";
    f64 mu = 0.0;
    f64 radius = 0.0;
    i32 degree = 0;
    i32 order = 0;
    string coefficients;
};

struct ScenarioPropagationConfig {
    bool translation = false;
    bool attitude = false;
};

struct ScenarioStateAttConfig {
    string type = "quaternion"; // axis_angle, euler_angles, dcm, crp, mrp
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
    string model;
    StateTr x_tr;
    ScenarioStateAttConfig x_att;
    ScenarioPropagationConfig propagation;
    ScenarioGravityConfig gravity;
};

struct ScenarioSatelliteConfig {
    string id;
    string name;
    StateTr x_tr;
    ScenarioStateAttConfig x_att;
    ScenarioPropagationConfig propagation;
    f64 mass = 0.0;
    vec3d inertia_diag = vec3d0;
};

struct ScenarioCovarianceConfig;
struct ScenarioInstrumentConfig;
struct ScenarioStationConfig;
struct ScenarioWorldStepperConfig;
struct ScenarioConfig {};

StatusCode load_scenario_json(const std::string& filepath, ScenarioConfig& out);
StatusCode validate_scenario_config(const ScenarioConfig& cfg);
