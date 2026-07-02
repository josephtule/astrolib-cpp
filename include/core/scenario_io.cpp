#pragma once

#include "core/estimation_common.hpp"
#include "core/state.hpp"
#include "nlohmann/json.hpp"
#include "nlohmann/json_fwd.hpp"
#include "util/typedefs.hpp"
#include "util/vecdefs.hpp"

using json = nlohmann::json;
using string = std::string;

struct ScenarioSchemaConfig {
    string name = "astrolib.scenario";
    i32 version = 1;
};
struct ScenarioMetadataConfig {
    string name = "earth_moon_sat_demo";
    i32 rng_seed = 12345;
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
    string attitude_quaternino = "passive_N_to_B";
    string station_local = "ENU";
};

struct ScenarioGravityProviderConfig {
    string id;
    string model;
    string filepath;
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
    string type = "quaternion";
    vec4d q = q_identity;
    vec3d axis = vec3d0;
    f64 angle = 0.0;
    vec3d w = vec3d0;
};

struct ScenarioCelesetialConfig {
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
    f64 mass = 0.0;
    vec3d inertia_diag = vec3d0;
};

inline StatusCode load_scenario_json() { return StatusCode::ok; }