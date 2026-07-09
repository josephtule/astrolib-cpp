#include "core/scenario_io.hpp"
#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/estimation_common.hpp"
#include "core/integrator.hpp"
#include "core/measurement.hpp"
#include "core/observation_type.hpp"
#include "core/orbital_elements.hpp"
#include "core/planets.hpp"
#include "core/state.hpp"
#include "core/time.hpp"
#include "core/transform.hpp"

#include "core/world.hpp"
#include "core/world_stepper.hpp"
#include "raylib.h"
#include "util/constants.hpp"
#include "util/math.hpp"
#include "util/tools.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include "nlohmann/json.hpp"

#include <cmath>
#include <filesystem>
#include <memory>

using json = nlohmann::json;

// TODO: organize all of this
// TODO: add alternate names for enums
// TODO: change units to default after conversion

// key queries

static StatusCode get_req_child(
    const json& object,
    const string& key,
    const json*& child,
    const string& path
) {
    if (!object.is_object()) {
        return StatusCode::invalid_input;
    }
    if (!object.contains(key)) {
        return StatusCode::invalid_input;
    }
    child = &object.at(key);

    return StatusCode::ok;
}

static StatusCode get_opt_child(
    const json& object,
    const string& key,
    const json*& child,
    bool& found,
    const string& path
) {
    found = false;
    child = nullptr;

    if (!object.is_object()) {
        return StatusCode::invalid_input;
    }
    if (!object.contains(key)) {
        return StatusCode::ok;
    }

    found = true;
    child = &object.at(key);

    return StatusCode::ok;
}

// type parsers

static StatusCode parse_req_bool(
    const json& object,
    const string& key,
    bool& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    if (!child->is_boolean()) return StatusCode::invalid_input;

    out = child->get<bool>();
    return StatusCode::ok;
}

static StatusCode parse_opt_bool(
    const json& object,
    const string& key,
    bool& found,
    bool& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_opt_child(object, key, child, found, path);
    if (status != StatusCode::ok) return status;

    if (!found) {
        // no default
        return StatusCode::ok;
    }

    if (!child->is_boolean()) return StatusCode::invalid_input;

    out = child->get<bool>();
    return StatusCode::ok;
}

static StatusCode parse_req_string(
    const json& object,
    const string& key,
    string& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    if (!child->is_string()) {
        return StatusCode::invalid_input;
    }

    out = child->get<string>();
    return StatusCode::ok;
}

static StatusCode parse_opt_string(
    const json& object,
    const string& key,
    bool& found,
    string& out,
    const string& path
) {
    found = false;
    const json* child = nullptr;

    StatusCode status = get_opt_child(object, key, child, found, path);
    if (status != StatusCode::ok) return status;
    if (!found) return StatusCode::ok;

    if (!child->is_string()) return StatusCode::invalid_input;

    out = child->get<string>();

    return StatusCode::ok;
}

static StatusCode parse_req_f64(
    const json& object,
    const string& key,
    f64& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    if (!child->is_number()) {
        return StatusCode::invalid_input;
    }

    f64 value = child->get<f64>();
    if (!std::isfinite(value)) return StatusCode::invalid_input;

    out = value;
    return StatusCode::ok;
}

static StatusCode parse_opt_f64(
    const json& object,
    const string& key,
    bool& found,
    f64& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_opt_child(object, key, child, found, path);
    if (status != StatusCode::ok) return status;
    if (!found) return StatusCode::ok;

    if (!child->is_number()) return StatusCode::invalid_input;

    f64 value = child->get<f64>();
    if (!std::isfinite(value)) return StatusCode::invalid_input;

    out = value;
    return StatusCode::ok;
}

static StatusCode parse_req_i32(
    const json& object,
    const string& key,
    i32& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    if (!child->is_number()) {
        return StatusCode::invalid_input;
    }

    out = child->get<i32>();
    return StatusCode::ok;
}

static StatusCode parse_opt_i32(
    const json& object,
    const string& key,
    bool& found,
    i32& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_opt_child(object, key, child, found, path);
    if (status != StatusCode::ok) return status;
    if (!found) return StatusCode::ok;

    if (!child->is_number()) return StatusCode::invalid_input;

    out = child->get<i32>();
    return StatusCode::ok;
}

static StatusCode parse_req_u32(
    const json& object,
    const string& key,
    u32& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    if (!child->is_number()) {
        return StatusCode::invalid_input;
    }

    out = child->get<u32>();
    return StatusCode::ok;
}

static StatusCode parse_opt_u32(
    const json& object,
    const string& key,
    bool& found,
    u32& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_opt_child(object, key, child, found, path);
    if (status != StatusCode::ok) return status;
    if (!found) return StatusCode::ok;

    if (!child->is_number()) return StatusCode::invalid_input;

    out = child->get<u32>();
    return StatusCode::ok;
}

static StatusCode parse_vec3d(const json& value, vec3d& out, const string& path) {
    if (!value.is_array() || value.size() != 3) {
        return StatusCode::invalid_input;
    }

    vec3d temp;
    for (i32 i = 0; i < 3; ++i) {
        if (!value.at(i).is_number()) {
            return StatusCode::invalid_input;
        }

        temp(i) = value[i].get<f64>();
        if (!std::isfinite(temp(i))) {
            return StatusCode::invalid_input;
        }
    }

    out = temp;
    return StatusCode::ok;
}

static StatusCode parse_vec4d(const json& value, vec4d& out, const string& path) {
    if (!value.is_array() || value.size() != 4) {
        return StatusCode::invalid_input;
    }

    vec4d temp;
    for (i32 i = 0; i < 4; ++i) {
        if (!value.at(i).is_number()) {
            return StatusCode::invalid_input;
        }

        temp(i) = value[i].get<f64>();
        if (!std::isfinite(temp(i))) {
            return StatusCode::invalid_input;
        }
    }
    out = temp;

    return StatusCode::ok;
}

static StatusCode parse_flat_mat3d(const json& value, mat3d& out, const string& path) {
    if (!value.is_array() || value.size() != 9) return StatusCode::invalid_input;

    mat3d temp;
    for (i32 i = 0; i < 9; ++i) {
        if (!value.at(i).is_number()) return StatusCode::invalid_input;
        temp(i / 3, i % 3) = value[i].get<f64>();
        if (!std::isfinite(temp(i / 3, i % 3))) return StatusCode::invalid_input;
    }

    out = temp;
    return StatusCode::ok;
}

static StatusCode parse_flat_matXd(
    const json& value,
    matXd& out,
    const i32 n,
    const i32 m,
    const string& path
) {
    if (!value.is_array()) return StatusCode::invalid_input;
    if (n <= 0 || m <= 0 || value.size() != n * m) return StatusCode::invalid_input;

    matXd temp(n, m);
    for (i32 i = 0; i < n * m; ++i) {
        if (!value.at(i).is_number()) return StatusCode::invalid_input;
        temp(i / m, i % m) = value[i].get<f64>();
        if (!std::isfinite(temp(i / m, i % m))) return StatusCode::invalid_input;
    }

    out = temp;
    return StatusCode::ok;
}

template <int N, int M>
static StatusCode parse_flat_matNMd(
    const json& value,
    mat<f64, N, M>& out,
    const string& path
) {
    if (!value.is_array()) return StatusCode::invalid_input;
    if (N <= 0 || M <= 0 || value.size() != N * M) return StatusCode::invalid_input;

    mat<f64, N, M> temp;
    for (i32 i = 0; i < N * M; ++i) {
        if (!value.at(i).is_number()) return StatusCode::invalid_input;
        temp(i / M, i % M) = value[i].get<f64>();
        if (!std::isfinite(temp(i / M, i % M))) return StatusCode::invalid_input;
    }

    out = temp;
    return StatusCode::ok;
}

static StatusCode parse_nested_matXd(const json& value, matXd& out, const string& path) {
    // array = [[row0], [row1], ...]
    if (!value.is_array()) return StatusCode::invalid_input;
    i32 n = value.size();
    if (n <= 0) return StatusCode::invalid_input;
    if (!value[0].is_array()) return StatusCode::invalid_input;
    i32 m = value[0].size();

    matXd temp(n, m);
    for (i32 i = 0; i < n; ++i) {
        if (!value[i].is_array()) return StatusCode::invalid_input;
        i32 m_current = value[i].size();
        if (m_current != m) return StatusCode::invalid_input; // check row size
        for (i32 j = 0; j < m; ++j) {
            if (!value.at(i).at(j).is_number()) return StatusCode::invalid_input;

            temp(i, j) = value[i][j].get<f64>();
            if (!std::isfinite(temp(i, j))) return StatusCode::invalid_input;
        }
    }

    out = temp;
    return StatusCode::ok;
}

static StatusCode parse_nested_mat3d(const json& value, mat3d& out, const string& path) {
    if (!value.is_array() || value.size() != 3) return StatusCode::invalid_input;

    mat3d temp;
    for (i32 i = 0; i < 3; ++i) {
        if (!value[i].is_array()) return StatusCode::invalid_input;
        i32 m_current = value[i].size();
        if (m_current != 3) return StatusCode::invalid_input;
        for (i32 j = 0; j < 3; ++j) {
            if (!value.at(i).at(j).is_number()) return StatusCode::invalid_input;

            temp(i, j) = value[i][j].get<f64>();
            if (!std::isfinite(temp(i, j))) return StatusCode::invalid_input;
        }
    }

    out = temp;
    return StatusCode::ok;
}

static StatusCode parse_vecXd(const json& value, vecXd& out, const string& path) {
    if (!value.is_array()) return StatusCode::invalid_input;

    i32 n = value.size();
    vecXd temp(n);
    for (i32 i = 0; i < n; ++i) {
        if (!value.at(i).is_number()) return StatusCode::invalid_input;
        temp(i) = value[i].get<f64>();
        if (!std::isfinite(temp(i))) return StatusCode::invalid_input;
    }

    out = temp;
    return StatusCode::ok;
}

static StatusCode parse_array_size(const json& value, i32& out, const string& path) {
    if (!value.is_array()) return StatusCode::invalid_input;

    out = value.size();
    return StatusCode::ok;
}

static StatusCode parse_array_size(
    const json& value,
    const string& key,
    i32& out,
    const string& path
) {
    if (!value.is_object()) return StatusCode::invalid_input;

    const json* array_child = nullptr;
    StatusCode status = get_req_child(value, key, array_child, path);
    if (status != StatusCode::ok) return status;

    status = parse_array_size(*array_child, out, path);
    if (status != StatusCode::ok) return status;

    return StatusCode::ok;
}

template <typename T, int N>
static StatusCode parse_array(
    const json& value,
    std::array<T, N>& out,
    const string& path
) {
    if (!value.is_array() || value.size() != N) return StatusCode::invalid_input;

    std::array<T, N> temp;
    for (i32 i = 0; i < N; ++i) {
        if (!value.at(i).is_number()) return StatusCode::invalid_input;
        temp[i] = value[i].get<T>();
        if (!std::isfinite(temp[i])) return StatusCode::invalid_input;
    }

    out = temp;
    return StatusCode::ok;
}

template <typename T, int N>
static StatusCode parse_req_array(
    const json& object,
    const std::string& key,
    std::array<T, N>& out,
    const std::string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_array<T, N>(*child, out, path + "." + key);
}

static StatusCode parse_req_flat_mat3d(
    const json& object,
    const string& key,
    mat3d& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_flat_mat3d(*child, out, path + "." + key);
}

static StatusCode parse_opt_flat_mat3d(
    const json& object,
    const string& key,
    bool& found,
    mat3d& out,
    const string& path
) {
    found = false;
    const json* child = nullptr;

    StatusCode status = get_opt_child(object, key, child, found, path);
    if (status != StatusCode::ok) return status;
    if (!found) return StatusCode::ok;

    return parse_flat_mat3d(*child, out, path + "." + key);
}

static StatusCode parse_req_flat_matXd(
    const json& object,
    const string& key,
    matXd& out,
    const i32 n,
    const i32 m,
    const string& path
) {
    const json* child = nullptr;

    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_flat_matXd(*child, out, n, m, path + "." + key);
}

template <int N, int M>
static StatusCode parse_req_flat_matNMd(
    const json& object,
    const string& key,
    mat<f64, N, M>& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_flat_matNMd<N, M>(*child, out, path + "." + key);
}

static StatusCode parse_req_nested_matXd(
    const json& object,
    const string& key,
    matXd& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_nested_matXd(*child, out, path + "." + key);
}

static StatusCode parse_req_nested_mat3d(
    const json& object,
    const string& key,
    mat3d& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_nested_mat3d(*child, out, path + "." + key);
}

static StatusCode parse_req_vec3d(
    const json& object,
    const std::string& key,
    vec3d& out,
    const std::string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_vec3d(*child, out, path + "." + key);
}

static StatusCode parse_opt_vec3d(
    const json& object,
    const string& key,
    bool& found,
    vec3d& out,
    const string& path
) {
    found = false;
    const json* child = nullptr;

    StatusCode status = get_opt_child(object, key, child, found, path);
    if (status != StatusCode::ok) return status;
    if (!found) return StatusCode::ok;

    return parse_vec3d(*child, out, path + "." + key);
}

static StatusCode parse_req_vec4d(
    const json& object,
    const string& key,
    vec4d& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_vec4d(*child, out, path + "." + key);
}

static StatusCode parse_req_vecXd(
    const json& object,
    const string& key,
    vecXd& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_vecXd(*child, out, path + "." + key);
}

static StatusCode parse_opt_vecXd(
    const json& object,
    const string& key,
    bool& found,
    vecXd& out,
    const string& path
) {
    found = false;
    const json* child = nullptr;

    StatusCode status = get_opt_child(object, key, child, found, path);
    if (status != StatusCode::ok) return status;
    if (!found) return StatusCode::ok;

    return parse_vecXd(*child, out, path + "." + key);
}

// enum parsers

static StatusCode parse_units_angle(const string& str, UAngle& out) {
    if (str == "radian")
        out = UAngle::radian;
    else if (str == "degree")
        out = UAngle::degree;
    else if (str == "arcminute")
        out = UAngle::arcminute;
    else if (str == "arcsecond")
        out = UAngle::arcsecond;
    else if (str == "milliarcsecond")
        out = UAngle::milliarcsecond;
    else
        return StatusCode::invalid_input;

    return StatusCode::ok;
}

static StatusCode parse_state_tr_type(const string& str, StateTrInputType& out) {
    if (str == "pos_vel")
        out = StateTrInputType::pos_vel;
    else if (str == "classical")
        out = StateTrInputType::classical;
    else
        return StatusCode::invalid_input;

    return StatusCode::ok;
}

static StatusCode parse_req_units_angle(
    const json& object,
    const string& key,
    UAngle& out,
    const string& path
) {
    string units_angle_str;
    StatusCode status = parse_req_string(object, key, units_angle_str, path);
    if (status != StatusCode::ok) return status;

    return parse_units_angle(units_angle_str, out);
}

static StatusCode parse_opt_units_angle(
    const json& object,
    const string& key,
    bool& found,
    UAngle& out,
    const string& path
) {
    string units_angle_str;
    StatusCode status = parse_opt_string(object, key, found, units_angle_str, path);
    if (status != StatusCode::ok) return status;
    if (!found) {
        out = UAngle::radian;
        return StatusCode::ok;
    }

    return parse_units_angle(units_angle_str, out);
}

static StatusCode parse_units_length(const string& str, ULength& out) {
    if (str == "kilometer")
        out = ULength::kilometer;
    else if (str == "millimeter")
        out = ULength::millimeter;
    else if (str == "centimeter")
        out = ULength::centimeter;
    else if (str == "meter")
        out = ULength::meter;
    else if (str == "inch")
        out = ULength::inch;
    else if (str == "foot")
        out = ULength::foot;
    else if (str == "mile")
        out = ULength::mile;
    else if (str == "au")
        out = ULength::au;
    else
        return StatusCode::invalid_input;

    return StatusCode::ok;
}

static StatusCode parse_req_units_length(
    const json& object,
    const string& key,
    ULength& out,
    const string& path
) {
    string units_length_str;
    StatusCode status = parse_req_string(object, key, units_length_str, path);
    if (status != StatusCode::ok) return status;

    return parse_units_length(units_length_str, out);
}

static StatusCode parse_opt_units_length(
    const json& object,
    const string& key,
    bool& found,
    ULength& out,
    const string& path
) {
    string units_length_str;
    StatusCode status = parse_opt_string(object, key, found, units_length_str, path);
    if (status != StatusCode::ok) return status;
    if (!found) {
        out = ULength::kilometer;
        return StatusCode::ok;
    }

    return parse_units_length(units_length_str, out);
}

static StatusCode parse_observation_type(string str, ObservationType& out) {
    // NOTE: update if adding more types
    if (str == "radec")
        out = ObservationType::radec;
    else if (str == "azel")
        out = ObservationType::azel;
    else if (str == "range")
        out = ObservationType::range;
    else if (str == "range_rate")
        out = ObservationType::range_rate;
    else if (str == "pos")
        out = ObservationType::pos;
    else if (str == "pos_vel")
        out = ObservationType::pos_vel;
    else if (str == "rel_pos")
        out = ObservationType::rel_pos;
    else if (str == "rel_pos_vel")
        out = ObservationType::rel_pos_vel;
    else
        return StatusCode::invalid_input;

    return StatusCode::ok;
}

static StatusCode parse_gravity_model(string str, GravityModel& out) {
    // NOTE: update if adding more types
    if (str == "pointmass")
        out = GravityModel::pointmass;
    else if (str == "zonal")
        out = GravityModel::zonal;
    else if (str == "spherical_harmonics")
        out = GravityModel::spherical_harmonics;
    else
        return StatusCode::invalid_input;

    return StatusCode::ok;
}

static StatusCode parse_celestial_attitude_model(
    string str,
    CelestialAttitudeModel& out
) {
    str = make_lower(str);

    if (str == "fixed")
        out = CelestialAttitudeModel::fixed;
    else if (str == "simple_spin" || str == "simple spin")
        out = CelestialAttitudeModel::simple_spin;
    else if (str == "provider")
        out = CelestialAttitudeModel::provider;
    else
        return StatusCode::invalid_input;

    return StatusCode::ok;
}

static StatusCode parse_opt_celestial_attitude_model(
    const json& object,
    const string& key,
    bool& found,
    CelestialAttitudeModel& out,
    const string& path
) {
    string model_str;
    StatusCode status = parse_opt_string(object, key, found, model_str, path);
    if (status != StatusCode::ok) return status;
    if (!found) return StatusCode::ok;

    return parse_celestial_attitude_model(model_str, out);
}

static StatusCode parse_req_gravity_model(
    const json& object,
    const string& key,
    GravityModel& out,
    const string& path
) {
    string type_str;
    StatusCode status = parse_req_string(object, key, type_str, path);
    if (status != StatusCode::ok) return status;

    return parse_gravity_model(type_str, out);
}

static StatusCode parse_opt_gravity_model(
    const json& object,
    const string& key,
    bool& found,
    GravityModel& out,
    const string& path
) {
    string model_str;
    StatusCode status = parse_opt_string(object, key, found, model_str, path);
    if (!found) {
        out = GravityModel::pointmass; // default
        return StatusCode::ok;
    }

    return parse_gravity_model(model_str, out);
}

static StatusCode parse_integrator_type(string str, IntegratorType& out) {
    // NOTE: update if adding more types
    if (str == "rk1")
        out = IntegratorType::rk1;
    else if (str == "rk2")
        out = IntegratorType::rk2;
    else if (str == "rk2_heun")
        out = IntegratorType::rk2_heun;
    else if (str == "rk2_ralston")
        out = IntegratorType::rk2_ralston;
    else if (str == "rk3")
        out = IntegratorType::rk3;
    else if (str == "rk4")
        out = IntegratorType::rk4;
    else
        return StatusCode::invalid_input;

    return StatusCode::ok;
}

static StatusCode parse_req_integrator_type(
    const json& object,
    const string& key,
    IntegratorType& out,
    const string& path
) {
    string integrator_type_str;
    StatusCode status = parse_req_string(object, key, integrator_type_str, path);
    if (status != StatusCode::ok) return status;

    return parse_integrator_type(integrator_type_str, out);
}

static StatusCode parse_opt_integrator_type(
    const json& object,
    const string& key,
    bool& found,
    IntegratorType& out,
    const string& path
) {
    string integrator_type_str;
    StatusCode status = parse_opt_string(object, key, found, integrator_type_str, path);
    if (status != StatusCode::ok) return status;
    if (!found) {
        out = IntegratorType::rk4;
        return StatusCode::ok;
    }

    return parse_integrator_type(integrator_type_str, out);
}

static StatusCode parse_attitude_type(string str, AttitudeType& out) {
    if (str == "quaternion")
        out = AttitudeType::quaternion;
    else if (str == "dcm")
        out = AttitudeType::dcm;
    else if (str == "axis_angle")
        out = AttitudeType::axis_angle;
    else if (str == "euler_angles")
        out = AttitudeType::euler_angles;
    else if (str == "crp")
        out = AttitudeType::crp;
    else if (str == "mrp")
        out = AttitudeType::mrp;
    else
        return StatusCode::invalid_input;

    return StatusCode::ok;
}

static StatusCode parse_req_attitude_type(
    const json& object,
    const string& key,
    AttitudeType& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    string type_str;
    StatusCode status = parse_req_string(object, key, type_str, path);
    if (status != StatusCode::ok) return status;

    return parse_attitude_type(type_str, out);
}

static StatusCode parse_opt_attitude_type(
    const json& object,
    const string& key,
    bool& found,
    AttitudeType& out,
    const string& path
) {
    string type_str;
    StatusCode status = parse_opt_string(object, key, found, type_str, path);
    if (status != StatusCode::ok) return status;
    if (!found) {
        out = AttitudeType::quaternion; // default
        return StatusCode::ok;
    }

    return parse_attitude_type(type_str, out);
}

static StatusCode parse_time_scale(const string& str, TimeScale& out) {
    // TODO: add other spellings
    if (str == "utc")
        out = TimeScale::utc;
    else if (str == "ut1" || str == "ut")
        out = TimeScale::ut1;
    else if (str == "tai")
        out = TimeScale::tai;
    else if (str == "tt")
        out = TimeScale::tt;
    else if (str == "tdb")
        out = TimeScale::tdb;
    else if (str == "gps")
        out = TimeScale::gps;
    else
        return StatusCode::invalid_input;

    return StatusCode::ok;
}

static StatusCode parse_req_time_scale(
    const json& object,
    const string& key,
    TimeScale& out,
    const string& path
) {
    string type_str;
    StatusCode status = parse_req_string(object, key, type_str, path);
    if (status != StatusCode::ok) return status;

    return parse_time_scale(type_str, out);
}

static StatusCode parse_opt_time_scale(
    const json& object,
    const string& key,
    bool& found,
    TimeScale& out,
    const string& path
) {
    string type_str;
    StatusCode status = parse_opt_string(object, key, found, type_str, path);
    if (status != StatusCode::ok) return status;
    if (!found) {
        out = TimeScale::utc;
        return StatusCode::ok;
    }

    return parse_time_scale(type_str, out);
}

static StatusCode parse_date_type(const string& str, DateType& out) {
    // TODO: add other spellings
    if (in_list<string>(
            str,
            {"cal", "calendar", "cal_time", "cal_t", "calendar_time", "calendar_t"}
        ))
        out = DateType::cal;
    else if (in_list<string>(str, {"julian date", "julian_date", "jd"}))
        out = DateType::jd;
    else if (
        in_list<string>(str, {"modified julian date", "modified_julian_date", "mjd"})
    )
        out = DateType::mjd;

    return StatusCode::ok;
}

static StatusCode parse_req_date_type(
    const json& object,
    const string& key,
    DateType& out,
    const string& path
) {
    string type_str;
    StatusCode status = parse_req_string(object, key, type_str, path);
    if (status != StatusCode::ok) return status;

    return parse_date_type(type_str, out);
}

static StatusCode parse_cal_style(const string& str, CalendarPrintStyle& out) {
    if (str == "separate")
        out = CalendarPrintStyle::separate;
    else if (str == "vector")
        out = CalendarPrintStyle::vector;
    else if (str == "string")
        out = CalendarPrintStyle::string;

    return StatusCode::ok;
}

static StatusCode parse_req_cal_style(
    const json& object,
    const string& key,
    CalendarPrintStyle& out,
    const string& path
) {
    string style_str;
    StatusCode status = parse_req_string(object, key, style_str, path);
    if (status != StatusCode::ok) return status;

    return parse_cal_style(style_str, out);
}
// config parsers

static StatusCode parse_state_tr_config(
    const json& object,
    ScenarioStateTrConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;
    const json* state_tr = nullptr;
    status = get_req_child(object, "state_tr", state_tr, path);
    if (status != StatusCode::ok) return status;

    bool found_ulength;
    status = parse_opt_units_length(
        *state_tr,
        "units_length",
        found_ulength,
        out.units_length,
        path
    );
    if (status != StatusCode::ok) return status;
    if (!found_ulength) out.units_length = ULength::kilometer; // default

    bool found_inputtype;
    string input_type_str = "";
    status = parse_opt_string(
        *state_tr,
        "input_type",
        found_inputtype,
        input_type_str,
        path
    );
    if (status != StatusCode::ok) return status;
    if (!found_inputtype) {
        out.input_type = StateTrInputType::pos_vel; // default
    } else {
        status = parse_state_tr_type(input_type_str, out.input_type);
    }

    switch (out.input_type) {
    case StateTrInputType::pos_vel: {
        status = parse_req_vec3d(*state_tr, "r", out.r, path);
        if (status != StatusCode::ok) return status;

        status = parse_req_vec3d(*state_tr, "v", out.v, path);
        if (status != StatusCode::ok) return status;
    } break;
    case StateTrInputType::classical: {
        status = parse_req_string(*state_tr, "center", out.central, path);
        if (status != StatusCode::ok)
            status = parse_req_string(*state_tr, "central", out.central, path);
        if (status != StatusCode::ok) return status;

        bool found_uangle;
        status = parse_opt_units_angle(
            *state_tr,
            "units_angle",
            found_uangle,
            out.units_angle,
            path
        );
        if (status != StatusCode::ok) return status;
        if (!found_uangle) out.units_angle = UAngle::degree; // assume degree input

        // TODO: add additional spellings
        status = parse_req_f64(*state_tr, "sma", out.coes.sma, path);
        if (status != StatusCode::ok) return status;
        status = parse_req_f64(*state_tr, "ecc", out.coes.ecc, path);
        if (status != StatusCode::ok) return status;
        status = parse_req_f64(*state_tr, "inc", out.coes.inc, path);
        if (status != StatusCode::ok) return status;
        status = parse_req_f64(*state_tr, "raan", out.coes.raan, path);
        if (status != StatusCode::ok) return status;
        status = parse_req_f64(*state_tr, "aop", out.coes.aop, path);
        if (status != StatusCode::ok) return status;
        status = parse_req_f64(*state_tr, "ta", out.coes.ta, path);
        if (status != StatusCode::ok) return status;

    } break;
    }

    return StatusCode::ok;
}

static StatusCode parse_covariance(
    const json& object,
    ScenarioCovarianceConfig& out,
    i32 dim,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    string type = "";
    StatusCode status = parse_req_string(object, "type", type, path);
    if (status != StatusCode::ok) return status;

    out.covariance = matXd::Zero(dim, dim);

    i32 input_size = -1;
    if (type == "diagonal" || type == "diag") {
        status = parse_array_size(object, "values", input_size, path);
        if (status != StatusCode::ok) return status;
        if (input_size != dim) return StatusCode::invalid_covariance;

        vecXd temp_diag;
        status = parse_req_vecXd(object, "values", temp_diag, path);
        if (status != StatusCode::ok) return status;

        out.covariance = temp_diag.asDiagonal();
    } else if (type == "flattened" || type == "flat") {
        status = parse_req_flat_matXd(object, "values", out.covariance, dim, dim, path);
        if (status != StatusCode::ok) return status;
    } else if (type == "matrix" || type == "mat" || type == "nested") {
        status = parse_array_size(object, "values", input_size, path);
        if (status != StatusCode::ok) return status;
        if (input_size != dim) return StatusCode::invalid_covariance;

        status = parse_req_nested_matXd(object, "values", out.covariance, path);
        if (status != StatusCode::ok) return status;
    } else {
        return StatusCode::invalid_input;
    }

    return StatusCode::ok;
}

// Config Parsing
// ----------------------------------------------------------------------------------

static StatusCode parse_instrument_config(
    const json& object,
    ScenarioInstrumentConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status = parse_req_string(object, "id", out.id, path);
    if (status != StatusCode::ok) return status;

    string type_str;
    status = parse_req_string(object, "type", type_str, path);
    if (status != StatusCode::ok) return status;
    status = parse_observation_type(type_str, out.type);
    if (status != StatusCode::ok) return status;
    i32 dim = measurement_dim(out.type);

    status = parse_req_bool(object, "enabled", out.enabled, path);
    if (status != StatusCode::ok) return status;

    const json* cov_child = nullptr;
    status = get_req_child(object, "covariance", cov_child, path);
    if (status != StatusCode::ok) return status;
    status = parse_covariance(*cov_child, out.covariance_cfg, dim, path);
    if (status != StatusCode::ok) return status;

    return StatusCode::ok;
}

static StatusCode parse_instrument_template(
    const json& object,
    ScenarioInstrumentConfig& out,
    const string& path
) {
    return parse_instrument_config(object, out, path);
}

static StatusCode parse_state_att_config(
    const json& object,
    ScenarioStateAttConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;

    bool found_x_att;
    const json* state_att = nullptr;
    status = get_opt_child(object, "state_att", state_att, found_x_att, path);
    if (status != StatusCode::ok) return status;
    if (!found_x_att) {
        out.input_type = AttitudeType::quaternion;
        out.q = q_identity;
        out.w = vec3d0;
        return StatusCode::ok;
    }

    bool found_uangle;
    status = parse_opt_units_angle(
        *state_att,
        "units_angle",
        found_uangle,
        out.units_angle,
        path
    );
    if (status != StatusCode::ok) return status;
    if (!found_uangle) out.units_angle = UAngle::radian;

    bool found_type = false;
    status = parse_opt_attitude_type(
        *state_att,
        "input_type",
        found_type,
        out.input_type,
        path
    );
    if (status != StatusCode::ok) return status;

    // TODO: FIXME: move this to scenario -> world building
    if (found_type) {
        switch (out.input_type) {
        case AttitudeType::quaternion: {
            status = parse_req_vec4d(*state_att, "q", out.q, path);
            if (status != StatusCode::ok) return status;
        } break;
        case AttitudeType::dcm: {
            status = parse_req_flat_matNMd(*state_att, "dcm", out.dcm, path);
            if (status != StatusCode::ok) return status;
            out.q = dcm_to_ep(out.dcm);
        } break;
        case AttitudeType::axis_angle: {
            status = parse_req_vec3d(*state_att, "axis", out.axis, path);
            if (status != StatusCode::ok) return status;

            status = parse_req_f64(*state_att, "angle", out.angle, path);
            if (status != StatusCode::ok) return status;
            out.axis.normalize();

            out.q = axis_angle_to_ep(out.axis, out.angle, out.units_angle);
        } break;
        case AttitudeType::euler_angles: {
            status = parse_req_vec3d(*state_att, "angles", out.angles, path);
            if (status != StatusCode::ok) return status;

            status = parse_req_array<i32, 3>(*state_att, "sequence", out.sequence, path);
            if (status != StatusCode::ok) return status;

            std::array<RotAxis, 3> rots;
            for (i32 i = 0; i < 3; ++i) {
                status = i32_to_rotaxis(out.sequence[i], rots[i]);
                if (status != StatusCode::ok) return status;
            }
            mat3d dcm = ea_to_dcm(out.angles, rots, out.units_angle);

            out.q = dcm_to_ep(dcm);
        } break;
        case AttitudeType::crp: {
            status = parse_req_vec3d(*state_att, "axis", out.axis, path);
            if (status != StatusCode::ok) return status;
            out.q = crp_to_ep(out.axis);
        } break;
        case AttitudeType::mrp: {
            status = parse_req_vec3d(*state_att, "axis", out.axis, path);
            if (status != StatusCode::ok) return status;
            out.q = mrp_to_ep(out.axis);
        } break;
        default: {
            return StatusCode::attitude_type_not_found;
        }
        }
        if (out.q.norm() <= tol12) return StatusCode::invalid_input;
        out.q.normalize();
    } else {
        out.q = q_identity;
    }

    bool found_w;
    status = parse_opt_vec3d(*state_att, "w", found_w, out.w, path);
    if (status != StatusCode::ok) return status;
    if (!found_w) {
        out.w = vec3d0;
    }

    if (out.units_angle != UAngle::radian) {
        for (i32 i = 0; i < 3; ++i) { // units/s -> radians/s
            // TODO: create overload for vectors
            out.w(i) = convert_angle(out.w(i), out.units_angle, UAngle::radian);
        }
    }

    return StatusCode::ok;
}

static StatusCode parse_gravity_provider_config(
    const json& object,
    ScenarioGravityProviderConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;

    status = parse_req_string(object, "id", out.id, path);
    if (status != StatusCode::ok) return status;

    status = parse_req_string(object, "format", out.format, path);
    if (status != StatusCode::ok) return status;

    status = parse_req_string(object, "path", out.filepath, path);
    if (status != StatusCode::ok) return status;

    bool found_lineskips = false;
    status = parse_opt_i32(object, "lineskips", found_lineskips, out.lineskips, path);
    if (status != StatusCode::ok) return status;
    if (!found_lineskips) {
        status
            = parse_opt_i32(object, "line_skips", found_lineskips, out.lineskips, path);
        if (status != StatusCode::ok) return status;
    }
    if (!found_lineskips) out.lineskips = 0; // default

    bool found_normalized = false;
    status = parse_opt_bool(object, "normalized", found_normalized, out.normalized, path);
    if (status != StatusCode::ok) return status;
    if (!found_normalized) {
        out.normalized = true; // default
    }

    return StatusCode::ok;
}

static StatusCode parse_gravity_config(
    const json& object,
    ScenarioGravityConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;
    const json* gravity = nullptr;
    status = get_req_child(object, "gravity", gravity, path);
    if (status != StatusCode::ok) return status;

    status = parse_req_gravity_model(*gravity, "model", out.model, path);
    if (status != StatusCode::ok) return status;

    bool found_mu;
    status = parse_opt_f64(*gravity, "mu", found_mu, out.mu, path);
    if (status != StatusCode::ok) return status;
    // use from celestial model provider if missing

    // may parse higher fidelity model but use lower fidelity model

    bool found_units_length = false;
    status = parse_opt_units_length(
        *gravity,
        "units_length",
        found_units_length,
        out.units_length,
        path
    );
    if (status != StatusCode::ok) return status;
    if (!found_units_length) out.units_length = ULength::kilometer;

    if (found_mu) {
        f64 length_scale = convert_length(1.0, out.units_length, ULength::kilometer);
        out.mu *= length_scale * length_scale * length_scale;
    }

    bool found_radius = false;

    status = parse_opt_f64(*gravity, "radius", found_radius, out.radius, path);
    if (status != StatusCode::ok) return status;
    if (found_radius) {
        out.radius = convert_length(out.radius, out.units_length, ULength::kilometer);
    }

    bool _ = false;
    if (out.model != GravityModel::pointmass) {
        status = parse_req_i32(*gravity, "degree", out.degree, path);
        if (status != StatusCode::ok) return status;

        status = parse_req_i32(*gravity, "order", out.order, path);
        if (status != StatusCode::ok) return status;
    } else {
        status = parse_opt_i32(*gravity, "degree", _, out.degree, path);
        if (status != StatusCode::ok) return status;

        status = parse_opt_i32(*gravity, "order", _, out.order, path);
        if (status != StatusCode::ok) return status;
    }

    bool found_provider = false;
    status = parse_opt_string(
        *gravity,
        "coefficients",
        found_provider,
        out.coefficients,
        path
    );
    if (status != StatusCode::ok) return status;

    bool found_J = false;
    status = parse_opt_vecXd(*gravity, "J", found_J, out.J, path);
    if (status != StatusCode::ok) return status;

    bool found_zonal_coefficients = false;
    if (!found_J) {
        status = parse_opt_vecXd(
            *gravity,
            "zonal_coefficients",
            found_zonal_coefficients,
            out.J,
            path
        );
        if (status != StatusCode::ok) return status;
        found_J = found_zonal_coefficients;
    }

    if (found_provider && found_J) return StatusCode::invalid_input;

    if (found_provider) {
        out.coefficient_source = GravityCoefficientSource::provider;
    } else if (found_J) {
        out.coefficient_source = GravityCoefficientSource::direct_zonal;
    } else {
        out.coefficient_source = GravityCoefficientSource::none;
    }

    return StatusCode::ok;
}

static StatusCode resolve_custom_celestial_shape_config(
    ScenarioCelestialModelConfig& out,
    const bool found_semimajor,
    const bool found_semiminor,
    const bool found_mean_radius
) {
    if (found_semimajor) {
        out.semimajor_axis
            = convert_length(out.semimajor_axis, out.units_length, ULength::kilometer);
        if (!finite_nonneg(out.semimajor_axis)) return StatusCode::invalid_input;
    }

    if (found_semiminor) {
        out.semiminor_axis
            = convert_length(out.semiminor_axis, out.units_length, ULength::kilometer);
        if (!finite_nonneg(out.semiminor_axis)) return StatusCode::invalid_input;
    }

    if (found_mean_radius) {
        out.mean_radius
            = convert_length(out.mean_radius, out.units_length, ULength::kilometer);
        if (!finite_nonneg(out.mean_radius)) return StatusCode::invalid_input;
    }

    bool has_semimajor = finite_pos(out.semimajor_axis);
    bool has_semiminor = finite_pos(out.semiminor_axis);
    bool has_mean_radius = finite_pos(out.mean_radius);

    f64 shape_fallback = 0.0;
    if (has_semimajor) {
        shape_fallback = out.semimajor_axis;
    } else if (has_semiminor) {
        shape_fallback = out.semiminor_axis;
    } else if (has_mean_radius) {
        shape_fallback = out.mean_radius;
    } else {
        return StatusCode::invalid_input;
    }

    if (!has_semimajor) out.semimajor_axis = shape_fallback;
    if (!has_semiminor) out.semiminor_axis = shape_fallback;

    if (out.semimajor_axis <= 0.0 || out.semiminor_axis <= 0.0)
        return StatusCode::invalid_input;
    if (out.semimajor_axis < out.semiminor_axis) return StatusCode::invalid_input;

    if (!has_mean_radius) {
        out.mean_radius = std::pow(
            out.semimajor_axis * out.semimajor_axis * out.semiminor_axis,
            1.0 / 3.0
        );
    }
    if (!finite_pos(out.mean_radius)) return StatusCode::invalid_input;

    return StatusCode::ok;
}

static StatusCode parse_celestial_model_config(
    const json& object,
    ScenarioCelestialModelConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;
    const json* model = nullptr;
    status = get_req_child(object, "model", model, path);
    if (status != StatusCode::ok) return status;

    bool found_id;
    status = parse_opt_string(*model, "id", found_id, out.id, path);
    if (status != StatusCode::ok) return status;
    if (!found_id) {
        out.id = "custom";
    }

    if (out.id == "custom") {
        status = parse_gravity_config(*model, out.gravity_model, path);
        if (status != StatusCode::ok) return status;

        bool found_units_length = false;
        status = parse_opt_units_length(
            *model,
            "units_length",
            found_units_length,
            out.units_length,
            path
        );
        if (status != StatusCode::ok) return status;
        if (!found_units_length) out.units_length = ULength::kilometer;

        bool found_semimajor = false;
        status = parse_opt_f64(
            *model,
            "semimajor_axis",
            found_semimajor,
            out.semimajor_axis,
            path
        );
        if (status != StatusCode::ok) return status;

        bool found_semiminor = false;
        status = parse_opt_f64(
            *model,
            "semiminor_axis",
            found_semiminor,
            out.semiminor_axis,
            path
        );
        if (status != StatusCode::ok) return status;

        bool found_mean_radius = false;
        status = parse_opt_f64(
            *model,
            "mean_radius",
            found_mean_radius,
            out.mean_radius,
            path
        );
        if (status != StatusCode::ok) return status;

        status = resolve_custom_celestial_shape_config(
            out,
            found_semimajor,
            found_semiminor,
            found_mean_radius
        );
        if (status != StatusCode::ok) return status;

        bool found_eccentricity = false;
        status = parse_opt_f64(
            *model,
            "eccentricity",
            found_eccentricity,
            out.eccentricity,
            path
        );
        if (status != StatusCode::ok) return status;
        if (!found_eccentricity) {
            out.eccentricity
                = ecc_from_semiaxes<f64>(out.semimajor_axis, out.semiminor_axis);
        }

        bool found_flattening = false;
        status
            = parse_opt_f64(*model, "flattening", found_flattening, out.flattening, path);
        if (status != StatusCode::ok) return status;
        if (!found_flattening) {
            out.flattening
                = (out.semimajor_axis - out.semiminor_axis) / out.semimajor_axis;
        }
    } else {
        status = parse_gravity_config(*model, out.gravity_model, path);
        if (status != StatusCode::ok) return status;

        // let celestial model provider handle the rest
    }

    return StatusCode::ok;
}

static StatusCode parse_propagation_config(
    const json& object,
    ScenarioPropagationConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;
    const json* propagation = nullptr;
    status = get_req_child(object, "propagation", propagation, path);
    if (status != StatusCode::ok) return status;

    bool found_tr = false;
    status = parse_opt_bool(*propagation, "translation", found_tr, out.translation, path);
    if (status != StatusCode::ok) return status;
    if (!found_tr) out.translation = true;

    bool found_att = false;
    status = parse_opt_bool(*propagation, "attitude", found_att, out.attitude, path);
    if (status != StatusCode::ok) return status;
    if (!found_att) out.attitude = false;

    return StatusCode::ok;
}

static StatusCode parse_celestial_config(
    const json& object,
    ScenarioCelestialConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;

    status = parse_req_string(object, "id", out.id, path);
    if (status != StatusCode::ok) return status;

    bool _;

    status = parse_opt_string(object, "name", _, out.name, path);
    if (status != StatusCode::ok) return status;

    status = parse_celestial_model_config(object, out.model, path);
    if (status != StatusCode::ok) return status;

    status = parse_state_tr_config(object, out.x_tr, path);
    if (status != StatusCode::ok) return status;

    status = parse_state_att_config(object, out.x_att, path);
    if (status != StatusCode::ok) return status;

    status = parse_opt_celestial_attitude_model(
        object,
        "attitude_model",
        out.has_attitude_model,
        out.attitude_model,
        path
    );
    if (status != StatusCode::ok) return status;

    status = parse_propagation_config(object, out.propagation, path);
    if (status != StatusCode::ok) return status;

    return StatusCode::ok;
}

static StatusCode parse_opt_mass_properties_config(
    const json& object,
    bool& found,
    ScenarioMassPropertiesConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;
    StatusCode status;

    const json* mass_properties = nullptr;
    status = get_opt_child(object, "mass_properties", mass_properties, found, path);
    if (status != StatusCode::ok) return status;
    if (!found) return StatusCode::ok;

    bool found_mass;
    status = parse_opt_f64(*mass_properties, "mass", found_mass, out.mass, path);
    if (status != StatusCode::ok) return status;
    if (!found_mass) out.mass = 0.0;

    bool found_inertia;
    const json* inertia_tensor = nullptr;
    status = get_opt_child(
        *mass_properties,
        "inertia_tensor",
        inertia_tensor,
        found_inertia,
        path
    );
    if (status != StatusCode::ok) return status;

    if (!found_inertia && !found_mass) {
        // empty mass properties
        found = false;
        return StatusCode::ok;
    }

    if (found_inertia) {
        bool found_type;
        status = parse_opt_string(*inertia_tensor, "type", found_type, out.type, path);
        if (status != StatusCode::ok) return status;
        if (found_type) {
            bool found_PA;
            status = parse_opt_bool(
                *inertia_tensor,
                "principle_axes",
                found_PA,
                out.principle_axes,
                path
            );
            if (!found_PA) out.principle_axes = true;

            if (out.type == "diagonal" || out.type == "diag") {
                vec3d diag;
                status = parse_req_vec3d(*inertia_tensor, "values", diag, path);
                if (status != StatusCode::ok) return status;

                out.inertia = diag.asDiagonal();
            } else if (out.type == "flattened" || out.type == "flat") {
                status
                    = parse_req_flat_mat3d(*inertia_tensor, "values", out.inertia, path);
                if (status != StatusCode::ok) return status;
            } else if (
                out.type == "matrix" || out.type == "mat" || out.type == "nested"
            ) {
                status = parse_req_nested_mat3d(
                    *inertia_tensor,
                    "values",
                    out.inertia,
                    path
                );
                if (status != StatusCode::ok) return status;
            } else {
                return StatusCode::invalid_input;
            }
        } else {
            found = found_type;
        }
    }

    if (!found_inertia && found_mass) {
        // default to cube 1x1x1m cube with mass 1kg
        out.principle_axes = true;
        f64 length = 1.0 / 1000;
        out.inertia = 1.0 / 6.0 * out.mass * length * length * mat3d1;
    }

    status
        = parse_opt_vec3d(*mass_properties, "offset", out.offset, out.offset_body, path);
    if (status != StatusCode::ok) return status;

    return StatusCode::ok;
}

static StatusCode parse_satellite_config(
    const json& object,
    ScenarioSatelliteConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;

    status = parse_req_string(object, "id", out.id, path);
    if (status != StatusCode::ok) return status;

    bool _;

    status = parse_opt_string(object, "name", _, out.name, path);
    if (status != StatusCode::ok) return status;

    status = parse_state_tr_config(object, out.x_tr, path);
    if (status != StatusCode::ok) return status;

    status = parse_state_att_config(object, out.x_att, path);
    if (status != StatusCode::ok) return status;

    status = parse_propagation_config(object, out.propagation, path);
    if (status != StatusCode::ok) return status;

    bool found_mass_properties;
    status = parse_opt_mass_properties_config(
        object,
        found_mass_properties,
        out.mass_properties,
        path
    );
    if (!found_mass_properties)
        out.propagation.attitude = false; // cannot propagate without mass properties

    return StatusCode::ok;
}

static StatusCode parse_station_config(
    const json& object,
    ScenarioStationConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;

    bool found_ulength;
    status = parse_opt_units_length(
        object,
        "units_length",
        found_ulength,
        out.units_length,
        path
    );
    if (status != StatusCode::ok) return status;
    if (!found_ulength) out.units_length = ULength::kilometer;

    bool found_uangle;
    status = parse_opt_units_angle(
        object,
        "units_angle",
        found_uangle,
        out.units_angle,
        path
    );
    if (status != StatusCode::ok) return status;
    if (!found_uangle) out.units_angle = UAngle::degree;

    status = parse_req_string(object, "id", out.id, path);
    if (status != StatusCode::ok) return status;

    bool _;

    status = parse_opt_string(object, "name", _, out.name, path);
    if (status != StatusCode::ok) return status;

    bool found_anchored;
    status = parse_opt_bool(object, "anchored", found_anchored, out.anchored, path);
    if (status != StatusCode::ok) return status;
    if (!found_anchored) out.anchored = true;

    if (out.anchored) {
        status = parse_req_string(object, "anchor", out.anchor, path);
        if (status != StatusCode::ok) return status;

        bool found_coord;
        status = parse_opt_string(
            object,
            "coordinate_type",
            found_coord,
            out.coordinate_type,
            path
        );
        if (status != StatusCode::ok) return status;
        if (!found_coord) out.coordinate_type = "detic_llh"; // default
        // TODO: create enum for this?
        if (out.coordinate_type == "detic_llh") {
            status = parse_req_vec3d(object, "llh", out.llh_BCBF, path);
            if (status != StatusCode::ok) return status;
        } else if (out.coordinate_type == "body_fixed") {
            status = parse_req_vec3d(object, "r_body", out.r_body, path);
            if (status != StatusCode::ok) return status;
        }

        bool found_frame;
        status
            = parse_opt_string(object, "local_frame", found_frame, out.local_frame, path);
        if (status != StatusCode::ok) return status;
        if (!found_frame) out.local_frame = "ENU";

        out.propagation.translation = true;
        out.propagation.attitude = true;
    } else {
        status = parse_state_tr_config(object, out.x_tr, path);
        if (status != StatusCode::ok) return status;

        status = parse_state_att_config(object, out.x_att, path);
        if (status != StatusCode::ok) return status;

        status = parse_propagation_config(object, out.propagation, path);
        if (status != StatusCode::ok) return status;
    }

    bool found_mass;
    status
        = parse_opt_mass_properties_config(object, found_mass, out.mass_properties, path);
    if (status != StatusCode::ok) return status;
    // no default

    bool found_instruments;
    const json* child = nullptr;
    status = get_opt_child(object, "instruments", child, found_instruments, path);
    if (status != StatusCode::ok) return status;
    if (found_instruments) {
        i32 num_instruments = 0;
        status = parse_array_size(object, "instruments", num_instruments, path);
        if (status != StatusCode::ok) return status;

        out.instruments.reserve(num_instruments);
        for (i32 i = 0; i < num_instruments; ++i) {
            ScenarioInstrumentConfig temp_instrument;
            status = parse_instrument_config(child->at(i), temp_instrument, path);
            if (status != StatusCode::ok) return status;

            out.instruments.push_back(temp_instrument);
        }
    }

    return StatusCode::ok;
}

static StatusCode parse_scenario_schema(
    const json& object,
    ScenarioSchemaConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;
    const json* child = nullptr;
    status = get_req_child(object, "schema", child, path);
    if (status != StatusCode::ok) return status;

    status = parse_req_string(*child, "name", out.name, path);
    if (status != StatusCode::ok) return status;

    status = parse_req_i32(*child, "version", out.version, path);
    if (status != StatusCode::ok) return status;

    return StatusCode::ok;
}

static StatusCode parse_metadata_config(
    const json& object,
    ScenarioMetadataConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;
    const json* child = nullptr;
    status = get_req_child(object, "metadata", child, path);
    if (status != StatusCode::ok) return status;

    status = parse_req_string(*child, "name", out.name, path);
    if (status != StatusCode::ok) return status;

    bool found_rng;
    status = parse_opt_u32(*child, "rng_seed", found_rng, out.rng_seed, path);
    if (status != StatusCode::ok) return status;
    if (!found_rng) out.rng_seed = 12345;

    return StatusCode::ok;
}

static StatusCode parse_opt_world_stepper_config(
    const json& object,
    bool& found,
    ScenarioWorldStepperConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;
    const json* child = nullptr;
    status = get_opt_child(object, "world_stepper", child, found, path);
    if (status != StatusCode::ok) return status;
    if (!found) {
        out = ScenarioWorldStepperConfig{}; // defaults
        return StatusCode::ok;
    }

    bool found_paused;
    status = parse_opt_bool(*child, "paused", found_paused, out.paused, path);
    if (status != StatusCode::ok) return status;
    if (!found_paused) {
        out.paused = false;
    }

    bool found_int_tr;
    status = parse_opt_integrator_type(
        *child,
        "translation_integrator",
        found_int_tr,
        out.integrator_tr,
        path
    );
    if (status != StatusCode::ok) return status;

    bool found_int_att;
    status = parse_opt_integrator_type(
        *child,
        "attitude_integrator",
        found_int_att,
        out.integrator_att,
        path
    );
    if (status != StatusCode::ok) return status;

    bool temp_found;

    status = parse_opt_u32(*child, "substeps", temp_found, out.substeps, path);
    if (status != StatusCode::ok) return status;
    if (!temp_found) out.substeps = 1;

    status = parse_opt_u32(*child, "ticks", temp_found, out.ticks, path);
    if (status != StatusCode::ok) return status;
    if (!temp_found) out.ticks = 1;

    status = parse_opt_f64(*child, "dt_scale", temp_found, out.dt_scale, path);
    if (status != StatusCode::ok) return status;
    if (!temp_found) out.dt_scale = 1.0;

    return StatusCode::ok;
}

static StatusCode parse_opt_time_config(
    const json& object,
    bool& found,
    ScenarioTimeConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;
    const json* time = nullptr;
    status = get_opt_child(object, "time", time, found, path);
    if (status != StatusCode::ok) return status;
    if (!found) {
        out = ScenarioTimeConfig{}; // defaults
        return StatusCode::ok;
    }

    status = parse_opt_f64(*time, "t0", found, out.t0, path);
    if (status != StatusCode::ok) return status;
    if (!found) out.t0 = 0.0; // default

    const json* date = nullptr;
    status = get_opt_child(*time, "date", date, found, path);
    if (!found) out.time_scale = TimeScale::utc;
    if (!found) {
        out.time_scale = TimeScale::utc;
        out.date_type = DateType::jd;
        out.jd = JulianDate{};
        return StatusCode::ok;
    }

    status = parse_opt_time_scale(*date, "time_scale", found, out.time_scale, path);
    if (status != StatusCode::ok) return status;
    if (!found) out.time_scale = TimeScale::utc;

    status = parse_req_date_type(*date, "input_type", out.date_type, path);
    if (status != StatusCode::ok) return status;

    switch (out.date_type) {
    case DateType::cal: {
        const json* cal = nullptr;
        status = get_opt_child(*date, "cal", cal, found, path);
        if (status != StatusCode::ok) return status;
        if (!found) {
            out.cal = CalendarTime{};
        } else {
            status = parse_req_cal_style(*cal, "input_type", out.cal_style, path);
            if (status != StatusCode::ok) return status;

            switch (out.cal_style) {

            case CalendarPrintStyle::separate_vertical: [[fallthrough]];
            case CalendarPrintStyle::separate: {
                status = parse_req_i32(*cal, "year", out.cal.year, path);
                if (status != StatusCode::ok) return status;
                status = parse_req_i32(*cal, "month", out.cal.month, path);
                if (status != StatusCode::ok) return status;
                status = parse_req_i32(*cal, "day", out.cal.day, path);
                if (status != StatusCode::ok) return status;
                status = parse_req_i32(*cal, "hour", out.cal.hour, path);
                if (status != StatusCode::ok) return status;
                status = parse_req_i32(*cal, "minute", out.cal.minute, path);
                if (status != StatusCode::ok) return status;
                status = parse_req_f64(*cal, "second", out.cal.second, path);
                if (status != StatusCode::ok) return status;
            } break;
            case CalendarPrintStyle::string: {
                string cal_string;
                status = parse_req_string(*cal, "string", cal_string, path);
                if (status != StatusCode::ok) return status;

                bool read_ok = parse_cal_str(cal_string, out.cal);
                if (!read_ok) return StatusCode::invalid_input;
            } break;
            case CalendarPrintStyle::vector: {
                vecXd cal_vec;
                status = parse_req_vecXd(*cal, "vector", cal_vec, path);
                if (status != StatusCode::ok) return status;

                out.cal.year = cal_vec(0);
                out.cal.month = cal_vec(1);
                out.cal.day = cal_vec(2);
                out.cal.hour = cal_vec(3);
                out.cal.minute = cal_vec(4);
                out.cal.second = cal_vec(5);
            } break;
            }
        }
    } break;
    case DateType::jd: {
        const json* jd = nullptr;
        status = get_opt_child(*date, "jd", jd, found, path);
        if (status != StatusCode::ok) return status;
        if (!found) {
            out.jd = JulianDate{};
        } else {
            status = parse_req_f64(*jd, "day", out.jd.day, path);
            if (status != StatusCode::ok) return status;
            status = parse_req_f64(*jd, "frac", out.jd.frac, path);
        }
    } break;
    case DateType::mjd: {
        const json* mjd = nullptr;
        status = get_opt_child(*date, "mjd", mjd, found, path);
        if (status != StatusCode::ok) return status;
        if (!found) {
            out.mjd = ModifiedJulianDate{};
        } else {
            status = parse_req_f64(*mjd, "day", out.mjd.day, path);
            if (status != StatusCode::ok) return status;
            status = parse_req_f64(*mjd, "frac", out.mjd.frac, path);
        }
    } break;
    }

    return StatusCode::ok;
}

static StatusCode parse_scenario_config(
    const json& object,
    ScenarioConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StatusCode status;
    ScenarioConfig temp;

    status = parse_scenario_schema(object, temp.schema, path);
    if (status != StatusCode::ok) return status;

    status = parse_metadata_config(object, temp.metadata, path);
    if (status != StatusCode::ok) return status;

    bool found_time;
    status = parse_opt_time_config(object, found_time, temp.time, path);
    if (status != StatusCode::ok) return status;

    const json* providers = nullptr;
    bool found_providers;
    status = get_opt_child(object, "providers", providers, found_providers, path);
    if (status != StatusCode::ok) return status;
    if (found_providers) {
        const json* gravity_provider = nullptr;
        bool found_gravity;
        status
            = get_opt_child(*providers, "gravity", gravity_provider, found_gravity, path);
        if (status != StatusCode::ok) return status;
        if (found_gravity) {
            i32 num_gravity_providers = 0;
            status = parse_array_size(*gravity_provider, num_gravity_providers, path);
            if (status != StatusCode::ok) return status;
            for (i32 i = 0; i < num_gravity_providers; ++i) {
                ScenarioGravityProviderConfig grav_prov;
                status = parse_gravity_provider_config(
                    gravity_provider->at(i),
                    grav_prov,
                    path
                );
                if (status != StatusCode::ok) return status;

                temp.gravity_providers.push_back(grav_prov);
            }
        }
    }

    const json* celestials = nullptr;
    bool found_cel;
    status = get_opt_child(object, "celestials", celestials, found_cel, path);
    if (status != StatusCode::ok) return status;
    if (found_cel) {
        i32 num_celestials = 0;
        status = parse_array_size(*celestials, num_celestials, path);
        if (status != StatusCode::ok) return status;
        for (i32 i = 0; i < num_celestials; ++i) {
            ScenarioCelestialConfig cel;
            status = parse_celestial_config(celestials->at(i), cel, path);
            if (status != StatusCode::ok) return status;
            if (cel.name.empty()) cel.name = cel.id;
            temp.celestials.push_back(cel);
        }
    }

    const json* satellites = nullptr;
    bool found_sat;
    status = get_opt_child(object, "satellites", satellites, found_sat, path);
    if (status != StatusCode::ok) return status;
    if (found_sat) {
        i32 num_satellites = 0;
        status = parse_array_size(*satellites, num_satellites, path);
        if (status != StatusCode::ok) return status;
        for (i32 i = 0; i < num_satellites; ++i) {
            ScenarioSatelliteConfig sat;
            status = parse_satellite_config(satellites->at(i), sat, path);
            if (status != StatusCode::ok) return status;
            if (sat.name.empty()) sat.name = sat.id;
            temp.satellites.push_back(sat);
        }
    }

    const json* stations = nullptr;
    bool found_stat;
    status = get_opt_child(object, "stations", stations, found_stat, path);
    if (status != StatusCode::ok) return status;
    if (found_stat) {
        i32 num_stations = 0;
        status = parse_array_size(*stations, num_stations, path);
        if (status != StatusCode::ok) return status;
        for (i32 i = 0; i < num_stations; ++i) {
            ScenarioStationConfig stat;
            status = parse_station_config(stations->at(i), stat, path);
            if (status != StatusCode::ok) return status;
            if (stat.name.empty()) stat.name = stat.id;
            temp.stations.push_back(stat);
        }
    }

    bool found_world_stepper;
    status = parse_opt_world_stepper_config(
        object,
        found_world_stepper,
        temp.world_stepper,
        path
    );
    if (status != StatusCode::ok) return status;

    // TODO: add missing

    out = temp;

    return StatusCode::ok;
}

// scenario loading and validation
// -----------------------------------------------------------------

// TODO: move some of these functions to a validation api
static bool string_empty(const string& s) { return s.empty(); }

static bool finite_state_att(const ScenarioStateAttConfig& x) {
    return finite_norm_nonzero(x.q, tol12) && finite_vec(x.w);
}

// id and lookups
static bool insert_unique_id(uset<string>& ids, const string& id) {
    if (string_empty(id)) return false;
    return ids.insert(id).second;
}

static const ScenarioGravityProviderConfig* find_gravity_provider_config(
    const ScenarioConfig& cfg,
    const string& id
) {
    for (const auto& provider : cfg.gravity_providers) {
        if (provider.id == id) return &provider;
    }

    return nullptr;
}
static const ScenarioCelestialConfig* find_celestial_config(
    const ScenarioConfig& cfg,
    const string& id
) {
    for (const auto& cel : cfg.celestials) {
        if (cel.id == id) return &cel;
    }

    return nullptr;
}
static const ScenarioSatelliteConfig* find_satellite_config(
    const ScenarioConfig& cfg,
    const string& id
) {
    for (const auto& sat : cfg.satellites) {
        if (sat.id == id) return &sat;
    }

    return nullptr;
}
static const ScenarioStationConfig* find_station_config(
    const ScenarioConfig& cfg,
    const string& id
) {
    for (const auto& stat : cfg.stations) {
        if (stat.id == id) return &stat;
    }

    return nullptr;
}

static bool find_entity_id(
    const umap<string, EntityId>& ids,
    const string& name,
    EntityId& out
) {
    auto it = ids.find(name);
    if (it == ids.end()) {
        out = kInvalidEntityId;
        return false;
    }

    out = it->second;
    return true;
}

// validation

static StatusCode validate_time_config(const ScenarioTimeConfig& time) {
    switch (time.date_type) {
    case DateType::cal: {
        if (!validate_cal(time.cal)) return StatusCode::invalid_input;
    } break;
    case DateType::jd: {
        if (!finite_nonneg(time.jd.frac)) return StatusCode::invalid_input;
    } break;
    case DateType::mjd: {
        if (!finite_nonneg(time.mjd.frac)) return StatusCode::invalid_input;
    } break;
    }

    return StatusCode::ok;
}

static StatusCode validate_state_tr_config(
    const ScenarioConfig& cfg,
    const ScenarioStateTrConfig& x
) {
    switch (x.input_type) {
    case StateTrInputType::pos_vel: {
        if (!finite_vec(x.r) || !finite_vec(x.v)) return StatusCode::invalid_state;
    } break;
    case StateTrInputType::classical: {
        if (string_empty(x.central)) return StatusCode::missing_reference;

        const auto* central = find_celestial_config(cfg, x.central);
        if (!central) return StatusCode::missing_reference;
        if (central->model.id == "custom"
            && !finite_pos(central->model.gravity_model.mu)) {
            return StatusCode::invalid_state;
        }

        if (!finite_pos(x.coes.sma)) {
            return StatusCode::invalid_state;
        }
        if (!finite_inrange(x.coes.ecc, 0.0, 1.0, true, false)) {
            return StatusCode::invalid_state;
        }
        if (!std::isfinite(x.coes.inc) || !std::isfinite(x.coes.raan)
            || !std::isfinite(x.coes.aop) || !std::isfinite(x.coes.ta)) {
            return StatusCode::invalid_state;
        }
    } break;
    default: return StatusCode::unsupported_type;
    }

    return StatusCode::ok;
}

static bool validate_covariance(const matXd& cov, const i32 dim, f64 tol = tol12) {
    if (dim <= 0) return false;
    if (cov.rows() != dim || cov.cols() != dim) return false;

    if (!finite_nonempty_mat(cov)) return false;
    for (i32 i = 0; i < dim; ++i) {
        if (!finite_nonneg(cov(i, i))) return false;
        for (i32 j = i + 1; j < dim; ++j) {
            if (std::abs(cov(i, j) - cov(j, i)) >= tol) return false;
        }
    }

    return true;
}

static StatusCode validate_covariance_config(
    const ScenarioCovarianceConfig& cov,
    const i32 dim,
    f64 tol = tol12
) {
    if (dim <= 0) return StatusCode::unsupported_type;

    if (!validate_covariance(cov.covariance, dim, tol))
        return StatusCode::invalid_covariance;

    return StatusCode::ok;
}

static bool validate_inertia(const mat3d& inertia, f64 tol = tol12) {
    for (i32 i = 0; i < 3; ++i) {
        if (!finite_pos(inertia(i, i))) return false;
    }
    if (!finite_nonzero(inertia.determinant(), tol)) return false;

    return true;
}

static StatusCode validate_mass_properties(
    const ScenarioMassPropertiesConfig& mp,
    const bool required,
    f64 tol = tol12
) {
    if (!required && mp.mass == 0.0) return StatusCode::ok;
    if (!finite_pos(mp.mass)) return StatusCode::invalid_mass_properties;
    if (!finite_mat(mp.inertia)) return StatusCode::invalid_mass_properties;

    if (!validate_inertia(mp.inertia, tol)) return StatusCode::invalid_mass_properties;

    return StatusCode::ok;
}

static bool supported_gravity_provider_format(
    const ScenarioGravityProviderConfig& provider,
    const ScenarioGravityConfig& gravity
) {
    const string format = make_lower(provider.format);
    return in_list<string>(format, {"egm", "sha", "nasa_sha", "nasa", "gfc", "icgem"});
}

static bool validate_direct_zonal_coefficients(const ScenarioGravityConfig& gravity) {
    if (!finite_nonempty_vec(gravity.J)) return false;
    if (gravity.degree < 0) return false;
    if (gravity.J.size() <= gravity.degree) return false;

    return true;
}

static StatusCode validate_gravity_config(
    const ScenarioConfig& cfg,
    const ScenarioGravityConfig& gravity
) {
    if (!finite_nonneg(gravity.mu)) return StatusCode::invalid_input;
    if (!finite_nonneg(gravity.radius)) return StatusCode::invalid_input;
    if (gravity.degree < 0 || gravity.order < 0) return StatusCode::invalid_input;

    if (gravity.model == GravityModel::pointmass) {
        if (gravity.coefficient_source == GravityCoefficientSource::none)
            return StatusCode::ok;
        if (gravity.coefficient_source != GravityCoefficientSource::provider) {
            return StatusCode::validation_failed;
        }

        const auto* provider = find_gravity_provider_config(cfg, gravity.coefficients);
        if (!provider) return StatusCode::gravity_model_not_found;
        if (!supported_gravity_provider_format(*provider, gravity)) {
            return StatusCode::unsupported_method;
        }

        return StatusCode::ok;
    }

    if (gravity.model == GravityModel::zonal
        && gravity.coefficient_source == GravityCoefficientSource::direct_zonal) {
        if (!validate_direct_zonal_coefficients(gravity))
            return StatusCode::validation_failed;
        return StatusCode::ok;
    }

    if (gravity.coefficient_source != GravityCoefficientSource::provider) {
        return StatusCode::gravity_model_not_found;
    }

    const auto* provider = find_gravity_provider_config(cfg, gravity.coefficients);
    if (!provider) return StatusCode::gravity_model_not_found;

    if (!supported_gravity_provider_format(*provider, gravity)) {
        return StatusCode::unsupported_method;
    }

    return StatusCode::ok;
}

static bool supported_celestial_models(const string& model_id) {
    string model = make_lower(model_id);
    return in_list<string>(model, {"custom", "wgs84", "iau_moon"});
}
static StatusCode validate_celestial_config(
    const ScenarioConfig& cfg,
    const ScenarioCelestialConfig& cel
) {
    if (string_empty(cel.id)) return StatusCode::invalid_input;
    StatusCode status = validate_state_tr_config(cfg, cel.x_tr);
    if (status != StatusCode::ok) return status;
    if (!finite_state_att(cel.x_att)) return StatusCode::invalid_attitude_state;
    if (string_empty(cel.model.id)) return StatusCode::invalid_input;

    status = validate_gravity_config(cfg, cel.model.gravity_model);
    if (status != StatusCode::ok) return status;

    if (!supported_celestial_models(cel.model.id))
        return StatusCode::celestial_model_not_found;
    if (cel.model.semimajor_axis != 0.0 || cel.model.semiminor_axis != 0.0) {
        if (!finite_pos(cel.model.semimajor_axis) || !finite_pos(cel.model.semiminor_axis)
            || cel.model.semimajor_axis < cel.model.semiminor_axis) {
            return StatusCode::invalid_shape;
        }
    }
    if (!finite_nonneg(cel.model.mean_radius)) {
        return StatusCode::invalid_shape;
    }
    if (!finite_inrange(cel.model.eccentricity, 0.0, 1.0, true, false)) {
        return StatusCode::invalid_shape;
    }
    if (!finite_nonneg(cel.model.flattening)) {
        return StatusCode::invalid_shape;
    }

    return StatusCode::ok;
}

static StatusCode validate_satellite_config(
    const ScenarioConfig& cfg,
    const ScenarioSatelliteConfig& sat
) {
    if (string_empty(sat.id)) return StatusCode::invalid_input;
    StatusCode status = validate_state_tr_config(cfg, sat.x_tr);
    if (status != StatusCode::ok) return status;
    if (!finite_state_att(sat.x_att)) return StatusCode::invalid_attitude_state;

    return validate_mass_properties(sat.mass_properties, sat.propagation.attitude);
}

static StatusCode validate_station_config(
    const ScenarioConfig& cfg,
    const ScenarioStationConfig& stat
) {
    if (string_empty(stat.id)) return StatusCode::invalid_input;

    if (stat.anchored) {
        if (string_empty(stat.anchor)) return StatusCode::invalid_anchor;

        const auto* anchor = find_celestial_config(cfg, stat.anchor);
        if (!anchor) return StatusCode::invalid_anchor;

        if (stat.coordinate_type == "detic_llh") {
            if (!finite_vec(stat.llh_BCBF)) return StatusCode::invalid_anchor;
            if (!finite_inrange(stat.llh_BCBF(0), -90.0, 90.0, true, true)) {
                return StatusCode::invalid_anchor;
            }
        } else if (stat.coordinate_type == "body_fixed") {
            if (!finite_norm_nonzero(stat.r_body, tol12))
                return StatusCode::invalid_anchor;
        } else {
            return StatusCode::unsupported_method;
        }

        // TODO: add other frames when available
        if (stat.local_frame != "ENU") return StatusCode::unsupported_method;
    } else {
        StatusCode status = validate_state_tr_config(cfg, stat.x_tr);
        if (status != StatusCode::ok) return status;
        if (!finite_state_att(stat.x_att)) return StatusCode::invalid_attitude_state;

        status
            = validate_mass_properties(stat.mass_properties, stat.propagation.attitude);
        if (status != StatusCode::ok) return status;
    }

    uset<string> instrument_ids;
    for (const ScenarioInstrumentConfig& instrument : stat.instruments) {
        if (!insert_unique_id(instrument_ids, instrument.id)) {
            return StatusCode::duplicate_id;
        }

        i32 dim = measurement_dim(instrument.type);
        if (dim <= 0) return StatusCode::unsupported_type;

        StatusCode status = validate_covariance_config(instrument.covariance_cfg, dim);
        if (status != StatusCode::ok) return status;
    }

    return StatusCode::ok;
}

StatusCode validate_scenario_config(const ScenarioConfig& cfg) {
    if (cfg.schema.name != "astrolib.scenario") return StatusCode::invalid_input;
    if (cfg.schema.version != 1) return StatusCode::unsupported_method;
    if (string_empty(cfg.metadata.name)) return StatusCode::invalid_input;

    StatusCode status = validate_time_config(cfg.time);
    if (status != StatusCode::ok) return status;

    uset<string> global_ids;
    uset<string> gravity_provider_ids;
    for (const ScenarioGravityProviderConfig& provider : cfg.gravity_providers) {
        if (!insert_unique_id(global_ids, provider.id)) return StatusCode::duplicate_id;
        if (!insert_unique_id(gravity_provider_ids, provider.id)) {
            return StatusCode::duplicate_id;
        }
        if (string_empty(provider.format)) return StatusCode::invalid_input;
        if (string_empty(provider.filepath)) return StatusCode::file_not_found;
        if (provider.lineskips < 0) return StatusCode::invalid_input;
        // TODO: check if file exists
    }

    uset<string> celestial_ids;
    for (const ScenarioCelestialConfig& cel : cfg.celestials) {
        if (!insert_unique_id(global_ids, cel.id)) return StatusCode::duplicate_id;
        if (!insert_unique_id(celestial_ids, cel.id)) return StatusCode::duplicate_id;

        StatusCode status = validate_celestial_config(cfg, cel);
        if (status != StatusCode::ok) return status;
    }

    for (const ScenarioSatelliteConfig& sat : cfg.satellites) {
        if (!insert_unique_id(global_ids, sat.id)) return StatusCode::duplicate_id;

        StatusCode status = validate_satellite_config(cfg, sat);
        if (status != StatusCode::ok) return status;
    }

    for (const ScenarioStationConfig& stat : cfg.stations) {
        if (!insert_unique_id(global_ids, stat.id)) return StatusCode::duplicate_id;

        StatusCode status = validate_station_config(cfg, stat);
        if (status != StatusCode::ok) return status;
    }

    if (global_ids.empty()) return StatusCode::invalid_input;
    if (cfg.world_stepper.substeps < 1) return StatusCode::invalid_input;
    if (cfg.world_stepper.ticks < 1) return StatusCode::invalid_input;
    if (!finite_pos(cfg.world_stepper.dt_scale)) {
        return StatusCode::invalid_input;
    }

    return StatusCode::ok;
}

StatusCode load_scenario_json(const std::string& filepath, ScenarioConfig& out) {
    std::ifstream file(filepath);
    if (!file) {
        return StatusCode::file_open_failed;
    }

    json root;
    try {
        // file >> root;
        root = json::parse(
            file,
            nullptr, // parser callback
            true,    // allow exceptions
            true     // ignore comments
        );
    } catch (const json::parse_error&) {
        return StatusCode::parse_failed;
    } catch (const json::type_error&) {
        return StatusCode::invalid_input;
    }

    StatusCode status;

    ScenarioConfig temp;
    status = parse_scenario_config(root, temp, "root");
    if (status != StatusCode::ok) return status;

    status = validate_scenario_config(temp);
    if (status != StatusCode::ok) return status;

    out = temp;

    return StatusCode::ok;
}

// scenario to world building

static StatusCode apply_state_att_config(
    const ScenarioStateAttConfig& cfg,
    StateAtt& x_att
) {
    StatusCode status;

    switch (cfg.input_type) {
    case AttitudeType::quaternion: {
        x_att.q = cfg.q;
    } break;
    case AttitudeType::dcm: {
        x_att.q = dcm_to_ep(cfg.dcm);
    } break;
    case AttitudeType::axis_angle: {
        x_att.q = axis_angle_to_ep(cfg.axis, cfg.angle, cfg.units_angle);
    } break;
    case AttitudeType::euler_angles: {
        std::array<RotAxis, 3> rots;
        for (i32 i = 0; i < 3; ++i) {
            status = i32_to_rotaxis(cfg.sequence[i], rots[i]);
            if (status != StatusCode::ok) return status;
        }
        mat3d dcm = ea_to_dcm(cfg.angles, rots, cfg.units_angle);

        x_att.q = dcm_to_ep(dcm);
    } break;
    case AttitudeType::crp: {
        x_att.q = crp_to_ep(cfg.axis);
    } break;
    case AttitudeType::mrp: {
        x_att.q = mrp_to_ep(cfg.axis);
    } break;
    default: {
        return StatusCode::attitude_type_not_found;
    }
    }
    if (x_att.q.norm() <= tol12) return StatusCode::invalid_attitude_state;
    x_att.q.normalize();
    x_att.w = cfg.w;

    return StatusCode::ok;
}

static vec3d convert_vec3_length(
    const vec3d& x,
    const ULength length_in,
    const ULength length_out
) {
    if (length_in == length_out) return x;

    vec3d out = x;
    for (i32 i = 0; i < 3; ++i) {
        out(i) = convert_length(out(i), length_in, length_out);
    }

    return out;
}

static vec3d convert_llh_units(
    const vec3d& llh,
    const UAngle angle_in,
    const ULength length_in
) {
    vec3d out = llh;
    out(0) = convert_angle(out(0), angle_in, UAngle::radian);
    out(1) = convert_angle(out(1), angle_in, UAngle::radian);
    out(2) = convert_length(out(2), length_in, ULength::kilometer);

    return out;
}

static StatusCode apply_state_tr_config(
    const ScenarioStateTrConfig& cfg,
    StateTr& x_tr,
    const f64 mu
) {
    StatusCode code;

    switch (cfg.input_type) {
    case StateTrInputType::pos_vel: {
        x_tr.r = convert_vec3_length(cfg.r, cfg.units_length, ULength::kilometer);
        x_tr.v = convert_vec3_length(cfg.v, cfg.units_length, ULength::kilometer);
    } break;
    case StateTrInputType::classical: {
        OEClassical coes = cfg.coes;
        coes.sma = convert_length(coes.sma, cfg.units_length, ULength::kilometer);
        x_tr = classical_to_rv(coes, mu, cfg.units_angle);
        if (!finite_vec(x_tr.r) || !finite_vec(x_tr.v)) return StatusCode::invalid_state;
    } break;
    }

    return StatusCode::ok;
}

static string resolve_project_path(const string& filepath) {
    // TODO: this is temporary, replace with path macros later
    if (filepath.empty()) return filepath;
    if (pwd.empty()) return filepath;

    if (filepath.front() == '/') {
        std::filesystem::path abs_path(filepath);
        if (std::filesystem::exists(abs_path)) return abs_path.string();

        string stripped = filepath.substr(1);
        if (stripped.empty()) return filepath;

        bool root_has_slash = pwd.back() == '/';
        string project_path = root_has_slash ? pwd + stripped : pwd + "/" + stripped;
        if (std::filesystem::exists(project_path)) return project_path;

        return filepath;
    }

    bool root_has_slash = pwd.back() == '/';
    bool path_has_slash = filepath.front() == '/';
    if (root_has_slash && path_has_slash) return pwd + filepath.substr(1);
    if (!root_has_slash && !path_has_slash) return pwd + "/" + filepath;

    return pwd + filepath;
}

static StatusCode apply_celestial_config(
    const ScenarioConfig& scenario,
    const ScenarioCelestialConfig& cfg,
    Celestial& cel
) {
    StatusCode status;
    Celestial temp;

    const auto& model = cfg.model;
    const auto& gravity = model.gravity_model;

    if (model.id != "custom") {
        if (!supported_celestial_models(model.id))
            return StatusCode::celestial_model_not_found;
        if (model.id == "wgs84") {
            temp = wgs84(ULength::kilometer);
        } else if (model.id == "iau_moon") {
            temp = iau_moon(ULength::kilometer);
        }
    } else {
        temp.gravity_model = model.gravity_model.model;
        temp.semimajor_axis = model.semimajor_axis;
        temp.semiminor_axis = model.semiminor_axis;
        temp.mean_radius = model.mean_radius;
        temp.eccentricity = model.eccentricity;
        temp.flattening = model.flattening;
    }

    temp.name = cfg.name; // override

    // gravity
    // fallback priority for reference radius (used for gravity)
    // not based on celestial model
    if (gravity.model != GravityModel::pointmass) {
        if (finite_pos(gravity.radius)) {
            temp.ref_radius = gravity.radius;
        } else if (finite_pos(temp.semimajor_axis)) {
            temp.ref_radius = temp.semimajor_axis;
        } else if (finite_pos(temp.mean_radius)) {
            temp.ref_radius = temp.mean_radius;
        } else {
            return StatusCode::invalid_input;
        }
    }
    // same as above, config > provided model > fail
    if (finite_pos(gravity.mu)) {
        temp.mu = gravity.mu;
    } else if (finite_pos(temp.mu)) {
        // cel.mu = cel.mu;
    } else {
        return StatusCode::invalid_input;
    }
    temp.gravity_model = gravity.model;
    temp.degree = gravity.degree;
    temp.order = gravity.order;
    // gravity.coefficients);
    switch (gravity.coefficient_source) {
    case GravityCoefficientSource::none: break;
    case GravityCoefficientSource::provider: {
        const auto* gravity_provider
            = find_gravity_provider_config(scenario, gravity.coefficients);
        if (gravity_provider) {
            string filepath = resolve_project_path(gravity_provider->filepath);
            if (in_list<string>(gravity_provider->format, {"gfc", "icgem"})) {
                bool read_ok = read_gfc(
                    filepath,
                    temp.C,
                    temp.S,
                    temp.degree,
                    temp.order,
                    gravity_provider->lineskips
                );
                if (!read_ok) return StatusCode::file_open_failed;
            } else if (
                in_list<string>(
                    gravity_provider->format,
                    {"egm", "sha", "nasa", "nasa_sha"}
                )
            ) {
                bool read_ok = read_sphh_coefs(
                    filepath,
                    temp.C,
                    temp.S,
                    temp.degree,
                    temp.order,
                    gravity_provider->lineskips
                );
                if (!read_ok) return StatusCode::file_open_failed;
            } else {
                return StatusCode::gravity_model_not_found;
            }
        } else {
            return StatusCode::gravity_model_not_found;
        }
    } break;
    case GravityCoefficientSource::direct_zonal: {
        auto n = std::min(temp.J.size(), gravity.J.size());
        for (i32 i = 0; i < n; ++i) {
            temp.J(i) = gravity.J(i);
        }
    } break;
    case GravityCoefficientSource::direct_spherical_harmonics:
        break;
        // TODO: do direct C an S loading
    }

    status = apply_state_att_config(cfg.x_att, temp.x_att);
    if (status != StatusCode::ok) return status;
    if (cfg.has_attitude_model) {
        temp.attitude_model = cfg.attitude_model;
    }
    temp.set_spin_rate(temp.x_att.w.norm());

    // build state_tr
    status = apply_state_tr_config(cfg.x_tr, temp.x_tr, temp.mu);
    if (status != StatusCode::ok) return status;

    temp.propagate_att = cfg.propagation.attitude;
    temp.propagate_tr = cfg.propagation.translation;

    cel = temp;
    return StatusCode::ok;
}

static StatusCode apply_mass_properties_config(
    const ScenarioMassPropertiesConfig& cfg,
    MassProperties& mp
) {
    mp.mass = cfg.mass;
    mp.I = cfg.inertia;

    mp.principal_axes = cfg.principle_axes;

    if (cfg.offset) {
        mp.I = inertia_PAT(mp.I, mp.mass, mp.offset_body);
        if (!mp.I.isDiagonal()) mp.principal_axes = false;
        mp.offset_body = vec3d0;
    }

    mp.I_inv = mp.I.inverse();
    if (!finite_mat(mp.I_inv)) return StatusCode::matrix_invert_failed;

    mp.active = true;

    return StatusCode::ok;
}

static StatusCode resolve_central_celestial(
    const World& world,
    const umap<string, EntityId>& celestial_ids,
    const string& central_name,
    const Celestial*& central
) {
    central = nullptr;

    if (central_name.empty()) {
        return StatusCode::ok;
    }

    auto it = celestial_ids.find(central_name);
    if (it == celestial_ids.end()) {
        return StatusCode::missing_reference;
    }

    central = world.celestial(it->second);
    if (central == nullptr) {
        return StatusCode::missing_reference;
    }

    if (!finite_pos(central->mu)) {
        return StatusCode::invalid_state;
    }
    if (!finite_vec(central->x_tr.r) || !finite_vec(central->x_tr.v))
        return StatusCode::invalid_state;

    return StatusCode::ok;
}

static StatusCode apply_satellite_config(
    const ScenarioConfig& scenario,
    const ScenarioSatelliteConfig& cfg,
    const umap<string, EntityId>& cel_ids,
    const World& world,
    Satellite& sat
) {
    StatusCode status;
    Satellite temp;

    const auto& mp = cfg.mass_properties;

    temp.name = cfg.name;

    status = apply_mass_properties_config(cfg.mass_properties, temp.mass_properties);
    if (status != StatusCode::ok) return status;

    status = apply_state_att_config(cfg.x_att, temp.x_att);
    if (status != StatusCode::ok) return status;

    const Celestial* cel;
    status = resolve_central_celestial(world, cel_ids, cfg.x_tr.central, cel);
    if (status != StatusCode::ok) return status;
    if (cfg.x_tr.central.empty()) {
        status = apply_state_tr_config(cfg.x_tr, temp.x_tr, 0.0);
        if (status != StatusCode::ok) return status;
    } else {
        status = apply_state_tr_config(cfg.x_tr, temp.x_tr, cel->mu);
        if (status != StatusCode::ok) return status;
        if (cfg.x_tr.input_type == StateTrInputType::classical) {
            temp.x_tr += cel->x_tr;
        }
    }

    temp.propagate_att = cfg.propagation.attitude;
    temp.propagate_tr = cfg.propagation.translation;

    sat = temp;
    return StatusCode::ok;
}

static StatusCode apply_instrument_config(
    const ScenarioInstrumentConfig& cfg,
    StationInstrument& instrument
) {
    instrument.name = cfg.id;
    instrument.type = cfg.type;
    instrument.R = cfg.covariance_cfg.covariance;
    instrument.enabled = cfg.enabled;

    return StatusCode::ok;
}

static StatusCode apply_station_config(
    const ScenarioConfig& scenario,
    const ScenarioStationConfig& cfg,
    const umap<string, EntityId>& cel_ids,
    const World& world,
    Station& stat
) {
    StatusCode status;
    Station temp;

    temp.name = cfg.name;

    temp.anchored = cfg.anchored;
    if (cfg.anchored) {
        auto it = cel_ids.find(cfg.anchor);
        if (it == cel_ids.end()) return StatusCode::invalid_anchor;
        temp.anchor_id = it->second;

        if (cfg.coordinate_type == "detic_llh") {
            temp.llh_BCBF
                = convert_llh_units(cfg.llh_BCBF, cfg.units_angle, cfg.units_length);
        } else if (cfg.coordinate_type == "body_fixed") {
            temp.r_body_BCBF
                = convert_vec3_length(cfg.r_body, cfg.units_length, ULength::kilometer);
            temp.llh_BCBF = vec3d0;
        } else {
            return StatusCode::unsupported_method;
        }

        temp.propagate_att = false;
        temp.propagate_tr = false;
    } else {
        status = apply_mass_properties_config(cfg.mass_properties, temp.mass_properties);
        if (status != StatusCode::ok) return status;

        status = apply_state_att_config(cfg.x_att, temp.x_att);
        if (status != StatusCode::ok) return status;

        const Celestial* cel;
        status = resolve_central_celestial(world, cel_ids, cfg.x_tr.central, cel);
        if (status != StatusCode::ok) return status;
        if (cfg.x_tr.central.empty()) {
            status = apply_state_tr_config(cfg.x_tr, temp.x_tr, 0.0);
            if (status != StatusCode::ok) return status;
        } else {
            status = apply_state_tr_config(cfg.x_tr, temp.x_tr, cel->mu);
            if (status != StatusCode::ok) return status;
            if (cfg.x_tr.input_type == StateTrInputType::classical) {
                temp.x_tr += cel->x_tr;
            }
        }

        temp.propagate_att = cfg.propagation.attitude;
        temp.propagate_tr = cfg.propagation.translation;
    }

    for (const auto& instrument_cfg : cfg.instruments) {
        StationInstrument instrument;
        status = apply_instrument_config(instrument_cfg, instrument);
        if (status != StatusCode::ok) return status;

        InstrumentId instrument_id;
        status = add_station_instrument(temp, instrument, instrument_id);
        if (status != StatusCode::ok) return status;
    }

    stat = temp;
    return StatusCode::ok;
}

StatusCode build_world_from_scenario_config(
    const ScenarioConfig& cfg,
    World& world,
    ScenarioBuildResult& result,
    WorldStepperConfig& stepper
) {
    StatusCode status;
    result = ScenarioBuildResult{};

    stepper.integrator_att = cfg.world_stepper.integrator_att;
    stepper.integrator_tr = cfg.world_stepper.integrator_tr;
    stepper.dt_scale = cfg.world_stepper.dt_scale;
    stepper.substeps = cfg.world_stepper.substeps;
    stepper.ticks = cfg.world_stepper.ticks;
    stepper.paused = cfg.world_stepper.paused;

    for (const auto& cel_cfg : cfg.celestials) {
        auto cel = std::make_unique<Celestial>();
        status = apply_celestial_config(cfg, cel_cfg, *cel);
        if (status != StatusCode::ok) return status;

        EntityId id = world.insert_celestial(std::move(cel));
        result.celestial_ids.insert({cel_cfg.id, id});
        result.body_ids.insert({cel_cfg.id, id});
    }

    for (const auto& sat_cfg : cfg.satellites) {
        auto sat = std::make_unique<Satellite>();
        status = apply_satellite_config(cfg, sat_cfg, result.celestial_ids, world, *sat);
        if (status != StatusCode::ok) return status;

        EntityId id = world.insert_satellite(std::move(sat));
        result.satellite_ids.insert({sat_cfg.id, id});
        result.body_ids.insert({sat_cfg.id, id});
    }

    for (const auto& stat_cfg : cfg.stations) {
        auto stat = std::make_unique<Station>();
        status = apply_station_config(cfg, stat_cfg, result.celestial_ids, world, *stat);
        if (status != StatusCode::ok) return status;

        EntityId id = world.insert_station(std::move(stat));
        if (stat_cfg.anchored) {
            auto* stat = world.station(id);
            if (stat == nullptr) return StatusCode::body_not_found;
            auto* cel = world.celestial(stat->anchor_id);
            if (cel == nullptr) return StatusCode::invalid_anchor;

            if (stat_cfg.coordinate_type == "detic_llh") {
                bool set_ok = world.set_stat_anchor_detic(
                    id,
                    stat->anchor_id,
                    stat->llh_BCBF,
                    UAngle::radian
                );
                if (!set_ok) return StatusCode::invalid_anchor;
            } else if (stat_cfg.coordinate_type == "body_fixed") {
                stat->llh_BCBF = bcbf_to_detic(stat->r_body_BCBF, *cel, UAngle::radian);
            } else {
                return StatusCode::unsupported_method;
            }
        }
        result.station_ids.insert({stat_cfg.id, id});
        result.body_ids.insert({stat_cfg.id, id});
    }

    return StatusCode::ok;
}

StatusCode build_scenario_config_from_world(
    ScenarioConfig& cfg,
    const World& world,
    ScenarioBuildResult& result,
    const WorldStepperConfig& stepper
) {
    cfg = ScenarioConfig{};
    ScenarioWorldStepperConfig& cfg_stepper = cfg.world_stepper;

    // TODO: some field are not filled out in the cfg, fill out later
    // TODO: add options later

    cfg_stepper.integrator_tr = stepper.integrator_tr;
    cfg_stepper.integrator_att = stepper.integrator_att;
    cfg_stepper.paused = stepper.paused;
    cfg_stepper.substeps = stepper.substeps;
    cfg_stepper.ticks = stepper.ticks;
    cfg_stepper.dt_scale = stepper.dt_scale;

    cfg.time.t0 = world.t_sim();
    // world.

    svec<ScenarioCelestialConfig>& celestials = cfg.celestials;
    for (const EntityId id : world.celestial_ids()) {
        const Celestial* cel = world.celestial(id);
        if (cel == nullptr) return StatusCode::body_not_found;

        ScenarioCelestialConfig out;
        out.id = cel->name;

        out.attitude_model = cel->attitude_model;

        out.x_tr.input_type = StateTrInputType::pos_vel;
        out.x_tr.r = cel->x_tr.r;
        out.x_tr.v = cel->x_tr.v;

        out.x_att.input_type = AttitudeType::quaternion;
        out.x_att.q = cel->x_att.q;
        out.x_att.w = cel->x_att.w;

        out.propagation.translation = cel->propagate_tr;
        out.propagation.attitude = cel->propagate_att;

        out.model.id = "custom";
        out.model.semimajor_axis = cel->semimajor_axis;
        out.model.semiminor_axis = cel->semiminor_axis;
        out.model.mean_radius = cel->mean_radius;
        out.model.eccentricity = cel->eccentricity;
        out.model.flattening = cel->flattening;

        ScenarioGravityConfig& gravity = out.model.gravity_model;
        gravity.model = cel->gravity_model;
        gravity.mu = cel->mu;
        gravity.radius = cel->ref_radius;
        gravity.degree = cel->degree;
        gravity.order = cel->order;
        gravity.J = cel->J;
        gravity.C = cel->C;
        gravity.S = cel->S;
        switch (gravity.model) {
        case GravityModel::pointmass:
            gravity.coefficient_source = GravityCoefficientSource::none;
        case GravityModel::zonal:
            gravity.coefficient_source = GravityCoefficientSource::direct_zonal;
        case GravityModel::spherical_harmonics:
            gravity.coefficient_source
                = GravityCoefficientSource::direct_spherical_harmonics;
        }

        celestials.push_back(out);
        // TODO: figure out how i should do this, maybe create new json file with coefs?
    }

    svec<ScenarioSatelliteConfig>& satellites = cfg.satellites;
    for (const EntityId id : world.satellite_ids()) {
        const Satellite* sat = world.satellite(id);
        if (sat == nullptr) return StatusCode::body_not_found;

        ScenarioSatelliteConfig out;
        out.id = sat->name;

        out.propagation.translation = sat->propagate_tr;
        out.propagation.attitude = sat->propagate_att;

        out.x_tr.input_type = StateTrInputType::pos_vel;
        out.x_tr.r = sat->x_tr.r;
        out.x_tr.v = sat->x_tr.v;

        out.x_att.input_type = AttitudeType::quaternion;
        out.x_att.q = sat->x_att.q;
        out.x_att.w = sat->x_att.w;

        out.mass_properties.mass = sat->mass_properties.mass;
        out.mass_properties.inertia = sat->mass_properties.I;
        out.mass_properties.principle_axes = sat->mass_properties.principal_axes;
        out.mass_properties.offset_body = sat->mass_properties.offset_body;

        satellites.push_back(out);
    }

    svec<ScenarioStationConfig>& stations = cfg.stations;
    for (const EntityId id : world.station_ids()) {
        const Station* stat = world.station(id);
        if (stat == nullptr) return StatusCode::body_not_found;

        ScenarioStationConfig out;
        out.id = stat->name;

        out.anchored = stat->anchored;
        if (out.anchored) {
            const Celestial* cel = world.celestial(stat->anchor_id);
            if (cel == nullptr) return StatusCode::invalid_anchor;
            out.anchor = cel->name;
            out.coordinate_type = "detic_llh";
            out.llh_BCBF = stat->llh_BCBF;
            out.local_frame = "ENU";
            out.units_angle = UAngle::radian;
        } else {
            out.propagation.translation = stat->propagate_tr;
            out.propagation.attitude = stat->propagate_att;

            out.x_tr.input_type = StateTrInputType::pos_vel;
            out.x_tr.r = stat->x_tr.r;
            out.x_tr.v = stat->x_tr.v;

            out.x_att.input_type = AttitudeType::quaternion;
            out.x_att.q = stat->x_att.q;
            out.x_att.w = stat->x_att.w;

            out.mass_properties.mass = stat->mass_properties.mass;
            out.mass_properties.inertia = stat->mass_properties.I;
            out.mass_properties.principle_axes = stat->mass_properties.principal_axes;
            out.mass_properties.offset_body = stat->mass_properties.offset_body;
        }

        for (const auto& [id, instr] : stat->instruments) {
            ScenarioInstrumentConfig out;

            out.id = instr.name;
            out.type = instr.type;
            out.covariance_cfg.covariance = instr.R;
            out.covariance_cfg.type = "matrix";
            out.enabled = instr.enabled;
        }
    }

    return StatusCode::ok;
}