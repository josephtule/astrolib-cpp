#include "core/scenario_io.hpp"
#include "core/body.hpp"
#include "core/estimation_common.hpp"
#include "core/integrator.hpp"
#include "core/measurement.hpp"
#include "core/observation_type.hpp"
#include "core/state.hpp"
#include "core/transform.hpp"

#include "util/constants.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include "nlohmann/json.hpp"

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

static StatusCode parse_units_angle(const string& str, UAngle& out) {
    if (str == "radian") {
        out = UAngle::radian;
    } else if (str == "degree") {
        out = UAngle::degree;
    } else if (str == "arcminute") {
        out = UAngle::arcminute;
    } else if (str == "arcsecond") {
        out = UAngle::arcsecond;
    } else if (str == "milliarcsecond") {
        out = UAngle::milliarcsecond;
    } else {
        return StatusCode::invalid_input;
    }

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
    if (!isfinite(value)) return StatusCode::invalid_input;

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
    if (!isfinite(value)) return StatusCode::invalid_input;

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
    if (!found) {
        out = AttitudeType::quaternion; // default
        return StatusCode::ok;
    }

    return parse_attitude_type(type_str, out);
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
        if (!isfinite(temp(i))) {
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
        if (!isfinite(temp(i))) {
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
        if (!isfinite(temp(i / 3, i % 3))) return StatusCode::invalid_input;
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
        if (!isfinite(temp(i / m, i % m))) return StatusCode::invalid_input;
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
        if (!isfinite(temp(i / M, i % M))) return StatusCode::invalid_input;
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
            if (!isfinite(temp(i, j))) return StatusCode::invalid_input;
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
            if (!isfinite(temp(i, j))) return StatusCode::invalid_input;
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
        if (!isfinite(temp(i))) return StatusCode::invalid_input;
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
        if (!isfinite(temp[i])) return StatusCode::invalid_input;
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

    status = parse_req_vec3d(*state_tr, "r", out.r, path);
    if (status != StatusCode::ok) return status;

    status = parse_req_vec3d(*state_tr, "v", out.v, path);
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

    if (out.units_length != ULength::kilometer) {
        for (i32 i = 0; i < 3; ++i) {
            out.r(i) = convert_length(out.r(i), out.units_length, ULength::kilometer);
            out.v(i) = convert_length(out.v(i), out.units_length, ULength::kilometer);
        }
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
    const json* state_att = nullptr;
    status = get_req_child(object, "state_att", state_att, path);
    if (status != StatusCode::ok) return status;

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
    status = parse_opt_attitude_type(*state_att, "type", found_type, out.type, path);
    if (status != StatusCode::ok) return status;

    if (found_type) {
        switch (out.type) {
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
            string temp_units_angle = "";
            status = parse_req_string(*state_att, "units_angle", temp_units_angle, path);
            if (status != StatusCode::ok) return status;
            status = parse_units_angle(temp_units_angle, out.units_angle);
            if (status != StatusCode::ok) return status;

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

    status = parse_req_string(object, "type", out.type, path);
    if (status != StatusCode::ok) return status;

    status = parse_req_string(object, "path", out.filepath, path);
    if (status != StatusCode::ok) return status;

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

    bool _ = false;

    status = parse_opt_f64(*gravity, "radius", _, out.radius, path);
    if (status != StatusCode::ok) return status;
    // TODO: use from celestial model or determine where to get reference radius (file or
    // computed)

    if (out.model == GravityModel::pointmass) {
        status = parse_opt_i32(*gravity, "degree", _, out.degree, path);
        if (status != StatusCode::ok) return status;

        status = parse_opt_i32(*gravity, "order", _, out.order, path);
        if (status != StatusCode::ok) return status;

        status = parse_opt_string(*gravity, "coefficients", _, out.coefficients, path);
        if (status != StatusCode::ok) return status;
    } else {
        status = parse_req_i32(*gravity, "degree", out.degree, path);
        if (status != StatusCode::ok) return status;

        status = parse_req_i32(*gravity, "order", out.order, path);
        if (status != StatusCode::ok) return status;

        status = parse_req_string(*gravity, "coefficients", out.coefficients, path);
        if (status != StatusCode::ok) return status;
    }

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

        status = parse_gravity_config(*model, out.gravity_model, path);
        if (status != StatusCode::ok) return status;

        status = parse_req_f64(*model, "semimajor_axis", out.semimajor_axis, path);
        if (status != StatusCode::ok) return status;

        status = parse_req_f64(*model, "semiminor_axis", out.semiminor_axis, path);
        if (status != StatusCode::ok) return status;

        bool found;
        status = parse_opt_f64(*model, "mean_radius", found, out.mean_radius, path);
        if (status != StatusCode::ok) return status;
        if (!found) {
            out.mean_radius = std::pow(
                out.semimajor_axis * out.semimajor_axis * out.semiminor_axis,
                1.0 / 3.0
            );
        }

        status = parse_opt_f64(*model, "eccentricity", found, out.eccentricity, path);
        if (status != StatusCode::ok) return status;
        if (!found) {
            out.eccentricity = std::sqrt(
                1.0
                - out.semiminor_axis * out.semiminor_axis
                      / (out.semimajor_axis * out.semimajor_axis)
            );
        }

        status = parse_opt_f64(*model, "flattening", found, out.flattening, path);
        if (status != StatusCode::ok) return status;
        if (!found) {
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

            if (out.units_angle != UAngle::degree) {
                out.llh_BCBF(0)
                    = convert_angle(out.llh_BCBF(0), out.units_angle, UAngle::degree);
                out.llh_BCBF(1)
                    = convert_angle(out.llh_BCBF(1), out.units_angle, UAngle::degree);
            }
            if (out.units_length != ULength::kilometer) {
                out.llh_BCBF(
                    2
                ) = convert_length(out.llh_BCBF(2), out.units_length, ULength::kilometer);
            }
        } else if (out.coordinate_type == "body_fixed") {
            status = parse_req_vec3d(object, "r_body", out.r_body, path);
            if (status != StatusCode::ok) return status;

            if (out.units_length != ULength::kilometer) {
                for (i32 i = 0; i < 3; ++i) {
                    out.r_body(i) = convert_length(
                        out.r_body(i),
                        out.units_length,
                        ULength::kilometer
                    );
                }
            }
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

    status = parse_opt_f64(*child, "time_scale", temp_found, out.time_scale, path);
    if (status != StatusCode::ok) return status;
    if (!temp_found) out.time_scale = 1.0;

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
// Scenario Loading
// --------------------------------------------------------------------------------

StatusCode validate_scenario_config(const ScenarioConfig& cfg) {
    // TODO: complete this
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
