#include "core/scenario_io.hpp"
#include "core/estimation_common.hpp"

#include "core/state.hpp"
#include "core/transform.hpp"
#include "nlohmann/json.hpp"
#include "util/constants.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

using json = nlohmann::json;

// TODO: add json specific status codes later

static StatusCode get_required_child(
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

static StatusCode get_optional_child(
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

static StatusCode parse_required_string(
    const json& object,
    const string& key,
    string& out,
    const string& path
) {
    const json* child = nullptr;

    StatusCode status = get_required_child(object, key, child, path);
    if (status != StatusCode::ok) {
        return status;
    }

    if (!child->is_string()) {
        return StatusCode::invalid_input;
    }

    out = child->get<string>();

    return StatusCode::ok;
}

static StatusCode parse_optional_string(
    const json& object,
    const string& key,
    bool found,
    string& out,
    const string& path
) {
    found = false;
    const json* child = nullptr;

    StatusCode status = get_optional_child(object, key, child, found, path);
    if (status != StatusCode::ok) {
        return status;
    }

    if (!child->is_string()) {
        return StatusCode::invalid_input;
    }

    out = child->get<string>();

    return StatusCode::ok;
}

static StatusCode parse_required_f64(
    const json& object,
    const string& key,
    f64& out,
    const string& path
) {
    const json* child = nullptr;

    StatusCode status = get_required_child(object, key, child, path);
    if (status != StatusCode::ok) {
        return status;
    }

    if (!child->is_number()) {
        return StatusCode::invalid_input;
    }

    f64 value = child->get<f64>();
    if (!isfinite(value)) {
        return StatusCode::invalid_input;
    }
    out = value;

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

static StatusCode parse_mat3d(const json& value, mat3d& out, const string& path) {
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

template <int N, int M>
static StatusCode parse_matNMd(
    const json& value,
    mat<f64, N, M>& out,
    const string& path
) {
    if (!value.is_array() || value.size() != N * M) return StatusCode::invalid_input;

    mat<f64, N, M> temp;
    for (i32 i = 0; i < N * M; ++i) {
        if (!value.at(i).is_number()) return StatusCode::invalid_input;
        temp(i / M, i % N) = value[i].get<f64>();
        if (!isfinite(temp(i / M, i % N))) return StatusCode::invalid_input;
    }
    out = temp;

    return StatusCode::ok;
}

static StatusCode parse_vecXd(const json& value, vecXd& out, const string& path) {
    if (!value.is_array()) return StatusCode::invalid_input;

    vecXd temp;
    i32 n = value.size();
    for (i32 i = 0; i < n; ++i) {
        if (!value.at(i).is_number()) return StatusCode::invalid_input;
        temp(i) = value[i].get<f64>();
        if (!isfinite(temp(i))) return StatusCode::invalid_input;
    }

    out = temp;

    return StatusCode::ok;
}

template <typename T, int N>
static StatusCode parse_array(
    const json& value,
    std::array<T, N>& out,
    const string& path
) {
    if (!value.is_array()) return StatusCode::invalid_input;

    std::array<T, N> temp;
    for (i32 i = 0; i < N; ++i) {
        if (!value.at(i).is_number()) return StatusCode::invalid_input;
        temp[i] = value[i].get<f64>();
        if (!isfinite(temp[i])) return StatusCode::invalid_input;
    }

    out = temp;
    return StatusCode::ok;
}

template <typename T, int N>
static StatusCode parse_required_array(
    const json& object,
    const std::string& key,
    std::array<T, N>& out,
    const std::string& path
) {
    const json* child = nullptr;
    StatusCode status = get_required_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_array<T, N>(*child, out, path + "." + key);
}

static StatusCode parse_required_mat3d(
    const json& object,
    const string& key,
    mat3d& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_required_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_mat3d(*child, out, path + "." + key);
}

template <int N, int M>
static StatusCode parse_required_matNMd(
    const json& object,
    const string& key,
    mat<f64, N, M>& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_required_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_matNMd<N, M>(*child, out, path + "." + key);
}

static StatusCode parse_required_vec3d(
    const json& object,
    const std::string& key,
    vec3d& out,
    const std::string& path
) {
    const json* child = nullptr;
    StatusCode status = get_required_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_vec3d(*child, out, path + "." + key);
}

static StatusCode parse_optional_vec3d(
    const json& object,
    const string& key,
    bool found,
    vec3d& out,
    const string& path
) {
    found = false;
    const json* child = nullptr;

    StatusCode status = get_optional_child(object, key, child, found, path);
    if (status != StatusCode::ok) {
        return status;
    }

    return parse_vec3d(*child, out, path + "." + key);
}

static StatusCode parse_required_vec4d(
    const json& object,
    const string& key,
    vec4d& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_required_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_vec4d(*child, out, path + "." + key);
}

static StatusCode parse_required_vecXd(
    const json& object,
    const string& key,
    vecXd& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_required_child(object, key, child, path);
    if (status != StatusCode::ok) return status;

    return parse_vecXd(*child, out, path + "." + key);
}

static StatusCode parse_state_tr(const json& object, StateTr& out, const string& path) {
    if (!object.is_object()) return StatusCode::invalid_input;

    StateTr temp;
    StatusCode status = parse_required_vec3d(object, "r", temp.r, path);
    if (status != StatusCode::ok) return status;

    status = parse_required_vec3d(object, "v", temp.v, path);
    if (status != StatusCode::ok) return status;

    out = temp;

    return StatusCode::ok;
}

static StatusCode parse_state_att(
    const json& object,
    ScenarioStateAttConfig& out,
    const string& path
) {
    if (!object.is_object()) return StatusCode::invalid_input;

    ScenarioStateAttConfig temp;
    StatusCode status = StatusCode::ok;
    if (temp.type == "quaternion") {
        status = parse_required_vec4d(object, "q", temp.q, path);
        if (status != StatusCode::ok) return status;

        if (temp.q.norm() <= tol12) return StatusCode::invalid_input;
        temp.q.normalize();
    } else if (temp.type == "axis_angle") {
        string temp_units_angle = "";
        status = parse_required_string(object, "units_angle", temp_units_angle, path);
        if (status != StatusCode::ok) return status;
        status = string_to_uangle(temp_units_angle, temp.units_angle);
        if (status != StatusCode::ok) return status;

        status = parse_required_vec3d(object, "axis", temp.axis, path);
        if (status != StatusCode::ok) return status;

        status = parse_required_f64(object, "angle", temp.angle, path);
        if (status != StatusCode::ok) return status;
        temp.axis.normalize();

        temp.q = quat_from_axis_angle_passive(temp.axis, temp.angle, temp.units_angle);
    } else if (temp.type == "euler_angles") {
        string temp_units_angle = "";
        status = parse_required_string(object, "units_angle", temp_units_angle, path);
        if (status != StatusCode::ok) return status;
        status = string_to_uangle(temp_units_angle, temp.units_angle);
        if (status != StatusCode::ok) return status;

        status = parse_required_vec3d(object, "angles", temp.angles, path);
        if (status != StatusCode::ok) return status;

        status = parse_required_array<i32, 3>(object, "sequence", temp.sequence, path);
        if (status != StatusCode::ok) return status;

        std::array<RotAxis, 3> rots;
        for (i32 i = 0; i < 3; ++i) {
            status = i32_to_rotaxis(temp.sequence[i], rots[i]);
            if (status != StatusCode::ok) return status;
        }
        mat3d dcm = ea_to_dcm(temp.angles, rots, temp.units_angle);

        temp.q = dcm_to_ep(dcm);
    } else if (temp.type == "dcm") {
        status = parse_required_matNMd(object, "dcm", temp.dcm, path);
        if (status != StatusCode::ok) return status;
        temp.q = dcm_to_ep(temp.dcm);
    } else if (temp.type == "crp") {
        status = parse_required_vec3d(object, "axis", temp.axis, path);
        if (status != StatusCode::ok) return status;
        temp.q = crp_to_ep(temp.axis);
    } else if (temp.type == "mrp") {
        status = parse_required_vec3d(object, "axis", temp.axis, path);
        if (status != StatusCode::ok) return status;
        temp.q = mrp_to_dcm(temp.axis);
    } else {
        return StatusCode::invalid_input;
    }

    bool w_found;
    status = parse_optional_vec3d(object, "w", w_found, temp.w, path);
    if (status != StatusCode::ok) return status;

    out = temp;

    return StatusCode::ok;
}

StatusCode load_scenario_json(const std::string& filepath, ScenarioConfig& out) {
    return StatusCode::ok;
}

StatusCode validate_scenario_config(const ScenarioConfig& cfg) { return StatusCode::ok; }
