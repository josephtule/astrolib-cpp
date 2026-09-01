// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/status.hpp"
#include "core/time.hpp"

#include "util/typedefs.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

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

// shared enum parsers

static StatusCode parse_units_angle(const string& str, UAngle& out) {
    UAngle temp;
    if (str == uangle_str(UAngle::radian))
        temp = UAngle::radian;
    else if (str == uangle_str(UAngle::degree))
        temp = UAngle::degree;
    else if (str == uangle_str(UAngle::arcminute))
        temp = UAngle::arcminute;
    else if (str == uangle_str(UAngle::arcsecond) || str == "arcseconds")
        temp = UAngle::arcsecond;
    else if (str == uangle_str(UAngle::milliarcsecond))
        temp = UAngle::milliarcsecond;
    else
        return StatusCode::unsupported_type;

    out = temp;
    return StatusCode::ok;
}

static StatusCode parse_req_units_angle(
    const json& object,
    const string& key,
    UAngle& out,
    const string& path
) {
    string value;
    StatusCode status = parse_req_string(object, key, value, path);
    if (status != StatusCode::ok) return status;

    return parse_units_angle(value, out);
}

static StatusCode parse_opt_units_angle(
    const json& object,
    const string& key,
    bool& found,
    UAngle& out,
    const string& path
) {
    string value;
    StatusCode status = parse_opt_string(object, key, found, value, path);
    if (status != StatusCode::ok || !found) return status;

    return parse_units_angle(value, out);
}

static StatusCode parse_units_length(const string& str, ULength& out) {
    ULength temp;
    if (str == ulength_str(ULength::nanometer))
        temp = ULength::nanometer;
    else if (str == ulength_str(ULength::millimeter))
        temp = ULength::millimeter;
    else if (str == ulength_str(ULength::centimeter))
        temp = ULength::centimeter;
    else if (str == ulength_str(ULength::meter))
        temp = ULength::meter;
    else if (str == ulength_str(ULength::kilometer))
        temp = ULength::kilometer;
    else if (str == ulength_str(ULength::inch))
        temp = ULength::inch;
    else if (str == ulength_str(ULength::foot))
        temp = ULength::foot;
    else if (str == ulength_str(ULength::mile))
        temp = ULength::mile;
    else if (str == ulength_str(ULength::au))
        temp = ULength::au;
    else
        return StatusCode::unsupported_type;

    out = temp;
    return StatusCode::ok;
}

static StatusCode parse_req_units_length(
    const json& object,
    const string& key,
    ULength& out,
    const string& path
) {
    string value;
    StatusCode status = parse_req_string(object, key, value, path);
    if (status != StatusCode::ok) return status;

    return parse_units_length(value, out);
}

static StatusCode parse_opt_units_length(
    const json& object,
    const string& key,
    bool& found,
    ULength& out,
    const string& path
) {
    string value;
    StatusCode status = parse_opt_string(object, key, found, value, path);
    if (status != StatusCode::ok || !found) return status;

    return parse_units_length(value, out);
}

static StatusCode parse_units_time(const string& str, UTime& out) {
    UTime temp;
    if (str == utime_str(UTime::year))
        temp = UTime::year;
    else if (str == utime_str(UTime::month))
        temp = UTime::month;
    else if (str == utime_str(UTime::day))
        temp = UTime::day;
    else if (str == utime_str(UTime::hour))
        temp = UTime::hour;
    else if (str == utime_str(UTime::minute))
        temp = UTime::minute;
    else if (str == utime_str(UTime::second))
        temp = UTime::second;
    else if (str == utime_str(UTime::millisecond))
        temp = UTime::millisecond;
    else if (str == utime_str(UTime::microsecond))
        temp = UTime::microsecond;
    else if (str == utime_str(UTime::nanosecond))
        temp = UTime::nanosecond;
    else
        return StatusCode::unsupported_type;

    out = temp;
    return StatusCode::ok;
}

static StatusCode parse_req_units_time(
    const json& object,
    const string& key,
    UTime& out,
    const string& path
) {
    string value;
    StatusCode status = parse_req_string(object, key, value, path);
    if (status != StatusCode::ok) return status;

    return parse_units_time(value, out);
}

static StatusCode parse_opt_units_time(
    const json& object,
    const string& key,
    bool& found,
    UTime& out,
    const string& path
) {
    string value;
    StatusCode status = parse_opt_string(object, key, found, value, path);
    if (status != StatusCode::ok || !found) return status;

    return parse_units_time(value, out);
}

static StatusCode parse_time_scale(const string& str, TimeScale& out) {
    TimeScale temp;
    if (str == time_scale_str(TimeScale::utc))
        temp = TimeScale::utc;
    else if (str == time_scale_str(TimeScale::ut1) || str == "ut")
        temp = TimeScale::ut1;
    else if (str == time_scale_str(TimeScale::tai))
        temp = TimeScale::tai;
    else if (str == time_scale_str(TimeScale::tt))
        temp = TimeScale::tt;
    else if (str == time_scale_str(TimeScale::tdb))
        temp = TimeScale::tdb;
    else if (str == time_scale_str(TimeScale::gps))
        temp = TimeScale::gps;
    else
        return StatusCode::unsupported_type;

    out = temp;
    return StatusCode::ok;
}

static StatusCode parse_req_time_scale(
    const json& object,
    const string& key,
    TimeScale& out,
    const string& path
) {
    string value;
    StatusCode status = parse_req_string(object, key, value, path);
    if (status != StatusCode::ok) return status;

    return parse_time_scale(value, out);
}

static StatusCode parse_opt_time_scale(
    const json& object,
    const string& key,
    bool& found,
    TimeScale& out,
    const string& path
) {
    string value;
    StatusCode status = parse_opt_string(object, key, found, value, path);
    if (status != StatusCode::ok || !found) return status;

    return parse_time_scale(value, out);
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

static StatusCode parse_opt_i64(
    const json& object,
    const string& key,
    bool& found,
    i64& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_opt_child(object, key, child, found, path);
    if (status != StatusCode::ok) return status;
    if (!found) return StatusCode::ok;

    if (!child->is_number_integer()) return StatusCode::invalid_input;

    out = child->get<i64>();
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

static StatusCode parse_opt_nested_matXd(
    const json& object,
    const string& key,
    bool& found,
    matXd& out,
    const string& path
) {
    found = false;
    const json* child = nullptr;

    StatusCode status = get_opt_child(object, key, child, found, path);
    if (status != StatusCode::ok) return status;
    if (!found) return StatusCode::ok;

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
