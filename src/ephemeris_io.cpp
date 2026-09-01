// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/ephemeris_io.hpp"
#include "core/ephemeris.hpp"
#include "core/state.hpp"
#include "core/status.hpp"
#include "core/time.hpp"
#include "core/transform.hpp"
#include "io_common.hpp"
#include "util/math.hpp"
#include "util/tools.hpp"
#include "util/units.hpp"

#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>

#include <filesystem>
#include <utility>

struct EphemerisSchemaConfig {
    string name = "astrolib.ephemeris";
    i32 version = 1;
};

struct CartesianEphemerisManifestConfig {
    EphemerisSchemaConfig schema{};
    EphemerisCSVLayout csv{};
    CartesianEphemerisMetadata metadata{};
    bool has_velocity = true;
};

struct OrientationEphemerisManifestConfig {
    EphemerisSchemaConfig schema{};
    EphemerisCSVLayout csv{};
    OrientationEphemerisMetadata metadata{};
    bool has_angular_velocity = true;
};

// Enum strings and parsers

static string quaternion_component_order_str(QuaternionComponentOrder order) {
    switch (order) {
    case QuaternionComponentOrder::vector_scalar: return "vector_scalar";
    case QuaternionComponentOrder::scalar_vector: return "scalar_vector";
    }

    return "unknown";
}

static StatusCode parse_quaternion_component_order(
    const string& str,
    QuaternionComponentOrder& out
) {
    QuaternionComponentOrder temp;
    if (str == quaternion_component_order_str(QuaternionComponentOrder::vector_scalar))
        temp = QuaternionComponentOrder::vector_scalar;
    else if (
        str == quaternion_component_order_str(QuaternionComponentOrder::scalar_vector)
    )
        temp = QuaternionComponentOrder::scalar_vector;
    else
        return StatusCode::unsupported_type;

    out = temp;
    return StatusCode::ok;
}

static string rotation_convention_str(RotationConvention convention) {
    switch (convention) {
    case RotationConvention::passive: return "passive";
    case RotationConvention::active: return "active";
    }

    return "unknown";
}

static StatusCode parse_rotation_convention(const string& str, RotationConvention& out) {
    RotationConvention temp;
    if (str == rotation_convention_str(RotationConvention::passive))
        temp = RotationConvention::passive;
    else if (str == rotation_convention_str(RotationConvention::active))
        temp = RotationConvention::active;
    else
        return StatusCode::unsupported_type;

    out = temp;
    return StatusCode::ok;
}

static string angular_velocity_expression_frame_str(
    AngularVelocityExpressionFrame frame
) {
    switch (frame) {
    case AngularVelocityExpressionFrame::source: return "source";
    case AngularVelocityExpressionFrame::target: return "target";
    }

    return "unknown";
}

static StatusCode parse_angular_velocity_expression_frame(
    const string& str,
    AngularVelocityExpressionFrame& out
) {
    AngularVelocityExpressionFrame temp;
    if (str
        == angular_velocity_expression_frame_str(AngularVelocityExpressionFrame::source))
        temp = AngularVelocityExpressionFrame::source;
    else if (
        str
        == angular_velocity_expression_frame_str(AngularVelocityExpressionFrame::target)
    )
        temp = AngularVelocityExpressionFrame::target;
    else
        return StatusCode::unsupported_type;

    out = temp;
    return StatusCode::ok;
}

static string angular_velocity_direction_str(AngularVelocityDirection direction) {
    switch (direction) {
    case AngularVelocityDirection::target_relative_source:
        return "target_relative_source";
    case AngularVelocityDirection::source_relative_target:
        return "source_relative_target";
    }

    return "unknown";
}

static StatusCode parse_angular_velocity_direction(
    const string& str,
    AngularVelocityDirection& out
) {
    AngularVelocityDirection temp;
    if (str
        == angular_velocity_direction_str(
            AngularVelocityDirection::target_relative_source
        ))
        temp = AngularVelocityDirection::target_relative_source;
    else if (
        str
        == angular_velocity_direction_str(
            AngularVelocityDirection::source_relative_target
        )
    )
        temp = AngularVelocityDirection::source_relative_target;
    else
        return StatusCode::unsupported_type;

    out = temp;
    return StatusCode::ok;
}

// Ephemeris parsers

static StatusCode ephemeris_metadata_status(StatusCode status) {
    if (status == StatusCode::ok || status == StatusCode::unsupported_type) return status;
    return StatusCode::invalid_ephemeris_metadata;
}

static StatusCode parse_ephemeris_schema(
    const json& object,
    EphemerisSchemaConfig& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, "schema", child, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    EphemerisSchemaConfig temp;
    status = parse_req_string(*child, "name", temp.name, path + ".schema");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    status = parse_req_i32(*child, "version", temp.version, path + ".schema");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    if (temp.name != "astrolib.ephemeris") return StatusCode::invalid_ephemeris_metadata;
    if (temp.version != 1) return StatusCode::unsupported_type;

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode parse_ephemeris_epoch_metadata(
    const json& object,
    EphemerisEpochMetadata& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, "epoch", child, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    EphemerisEpochMetadata temp;
    status = parse_req_f64(*child, "day", temp.ref_epoch.day, path + ".epoch");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    status = parse_req_f64(*child, "frac", temp.ref_epoch.frac, path + ".epoch");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    if (temp.ref_epoch.frac < 0.0 || temp.ref_epoch.frac >= 1.0)
        return StatusCode::invalid_ephemeris_metadata;

    string value;
    status = parse_req_string(*child, "time_scale", value, path + ".epoch");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_time_scale(value, temp.time_scale);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    status = parse_req_string(*child, "offset_unit", value, path + ".epoch");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_units_time(value, temp.offset_unit);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode parse_ephemeris_source(
    const json& object,
    EphemerisSource& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, "source", child, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    EphemerisSource temp;
    status = parse_req_string(*child, "source_type", temp.source_type, path + ".source");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    status = parse_req_string(*child, "source_name", temp.source_name, path + ".source");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    if (temp.source_type.empty() || temp.source_name.empty())
        return StatusCode::invalid_ephemeris_metadata;

    bool found = false;
    status = parse_opt_string(
        *child,
        "source_path",
        found,
        temp.source_path,
        path + ".source"
    );
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    status = parse_opt_string(
        *child,
        "description",
        found,
        temp.description,
        path + ".source"
    );
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode parse_ephemeris_csv_layout(
    const json& object,
    EphemerisCSVLayout& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, "data", child, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    EphemerisCSVLayout temp;
    status = parse_req_string(*child, "filepath", temp.filepath, path + ".data");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    status = parse_req_string(*child, "delimiter", temp.delimiter, path + ".data");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    status = parse_req_bool(*child, "has_header", temp.has_header, path + ".data");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    bool found_header_lines = false;
    status = parse_opt_i32(
        *child,
        "header_lines",
        found_header_lines,
        temp.header_lines,
        path + ".data"
    );
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    if (temp.has_header) {
        if (temp.header_lines < 1) return StatusCode::invalid_ephemeris_metadata;
    } else {
        if (found_header_lines && temp.header_lines != 0)
            return StatusCode::invalid_ephemeris_metadata;
        temp.header_lines = 0;
    }

    if (temp.filepath.empty() || temp.delimiter.size() != 1)
        return StatusCode::invalid_ephemeris_metadata;

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode parse_cartesian_ephemeris_frame_metadata(
    const json& object,
    CartesianEphemerisFrameMetadata& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, "frame", child, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    CartesianEphemerisFrameMetadata temp;
    status = parse_req_string(*child, "object", temp.object, path + ".frame");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_req_string(*child, "center", temp.center, path + ".frame");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_req_string(*child, "frame", temp.frame, path + ".frame");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    if (temp.object.empty() || temp.center.empty() || temp.frame.empty())
        return StatusCode::invalid_ephemeris_metadata;

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode parse_cartesian_ephemeris_unit_metadata(
    const json& object,
    CartesianEphemerisUnitMetadata& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, "units", child, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    CartesianEphemerisUnitMetadata temp;
    string value;
    status = parse_req_string(*child, "length", value, path + ".units");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_units_length(value, temp.length);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    status = parse_req_string(*child, "time", value, path + ".units");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_units_time(value, temp.time);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode parse_orientation_ephemeris_frame_metadata(
    const json& object,
    OrientationEphemerisFrameMetadata& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, "frame", child, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    OrientationEphemerisFrameMetadata temp;
    status = parse_req_string(*child, "object", temp.object, path + ".frame");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_req_string(*child, "source_frame", temp.source_frame, path + ".frame");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_req_string(*child, "target_frame", temp.target_frame, path + ".frame");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    if (temp.object.empty() || temp.source_frame.empty() || temp.target_frame.empty())
        return StatusCode::invalid_ephemeris_metadata;

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode parse_orientation_ephemeris_unit_metadata(
    const json& object,
    OrientationEphemerisUnitMetadata& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, "units", child, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    OrientationEphemerisUnitMetadata temp;
    string value;
    status = parse_req_string(*child, "angular_velocity_angle", value, path + ".units");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_units_angle(value, temp.angular_velocity_angle);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    status = parse_req_string(*child, "angular_velocity_time", value, path + ".units");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_units_time(value, temp.angular_velocity_time);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode parse_orientation_ephemeris_convention_metadata(
    const json& object,
    OrientationEphemerisConventionMetadata& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, "convention", child, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    OrientationEphemerisConventionMetadata temp;
    string value;
    status = parse_req_string(*child, "quaternion_order", value, path + ".convention");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_quaternion_component_order(value, temp.quaternion_order);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    status = parse_req_string(*child, "rotation", value, path + ".convention");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_rotation_convention(value, temp.rotation);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    status
        = parse_req_string(*child, "angular_velocity_frame", value, path + ".convention");
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_angular_velocity_expression_frame(value, temp.angular_velocity_frame);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    status = parse_req_string(
        *child,
        "angular_velocity_direction",
        value,
        path + ".convention"
    );
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    status = parse_angular_velocity_direction(value, temp.angular_velocity_direction);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode parse_cartesian_ephemeris_metadata(
    const json& object,
    CartesianEphemerisMetadata& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, "metadata", child, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    CartesianEphemerisMetadata temp;
    status = parse_ephemeris_epoch_metadata(*child, temp.epoch, path + ".metadata");
    if (status != StatusCode::ok) return status;
    status = parse_cartesian_ephemeris_frame_metadata(
        *child,
        temp.frame,
        path + ".metadata"
    );
    if (status != StatusCode::ok) return status;
    status
        = parse_cartesian_ephemeris_unit_metadata(*child, temp.units, path + ".metadata");
    if (status != StatusCode::ok) return status;
    status = parse_ephemeris_source(*child, temp.source, path + ".metadata");
    if (status != StatusCode::ok) return status;

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode parse_orientation_ephemeris_metadata(
    const json& object,
    OrientationEphemerisMetadata& out,
    const string& path
) {
    const json* child = nullptr;
    StatusCode status = get_req_child(object, "metadata", child, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    OrientationEphemerisMetadata temp;
    status = parse_ephemeris_epoch_metadata(*child, temp.epoch, path + ".metadata");
    if (status != StatusCode::ok) return status;
    status = parse_orientation_ephemeris_frame_metadata(
        *child,
        temp.frame,
        path + ".metadata"
    );
    if (status != StatusCode::ok) return status;
    status = parse_orientation_ephemeris_unit_metadata(
        *child,
        temp.units,
        path + ".metadata"
    );
    if (status != StatusCode::ok) return status;
    status = parse_orientation_ephemeris_convention_metadata(
        *child,
        temp.convention,
        path + ".metadata"
    );
    if (status != StatusCode::ok) return status;
    status = parse_ephemeris_source(*child, temp.source, path + ".metadata");
    if (status != StatusCode::ok) return status;

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode parse_cartesian_ephemeris_manifest(
    const json& root,
    CartesianEphemerisManifestConfig& out,
    const string& path
) {
    if (!root.is_object()) return StatusCode::invalid_ephemeris_metadata;

    CartesianEphemerisManifestConfig temp;
    StatusCode status = parse_ephemeris_schema(root, temp.schema, path);
    if (status != StatusCode::ok) return status;

    string domain;
    status = parse_req_string(root, "domain", domain, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    if (domain != "cartesian") return StatusCode::invalid_ephemeris_metadata;

    status = parse_ephemeris_csv_layout(root, temp.csv, path);
    if (status != StatusCode::ok) return status;
    status = parse_cartesian_ephemeris_metadata(root, temp.metadata, path);
    if (status != StatusCode::ok) return status;
    status = parse_req_bool(root, "has_velocity", temp.has_velocity, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode parse_orientation_ephemeris_manifest(
    const json& root,
    OrientationEphemerisManifestConfig& out,
    const string& path
) {
    if (!root.is_object()) return StatusCode::invalid_ephemeris_metadata;

    OrientationEphemerisManifestConfig temp;
    StatusCode status = parse_ephemeris_schema(root, temp.schema, path);
    if (status != StatusCode::ok) return status;

    string domain;
    status = parse_req_string(root, "domain", domain, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);
    if (domain != "orientation") return StatusCode::invalid_ephemeris_metadata;

    status = parse_ephemeris_csv_layout(root, temp.csv, path);
    if (status != StatusCode::ok) return status;
    status = parse_orientation_ephemeris_metadata(root, temp.metadata, path);
    if (status != StatusCode::ok) return status;
    status
        = parse_req_bool(root, "has_angular_velocity", temp.has_angular_velocity, path);
    if (status != StatusCode::ok) return ephemeris_metadata_status(status);

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode resolve_ephemeris_data_path(
    const string& manifest_filepath,
    const string& data_filepath,
    std::filesystem::path& out
) {
    // resolves manifest and data filepaths to absolute filepaths
    // example:
    // manifests stored in root/data/ephemerides/moon_manifest.json
    // data stored in root/data/tables/moon_data.csv

    if (manifest_filepath.empty() || data_filepath.empty())
        return StatusCode::invalid_input;

    const std::filesystem::path manifest_path{manifest_filepath};
    const std::filesystem::path data_path{data_filepath};
    std::filesystem::path temp = data_path;

    // removes redundant path components
    if (data_path.is_relative()) temp = manifest_path.parent_path() / data_path;
    temp = temp.lexically_normal();

    std::error_code ec;
    if (!std::filesystem::is_regular_file(temp, ec) || ec)
        return StatusCode::file_not_found;

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode parse_ephemeris_csv_row(
    const string& line,
    const string& delimiter,
    i32 expected_columns,
    svec<f64>& out
) {
    if (delimiter.size() != 1 || expected_columns <= 0) {
        return StatusCode::invalid_input;
    }

    svec<f64> temp;
    temp.reserve(expected_columns);

    string::size_type start = 0;
    while (true) {
        const string::size_type end = line.find(delimiter, start);

        string token = (end == string::npos) ? line.substr(start)
                                             : line.substr(start, end - start);

        token = trim(token);
        // empty field
        if (token.empty()) {
            return StatusCode::invalid_input;
        }

        try {
            // allow fortran scientific notation format
            for (char& c : token) {
                if (c == 'd' || c == 'D') c = 'e';
            }

            std::size_t consumed = 0;
            const f64 value = std::stod(token, &consumed);

            if (consumed != token.size() || !std::isfinite(value))
                return StatusCode::invalid_input;

            temp.push_back(value);
        } catch (const std::invalid_argument&) {
            return StatusCode::invalid_input;
        } catch (const std::out_of_range&) {
            return StatusCode::invalid_input;
        }

        if (static_cast<i32>(temp.size()) > expected_columns) {
            return StatusCode::invalid_input;
        }

        if (end == string::npos) {
            break;
        }

        start = end + delimiter.size();
    }

    if (static_cast<i32>(temp.size()) != expected_columns) {
        return StatusCode::invalid_input;
    }

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode read_cartesian_ephemeris_csv(
    const std::filesystem::path& filepath,
    const CartesianEphemerisManifestConfig& manifest,
    CartesianEphemerisTable& out
) {
    /* columns
    0: dt
    1: rx
    2: ry
    3: rz
    4: vx, when enabled
    5: vy, when enabled
    6: vz, when enabled
    */
    // NOTE: these are not yet canonicalized to the preferred convention

    CartesianEphemerisTable temp{};
    temp.metadata = manifest.metadata;
    temp.has_velocity = manifest.has_velocity;

    const i32 expected_columns = manifest.has_velocity ? 7 : 4;

    std::ifstream file(filepath);

    if (!file.is_open()) {
        return StatusCode::file_open_failed;
    }

    std::string line;
    StatusCode status;
    svec<f64> values(expected_columns);

    if (manifest.csv.has_header) {
        for (i32 i = 0; i < manifest.csv.header_lines; ++i) {
            if (!std::getline(file, line)) return StatusCode::parse_failed;
        }
    }

    while (std::getline(file, line)) {
        // process one row
        status = parse_ephemeris_csv_row(
            line,
            manifest.csv.delimiter,
            expected_columns,
            values
        );
        if (status != StatusCode::ok) return StatusCode::parse_failed;

        temp.dt.push_back(values[0]);

        StateTr x{.r = vec3d{values[1], values[2], values[3]}};
        if (manifest.has_velocity) {
            x.v = vec3d{values[4], values[5], values[6]};
        }
        temp.states.push_back(x);
    }

    if (file.bad()) return StatusCode::parse_failed;
    if (temp.states.empty()) return StatusCode::empty_ephemeris;

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode read_orientation_ephemeris_csv(
    const std::filesystem::path& filepath,
    const OrientationEphemerisManifestConfig& manifest,
    OrientationEphemerisTable& out
) {
    /* columns
    0: dt
    1-4: raw quaternion components
    5-7: angular velocity, when enabled
    */
    // NOTE: these are not yet canonicalized to the preferred convention

    OrientationEphemerisTable temp{};
    temp.metadata = manifest.metadata;
    temp.has_angular_velocity = manifest.has_angular_velocity;

    const i32 expected_columns = manifest.has_angular_velocity ? 8 : 5;

    std::ifstream file(filepath);

    if (!file.is_open()) {
        return StatusCode::file_open_failed;
    }

    std::string line;
    StatusCode status;
    svec<f64> values(expected_columns);

    if (manifest.csv.has_header) {
        for (i32 i = 0; i < manifest.csv.header_lines; ++i) {
            if (!std::getline(file, line)) return StatusCode::parse_failed;
        }
    }

    while (std::getline(file, line)) {
        // process one row
        status = parse_ephemeris_csv_row(
            line,
            manifest.csv.delimiter,
            expected_columns,
            values
        );
        if (status != StatusCode::ok) return StatusCode::parse_failed;

        temp.dt.push_back(values[0]);

        StateAtt x{.q = vec4d{values[1], values[2], values[3], values[4]}};
        if (manifest.has_angular_velocity) {
            x.w = vec3d{values[5], values[6], values[7]};
        }
        temp.states.push_back(x);
    }

    if (file.bad()) return StatusCode::parse_failed;
    if (temp.states.empty()) return StatusCode::empty_ephemeris;

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode convert_cartesian_ephemeris_canon(
    const CartesianEphemerisTable& in,
    CartesianEphemerisTable& out
) {
    CartesianEphemerisTable temp = in;

    f64 scale_pos = convert_length(1.0, temp.metadata.units.length, ULength::kilometer);
    f64 scale_time_epoch;
    if (!convert_time(
            1.0,
            temp.metadata.epoch.offset_unit,
            UTime::second,
            scale_time_epoch
        ))
        return StatusCode::unsupported_type;
    f64 scale_time_vel;
    if (!convert_time(1.0, temp.metadata.units.time, UTime::second, scale_time_vel))
        return StatusCode::unsupported_type;
    f64 scale_vel = scale_pos / scale_time_vel;

    for (i32 i = 0; i < temp.dt.size(); ++i) {
        temp.dt[i] *= scale_time_epoch;
        temp.states[i].r *= scale_pos;
        if (temp.has_velocity) temp.states[i].v *= scale_vel;
    }

    temp.metadata.epoch.offset_unit = UTime::second;
    temp.metadata.units.length = ULength::kilometer;
    temp.metadata.units.time = UTime::second;

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode convert_orientation_ephemeris_canon(
    const OrientationEphemerisTable& in,
    OrientationEphemerisTable& out,
    f64 q_tol = tol12
) {
    if (!finite_pos(q_tol)) return StatusCode::invalid_input;

    OrientationEphemerisTable temp = in;

    f64 scale_angle
        = convert_angle(1.0, temp.metadata.units.angular_velocity_angle, UAngle::radian);
    f64 scale_time_epoch;
    if (!convert_time(
            1.0,
            temp.metadata.epoch.offset_unit,
            UTime::second,
            scale_time_epoch
        ))
        return StatusCode::unsupported_type;
    f64 scale_time_ang_vel;
    if (!convert_time(
            1.0,
            temp.metadata.units.angular_velocity_time,
            UTime::second,
            scale_time_ang_vel
        ))
        return StatusCode::unsupported_type;
    f64 scale_ang_vel = scale_angle / scale_time_ang_vel;

    auto& convention = temp.metadata.convention;
    for (i32 i = 0; i < temp.dt.size(); ++i) {
        temp.dt[i] *= scale_time_epoch;

        vec4d& q = temp.states[i].q;
        if (convention.quaternion_order == QuaternionComponentOrder::scalar_vector) {
            vec4d q_temp = q;
            q(3) = q(0);
            q.segment(0, 3) = q_temp.segment(1, 3);
        }
        if (convention.rotation == RotationConvention::active) {
            q = ep_conj(q);
        }
        if (!valid_quaternion(q, q_tol)) return StatusCode::invalid_att_state;
        normalize_quaternion_inplace(q, q_tol);

        if (temp.has_angular_velocity) {
            vec3d& w = temp.states[i].w;
            w *= scale_ang_vel;
            if (convention.angular_velocity_direction
                == AngularVelocityDirection::source_relative_target) {
                w = -w;
            }
            if (convention.angular_velocity_frame
                == AngularVelocityExpressionFrame::source) {
                w = ep_rotate_fast_passive(q, w);
            }
        }
    }

    temp.metadata.epoch.offset_unit = UTime::second;
    temp.metadata.units.angular_velocity_angle = UAngle::radian;
    temp.metadata.units.angular_velocity_time = UTime::second;
    temp.metadata.convention.quaternion_order = QuaternionComponentOrder::vector_scalar;
    temp.metadata.convention.rotation = RotationConvention::passive;
    temp.metadata.convention.angular_velocity_frame
        = AngularVelocityExpressionFrame::target;
    temp.metadata.convention.angular_velocity_direction
        = AngularVelocityDirection::target_relative_source;

    StatusCode status = canonicalize_orientation_ephemeris_samples(temp, q_tol);
    if (status != StatusCode::ok) return status;

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode read_ephemeris_manifest_json(
    const string& manifest_filepath,
    json& out
) {
    if (manifest_filepath.empty()) return StatusCode::invalid_input;

    std::ifstream file(manifest_filepath);
    if (!file.is_open()) return StatusCode::file_open_failed;

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

    out = root;

    return StatusCode::ok;
}

StatusCode load_cartesian_ephemeris(
    const string& manifest_filepath,
    CartesianEphemerisTable& out
) {
    json root;
    CartesianEphemerisManifestConfig manifest;
    std::filesystem::path data_filepath;
    CartesianEphemerisTable raw;
    CartesianEphemerisTable canonical;

    StatusCode status;

    status = read_ephemeris_manifest_json(manifest_filepath, root);
    if (status != StatusCode::ok) return status;

    status = parse_cartesian_ephemeris_manifest(root, manifest, "root");
    if (status != StatusCode::ok) return status;

    status = resolve_ephemeris_data_path(
        manifest_filepath,
        manifest.csv.filepath,
        data_filepath
    );
    if (status != StatusCode::ok) return status;

    status = read_cartesian_ephemeris_csv(data_filepath, manifest, raw);
    if (status != StatusCode::ok) return status;

    status = convert_cartesian_ephemeris_canon(raw, canonical);
    if (status != StatusCode::ok) return status;

    status = validate_cartesian_ephemeris_table(canonical);
    if (status != StatusCode::ok) return status;

    out = std::move(canonical);
    return StatusCode::ok;
}

StatusCode load_orientation_ephemeris(
    const string& manifest_filepath,
    OrientationEphemerisTable& out,
    f64 tol
) {
    json root;
    OrientationEphemerisManifestConfig manifest;
    std::filesystem::path data_filepath;
    OrientationEphemerisTable raw;
    OrientationEphemerisTable canonical;

    StatusCode status;

    status = read_ephemeris_manifest_json(manifest_filepath, root);
    if (status != StatusCode::ok) return status;

    status = parse_orientation_ephemeris_manifest(root, manifest, "root");
    if (status != StatusCode::ok) return status;

    status = resolve_ephemeris_data_path(
        manifest_filepath,
        manifest.csv.filepath,
        data_filepath
    );
    if (status != StatusCode::ok) return status;

    status = read_orientation_ephemeris_csv(data_filepath, manifest, raw);
    if (status != StatusCode::ok) return status;

    status = convert_orientation_ephemeris_canon(raw, canonical, tol);
    if (status != StatusCode::ok) return status;

    status = validate_orientation_ephemeris_table(canonical, tol);
    if (status != StatusCode::ok) return status;

    out = std::move(canonical);
    return StatusCode::ok;
}

static json json_from_ephemeris_schema() {
    return json{{"name", "astrolib.ephemeris"}, {"version", 1}};
}

static json json_from_ephemeris_epoch_metadata(const EphemerisEpochMetadata& epoch) {
    return json{
        {"day", epoch.ref_epoch.day},
        {"frac", epoch.ref_epoch.frac},
        {"time_scale", time_scale_str(epoch.time_scale)},
        {"offset_unit", utime_str(epoch.offset_unit)}
    };
}

static json json_from_ephemeris_source(const EphemerisSource& source) {
    return json{
        {"source_type", source.source_type},
        {"source_name", source.source_name},
        {"source_path", source.source_path},
        {"description", source.description}
    };
}

static json json_from_ephemeris_csv_layout(const EphemerisCSVLayout& layout) {
    return json{
        {"filepath", layout.filepath},
        {"delimiter", layout.delimiter},
        {"has_header", layout.has_header},
        {"header_lines", layout.header_lines}
    };
}

static json json_from_cartesian_ephemeris_metadata(
    const CartesianEphemerisMetadata& metadata
) {
    json frame{
        {"object", metadata.frame.object},
        {"center", metadata.frame.center},
        {"frame", metadata.frame.frame}
    };
    json units{
        {"length", ulength_str(metadata.units.length)},
        {"time", utime_str(metadata.units.time)}
    };

    return json{
        {"epoch", json_from_ephemeris_epoch_metadata(metadata.epoch)},
        {"frame", std::move(frame)},
        {"units", std::move(units)},
        {"source", json_from_ephemeris_source(metadata.source)}
    };
}

static json json_from_orientation_ephemeris_metadata(
    const OrientationEphemerisMetadata& metadata
) {
    json frame{
        {"object", metadata.frame.object},
        {"source_frame", metadata.frame.source_frame},
        {"target_frame", metadata.frame.target_frame}
    };
    json units{
        {"angular_velocity_angle", uangle_str(metadata.units.angular_velocity_angle)},
        {"angular_velocity_time", utime_str(metadata.units.angular_velocity_time)}
    };
    json convention{
        {"quaternion_order", quaternion_component_order_str(
             metadata.convention.quaternion_order
         )},
        {"rotation", rotation_convention_str(metadata.convention.rotation)},
        {"angular_velocity_frame", angular_velocity_expression_frame_str(
             metadata.convention.angular_velocity_frame
         )},
        {"angular_velocity_direction", angular_velocity_direction_str(
             metadata.convention.angular_velocity_direction
         )}
    };

    return json{
        {"epoch", json_from_ephemeris_epoch_metadata(metadata.epoch)},
        {"frame", std::move(frame)},
        {"units", std::move(units)},
        {"convention", std::move(convention)},
        {"source", json_from_ephemeris_source(metadata.source)}
    };
}

static StatusCode normalize_ephemeris_write_layout(
    const EphemerisWriteOptions& opts,
    EphemerisCSVLayout& out
) {
    if (opts.precision <= 0 || opts.csv.filepath.empty()
        || opts.csv.delimiter.size() != 1)
        return StatusCode::invalid_input;

    EphemerisCSVLayout temp = opts.csv;
    temp.header_lines = temp.has_header ? 1 : 0;

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode resolve_ephemeris_output_path(
    const string& manifest_filepath,
    const string& data_filepath,
    std::filesystem::path& out
) {
    if (manifest_filepath.empty() || data_filepath.empty())
        return StatusCode::invalid_input;

    const std::filesystem::path manifest_path{manifest_filepath};
    const std::filesystem::path data_path{data_filepath};
    std::filesystem::path temp = data_path;

    if (data_path.is_relative())
        temp = manifest_path.parent_path() / data_path;
    temp = temp.lexically_normal();

    out = std::move(temp);
    return StatusCode::ok;
}

struct EphemerisSavePaths {
    std::filesystem::path manifest;
    std::filesystem::path data;
    std::filesystem::path manifest_temp;
    std::filesystem::path data_temp;
    std::filesystem::path manifest_backup;
    std::filesystem::path data_backup;
};

static StatusCode make_ephemeris_sidecar_path(
    const std::filesystem::path& destination,
    const string& label,
    std::filesystem::path& out
) {
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();

    for (i32 i = 0; i < 100; ++i) {
        std::filesystem::path temp = destination;
        temp += ".astrolib-" + label + "-" + std::to_string(stamp) + "-"
                + std::to_string(i);

        std::error_code ec;
        const bool exists = std::filesystem::exists(temp, ec);
        if (ec) return StatusCode::file_publish_failed;
        if (!exists) {
            out = std::move(temp);
            return StatusCode::ok;
        }
    }

    return StatusCode::file_publish_failed;
}

static StatusCode normalize_ephemeris_save_destination(
    const std::filesystem::path& filepath,
    std::filesystem::path& out
) {
    if (filepath.empty()) return StatusCode::invalid_input;

    std::error_code ec;
    std::filesystem::path temp = std::filesystem::absolute(filepath, ec);
    if (ec) return StatusCode::file_publish_failed;
    temp = temp.lexically_normal();

    const std::filesystem::path parent = temp.parent_path();
    if (!std::filesystem::is_directory(parent, ec) || ec)
        return StatusCode::file_not_found;

    temp = std::filesystem::weakly_canonical(temp, ec);
    if (ec) return StatusCode::file_publish_failed;

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode validate_ephemeris_save_destination(
    const std::filesystem::path& filepath,
    bool overwrite
) {
    std::error_code ec;
    const bool exists = std::filesystem::exists(filepath, ec);
    if (ec) return StatusCode::file_publish_failed;
    if (!exists) return StatusCode::ok;

    if (!std::filesystem::is_regular_file(filepath, ec) || ec)
        return StatusCode::invalid_input;
    if (!overwrite) return StatusCode::file_already_exists;

    return StatusCode::ok;
}

static StatusCode prepare_ephemeris_save_paths(
    const std::filesystem::path& manifest_filepath,
    const std::filesystem::path& data_filepath,
    bool overwrite,
    EphemerisSavePaths& out
) {
    EphemerisSavePaths temp;
    StatusCode status = normalize_ephemeris_save_destination(
        manifest_filepath,
        temp.manifest
    );
    if (status != StatusCode::ok) return status;

    status = normalize_ephemeris_save_destination(data_filepath, temp.data);
    if (status != StatusCode::ok) return status;
    if (temp.manifest == temp.data) return StatusCode::invalid_input;

    status = validate_ephemeris_save_destination(temp.manifest, overwrite);
    if (status != StatusCode::ok) return status;
    status = validate_ephemeris_save_destination(temp.data, overwrite);
    if (status != StatusCode::ok) return status;

    status = make_ephemeris_sidecar_path(temp.manifest, "tmp", temp.manifest_temp);
    if (status != StatusCode::ok) return status;
    status = make_ephemeris_sidecar_path(temp.data, "tmp", temp.data_temp);
    if (status != StatusCode::ok) return status;
    status = make_ephemeris_sidecar_path(
        temp.manifest,
        "backup",
        temp.manifest_backup
    );
    if (status != StatusCode::ok) return status;
    status = make_ephemeris_sidecar_path(temp.data, "backup", temp.data_backup);
    if (status != StatusCode::ok) return status;

    out = std::move(temp);
    return StatusCode::ok;
}

static void cleanup_ephemeris_temp_files(const EphemerisSavePaths& paths) {
    std::error_code ec;
    std::filesystem::remove(paths.data_temp, ec);
    ec.clear();
    std::filesystem::remove(paths.manifest_temp, ec);
}

static bool ephemeris_regular_file_exists(
    const std::filesystem::path& filepath,
    bool& exists
) {
    std::error_code ec;
    exists = std::filesystem::exists(filepath, ec);
    if (ec) return false;
    if (!exists) return true;

    return std::filesystem::is_regular_file(filepath, ec) && !ec;
}

static StatusCode publish_ephemeris_files(
    const EphemerisSavePaths& paths,
    bool overwrite
) {
    bool data_exists = false;
    bool manifest_exists = false;
    if (!ephemeris_regular_file_exists(paths.data, data_exists)
        || !ephemeris_regular_file_exists(paths.manifest, manifest_exists)) {
        cleanup_ephemeris_temp_files(paths);
        return StatusCode::file_publish_failed;
    }
    if (!overwrite && (data_exists || manifest_exists)) {
        cleanup_ephemeris_temp_files(paths);
        return StatusCode::file_already_exists;
    }

    bool data_backed_up = false;
    bool manifest_backed_up = false;
    bool data_published = false;
    bool manifest_published = false;
    std::error_code ec;

    auto rollback = [&]() {
        if (manifest_published) {
            std::filesystem::remove(paths.manifest, ec);
            ec.clear();
        }
        if (data_published) {
            std::filesystem::remove(paths.data, ec);
            ec.clear();
        }
        if (manifest_backed_up) {
            std::filesystem::rename(paths.manifest_backup, paths.manifest, ec);
            ec.clear();
        }
        if (data_backed_up) {
            std::filesystem::rename(paths.data_backup, paths.data, ec);
            ec.clear();
        }
        cleanup_ephemeris_temp_files(paths);
    };

    if (data_exists) {
        std::filesystem::rename(paths.data, paths.data_backup, ec);
        if (ec) {
            rollback();
            return StatusCode::file_publish_failed;
        }
        data_backed_up = true;
    }
    if (manifest_exists) {
        std::filesystem::rename(paths.manifest, paths.manifest_backup, ec);
        if (ec) {
            rollback();
            return StatusCode::file_publish_failed;
        }
        manifest_backed_up = true;
    }

    std::filesystem::rename(paths.data_temp, paths.data, ec);
    if (ec) {
        rollback();
        return StatusCode::file_publish_failed;
    }
    data_published = true;

    std::filesystem::rename(paths.manifest_temp, paths.manifest, ec);
    if (ec) {
        rollback();
        return StatusCode::file_publish_failed;
    }
    manifest_published = true;

    if (data_backed_up) {
        std::filesystem::remove(paths.data_backup, ec);
        if (ec) return StatusCode::file_publish_failed;
    }
    if (manifest_backed_up) {
        std::filesystem::remove(paths.manifest_backup, ec);
        if (ec) return StatusCode::file_publish_failed;
    }

    return StatusCode::ok;
}

static EphemerisCSVLayout ephemeris_manifest_csv_layout(
    const string& manifest_filepath,
    const std::filesystem::path& data_filepath,
    const EphemerisCSVLayout& layout
) {
    EphemerisCSVLayout temp = layout;
    const std::filesystem::path configured_path{layout.filepath};

    if (configured_path.is_relative()) {
        temp.filepath = configured_path.lexically_normal().generic_string();
        return temp;
    }

    std::error_code ec;
    std::filesystem::path manifest_absolute
        = std::filesystem::absolute(manifest_filepath, ec);
    if (!ec) {
        std::filesystem::path data_absolute = std::filesystem::absolute(data_filepath, ec);
        if (!ec) {
            const std::filesystem::path relative = data_absolute.lexically_normal()
                                                       .lexically_relative(
                                                           manifest_absolute.parent_path()
                                                               .lexically_normal()
                                                       );
            const auto first = relative.begin();
            if (!relative.empty() && first != relative.end() && *first != "..") {
                temp.filepath = relative.generic_string();
                return temp;
            }
        }
    }

    temp.filepath = data_filepath.lexically_normal().generic_string();
    return temp;
}

static StatusCode write_cartesian_ephemeris_csv(
    const std::filesystem::path& filepath,
    const CartesianEphemerisTable& table,
    const EphemerisCSVLayout& layout,
    i32 precision
) {
    std::ofstream file(filepath);
    if (!file.is_open()) return StatusCode::file_open_failed;

    const char delimiter = layout.delimiter.front();
    file << std::setprecision(precision);

    if (layout.has_header) {
        file << "dt" << delimiter << "rx" << delimiter << "ry" << delimiter << "rz";
        if (table.has_velocity)
            file << delimiter << "vx" << delimiter << "vy" << delimiter << "vz";
        file << '\n';
    }

    for (i32 i = 0; i < table.dt.size(); ++i) {
        const StateTr& state = table.states[i];
        file << table.dt[i] << delimiter << state.r(0) << delimiter << state.r(1)
             << delimiter << state.r(2);
        if (table.has_velocity)
            file << delimiter << state.v(0) << delimiter << state.v(1) << delimiter
                 << state.v(2);
        file << '\n';
    }

    if (!file.good()) return StatusCode::file_write_failed;
    file.close();
    if (file.fail()) return StatusCode::file_close_failed;

    return StatusCode::ok;
}

static StatusCode write_orientation_ephemeris_csv(
    const std::filesystem::path& filepath,
    const OrientationEphemerisTable& table,
    const EphemerisCSVLayout& layout,
    i32 precision
) {
    std::ofstream file(filepath);
    if (!file.is_open()) return StatusCode::file_open_failed;

    const char delimiter = layout.delimiter.front();
    file << std::setprecision(precision);

    if (layout.has_header) {
        file << "dt" << delimiter << "qx" << delimiter << "qy" << delimiter << "qz"
             << delimiter << "qw";
        if (table.has_angular_velocity)
            file << delimiter << "wx" << delimiter << "wy" << delimiter << "wz";
        file << '\n';
    }

    for (i32 i = 0; i < table.dt.size(); ++i) {
        const StateAtt& state = table.states[i];
        file << table.dt[i] << delimiter << state.q(0) << delimiter << state.q(1)
             << delimiter << state.q(2) << delimiter << state.q(3);
        if (table.has_angular_velocity)
            file << delimiter << state.w(0) << delimiter << state.w(1) << delimiter
                 << state.w(2);
        file << '\n';
    }

    if (!file.good()) return StatusCode::file_write_failed;
    file.close();
    if (file.fail()) return StatusCode::file_close_failed;

    return StatusCode::ok;
}

static json json_from_cartesian_ephemeris_manifest(
    const CartesianEphemerisTable& table,
    const EphemerisCSVLayout& layout
) {
    return json{
        {"schema", json_from_ephemeris_schema()},
        {"domain", "cartesian"},
        {"data", json_from_ephemeris_csv_layout(layout)},
        {"metadata", json_from_cartesian_ephemeris_metadata(table.metadata)},
        {"has_velocity", table.has_velocity}
    };
}

static json json_from_orientation_ephemeris_manifest(
    const OrientationEphemerisTable& table,
    const EphemerisCSVLayout& layout
) {
    return json{
        {"schema", json_from_ephemeris_schema()},
        {"domain", "orientation"},
        {"data", json_from_ephemeris_csv_layout(layout)},
        {"metadata", json_from_orientation_ephemeris_metadata(table.metadata)},
        {"has_angular_velocity", table.has_angular_velocity}
    };
}

static StatusCode write_ephemeris_manifest_json(
    const std::filesystem::path& filepath,
    const json& root
) {
    std::ofstream file(filepath);
    if (!file.is_open()) return StatusCode::file_open_failed;

    file << std::setw(4) << root << '\n';
    if (!file.good()) return StatusCode::file_write_failed;
    file.close();
    if (file.fail()) return StatusCode::file_close_failed;

    return StatusCode::ok;
}

StatusCode save_cartesian_ephemeris(
    const string& manifest_filepath,
    const CartesianEphemerisTable& in,
    const EphemerisWriteOptions& opts
) {
    StatusCode status = validate_cartesian_ephemeris_table(in);
    if (status != StatusCode::ok) return status;

    EphemerisCSVLayout layout;
    status = normalize_ephemeris_write_layout(opts, layout);
    if (status != StatusCode::ok) return status;

    std::filesystem::path data_filepath;
    status = resolve_ephemeris_output_path(
        manifest_filepath,
        layout.filepath,
        data_filepath
    );
    if (status != StatusCode::ok) return status;

    EphemerisSavePaths paths;
    status = prepare_ephemeris_save_paths(
        manifest_filepath,
        data_filepath,
        opts.overwrite,
        paths
    );
    if (status != StatusCode::ok) return status;

    const EphemerisCSVLayout manifest_layout = ephemeris_manifest_csv_layout(
        paths.manifest.string(),
        paths.data,
        layout
    );
    const json root = json_from_cartesian_ephemeris_manifest(in, manifest_layout);

    status = write_cartesian_ephemeris_csv(
        paths.data_temp,
        in,
        layout,
        opts.precision
    );
    if (status != StatusCode::ok) {
        cleanup_ephemeris_temp_files(paths);
        return status;
    }

    status = write_ephemeris_manifest_json(paths.manifest_temp, root);
    if (status != StatusCode::ok) {
        cleanup_ephemeris_temp_files(paths);
        return status;
    }

    return publish_ephemeris_files(paths, opts.overwrite);
}

StatusCode save_orientation_ephemeris(
    const string& manifest_filepath,
    const OrientationEphemerisTable& in,
    const EphemerisWriteOptions& opts
) {
    StatusCode status = validate_orientation_ephemeris_table(in);
    if (status != StatusCode::ok) return status;

    EphemerisCSVLayout layout;
    status = normalize_ephemeris_write_layout(opts, layout);
    if (status != StatusCode::ok) return status;

    std::filesystem::path data_filepath;
    status = resolve_ephemeris_output_path(
        manifest_filepath,
        layout.filepath,
        data_filepath
    );
    if (status != StatusCode::ok) return status;

    EphemerisSavePaths paths;
    status = prepare_ephemeris_save_paths(
        manifest_filepath,
        data_filepath,
        opts.overwrite,
        paths
    );
    if (status != StatusCode::ok) return status;

    const EphemerisCSVLayout manifest_layout = ephemeris_manifest_csv_layout(
        paths.manifest.string(),
        paths.data,
        layout
    );
    const json root = json_from_orientation_ephemeris_manifest(in, manifest_layout);

    status = write_orientation_ephemeris_csv(
        paths.data_temp,
        in,
        layout,
        opts.precision
    );
    if (status != StatusCode::ok) {
        cleanup_ephemeris_temp_files(paths);
        return status;
    }

    status = write_ephemeris_manifest_json(paths.manifest_temp, root);
    if (status != StatusCode::ok) {
        cleanup_ephemeris_temp_files(paths);
        return status;
    }

    return publish_ephemeris_files(paths, opts.overwrite);
}
