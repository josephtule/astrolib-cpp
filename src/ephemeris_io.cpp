// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/ephemeris_io.hpp"
#include "core/status.hpp"
#include "core/time.hpp"
#include "io_common.hpp"

#include <nlohmann/json.hpp>

struct EphemerisSchemaConfig {
    string name = "astrolib.ephemeris";
    i32 version = 1;
};

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
    else if (str == quaternion_component_order_str(QuaternionComponentOrder::scalar_vector))
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

static StatusCode parse_rotation_convention(
    const string& str,
    RotationConvention& out
) {
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
    if (str == angular_velocity_expression_frame_str(AngularVelocityExpressionFrame::source))
        temp = AngularVelocityExpressionFrame::source;
    else if (
        str == angular_velocity_expression_frame_str(AngularVelocityExpressionFrame::target)
    )
        temp = AngularVelocityExpressionFrame::target;
    else
        return StatusCode::unsupported_type;

    out = temp;
    return StatusCode::ok;
}

static string angular_velocity_direction_str(AngularVelocityDirection direction) {
    switch (direction) {
    case AngularVelocityDirection::target_relative_source: return "target_relative_source";
    case AngularVelocityDirection::source_relative_target: return "source_relative_target";
    }

    return "unknown";
}

static StatusCode parse_angular_velocity_direction(
    const string& str,
    AngularVelocityDirection& out
) {
    AngularVelocityDirection temp;
    if (str == angular_velocity_direction_str(AngularVelocityDirection::target_relative_source))
        temp = AngularVelocityDirection::target_relative_source;
    else if (
        str == angular_velocity_direction_str(AngularVelocityDirection::source_relative_target)
    )
        temp = AngularVelocityDirection::source_relative_target;
    else
        return StatusCode::unsupported_type;

    out = temp;
    return StatusCode::ok;
}
