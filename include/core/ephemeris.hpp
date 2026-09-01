// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/state.hpp"
#include "core/status.hpp"
#include "core/time.hpp"
#include "util/constants.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"

// Time ------------------------------------------------------------------------

struct EphemerisEpochMetadata {
    JulianDate ref_epoch{};
    TimeScale time_scale = TimeScale::tdb;
    UTime offset_unit = UTime::second;
};

// Cartesian/Position ----------------------------------------------------------

struct CartesianEphemerisFrameMetadata {
    string object;
    string center;
    string frame;
    /* Example
    object = "MOON";
    center = "EARTH";
    frame = "J2000";
    */
};

struct CartesianEphemerisUnitMetadata {
    ULength length = ULength::kilometer;
    UTime time = UTime::second;
};

struct EphemerisSource {
    string source_type;
    // TODO: "native_csv","ccsds_oem","cspice_spk","generated","world_history"
    string source_name;
    string source_path;
    string description;
};

struct CartesianEphemerisMetadata {
    EphemerisEpochMetadata epoch;
    CartesianEphemerisFrameMetadata frame;
    CartesianEphemerisUnitMetadata units;
    EphemerisSource source;
};

struct CartesianEphemerisTable {
    CartesianEphemerisMetadata metadata{};
    svec<f64> dt;
    svec<StateTr> states;
    bool has_velocity = true;
};

StatusCode validate_cartesian_ephemeris_table(const CartesianEphemerisTable& table);

// Orientation/Attitude --------------------------------------------------------
enum struct QuaternionComponentOrder { vector_scalar, scalar_vector };

enum struct RotationConvention {
    passive, // transforms coordinates from the source frame to the target frame
    active   // rotates a vector while keeping its coordinate frame fixed
};

enum struct AngularVelocityExpressionFrame {
    source, // components are expressed in the source frame, primary: simulation frame
    target  // components are expressed in the target frame, primary: body frame
};

enum struct AngularVelocityDirection {
    target_relative_source, // angular velocity of the target frame relative to the source
                            // frame
    source_relative_target  // angular velocity of the source frame relative to the target
                            // frame
};

struct OrientationEphemerisFrameMetadata {
    string object;
    string source_frame;
    string target_frame;
    /* Example
    object       = "MARS"
    source_frame = "J2000"
    target_frame = "IAU_MARS"

    or

    object       = "SPACECRAFT"
    source_frame = "J2000"
    target_frame = "SPACECRAFT_BODY"
    */
};

struct OrientationEphemerisUnitMetadata {
    UAngle angular_velocity_angle = UAngle::radian;
    UTime angular_velocity_time = UTime::second;
};

struct OrientationEphemerisConventionMetadata {
    QuaternionComponentOrder quaternion_order = QuaternionComponentOrder::vector_scalar;
    RotationConvention rotation = RotationConvention::passive;
    AngularVelocityExpressionFrame angular_velocity_frame
        = AngularVelocityExpressionFrame::target;
    AngularVelocityDirection angular_velocity_direction
        = AngularVelocityDirection::target_relative_source;
};

struct OrientationEphemerisMetadata {
    EphemerisEpochMetadata epoch{};
    OrientationEphemerisFrameMetadata frame{};
    OrientationEphemerisUnitMetadata units{};
    OrientationEphemerisConventionMetadata convention{};
    EphemerisSource source{};
};

struct OrientationEphemerisTable {
    OrientationEphemerisMetadata metadata{};
    svec<f64> dt;
    svec<StateAtt> states;
    bool has_angular_velocity = true;
};

StatusCode validate_orientation_ephemeris_table(
    const OrientationEphemerisTable& table,
    f64 quaternion_tol = tol9
);


StatusCode canonicalize_orientation_ephemeris_samples(
    OrientationEphemerisTable& table,
    f64 tol = tol12
);