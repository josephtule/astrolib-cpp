// Copyright 2025-2026 Joseph Tu Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "util/typedefs.hpp"
#include <string>

enum struct ObservationType : i32 {
    radec,
    azel,
    range,
    range_rate,
    // in estimation frame (not observer relative)
    pos,
    pos_vel,
    // relative
    rel_pos,
    rel_pos_vel,
};
// TODO: add vel and rel_vel measurement types

inline std::string observation_type_str(ObservationType type) {
    switch (type) {
    case ObservationType::radec: return "Right-Ascension + Declination";
    case ObservationType::azel: return "Azimuth + Elevation";
    case ObservationType::range: return "Range";
    case ObservationType::range_rate: return "Range-Rate";
    case ObservationType::pos: return "Position";
    case ObservationType::pos_vel: return "Position + Velocity";
    case ObservationType::rel_pos: return "Relative Position";
    case ObservationType::rel_pos_vel: return "Relative Position + Velocity";
    }
    return "Unknown";
}

inline string observation_type_str_simple(ObservationType type) {
    switch (type) {
    case ObservationType::radec: return "radec";
    case ObservationType::azel: return "azel";
    case ObservationType::range: return "range";
    case ObservationType::range_rate: return "range_rate";
    case ObservationType::pos: return "pos";
    case ObservationType::pos_vel: return "pos_vel";
    case ObservationType::rel_pos: return "rel_pos";
    case ObservationType::rel_pos_vel: return "rel_pos_vel";
    }
    return "unknown";
}