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
    default: return "Unknown";
    }
}