#pragma once

#include "core/estimation_common.hpp"
#include "util/constants.hpp"
#include "util/typedefs.hpp"

enum struct UAngle : i32 { radian, degree, arcminute, arcsecond, milliarcsecond };

inline StatusCode string_to_uangle(std::string str, UAngle& out) {
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

enum struct ULength : i32 {
    nanometer,
    millimeter,
    centimeter,
    meter,
    kilometer,
    inch,
    foot,
    mile,
    au,
};

enum struct UTime : i32 {
    year,
    month,
    day,
    hour,
    minute,
    second,
    millisecond,
    microsecond,
    nanosecond,
};

enum struct UMass : i32 {
    MICROGRAM,
    MILLIGRAM,
    GRAM,
    KILOGRAM,
    POUNDMASS,
    OUNCEMASS,
};

template <typename T>
T convert_angle(T val, UAngle uin, UAngle uout) {
    if (uin == uout) return val;

    // convert to radians
    switch (uin) {
    case UAngle::radian: break;
    case UAngle::milliarcsecond: val /= static_cast<T>(1000); [[fallthrough]];
    case UAngle::arcsecond: val /= static_cast<T>(60); [[fallthrough]];
    case UAngle::arcminute: val /= static_cast<T>(60); [[fallthrough]];
    case UAngle::degree: val *= deg_to_rad; break;
    }

    // convert to output units
    switch (uout) {
    case UAngle::radian: break;
    case UAngle::milliarcsecond: val *= static_cast<T>(1000); [[fallthrough]];
    case UAngle::arcsecond: val *= static_cast<T>(60); [[fallthrough]];
    case UAngle::arcminute: val *= static_cast<T>(60); [[fallthrough]];
    case UAngle::degree: val *= rad_to_deg; break;
    }

    return val;
}

template <typename T>
constexpr T length_factor(ULength u) {
    switch (u) {
    case ULength::kilometer: return static_cast<T>(1000.0);
    case ULength::meter: return static_cast<T>(1.0); // base units
    case ULength::centimeter: return static_cast<T>(0.01);
    case ULength::millimeter: return static_cast<T>(0.001);
    case ULength::nanometer: return static_cast<T>(1e-9);

    case ULength::inch: return static_cast<T>(0.0254);
    case ULength::foot: return static_cast<T>(0.3048);
    case ULength::mile: return static_cast<T>(1609.344);

    case ULength::au: return static_cast<T>(1.495978707e11);
    }
    return static_cast<T>(0); // unreachable
}

template <typename T>
T convert_length(T val, ULength uin, ULength uout) {
    return val * length_factor<T>(uin) / length_factor<T>(uout);
}