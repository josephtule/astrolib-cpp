// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "util/constants.hpp"
#include "util/typedefs.hpp"

enum struct UAngle : i32 { radian, degree, arcminute, arcsecond, milliarcsecond };
inline string uangle_str(UAngle type){
    switch (type) {
    case UAngle::radian: return "radian";
    case UAngle::degree: return "degree";
    case UAngle::arcminute: return "arcminute";
    case UAngle::arcsecond: return "arcsecond";
    case UAngle::milliarcsecond: return "milliarcsecond";
    }

    return "unknown";
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
inline string ulength_str(ULength type) {
    switch (type) {
    case ULength::nanometer: return "nanometer";
    case ULength::millimeter: return "millimeter";
    case ULength::centimeter: return "centimeter";
    case ULength::meter: return "meter";
    case ULength::kilometer: return "kilometer";
    case ULength::inch: return "inch";
    case ULength::foot: return "foot";
    case ULength::mile: return "mile";
    case ULength::au: return "au";
    }

    return "unknown";
}
inline constexpr ULength length_default = ULength::kilometer;

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
inline string utime_str(UTime type) {
    switch (type) {
    case UTime::year: return "year";
    case UTime::month: return "month";
    case UTime::day: return "day";
    case UTime::hour: return "hour";
    case UTime::minute: return "minute";
    case UTime::second: return "second";
    case UTime::millisecond: return "millisecond";
    case UTime::microsecond: return "microsecond";
    case UTime::nanosecond: return "nanosecond";
    }

    return "unknown";
}
inline constexpr UTime time_default = UTime::second;

template <typename T>
constexpr bool time_factor(UTime unit, T& factor) {
    // multiplicative scale factor to convert to seconds
    switch (unit) {
    case UTime::day: factor = static_cast<T>(86400.0); break;
    case UTime::hour: factor = static_cast<T>(3600.0); break;
    case UTime::minute: factor = static_cast<T>(60.0); break;
    case UTime::second: factor = static_cast<T>(1.0); break;
    case UTime::millisecond: factor = static_cast<T>(1e-3); break;
    case UTime::microsecond: factor = static_cast<T>(1e-6); break;
    case UTime::nanosecond: factor = static_cast<T>(1e-9); break;
    case UTime::year:
    case UTime::month: return false;
    }

    return true;
}

template <typename T>
bool convert_time(T val, UTime uin, UTime uout, T& converted) {
    if (uin == uout) {
        converted = val;
        return true;
    }

    T factor_in;
    T factor_out;
    if (!time_factor(uin, factor_in) || !time_factor(uout, factor_out)) return false;

    converted = val * factor_in / factor_out;
    return true;
}

enum struct UMass : i32 {
    MICROGRAM,
    MILLIGRAM,
    GRAM,
    KILOGRAM,
    POUNDMASS,
    OUNCEMASS,
};
inline constexpr UMass mass_default = UMass::KILOGRAM;

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
    // multiplicative scale factor to convert to meter
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
