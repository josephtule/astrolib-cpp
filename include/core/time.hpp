#pragma once

#include "util/typedefs.hpp"
#include <cmath>

enum struct TimeScale {
    utc, // coordinated universal time
    ut1, // universal time, UT1 = UT
    tai, // international atomic time
    tt,  // terrestrial time
    tdb, // dynamical barycentric time
    gps, // global positioning system time
};

struct CalendarTime {
    i32 year = 2000;
    i32 month = 1;
    i32 day = 1;
    i32 hour = 0;
    i32 minute = 0;
    f64 second = 0.0;
};
struct HMSTime {
    i32 hour = 0;
    i32 minute = 0;
    f64 second = 0.0;
};

struct JulianDate { // split for precision
    f64 day = 2451545.0;
    f64 frac = 0.0;
};

struct ModifiedJulianDate { // split for precision
    f64 day = 51544.0;
    f64 frac = 0.5;
};

struct TimeOffsets {
    f64 dut1 = 0.0; // UT1 - UTC [s]
    f64 dat = 0.0;  // TAI - UTC [s]
};

inline bool is_leap_year(i32 year) {
    if (year % 4 == 0 && year % 100 != 0 || year % 400 == 0) {
        return true;
    }
    return false;
}

const svec<i32> days_leap = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
const svec<i32> days_reg = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
inline i32 days_in_month(i32 year, i32 month) {
    if (month > 12 || month < 1) {
        return 0;
    }
    if (is_leap_year(year)) {
        return days_leap[month - 1];
    } else {
        return days_reg[month - 1];
    }
}

inline f64 seconds_of_day(i32 hour, i32 minute, f64 second) {
    return hour * 3600.0 + minute * 60.0 + second;
}

inline HMSTime frac_day_to_hms(f64 frac_day) {

    f64 sec_in_day = 24 * 3600; // assumes uniform 86400s day
    f64 total_sec = frac_day * sec_in_day;
    i32 hours = i32(total_sec / 3600);
    total_sec -= f64(hours) * 3600.0; // remaining seconds
    i32 minutes = i32(total_sec / 60);
    f64 seconds = total_sec - f64(minutes) * 60.0;

    return HMSTime{.hour = hours, .minute = minutes, .second = seconds};
}

inline JulianDate normalize_jd(JulianDate jd) {
    f64 whole = std::floor(jd.frac);
    jd.day += whole;
    jd.frac -= whole;
    return jd;
}
inline ModifiedJulianDate normalize_mjd(ModifiedJulianDate mjd) {
    f64 whole = std::floor(mjd.frac);
    mjd.day += whole;
    mjd.frac -= whole;
    return mjd;
}

inline f64 jd_to_scalar(const JulianDate& jd) { return jd.day + jd.frac; }
inline JulianDate jd_from_scalar(f64 jd_scalar) {
    JulianDate jd;
    jd.day = std::floor(jd_scalar);
    jd.frac = jd_scalar - jd.day;
    return jd;
}

inline f64 mjd_to_scalar(const ModifiedJulianDate& mjd) { return mjd.day + mjd.frac; }
inline ModifiedJulianDate mjd_from_scalar(f64 mjd_scalar) {
    ModifiedJulianDate mjd;
    mjd.day = std::floor(mjd_scalar);
    mjd.frac = mjd_scalar - mjd.day;
    return mjd;
}

inline ModifiedJulianDate jd_to_mjd(const JulianDate& jd) {
    return mjd_from_scalar(jd_to_scalar(jd) - 2400000.5);
}

inline JulianDate mjd_to_jd(const ModifiedJulianDate& mjd) {
    return jd_from_scalar(mjd_to_scalar(mjd) + 2400000.5);
}