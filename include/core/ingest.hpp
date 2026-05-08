#pragma once

#include "core/body.hpp"
#include "core/time.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

// TLE
struct TLEData {
    // computed values
    JulianDate jd_utc;
    f64 eccen_anom;

    // line 0
    std::string name;

    // line 1
    i32 sat_num;
    i32 launch_year;
    i32 launch_num;
    std::string launch_piece;
    i32 epoch_year;
    f64 epoch_day_frac;
    f64 d_mean_motion;
    f64 dd_mean_motion;
    f64 b_star;
    i32 ephemeris_type;
    i32 element_set_number;

    // line 2
    f64 inc;
    f64 raan;
    f64 ecc;
    f64 aop;
    f64 mean_anom;
    f64 mean_motion;
    i32 rev;

    UAngle units_angle;
};

bool read_TLE_single(
    const std::string& filename,
    Satellite& sat,
    JulianDate& jd,
    f64 mu,
    i32 millenium = 2000,
    i32 lineskips = 0,
    i32 skip_sats = 0,
    UAngle angle_out = UAngle::radian
);

bool read_TLE_multiple(
    const std::string& filename,
    svec<Satellite>& sats,
    svec<JulianDate>& jds,
    f64 mu,
    i32 millenium = 2000,
    i32 skip_sats = 0,
    i32 num_sats = 1,
    i32 lineskips = 0,
    UAngle angle_out = UAngle::radian
);

bool read_TLE_index(
    const std::string& filename,
    svec<Satellite>& sats,
    svec<JulianDate>& jds,
    const svec<i32>& idx, // may need to sort
    f64 mu,
    i32 millenium = 2000,
    i32 lineskips = 0,
    UAngle angle_out = UAngle::radian
);