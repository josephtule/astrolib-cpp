#pragma once

#include "core/body.hpp"
#include "core/time.hpp"
#include "util/units.hpp"

// TLE (Earth only)
struct TLEData {
    // computed values
    JulianDate jd_utc;
    f64 sma = 0.0;
    f64 eccen_anom = 0.0;
    f64 ta = 0.0;

    // line 0
    std::string name = "";

    // line 1
    i32 sat_num = 0;
    i32 launch_year = 0;
    i32 launch_num = 0;
    std::string launch_piece = "";
    i32 epoch_year = 0;
    f64 epoch_day_frac = 0.0;
    f64 d_mean_motion = 0.0;
    f64 dd_mean_motion = 0.0;
    f64 b_star = 0.0;
    i32 ephemeris_type = 0;
    i32 element_set_number = 0;

    // line 2
    f64 inc = 0.0;
    f64 raan = 0.0;
    f64 ecc = 0.0;
    f64 aop = 0.0;
    f64 mean_anom = 0.0;
    f64 mean_motion = 0.0;
    i32 rev = 0;

    UAngle units_angle = UAngle::radian;
};

bool read_TLE_single(
    const std::string& filename,
    TLEData& tle,
    f64 mu,
    i32 millenium,
    i32 lineskips,
    i32 skip_sats,
    UAngle angle_out = UAngle::radian
);

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

void sat_from_tle_data(Satellite &sat, const TLEData &tle, f64 mu);
Satellite sat_from_tle_data(const TLEData& tle, f64 mu);