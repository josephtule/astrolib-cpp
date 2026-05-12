#pragma once

#include "core/body.hpp"
#include "core/orbital_elements.hpp"
#include "core/time.hpp"
#include "util/tools.hpp"
#include "util/units.hpp"
#include <memory>

enum struct TLEStatus {
    ok,
    solver_failed,
    file_not_found,
    checksum_failed,
    non_numeric_entry,
    line_mismatch,
    satellite_mismatch,
    satellite_not_found,
    empty_file,
    empty_line,
    unknown_line,
    invalid_millennium,
    invalid_lineskips,
    invalid_mu,
    count_too_large,
}; // TODO: use these instead of boolean returns

struct TLEReadOptions {
    i32 millennium = 2000;
    i32 lineskips = 0;
    bool convert = false;
    f64 mu = 0.0;
    UAngle angle_out = UAngle::radian; // library uses radians internally
};

// TLE (Earth only)
struct TLEData {
    // computed values
    JulianDate jd_utc;
    f64 sma = 0.0;
    f64 n_rad_s = 0.0;
    f64 d_n_rad_s2 = 0.0;
    f64 dd_n_rad_s3 = 0.0;
    f64 eccen_anom = 0.0;
    f64 ta = 0.0;

    // line 0
    std::string name = "";

    // line 1
    std::string sat_id = "";
    i32 sat_num = 0;
    i32 launch_year = 0;
    i32 launch_num = 0;
    std::string launch_piece = "";
    i32 epoch_year = 0;
    f64 epoch_day_frac = 0.0;
    f64 d_mean_motion_rev_day2 = 0.0;
    f64 dd_mean_motion_rev_day3 = 0.0;
    f64 b_star = 0.0;
    i32 ephemeris_type = 0;
    i32 element_set_number = 0;

    // line 2
    f64 inc = 0.0;
    f64 raan = 0.0;
    f64 ecc = 0.0;
    f64 aop = 0.0;
    f64 mean_anom = 0.0;
    f64 mean_motion_rev_day = 0.0;
    i32 rev = 0;

    bool converted = false;
    UAngle units_angle = UAngle::degree; // raw angles are in degrees
};

TLEStatus convert_TLE(TLEData& tle, const TLEReadOptions& opts);

TLEStatus read_tle_data_single(
    const std::string& filename,
    TLEData& tle,
    const TLEReadOptions& opts = TLEReadOptions{}
);
TLEStatus read_tle_data_all(
    const std::string& filename,
    svec<TLEData>& tles,
    const TLEReadOptions& opts = TLEReadOptions{}
);
TLEStatus read_tle_data_count(
    const std::string& filename,
    svec<TLEData>& tles,
    i32 count,
    const TLEReadOptions& opts = TLEReadOptions{}
);
TLEStatus read_tle_data_single_satnum(
    const std::string& filename,
    TLEData& tle,
    i32 sat_num,
    const TLEReadOptions& opts = TLEReadOptions{}
);
TLEStatus read_tle_data_single_satid(
    const std::string& filename,
    TLEData& tle,
    const std::string& sat_id,
    const TLEReadOptions& opts = TLEReadOptions{}
);
TLEStatus read_tle_data_index(
    const std::string& filename,
    svec<TLEData>& tles,
    const svec<i32>& idx,
    const TLEReadOptions& opts = TLEReadOptions{}
);
TLEStatus read_tle_data_satnums(
    const std::string& filename,
    svec<TLEData>& tles,
    const svec<i32>& sat_nums,
    const TLEReadOptions& opts = TLEReadOptions{}
);
TLEStatus read_tle_data_satids(
    const std::string& filename,
    svec<TLEData>& tles,
    const svec<std::string>& sat_ids,
    const TLEReadOptions& opts = TLEReadOptions{}
);
TLEStatus sat_from_tle_data(
    Satellite& sat,
    const TLEData& tle,
    const TLEReadOptions& opts
);


inline OEClassical coe_from_tle(const TLEData& tle) {
    return OEClassical{
        .sma = tle.sma,
        .ecc = tle.ecc,
        .inc = tle.inc,
        .raan = tle.raan,
        .aop = tle.aop,
        .ta = tle.ta
    };
}

inline bool tle_checksum(const std::string& raw_line) {
    // get rid of returns
    std::string line = remove_returns(raw_line);

    if (line.size() < 69) return false;

    char checksum_char = line[68];
    if (checksum_char < '0' || checksum_char > '9') return false;

    int sum = 0;
    // tle line length always 69, last index is checksum, skip in sum
    for (std::size_t i = 0; i < 68; ++i) {
        char c = line[i];

        if (c >= '0' && c <= '9') {
            sum += c - '0';
        } else if (c == '-') {
            sum += 1;
        }
    }

    return (sum % 10) == (checksum_char - '0');
}


std::string tle_status_string(const TLEStatus& status) ;
