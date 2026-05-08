#include "core/ingest.hpp"
#include "core/body.hpp"
#include "core/orbital_elements.hpp"
#include "core/time.hpp"
#include "util/astrodynamics.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>

// TODO: use TLEData struct

bool read_line_one(std::string line, TLEData& tle) {
    try {
        tle.sat_num = std::stoi(line.substr(2, 5));
        tle.launch_year = std::stoi(line.substr(9, 2));
        tle.launch_num = std::stoi(line.substr(11, 3));
        tle.launch_piece = line.substr(14, 3);
        tle.epoch_year = std::stoi(line.substr(18, 2));
        tle.epoch_day_frac = std::stod(line.substr(20, 12));
        tle.d_mean_motion = std::stod(line.substr(33, 10));
        tle.dd_mean_motion = std::stod("0." + line.substr(44, 6))
                             * std::pow(10, std::stod(line.substr(50, 2)));
        tle.b_star = std::stod("0." + line.substr(54, 6))
                     * std::pow(10, std::stod(line.substr(59, 2)));

        tle.ephemeris_type = std::stoi(line.substr(62, 1));
        tle.element_set_number = std::stoi(line.substr(64, 4));
    } catch (const std::invalid_argument&) {
        return false;
    } catch (const std::out_of_range&) {
        return false;
    }

    return true;
}

bool read_line_two(std::string line, TLEData& tle) {
    try {
        i32 sat_num_2 = std::stoi(line.substr(2, 5));
        if (tle.sat_num != sat_num_2) return false;

        tle.inc = std::stod(line.substr(8, 8));
        tle.raan = std::stod(line.substr(17, 8));
        tle.ecc = std::stod("0." + line.substr(26, 7));
        tle.aop = std::stod(line.substr(34, 8));
        tle.mean_anom = std::stod(line.substr(43, 8));
        tle.mean_motion = std::stod(line.substr(52, 11));
        tle.rev = std::stoi(line.substr(63, 5));
    } catch (const std::invalid_argument&) {
        return false;
    } catch (const std::out_of_range&) {
        return false;
    }

    return true;
}

bool read_TLE_raw_single(
    const std::string& filename,
    TLEData& tle,
    f64 mu,
    i32 millennium,
    i32 lineskips,
    i32 skip_sats
) {
    if (millennium != 1900 && millennium != 2000) return false;

    // reads tle file with only one satellite entry
    std::ifstream file(filename);
    if (!file) return false;

    std::string line;
    for (i32 i = 0; i < lineskips && std::getline(file, line); ++i) {}

    bool found_sat = false;

    i32 sats_seen = 0;
    bool have_line1 = false;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        i32 linenum;

        // line 0 may not have 0
        if (std::isdigit(static_cast<unsigned char>(line[0]))) {
            linenum = std::stoi(line.substr(0, 1));
        } else {
            linenum = 0;
        }

        switch (linenum) {
        case 0: {
            if (sats_seen == skip_sats) {
                if (line.size() > 2 && line[0] == '0')
                    tle.name = line.substr(2);
                else
                    tle.name = line;
            }
        } break;
        case 1: {
            have_line1 = true;

            if (sats_seen == skip_sats) {
                if (!read_line_one(line, tle)) return false;
            }
        } break;
        case 2: {
            if (!have_line1) return false;

            if (sats_seen == skip_sats) {
                if (!read_line_two(line, tle)) return false;

                found_sat = true;
                break;
            }

            // sat complete once line 2 read
            ++sats_seen;
            have_line1 = false;
        } break;
        default: break;
        }

        if (found_sat) break;
    }
    if (!found_sat) return false;

    CalendarTime cal = doy_to_cal(tle.epoch_day_frac, millennium + tle.epoch_year);
    tle.jd_utc = cal_to_jd(cal);

    return true;
}

bool read_TLE_converted_single(
    const std::string& filename,
    TLEData& tle,
    f64 mu,
    i32 millennium,
    i32 lineskips,
    i32 skip_sats,
    UAngle angle_out
) {
    bool tle_ok
        = read_TLE_raw_single(filename, tle, mu, millennium, lineskips, skip_sats);
    if (!tle_ok) return false;

    bool tle_convert_ok = convert_TLE(tle, mu, angle_out);
    if (!tle_convert_ok) return false;

    return true;
}
bool read_TLE_converted_single(
    const std::string& filename,
    Satellite& sat,
    JulianDate& jd_utc,
    f64 mu,
    i32 millennium,
    i32 lineskips,
    i32 skip_sats,
    UAngle angle_out
) {
    TLEData tle;
    bool read_ok = read_TLE_converted_single(
        filename,
        tle,
        mu,
        millennium,
        lineskips,
        skip_sats,
        angle_out
    );
    if (!read_ok) return false;

    sat = sat_from_tle_data(tle, mu);
    jd_utc = tle.jd_utc;

    return true;
}

bool convert_TLE(TLEData& tle, f64 mu, UAngle angle_out) {
    // convert from rev/day^i to deg/s^i, assume 86400s day
    tle.mean_motion *= 360.0 / 86400.0;
    tle.d_mean_motion *= 360.0 / std::pow(86400.0, 2);
    tle.dd_mean_motion *= 360.0 / std::pow(86400.0, 3);

    bool eccen_anom_ok = mean_anom_to_eccen_anom(
        tle.mean_anom,
        tle.ecc,
        tle.eccen_anom,
        UAngle::degree,
        UAngle::radian
    ); // rad
    if (!eccen_anom_ok) return false;

    // true anomaly
    f64 ta = std::atan2(
        std::sqrt(1.0 - tle.ecc * tle.ecc) * std::sin(tle.eccen_anom),
        std::cos(tle.eccen_anom) - tle.ecc
    );
    f64 n = convert_angle(
        tle.mean_motion,
        UAngle::degree,
        UAngle::radian
    ); // mean motion in rad/s
    tle.sma = std::cbrt(mu / (n * n));

    // TLE data is in degrees
    UAngle angle_in = UAngle::degree;
    if (angle_out != UAngle::degree) {
        tle.inc = convert_angle(tle.inc, angle_in, angle_out);
        tle.raan = convert_angle(tle.raan, angle_in, angle_out);
        tle.aop = convert_angle(tle.aop, angle_in, angle_out);
        tle.mean_anom = convert_angle(tle.mean_anom, angle_in, angle_out);
        tle.mean_motion = convert_angle(tle.mean_motion, angle_in, angle_out);
        tle.d_mean_motion = convert_angle(tle.d_mean_motion, angle_in, angle_out);
        tle.dd_mean_motion = convert_angle(tle.dd_mean_motion, angle_in, angle_out);
    }
    tle.units_angle = angle_out;

    tle.converted = true;

    return true;
}

bool read_TLE_multiple(
    const std::string& filename,
    svec<std::unique_ptr<Satellite>>& sats,
    svec<JulianDate>& jds,
    f64 mu,
    i32 num_sats,
    i32 millennium,
    i32 skip_sats,
    i32 lineskips,
    UAngle angle_out
) {
    svec<i32> idx;
    for (i32 i = 0; i < num_sats; ++i) idx.push_back(i);
    return read_TLE_index(
        filename,
        sats,
        jds,
        idx,
        mu,
        millennium,
        lineskips,
        0,
        angle_out
    );
}

bool read_TLE_index(
    const std::string& filename,
    svec<std::unique_ptr<Satellite>>& sats,
    svec<JulianDate>& jds,
    const svec<i32>& idx,
    f64 mu,
    i32 millennium,
    i32 lineskips,
    i32 idx_start,
    UAngle angle_out
) {
    if (millennium != 1900 && millennium != 2000) return false;
    if (idx_start != 0 && idx_start != 1) return false;

    svec<i32> idx_sorted = idx;
    std::sort(idx_sorted.begin(), idx_sorted.end());

    std::ifstream file(filename);
    if (!file) return false;

    std::string line;
    for (i32 i = 0; i < lineskips && std::getline(file, line); ++i) {}

    bool found_sat = false;

    i32 sats_seen = 0;
    bool have_line1 = false;

    TLEData tle;
    i32 i = 0;
    i32 current_idx = idx_sorted[i] - idx_start;
    while (std::getline(file, line)) {

        if (line.empty()) continue;

        i32 linenum;

        // line 0 may not have 0
        if (std::isdigit(static_cast<unsigned char>(line[0]))) {
            linenum = std::stoi(line.substr(0, 1));
        } else {
            linenum = 0;
        }

        // TODO: maybe switch to using another getline
        switch (linenum) {
        case 0: {
            if (sats_seen == current_idx) {
                if (line.size() > 2 && line[0] == '0')
                    tle.name = line.substr(2);
                else
                    tle.name = line;
            }
        } break;
        case 1: {
            have_line1 = true;

            if (sats_seen == current_idx) {
                if (!read_line_one(line, tle)) return false;
            }
        } break;
        case 2: {
            if (!have_line1) return false;

            if (sats_seen == current_idx) {
                if (!read_line_two(line, tle)) return false;

                if (!convert_TLE(tle, mu)) return false;

                sats.emplace_back(
                    std::make_unique<Satellite>(sat_from_tle_data(tle, mu))
                );

                jds.push_back(tle.jd_utc);

                ++i;

                if (i >= idx_sorted.size()) {
                    found_sat = true;
                    break;
                }

                current_idx = idx_sorted[i] - idx_start;
            }

            ++sats_seen;
            have_line1 = false;
        } break;
        default: break;
        }
    }

    return true;
}

void sat_from_tle_data(Satellite& sat, const TLEData& tle, f64 mu) {
    OEClassical coe = coe_from_tle(tle);

    sat.name = tle.name;
    sat.x_tr = classical_to_rv(coe, mu, tle.units_angle);
}

Satellite sat_from_tle_data(const TLEData& tle, f64 mu) {
    Satellite sat;
    sat_from_tle_data(sat, tle, mu);

    return sat;
}
