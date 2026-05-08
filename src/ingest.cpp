#include "core/ingest.hpp"
#include "core/orbital_elements.hpp"
#include "core/time.hpp"
#include "util/astrodynamics.hpp"
#include "util/units.hpp"

#include <fstream>
#include <iostream>

// TODO: use TLEData struct

bool read_TLE_single(
    const std::string& filename,
    Satellite& sat,
    JulianDate& jd_utc,
    f64 mu,
    i32 millenium,
    i32 lineskips,
    i32 skip_sats,
    UAngle angle_out
) {
    if (millenium != 1900 && millenium != 2000) return false;

    // reads tle file with only one satellite entry
    std::ifstream file(filename);
    if (!file) return false;

    std::string line;
    for (i32 i = 0; i < lineskips && std::getline(file, line); ++i) {}

    // common tle fields
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
                    name = line.substr(2);
                else
                    name = line;
            }
        } break;

        case 1: {
            have_line1 = true;

            if (sats_seen == skip_sats) {
                sat_num = std::stoi(line.substr(2, 5));
                launch_year = std::stoi(line.substr(9, 2));
                launch_num = std::stoi(line.substr(11, 3));
                launch_piece = line.substr(14, 3);
                epoch_year = std::stoi(line.substr(18, 2));
                epoch_day_frac = std::stod(line.substr(20, 12));
                d_mean_motion = std::stod(line.substr(33, 10));

                dd_mean_motion = std::stod(line.substr(44, 6)) / std::pow(10, 5)
                                 * std::pow(10, std::stod(line.substr(50, 2)));

                b_star = std::stod(line.substr(54, 6)) / std::pow(10, 5)
                         * std::pow(10, std::stod(line.substr(59, 2)));

                ephemeris_type = std::stoi(line.substr(62, 1));
                element_set_number = std::stoi(line.substr(64, 4));
            }
        } break;

        case 2: {
            if (!have_line1) return false;

            if (sats_seen == skip_sats) {
                i32 sat_num_2 = std::stoi(line.substr(2, 5));
                if (sat_num != sat_num_2) return false;

                inc = std::stod(line.substr(8, 8));
                raan = std::stod(line.substr(17, 8));
                ecc = std::stod(line.substr(26, 7)) / std::pow(10, 7);
                aop = std::stod(line.substr(34, 8));
                mean_anom = std::stod(line.substr(43, 8));
                mean_motion = std::stod(line.substr(52, 11));
                rev = std::stoi(line.substr(63, 5));

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

    // convert from rev/day^i to deg/s^i, assume 86400s day
    mean_motion *= 360.0 / 86400.0;
    d_mean_motion *= 360.0 / std::pow(86400.0, 2);
    dd_mean_motion *= 360.0 / std::pow(86400.0, 3);

    f64 eccen_anom;
    bool eccen_anom_ok = mean_anom_to_eccen_anom(
        mean_anom,
        ecc,
        eccen_anom,
        UAngle::degree,
        UAngle::radian
    ); // rad
    if (!eccen_anom_ok) return false;

    // true anomaly
    f64 ta = std::atan2(
        std::sqrt(1.0 - ecc * ecc) * std::sin(eccen_anom),
        std::cos(eccen_anom) - ecc
    );
    f64 n = convert_angle(
        mean_motion,
        UAngle::degree,
        UAngle::radian
    ); // mean motion in rad/s
    double sma = std::cbrt(mu / (n * n));

    // TLE data is in degrees
    UAngle angle_in = UAngle::degree;
    if (angle_out != UAngle::degree) {
        inc = convert_angle(inc, angle_in, angle_out);
        raan = convert_angle(raan, angle_in, angle_out);
        aop = convert_angle(aop, angle_in, angle_out);
        mean_anom = convert_angle(mean_anom, angle_in, angle_out);
        mean_motion = convert_angle(mean_motion, angle_in, angle_out);
        d_mean_motion = convert_angle(d_mean_motion, angle_in, angle_out);
        dd_mean_motion = convert_angle(dd_mean_motion, angle_in, angle_out);
    }

    CalendarTime cal = doy_to_cal(epoch_day_frac, millenium + epoch_year);
    jd_utc = cal_to_jd(cal);

    OEClassical
        coe{.sma = sma, .ecc = ecc, .inc = inc, .raan = raan, .aop = aop, .ta = ta};

    sat.name = name;
    sat.x_tr = classical_to_rv(coe, mu, angle_out);

    return true;
}

bool read_TLE_multiple(
    const std::string& filename,
    svec<Satellite>& sats,
    svec<JulianDate>& jds,
    f64 mu,
    i32 millenium,
    i32 skip_sats,
    i32 num_sats,
    i32 lineskips,
    UAngle angle_out
) {
    if (millenium != 1900 || millenium != 2000) return false;

    // reads tle file with only one satellite entry
    std::ifstream file(filename);
    if (!file) return false;

    std::string line;
    for (i32 i = 0; i < lineskips && std::getline(file, line); ++i) {}
}
