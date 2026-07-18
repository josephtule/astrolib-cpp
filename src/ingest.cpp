// Copyright 2025-2026 Joseph Tu Le
// SPDX-License-Identifier: Apache-2.0

#include "core/ingest.hpp"
#include "core/astrodynamics.hpp"
#include "core/body.hpp"
#include "core/orbital_elements.hpp"
#include "core/time.hpp"
#include "util/tools.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

// TODO: use TLEData struct

bool parse_optional_tle_int(const std::string& field, i32& out) {
    // launch year/num optional
    std::string val = trim(field);
    if (val.empty()) {
        out = 0;
        return true;
    }
    if (!is_numeric(val)) return false;

    out = std::stoi(val);
    return true;
}

bool parse_tle_compact_exp(
    const std::string& frac_field,
    const std::string& exponent_field,
    f64& out
) {
    std::string frac = trim(frac_field);
    std::string exponent = trim(exponent_field);

    if (frac.empty()) {
        out = 0.0;
        return true;
    }
    if (exponent.empty() || !is_numeric(exponent)) return false;

    f64 sign = 1.0;
    if (frac.front() == '+' || frac.front() == '-') {
        sign = frac.front() == '-' ? -1.0 : 1.0;
        frac.erase(frac.begin());
    }
    if (frac.empty() || frac.find_first_not_of("0123456789") != std::string::npos) {
        return false;
    }

    out = sign * std::stod("0." + frac) * std::pow(10.0, std::stoi(exponent));
    return true;
}

TLEStatus read_line_one(const std::string& line, TLEData& tle) {
    if (!tle_checksum(line)) return TLEStatus::checksum_failed;

    std::string sat_num = line.substr(2, 5);
    std::string launch_year = line.substr(9, 2);
    std::string launch_num = line.substr(11, 3);
    std::string epoch_year = line.substr(18, 2);
    std::string epoch_day_frac = line.substr(20, 12);
    std::string d_mean_motion = line.substr(33, 10);
    std::string dd_frac = line.substr(44, 6);
    std::string dd_exponent = line.substr(50, 2);
    std::string bstar_frac = line.substr(53, 6);
    std::string bstar_exponent = line.substr(59, 2);
    std::string ephemeris_type = line.substr(62, 1);
    std::string element_set_num = line.substr(64, 4);

    std::string sat_id = trim(sat_num);

    if (sat_id.empty() || !is_numeric(epoch_year) || !is_numeric(epoch_day_frac)
        || !is_numeric(d_mean_motion) || !is_numeric(ephemeris_type)
        || !is_numeric(element_set_num)) {
        return TLEStatus::non_numeric_entry;
    }

    tle.sat_id = sat_id;
    tle.sat_num = is_numeric(sat_id) ? std::stoi(sat_id) : 0;
    if (!parse_optional_tle_int(launch_year, tle.launch_year)) {
        return TLEStatus::non_numeric_entry;
    }
    if (!parse_optional_tle_int(launch_num, tle.launch_num)) {
        return TLEStatus::non_numeric_entry;
    }
    tle.launch_piece = line.substr(14, 3);
    tle.epoch_year = std::stoi(epoch_year);
    tle.epoch_day_frac = std::stod(epoch_day_frac);
    tle.d_mean_motion_rev_day2 = std::stod(d_mean_motion);
    if (!parse_tle_compact_exp(dd_frac, dd_exponent, tle.dd_mean_motion_rev_day3)) {
        return TLEStatus::non_numeric_entry;
    }
    if (!parse_tle_compact_exp(bstar_frac, bstar_exponent, tle.b_star)) {
        return TLEStatus::non_numeric_entry;
    }
    tle.ephemeris_type = std::stoi(ephemeris_type);
    tle.element_set_number = std::stoi(element_set_num);

    return TLEStatus::ok;
}

TLEStatus read_line_two(const std::string& line, TLEData& tle) {
    if (!tle_checksum(line)) return TLEStatus::checksum_failed;

    std::string sat_id = trim(line.substr(2, 5));
    std::string inc = line.substr(8, 8);
    std::string raan = line.substr(17, 8);
    std::string ecc = line.substr(26, 7);
    std::string aop = line.substr(34, 8);
    std::string mean_anom = line.substr(43, 8);
    std::string mean_motion = line.substr(52, 11);
    std::string rev = line.substr(63, 5);

    if (sat_id.empty() || !is_numeric(inc) || !is_numeric(raan) || !is_numeric(ecc)
        || !is_numeric(aop) || !is_numeric(mean_anom) || !is_numeric(mean_motion)
        || !is_numeric(rev)) {
        return TLEStatus::non_numeric_entry;
    }

    if (tle.sat_id != sat_id) return TLEStatus::satellite_mismatch;

    tle.inc = std::stod(inc);
    tle.raan = std::stod(raan);
    tle.ecc = std::stod("0." + ecc);
    tle.aop = std::stod(aop);
    tle.mean_anom = std::stod(mean_anom);
    tle.mean_motion_rev_day = std::stod(mean_motion);
    tle.rev = std::stoi(rev);

    return TLEStatus::ok;
}

TLEStatus validate_tle_options(const TLEReadOptions& opts) {
    if (opts.millennium != 1900 && opts.millennium != 2000) {
        return TLEStatus::invalid_millennium;
    }
    if (opts.lineskips < 0) {
        return TLEStatus::invalid_lineskips;
    }
    if (opts.convert && opts.mu <= 0) {
        return TLEStatus::invalid_mu;
    }

    return TLEStatus::ok;
}

TLEStatus finalize_tle_data(TLEData& tle, const TLEReadOptions& opts) {
    CalendarTime cal = doy_to_cal(tle.epoch_day_frac, opts.millennium + tle.epoch_year);
    tle.jd_utc = cal_to_jd(cal);

    return TLEStatus::ok;
}

TLEStatus read_tle_raw_line(const std::string& line, TLEData& tle) {
    if (line.empty()) return TLEStatus::empty_line;

    i32 linenum;

    // line 0 may not have 0
    if (std::isdigit(static_cast<unsigned char>(line[0]))) {
        linenum = std::stoi(line.substr(0, 1));
    } else {
        linenum = 0;
    }

    switch (linenum) {
    case 0: {
        if (line.size() > 2 && line[0] == '0') {
            tle.name = remove_returns(line.substr(2));
        } else {
            tle.name = remove_returns(line);
        }
    } break;
    case 1: {
        TLEStatus line1_status = read_line_one(line, tle);
        if (line1_status != TLEStatus::ok) {
            return line1_status;
        }
    } break;
    case 2: {
        TLEStatus line2_status = read_line_two(line, tle);
        if (line2_status != TLEStatus::ok) {
            return line2_status;
        }
    } break;
    default: return TLEStatus::unknown_line;
    }

    return TLEStatus::ok;
}

TLEStatus convert_TLE(TLEData& tle, const TLEReadOptions& opts) {
    // keep raw rev/day values intact; use computed rad/s values for state creation
    tle.n_rad_s = tle.mean_motion_rev_day * twopi / 86400.0;
    tle.d_n_rad_s2 = tle.d_mean_motion_rev_day2 * twopi / std::pow(86400.0, 2);
    tle.dd_n_rad_s3 = tle.dd_mean_motion_rev_day3 * twopi / std::pow(86400.0, 3);

    bool eccen_anom_ok = mean_anom_to_eccen_anom(
        tle.mean_anom,
        tle.ecc,
        tle.eccen_anom,
        UAngle::degree,
        UAngle::radian
    ); // rad
    if (!eccen_anom_ok) return TLEStatus::solver_failed;

    // true anomaly
    tle.ta = std::atan2(
        std::sqrt(1.0 - tle.ecc * tle.ecc) * std::sin(tle.eccen_anom),
        std::cos(tle.eccen_anom) - tle.ecc
    );
    tle.sma = std::cbrt(opts.mu / (tle.n_rad_s * tle.n_rad_s));

    // TLE data is in degrees
    UAngle angle_in = UAngle::degree;
    if (opts.angle_out != UAngle::degree) {
        tle.inc = convert_angle(tle.inc, angle_in, opts.angle_out);
        tle.raan = convert_angle(tle.raan, angle_in, opts.angle_out);
        tle.aop = convert_angle(tle.aop, angle_in, opts.angle_out);
        tle.mean_anom = convert_angle(tle.mean_anom, angle_in, opts.angle_out);
    }
    tle.units_angle = opts.angle_out;

    tle.converted = true;

    return TLEStatus::ok;
}

TLEStatus read_tle_data_single(
    const std::string& filename,
    TLEData& tle,
    const TLEReadOptions& opts
) {
    TLEStatus opts_status = validate_tle_options(opts);
    if (opts_status != TLEStatus::ok) return opts_status;

    std::ifstream file(filename);
    if (!file) return TLEStatus::file_not_found;

    std::string line;
    for (i32 i = 0; i < opts.lineskips && std::getline(file, line); ++i) {}

    TLEStatus have_line1 = TLEStatus::line_mismatch;

    i32 linenum = -1;
    while (std::getline(file, line)) {
        if (std::isdigit(static_cast<unsigned char>(line[0]))) {
            linenum = std::stoi(line.substr(0, 1));
        } else {
            linenum = 0;
        }

        if (linenum == 2 && have_line1 != TLEStatus::ok) {
            return TLEStatus::line_mismatch;
        }

        TLEStatus read_line_status = read_tle_raw_line(line, tle);
        if (linenum == 1 && read_line_status == TLEStatus::ok) {
            have_line1 = TLEStatus::ok;
        }
        if (read_line_status != TLEStatus::ok) return read_line_status;

        if (linenum == 2) {
            TLEStatus finalize_status = finalize_tle_data(tle, opts);
            if (finalize_status != TLEStatus::ok) return finalize_status;

            if (opts.convert) {
                TLEStatus convert_ok = convert_TLE(tle, opts);
                if (convert_ok != TLEStatus::ok) return convert_ok;
            }
            break;
        }
    }
    if (linenum == -1) {
        return TLEStatus::empty_file;
    }
    if (linenum == 1) {
        return TLEStatus::line_mismatch;
    }

    return TLEStatus::ok;
}

TLEStatus read_tle_data_all(
    const std::string& filename,
    svec<TLEData>& tles,
    const TLEReadOptions& opts
) {
    TLEStatus opts_status = validate_tle_options(opts);
    if (opts_status != TLEStatus::ok) return opts_status;

    std::ifstream file(filename);
    if (!file) return TLEStatus::file_not_found;

    std::string line;
    for (i32 i = 0; i < opts.lineskips && std::getline(file, line); ++i) {}

    TLEStatus have_line1 = TLEStatus::line_mismatch;

    i32 linenum = -1;
    TLEData tle;
    while (std::getline(file, line)) {
        if (std::isdigit(static_cast<unsigned char>(line[0]))) {
            linenum = std::stoi(line.substr(0, 1));
        } else {
            linenum = 0;
        }

        if (linenum == 2 && have_line1 != TLEStatus::ok) {
            return TLEStatus::line_mismatch;
        }

        TLEStatus read_line_status = read_tle_raw_line(line, tle);
        if (linenum == 1 && read_line_status == TLEStatus::ok) {
            have_line1 = TLEStatus::ok;
        }
        if (read_line_status != TLEStatus::ok) return read_line_status;

        if (linenum == 2) {
            TLEStatus finalize_status = finalize_tle_data(tle, opts);
            if (finalize_status != TLEStatus::ok) return finalize_status;

            if (opts.convert) {
                TLEStatus convert_ok = convert_TLE(tle, opts);
                if (convert_ok != TLEStatus::ok) return convert_ok;
            }
            tles.push_back(tle);
            tle = TLEData{};
            have_line1 = TLEStatus::line_mismatch;
        }
    }
    if (linenum == -1) {
        return TLEStatus::empty_file;
    }
    if (linenum == 1) {
        return TLEStatus::line_mismatch;
    }

    return TLEStatus::ok;
}

TLEStatus read_tle_data_count(
    const std::string& filename,
    svec<TLEData>& tles,
    i32 count,
    const TLEReadOptions& opts
) {
    TLEStatus opts_status = validate_tle_options(opts);
    if (opts_status != TLEStatus::ok) return opts_status;

    std::ifstream file(filename);
    if (!file) return TLEStatus::file_not_found;

    std::string line;
    for (i32 i = 0; i < opts.lineskips && std::getline(file, line); ++i) {}

    TLEStatus have_line1 = TLEStatus::line_mismatch;

    i32 linenum = -1;
    i32 sats_seen = 0;
    TLEData tle;
    while (std::getline(file, line) && sats_seen < count) {
        if (std::isdigit(static_cast<unsigned char>(line[0]))) {
            linenum = std::stoi(line.substr(0, 1));
        } else {
            linenum = 0;
        }

        if (linenum == 2 && have_line1 != TLEStatus::ok) {
            return TLEStatus::line_mismatch;
        }

        TLEStatus read_line_status = read_tle_raw_line(line, tle);
        if (linenum == 1 && read_line_status == TLEStatus::ok) {
            have_line1 = TLEStatus::ok;
        }
        if (read_line_status != TLEStatus::ok) return read_line_status;

        if (linenum == 2) {
            TLEStatus finalize_status = finalize_tle_data(tle, opts);
            if (finalize_status != TLEStatus::ok) return finalize_status;

            if (opts.convert) {
                TLEStatus convert_ok = convert_TLE(tle, opts);
                if (convert_ok != TLEStatus::ok) return convert_ok;
            }
            tles.push_back(tle);
            tle = TLEData{};
            have_line1 = TLEStatus::line_mismatch;
            ++sats_seen;
        }
    }
    if (linenum == -1) {
        return TLEStatus::empty_file;
    }
    if (linenum == 1) {
        return TLEStatus::line_mismatch;
    }
    if (sats_seen < count) {
        return TLEStatus::count_too_large;
    }

    return TLEStatus::ok;
}

TLEStatus read_tle_data_single_satnum(
    const std::string& filename,
    TLEData& tle_out,
    i32 sat_num,
    const TLEReadOptions& opts
) {
    TLEStatus opts_status = validate_tle_options(opts);
    if (opts_status != TLEStatus::ok) return opts_status;
    if (sat_num <= 0) return TLEStatus::satellite_not_found;

    std::ifstream file(filename);
    if (!file) return TLEStatus::file_not_found;

    std::string line;
    for (i32 i = 0; i < opts.lineskips && std::getline(file, line); ++i) {}

    TLEStatus have_line1 = TLEStatus::line_mismatch;

    i32 linenum = -1;
    TLEData tle;
    while (std::getline(file, line)) {
        if (std::isdigit(static_cast<unsigned char>(line[0]))) {
            linenum = std::stoi(line.substr(0, 1));
        } else {
            linenum = 0;
        }

        if (linenum == 2 && have_line1 != TLEStatus::ok) {
            return TLEStatus::line_mismatch;
        }

        TLEStatus read_line_status = read_tle_raw_line(line, tle);
        if (linenum == 1 && read_line_status == TLEStatus::ok) {
            have_line1 = TLEStatus::ok;
        }
        if (read_line_status != TLEStatus::ok) return read_line_status;

        if (linenum == 2) {
            TLEStatus finalize_status = finalize_tle_data(tle, opts);
            if (finalize_status != TLEStatus::ok) return finalize_status;

            if (tle.sat_num == sat_num) {
                if (opts.convert) {
                    TLEStatus convert_ok = convert_TLE(tle, opts);
                    if (convert_ok != TLEStatus::ok) return convert_ok;
                }

                tle_out = tle;
                return TLEStatus::ok;
            }

            tle = TLEData{};
            have_line1 = TLEStatus::line_mismatch;
        }
    }
    if (linenum == -1) {
        return TLEStatus::empty_file;
    }
    if (linenum == 1) {
        return TLEStatus::line_mismatch;
    }

    return TLEStatus::satellite_not_found;
}

TLEStatus read_tle_data_single_satid(
    const std::string& filename,
    TLEData& tle_out,
    const std::string& sat_id,
    const TLEReadOptions& opts
) {
    TLEStatus opts_status = validate_tle_options(opts);
    if (opts_status != TLEStatus::ok) return opts_status;

    std::string target_id = trim(sat_id);
    if (target_id.empty()) return TLEStatus::satellite_not_found;

    std::ifstream file(filename);
    if (!file) return TLEStatus::file_not_found;

    std::string line;
    for (i32 i = 0; i < opts.lineskips && std::getline(file, line); ++i) {}

    TLEStatus have_line1 = TLEStatus::line_mismatch;

    i32 linenum = -1;
    TLEData tle;
    while (std::getline(file, line)) {
        if (std::isdigit(static_cast<unsigned char>(line[0]))) {
            linenum = std::stoi(line.substr(0, 1));
        } else {
            linenum = 0;
        }

        if (linenum == 2 && have_line1 != TLEStatus::ok) {
            return TLEStatus::line_mismatch;
        }

        TLEStatus read_line_status = read_tle_raw_line(line, tle);
        if (linenum == 1 && read_line_status == TLEStatus::ok) {
            have_line1 = TLEStatus::ok;
        }
        if (read_line_status != TLEStatus::ok) return read_line_status;

        if (linenum == 2) {
            TLEStatus finalize_status = finalize_tle_data(tle, opts);
            if (finalize_status != TLEStatus::ok) return finalize_status;

            if (tle.sat_id == target_id) {
                if (opts.convert) {
                    TLEStatus convert_ok = convert_TLE(tle, opts);
                    if (convert_ok != TLEStatus::ok) return convert_ok;
                }

                tle_out = tle;
                return TLEStatus::ok;
            }

            tle = TLEData{};
            have_line1 = TLEStatus::line_mismatch;
        }
    }
    if (linenum == -1) {
        return TLEStatus::empty_file;
    }
    if (linenum == 1) {
        return TLEStatus::line_mismatch;
    }

    return TLEStatus::satellite_not_found;
}

TLEStatus read_tle_data_index(
    const std::string& filename,
    svec<TLEData>& tles,
    const svec<i32>& idx,
    const TLEReadOptions& opts
) {
    TLEStatus opts_status = validate_tle_options(opts);
    if (opts_status != TLEStatus::ok) return opts_status;

    if (idx.empty()) return TLEStatus::ok;

    svec<i32> idx_sorted = idx;
    std::sort(idx_sorted.begin(), idx_sorted.end());
    idx_sorted.erase(std::unique(idx_sorted.begin(), idx_sorted.end()), idx_sorted.end());

    if (idx_sorted.front() < 0) return TLEStatus::satellite_not_found;

    std::ifstream file(filename);
    if (!file) return TLEStatus::file_not_found;

    std::string line;
    for (i32 i = 0; i < opts.lineskips && std::getline(file, line); ++i) {}

    TLEStatus have_line1 = TLEStatus::line_mismatch;

    i32 linenum = -1;
    i32 sats_seen = 0;
    i32 idx_i = 0;
    TLEData tle;
    while (std::getline(file, line)) {
        if (std::isdigit(static_cast<unsigned char>(line[0]))) {
            linenum = std::stoi(line.substr(0, 1));
        } else {
            linenum = 0;
        }

        if (linenum == 2 && have_line1 != TLEStatus::ok) {
            return TLEStatus::line_mismatch;
        }

        TLEStatus read_line_status = read_tle_raw_line(line, tle);
        if (linenum == 1 && read_line_status == TLEStatus::ok) {
            have_line1 = TLEStatus::ok;
        }
        if (read_line_status != TLEStatus::ok) return read_line_status;

        if (linenum == 2) {
            TLEStatus finalize_status = finalize_tle_data(tle, opts);
            if (finalize_status != TLEStatus::ok) return finalize_status;

            while (idx_i < static_cast<i32>(idx_sorted.size())
                   && idx_sorted[idx_i] < sats_seen) {
                ++idx_i;
            }

            while (idx_i < static_cast<i32>(idx_sorted.size())
                   && idx_sorted[idx_i] == sats_seen) {
                TLEData tle_selected = tle;
                if (opts.convert) {
                    TLEStatus convert_ok = convert_TLE(tle_selected, opts);
                    if (convert_ok != TLEStatus::ok) return convert_ok;
                }

                tles.push_back(tle_selected);
                ++idx_i;
            }

            if (idx_i >= static_cast<i32>(idx_sorted.size())) return TLEStatus::ok;

            tle = TLEData{};
            have_line1 = TLEStatus::line_mismatch;
            ++sats_seen;
        }
    }
    if (linenum == -1) {
        return TLEStatus::empty_file;
    }
    if (linenum == 1) {
        return TLEStatus::line_mismatch;
    }

    return TLEStatus::satellite_not_found;
}

TLEStatus read_tle_data_satnums(
    const std::string& filename,
    svec<TLEData>& tles,
    const svec<i32>& sat_nums,
    const TLEReadOptions& opts
) {
    TLEStatus opts_status = validate_tle_options(opts);
    if (opts_status != TLEStatus::ok) return opts_status;

    if (sat_nums.empty()) return TLEStatus::ok;

    svec<i32> sat_nums_sorted;
    for (i32 sat_num : sat_nums) {
        if (sat_num > 0) sat_nums_sorted.push_back(sat_num);
    }
    if (sat_nums_sorted.empty()) return TLEStatus::satellite_not_found;

    std::sort(sat_nums_sorted.begin(), sat_nums_sorted.end());
    sat_nums_sorted.erase(
        std::unique(sat_nums_sorted.begin(), sat_nums_sorted.end()),
        sat_nums_sorted.end()
    );

    std::ifstream file(filename);
    if (!file) return TLEStatus::file_not_found;

    std::string line;
    for (i32 i = 0; i < opts.lineskips && std::getline(file, line); ++i) {}

    TLEStatus have_line1 = TLEStatus::line_mismatch;

    i32 linenum = -1;
    TLEData tle;
    while (std::getline(file, line)) {
        if (std::isdigit(static_cast<unsigned char>(line[0]))) {
            linenum = std::stoi(line.substr(0, 1));
        } else {
            linenum = 0;
        }

        if (linenum == 2 && have_line1 != TLEStatus::ok) {
            return TLEStatus::line_mismatch;
        }

        TLEStatus read_line_status = read_tle_raw_line(line, tle);
        if (linenum == 1 && read_line_status == TLEStatus::ok) {
            have_line1 = TLEStatus::ok;
        }
        if (read_line_status != TLEStatus::ok) return read_line_status;

        if (linenum == 2) {
            TLEStatus finalize_status = finalize_tle_data(tle, opts);
            if (finalize_status != TLEStatus::ok) return finalize_status;

            auto sat_it = std::lower_bound(
                sat_nums_sorted.begin(),
                sat_nums_sorted.end(),
                tle.sat_num
            );

            if (sat_it != sat_nums_sorted.end() && *sat_it == tle.sat_num) {
                TLEData tle_selected = tle;
                if (opts.convert) {
                    TLEStatus convert_ok = convert_TLE(tle_selected, opts);
                    if (convert_ok != TLEStatus::ok) return convert_ok;
                }

                tles.push_back(tle_selected);
                sat_nums_sorted.erase(sat_it);
            }

            if (sat_nums_sorted.empty()) return TLEStatus::ok;

            tle = TLEData{};
            have_line1 = TLEStatus::line_mismatch;
        }
    }
    if (linenum == -1) {
        return TLEStatus::empty_file;
    }
    if (linenum == 1) {
        return TLEStatus::line_mismatch;
    }

    return TLEStatus::satellite_not_found;
}

TLEStatus read_tle_data_satids(
    const std::string& filename,
    svec<TLEData>& tles,
    const svec<std::string>& sat_ids,
    const TLEReadOptions& opts
) {
    TLEStatus opts_status = validate_tle_options(opts);
    if (opts_status != TLEStatus::ok) return opts_status;

    if (sat_ids.empty()) return TLEStatus::ok;

    svec<std::string> sat_ids_sorted;
    for (const std::string& sat_id : sat_ids) {
        std::string id = trim(sat_id);
        if (!id.empty()) sat_ids_sorted.push_back(id);
    }
    if (sat_ids_sorted.empty()) return TLEStatus::satellite_not_found;

    std::sort(sat_ids_sorted.begin(), sat_ids_sorted.end());
    sat_ids_sorted.erase(
        std::unique(sat_ids_sorted.begin(), sat_ids_sorted.end()),
        sat_ids_sorted.end()
    );

    std::ifstream file(filename);
    if (!file) return TLEStatus::file_not_found;

    std::string line;
    for (i32 i = 0; i < opts.lineskips && std::getline(file, line); ++i) {}

    TLEStatus have_line1 = TLEStatus::line_mismatch;

    i32 linenum = -1;
    TLEData tle;
    while (std::getline(file, line)) {
        if (std::isdigit(static_cast<unsigned char>(line[0]))) {
            linenum = std::stoi(line.substr(0, 1));
        } else {
            linenum = 0;
        }

        if (linenum == 2 && have_line1 != TLEStatus::ok) {
            return TLEStatus::line_mismatch;
        }

        TLEStatus read_line_status = read_tle_raw_line(line, tle);
        if (linenum == 1 && read_line_status == TLEStatus::ok) {
            have_line1 = TLEStatus::ok;
        }
        if (read_line_status != TLEStatus::ok) return read_line_status;

        if (linenum == 2) {
            TLEStatus finalize_status = finalize_tle_data(tle, opts);
            if (finalize_status != TLEStatus::ok) return finalize_status;

            auto sat_it = std::lower_bound(
                sat_ids_sorted.begin(),
                sat_ids_sorted.end(),
                tle.sat_id
            );

            if (sat_it != sat_ids_sorted.end() && *sat_it == tle.sat_id) {
                TLEData tle_selected = tle;
                if (opts.convert) {
                    TLEStatus convert_ok = convert_TLE(tle_selected, opts);
                    if (convert_ok != TLEStatus::ok) return convert_ok;
                }

                tles.push_back(tle_selected);
                sat_ids_sorted.erase(sat_it);
            }

            if (sat_ids_sorted.empty()) return TLEStatus::ok;

            tle = TLEData{};
            have_line1 = TLEStatus::line_mismatch;
        }
    }
    if (linenum == -1) {
        return TLEStatus::empty_file;
    }
    if (linenum == 1) {
        return TLEStatus::line_mismatch;
    }

    return TLEStatus::satellite_not_found;
}

TLEStatus sat_from_tle_data(
    Satellite& sat,
    const TLEData& tle,
    const TLEReadOptions& opts
) {
    TLEData tle_converted = tle;
    if (!tle.converted) {
        TLEStatus convert_ok = convert_TLE(tle_converted, opts);
        if (convert_ok != TLEStatus::ok) {
            return convert_ok;
        }
    }

    sat.name = tle.name + " (" + tle.sat_id + "-" + std::to_string(tle.launch_year)
               + std::to_string(tle.launch_num) + trim(tle.launch_piece) + ")";

    OEClassical coe = coe_from_tle(tle_converted);
    sat.x_tr = classical_to_rv(coe, opts.mu, tle_converted.units_angle);

    return TLEStatus::ok;
}

TLEStatus sats_from_tle_data(
    svec<uptr<Satellite>>& sats,
    const svec<TLEData>& tles,
    const TLEReadOptions& opts
) {
    sats.clear();
    sats.reserve(tles.size());

    for (const TLEData& tle : tles) {
        auto sat = std::make_unique<Satellite>();
        TLEStatus sat_status = sat_from_tle_data(*sat, tle, opts);
        if (sat_status != TLEStatus::ok) {
            sats.clear();
            return sat_status;
        }

        sats.emplace_back(std::move(sat));
    }

    return TLEStatus::ok;
}

std::string tle_status_string(const TLEStatus& status) {
    switch (status) {
    case TLEStatus::ok: return "ok";
    case TLEStatus::solver_failed: return "solver_failed";
    case TLEStatus::file_not_found: return "file_not_found";
    case TLEStatus::checksum_failed: return "checksum_failed";
    case TLEStatus::non_numeric_entry: return "non_numeric_entry";
    case TLEStatus::line_mismatch: return "line_mismatch";
    case TLEStatus::satellite_mismatch: return "satellite_mismatch";
    case TLEStatus::satellite_not_found: return "satellite_not_found";
    case TLEStatus::empty_file: return "empty_file";
    case TLEStatus::empty_line: return "empty_line";
    case TLEStatus::unknown_line: return "unknown_line";
    case TLEStatus::invalid_millennium: return "invalid_millennium";
    case TLEStatus::invalid_lineskips: return "invalid_lineskips";
    case TLEStatus::invalid_mu: return "invalid_mu";
    case TLEStatus::count_too_large: return "count_too_large";
    }

    return "unknown_status";
}
