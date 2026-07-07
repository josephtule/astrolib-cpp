#pragma once

#include "util/math.hpp"
#include "util/printing.hpp"
#include "util/tools.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include <chrono>
#include <cmath>
#include <print>
#include <string>

enum struct TimeScale {
    utc, // coordinated universal time
    ut1, // universal time, UT1 = UT
    tai, // international atomic time
    tt,  // terrestrial time
    tdb, // dynamical barycentric time
    gps, // global positioning system time
};
inline string time_scale_str(const TimeScale& scale) {
    switch (scale) {
    case TimeScale::utc: return "utc";
    case TimeScale::ut1: return "ut1";
    case TimeScale::tai: return "tai";
    case TimeScale::tt: return "tt";
    case TimeScale::tdb: return "tdb";
    case TimeScale::gps: return "gps";
    }

    return "unknown";
}

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

struct JulianDate {      // split for precision
    f64 day = 2451545.0; // J2000
    f64 frac = 0.0;
};

struct ModifiedJulianDate { // split for precision
    f64 day = 51544.0;
    f64 frac = 0.5;
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

inline f64 hms_to_frac_day(HMSTime hms) {
    f64 sec_in_day = 24 * 3600; // assumes uniform 86400s day
    f64 total_sec = seconds_of_day(hms.hour, hms.minute, hms.second);
    return total_sec / sec_in_day;
}

inline HMSTime hms_from_cal(const CalendarTime& cal) {
    return HMSTime{.hour = cal.hour, .minute = cal.minute, .second = cal.second};
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

inline JulianDate cal_to_jd(const CalendarTime& cal) {
    i32 y = cal.year;
    i32 m = cal.month;

    if (m <= 2) {
        y -= 1;
        m += 12;
    }

    i32 A = y / 100;
    i32 B = 2 - A + A / 4;
    // NOTE: this ignores Julian-Gregorian crossover

    f64 day_frac = seconds_of_day(cal.hour, cal.minute, cal.second) / 86400.0;

    f64 jd_scalar = std::floor(365.25 * (y + 4716)) + std::floor(30.6001 * (m + 1))
                    + cal.day + day_frac + B - 1524.5;

    return jd_from_scalar(jd_scalar);
}

inline CalendarTime jd_to_cal(const JulianDate& jd) {
    // https://quasar.as.utexas.edu/BillInfo/JulianDatesG.html
    // NOTE: this ignores Julian-Gregorian crossover

    f64 jd_shift = jd_to_scalar(jd) + 0.5;

    i32 Z = static_cast<i32>(std::floor(jd_shift));
    f64 F = jd_shift - static_cast<f64>(Z);

    i32 alpha = static_cast<i32>(std::floor((Z - 1867216.25) / 36524.25));
    i32 A = Z + 1 + alpha - alpha / 4;

    i32 B = A + 1524;
    i32 C = static_cast<i32>(std::floor((B - 122.1) / 365.25));
    i32 D = static_cast<i32>(std::floor(365.25 * C));
    i32 E = static_cast<i32>(std::floor((B - D) / 30.6001));

    f64 day_real = B - D - std::floor(30.6001 * E) + F;

    CalendarTime cal;
    cal.day = static_cast<i32>(std::floor(day_real));

    if (E < 14) {
        cal.month = E - 1;
    } else {
        cal.month = E - 13;
    }

    if (cal.month > 2) {
        cal.year = C - 4716;
    } else {
        cal.year = C - 4715;
    }

    f64 frac_day = day_real - std::floor(day_real);
    HMSTime hms = frac_day_to_hms(frac_day);

    cal.hour = hms.hour;
    cal.minute = hms.minute;
    cal.second = hms.second;

    return cal;
}

struct TimeOffsets {
    // https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop
    // Standard EOP data files contain [UT1-UTC]
    // https://data.iana.org/time-zones/data/leap-seconds.list
    // contains leapseconds [TAI-UTC]
    f64 ut1_utc = 0.0;   // UT1 - UTC [s]
    f64 tai_utc = 0.0;   // TAI - UTC [s], leap seconds
    f64 tt_tai = 32.184; // TT - TAI [s]
    f64 tai_gps = 19.0;  // TAI - GPS [s]
};

inline f64 time_scale_convert(
    f64 t,
    TimeScale scale_in,
    TimeScale scale_out,
    TimeOffsets offsets = TimeOffsets{}
) {
    // convert input to TAI
    switch (scale_in) {
    case TimeScale::ut1: {
        t -= offsets.ut1_utc; // UT1 -> UTC
        t += offsets.tai_utc; // UTC -> TAI
        break;
    }
    case TimeScale::utc: {
        t += offsets.tai_utc;
        break;
    }
    case TimeScale::gps: {
        t += offsets.tai_gps;
        break;
    }
    case TimeScale::tt: {
        t -= offsets.tt_tai;
        break;
    }
    case TimeScale::tai: break;
    case TimeScale::tdb: break; // TODO: not sure what to do here
    default: break;
    }

    // convert TAI to output
    switch (scale_out) {
    case TimeScale::ut1: {
        t -= offsets.tai_utc; // TAI -> UTC
        t += offsets.ut1_utc; // UTC -> UT1
        break;
    }
    case TimeScale::utc: {
        t -= offsets.tai_utc;
        break;
    }
    case TimeScale::gps: {
        t -= offsets.tai_gps;
        break;
    }
    case TimeScale::tt: {
        t += offsets.tt_tai;
        break;
    }
    case TimeScale::tai: break;
    case TimeScale::tdb: break; // TODO: not sure what to do here
    default: break;
    }

    return t;
}

inline JulianDate jd_scale_convert(
    const JulianDate& jd,
    TimeScale scale_in,
    TimeScale scale_out,
    TimeOffsets offsets = TimeOffsets{}
) {
    f64 jd_new_scalar = jd_to_scalar(jd);
    f64 dt = time_scale_convert(0.0, scale_in, scale_out, offsets) / 86400.0;
    return jd_from_scalar(jd_new_scalar + dt);
}

inline CalendarTime cal_scale_convert(
    const CalendarTime& cal,
    TimeScale scale_in,
    TimeScale scale_out,
    TimeOffsets offsets = TimeOffsets{}
) {
    return jd_to_cal(jd_scale_convert(cal_to_jd(cal), scale_in, scale_out, offsets));
}

inline CalendarTime doy_to_cal(f64 doy_frac, i32 year) {
    i32 doy = static_cast<i32>(std::floor(doy_frac));
    f64 frac_day = doy_frac - static_cast<f64>(doy);
    i32 max_doy = is_leap_year(year) ? 366 : 365;
    if (doy < 1 || doy > max_doy) {
        return CalendarTime{};
    }

    i32 month = 1;
    while (true) {
        i32 dim = days_in_month(year, month);
        if (dim < doy) {
            month++;
            doy -= dim;
        } else {
            break;
        }
    }

    HMSTime hms = frac_day_to_hms(frac_day);

    CalendarTime cal{
        .year = year,
        .month = month,
        .day = doy,
        .hour = hms.hour,
        .minute = hms.minute,
        .second = hms.second
    };

    return cal;
}

inline f64 cal_to_doy(const CalendarTime& cal) {
    f64 doy = 0.0;
    if (cal.day < 1 || cal.day > days_in_month(cal.year, cal.month)) {
        return doy;
    }

    for (i32 month = 1; month < cal.month; ++month) {
        doy += static_cast<f64>(days_in_month(cal.year, month));
    }
    doy += cal.day + hms_to_frac_day(hms_from_cal(cal));

    return doy;
}

inline i32 get_jd_index(const JulianDate& jd, ecref<vecX<JulianDate>> jds) {
    i32 i = 0;
    i32 n = jds.size();

    f64 jd_scalar_curr = jd_to_scalar(jds(i));
    f64 jd_scalar_next = 0.0;
    if (i + 1 < n) {
        jd_scalar_next = jd_to_scalar(jds(i + 1));
    } else {
        return i;
    }

    f64 jd_scalar = jd_to_scalar(jd);
    while (jd_scalar >= jd_scalar_next) {
        jd_scalar_curr = jd_to_scalar(jds(i));
        if (i + 1 < n) {
            jd_scalar_next = jd_to_scalar(jds(i + 1));
        } else {
            return i;
        }
        ++i;
    }
    f64 delta_curr = std::abs(jd_scalar - jd_scalar_curr);
    f64 delta_next = std::abs(jd_scalar - jd_scalar_next);
    if (delta_curr < delta_next) {
        if (i == 0) {
            return 0;
        } else {
            return i - 1;
        }
    } else {
        return i;
    }

    return i;
}

struct DHMStime {
    i32 day = 0;
    i32 hour = 0;
    i32 minute = 0;
    f64 second = 0.0;
};
inline DHMStime sec_to_dhms(f64 s, f64 sec_in_day = 86400.0) {
    DHMStime out;

    if (!std::isfinite(s) || !std::isfinite(sec_in_day) || sec_in_day <= 0.0) {
        return out;
    }

    f64 day_f = std::floor(s / sec_in_day);
    out.day = static_cast<i32>(day_f);

    f64 rem = s - day_f * sec_in_day;

    out.hour = static_cast<i32>(std::floor(rem / 3600.0));
    rem -= static_cast<f64>(out.hour) * 3600.0;

    out.minute = static_cast<i32>(std::floor(rem / 60.0));
    rem -= static_cast<f64>(out.minute) * 60.0;

    out.second = rem;

    return out;
}

inline string cal_str(
    const CalendarTime& cal,
    i32 precision = 4,
    string date_separator = "-",
    string DT_separator = "T",
    string time_separator = ":",
    string filler = "0"
) {
    // "YYYY-MM-DDTHH:MM:SS.FFFF"
    i32 num_fill = 0;
    string zeros = "";

    i32 y_digits = count_digits(cal.year);
    num_fill = 4 - y_digits - (cal.year < 0 ? 1 : 0);
    string year
        = (y_digits < 4 ? repeat_char(filler, num_fill) : "") + std::to_string(cal.year);

    i32 m_digits = count_digits(cal.month);
    num_fill = 2 - m_digits - (cal.month < 0 ? 1 : 0);
    string month
        = (m_digits < 2 ? repeat_char(filler, num_fill) : "") + std::to_string(cal.month);

    i32 d_digits = count_digits(cal.day);
    num_fill = 2 - d_digits - (cal.day < 0 ? 1 : 0);
    string day
        = (d_digits < 2 ? repeat_char(filler, num_fill) : "") + std::to_string(cal.day);

    i32 h_digits = count_digits(cal.hour);
    num_fill = 2 - h_digits - (cal.day < 0 ? 1 : 0);
    string hour
        = (h_digits < 2 ? repeat_char(filler, num_fill) : "") + std::to_string(cal.hour);

    i32 min_digits = count_digits(cal.minute);
    num_fill = 2 - min_digits - (cal.minute < 0 ? 1 : 0);
    string minute = (min_digits < 2 ? repeat_char(filler, num_fill) : "")
                    + std::to_string(cal.minute);

    i32 s_int = static_cast<i32>(cal.second);
    f64 s_frac = cal.second - s_int;
    i32 s_digits = count_digits(s_int);
    num_fill = 2 - s_digits - (cal.second < 0.0 ? 1 : 0);
    string second = (s_digits < 2 ? repeat_char(filler, num_fill) : "")
                    + to_string_precision(cal.second, precision);

    string str = year + date_separator + month + date_separator + day + DT_separator
                 + hour + time_separator + minute + time_separator + second;
    return str;
}

enum struct CalendarPrintStyle { separate, separate_vertical, string, vector };
inline void print_cal(
    const CalendarTime& cal,
    CalendarPrintStyle style = CalendarPrintStyle::separate
) {
    switch (style) {
    case CalendarPrintStyle::separate: {
        std::println(
            "Year: {}, Month: {}, Day: {}, Hour: {}, Minute: {}, Second: {}",
            cal.year,
            cal.month,
            cal.day,
            cal.hour,
            cal.minute,
            cal.second
        );
    } break;
    case CalendarPrintStyle::separate_vertical: {
        std::println(
            "Year: {}\nMonth: {}\nDay: {}\nHour: {}\nMinute: {}\nSecond: {}",
            cal.year,
            cal.month,
            cal.day,
            cal.hour,
            cal.minute,
            cal.second
        );
    } break;
    case CalendarPrintStyle::vector: {
        vec6d cal_vec;
        cal_vec << cal.year, cal.month, cal.day, cal.hour, cal.minute, cal.second;
        std::println("{}", vec_string(cal_vec));
    } break;
    case CalendarPrintStyle::string: {
        std::println("{}", cal_str(cal));
    } break;
    }
}

inline CalendarTime normalize_cal(const CalendarTime& cal) {
    CalendarTime out = cal;

    if (out.second < 0.0 || out.second >= 60.0) {
        i32 dt_min = static_cast<i32>(std::floor(out.second / 60.0));

        out.second -= static_cast<f64>(dt_min) * 60.0;
        out.minute += dt_min;
    }
    if (out.second >= 60.0) {
        out.second -= 60.0;
        out.minute += 1;
    } else if (out.second < 0.0) {
        out.second += 60.0;
        out.minute -= 1;
    }

    if (out.minute < 0 || out.minute >= 60) {
        i32 dt_h = floor_div(out.minute, 60);

        out.minute -= dt_h * 60;
        out.hour += dt_h;
    }

    if (out.hour < 0 || out.hour >= 24) {
        i32 dt_d = floor_div(out.hour, 24);

        out.hour -= dt_d * 24;
        out.day += dt_d;
    }

    if (out.month < 1 || out.month > 12) {
        i32 month0 = out.month - 1;

        i32 dt_y = floor_div(month0, 12);
        out.month = floor_mod(month0, 12) + 1;
        out.year += dt_y;
    }

    while (out.day <= 0) {
        --out.month;

        if (out.month <= 0) {
            out.month = 12;
            --out.year;
        }

        out.day += days_in_month(out.year, out.month);
    }

    while (out.day > days_in_month(out.year, out.month)) {
        out.day -= days_in_month(out.year, out.month);

        ++out.month;

        if (out.month > 12) {
            out.month = 1;
            ++out.year;
        }
    }

    return out;
}

inline bool validate_hms(const HMSTime& hms) {
    i32 h = hms.hour;
    i32 m = hms.minute;
    f64 s = hms.second;

    if (h < 0 || h > 23) return false;
    if (m < 0 || m > 59) return false;
    if (s < 0.0) return false;

    return true;
}
inline bool validate_cal(const CalendarTime& cal) {
    i32 y = cal.year;
    i32 m = cal.month;
    i32 d = cal.day;

    if (m <= 0 || m > 12) return false;
    if (d <= 0 || d > days_in_month(y, m)) return false;

    return validate_hms(hms_from_cal(cal));
}

inline bool parse_cal_str(const string& str, CalendarTime& cal) {
    // "YYYY-MM-DDTHH:MM:SS.FFFF...[SCALE]", assume strictly is this format

    string y_str = str.substr(0, 4);
    if (!is_numeric(y_str, true)) return false;
    cal.year = std::stoi(y_str);

    return true;
}

inline void print_chrono(auto duration, UTime units) {
    switch (units) {
    case UTime::minute:
        std::println(
            "Elapsed Time: {} min",
            std::chrono::duration<f64>(duration).count() / 60.0
        );

        break;
    case UTime::second:
        std::println("Elapsed Time: {} s", std::chrono::duration<f64>(duration).count());
        break;
    case UTime::millisecond:
        std::println(
            "Elapsed Time: {} ms",
            std::chrono::duration<f64, std::milli>(duration).count()
        );
        break;
    case UTime::microsecond:
        std::println(
            "Elapsed Time: {} µs",
            std::chrono::duration<f64, std::micro>(duration).count()
        );
        break;
    case UTime::nanosecond:
        std::println(
            "Elapsed Time: {} ns",
            std::chrono::duration<f64, std::nano>(duration).count()
        );
        break;
    default: std::println("Wrong time units for duration"); break;
    }
}