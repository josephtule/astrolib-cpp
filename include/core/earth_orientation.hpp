#pragma once

#include "core/time.hpp"
#include "core/transform.hpp"
#include "util/constants.hpp"
#include "util/math.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

inline f64 gmst_from_jd(
    const JulianDate& jd,
    TimeScale scale_in = TimeScale::ut1,
    TimeOffsets offsets = TimeOffsets{},
    UAngle angle_out = UAngle::degree
) {
    JulianDate jd_ut1 = jd_scale_convert(jd, scale_in, TimeScale::ut1, offsets);
    f64 jd_ut1_scalar = jd_to_scalar(jd_ut1);

    f64 T = (jd_ut1_scalar - 2451545.0) / 36525; // Centuries from J2000
    f64 GMST = 67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * T + 0.093104 * T * T
               - 6.2e-6 * T * T * T;

    // Degrees
    f64 GMST_angle = std::remainder(GMST, 24.0 * 3600.0) / 240.0;
    if (GMST_angle < 0.0) {
        GMST_angle += 360.0;
    }

    if (angle_out != UAngle::degree) {
        GMST_angle = convert_angle(GMST_angle, UAngle::degree, angle_out);
    }

    return GMST_angle;
}

inline f64 lmst_from_jd(
    const JulianDate& jd,
    f64 longitude,
    TimeScale scale_in = TimeScale::ut1,
    TimeOffsets offsets = TimeOffsets{},
    UAngle angle_in = UAngle::degree,
    UAngle angle_out = UAngle::degree
) {
    if (angle_in != UAngle::degree) {
        longitude = convert_angle(longitude, angle_in, UAngle::degree);
    }

    f64 GMST
        = gmst_from_jd(jd, scale_in, offsets); // return degrees, convert in lmst scope
    f64 LMST = GMST + longitude;
    if (LMST < 0.0) {
        LMST += 360.0;
    }
    LMST = wrap_angle(LMST, 0.0, 360.0, UAngle::degree, UAngle::degree);

    if (angle_out != UAngle::degree) {
        LMST = convert_angle(LMST, UAngle::degree, angle_out);
    }

    return LMST;
}

/*
Precession: TT
Nutation: TT
Earth rotation (sidereal rotation / ERA): UT1
Polar motion: UT1
Celestial Intermediate Pole (CIP) coordinates: TT
Earth Orientation Parameters (EOP, combined): mixed (UT1 for rotation, TT for
celestial, UTC for dissemination)
*/

// inline TimeOffsets offsets_from_EOP2(const JulianDate& jd, EOP,){}
enum struct EarthNutationModel {
    // https://hpiers.obspm.fr/eop-pc/models/nutations/
    IAU1980,
    IAU1996
};
enum struct EarthPrecessionModel {};
enum struct EarthPolarMotionModel {};
struct EarthOrientationParams {
    JulianDate jd;
    TimeOffsets offsets;
    EarthNutationModel nutation_model;
    EarthPrecessionModel precession_model;
    EarthPolarMotionModel polar_motion_model;
};
inline std::unique_ptr<EarthOrientationParams> load_EOP() {
    auto eop = std::make_unique<EarthOrientationParams>();

    return eop;
}
// load_nutation_model
// load_precession_model
// load polar motion model

struct EarthNutationParams {
    // Common
    vecXd lm, ls, F, D, Om, period;
    // IAU1980
    vecXd Psisin, t_sin, epscos, t_cos;
    // IAU1996
    // https://iers-conventions.obspm.fr/content/tn36.pdf
    vecXd lme, lv, le, lma, lj, lS, lu, ln, pa, Psicos, epssin;
};
inline bool load_nutation_model(
    std::string filename,
    EarthNutationParams& params,
    EarthNutationModel model,
    i32 precision = 1,
    i32 lineskips = 0
) {
    std::ifstream file(filename);
    if (!file) return false;
    std::string line;

    for (i32 i = 0; i < lineskips && std::getline(file, line); ++i) {}

    if (model == EarthNutationModel::IAU1980) {
        params.lm.resize(precision);
        params.ls.resize(precision);
        params.F.resize(precision);
        params.D.resize(precision);
        params.Om.resize(precision);
        params.period.resize(precision);
        params.Psisin.resize(precision);
        params.t_sin.resize(precision);
        params.epscos.resize(precision);
        params.t_cos.resize(precision);

        i32 i = 0;
        while (getline(file, line)) {
            if (i >= precision) break;
            std::istringstream iss(line);

            bool failed
                = !(iss >> params.lm(i) >> params.ls(i) >> params.F(i) >> params.D(i)
                    >> params.Om(i) >> params.period(i) >> params.Psisin(i)
                    >> params.t_sin(i) >> params.epscos(i) >> params.t_cos(i));
            if (failed) {
                continue;
            }

            ++i;
        }
    } else if (model == EarthNutationModel::IAU1996) {

        params.lm.resize(precision);
        params.ls.resize(precision);
        params.F.resize(precision);
        params.D.resize(precision);
        params.Om.resize(precision);
        params.lme.resize(precision);
        params.lv.resize(precision);
        params.le.resize(precision);
        params.lma.resize(precision);
        params.lj.resize(precision);
        params.lS.resize(precision);
        params.lu.resize(precision);
        params.ln.resize(precision);
        params.pa.resize(precision);
        params.period.resize(precision);
        params.Psisin.resize(precision);
        params.Psicos.resize(precision);
        params.epssin.resize(precision);
        params.epscos.resize(precision);

        std::istringstream iss(line);
        i32 i = 0;
        while (getline(file, line)) {
            if (i >= precision) break;
            std::istringstream iss(line);

            f64 temp;
            bool failed
                = !(iss >> temp >> params.lm(i) >> params.ls(i) >> params.F(i)
                    >> params.D(i) >> params.Om(i) >> params.lme(i) >> params.lv(i)
                    >> params.le(i) >> params.lma(i) >> params.lj(i) >> params.lS(i)
                    >> params.lu(i) >> params.ln(i) >> params.pa(i) >> params.period(i)
                    >> params.Psisin(i) >> params.Psicos(i) >> params.epssin(i)
                    >> params.epscos(i));

            if (failed) {
                continue;
            }

            ++i;
        }
    }
    return true;
}

inline mat3d rot_earth_precession(
    const JulianDate& jd,
    TimeScale scale_in = TimeScale::tt,
    TimeOffsets offsets = TimeOffsets{}
) {
    JulianDate jd_TT = jd;
    if (scale_in != TimeScale::tt) {
        JulianDate jd_TT = jd_scale_convert(jd, scale_in, TimeScale::tt);
    }
    f64 jd_TT_scalar = jd_to_scalar(jd_TT);

    f64 T = (jd_TT_scalar - 2451545.0) / 36525; // Centuries from J2000

    // Precession angles [arcseconds]
    f64 T2 = T * T;
    f64 T3 = T * T * T;
    f64 zeta = 2306.2181 * T + 0.30188 * T2 + 0.017998 * T3;
    f64 theta = 2004.3109 * T - 0.42665 * T2 - 0.041833 * T3;
    f64 z = zeta + 0.79280 * T2 + 0.000205 * T3;

    zeta /= 3600.0;
    theta /= 3600.0;
    z /= 3600.0;

    // Rotation matrix
    return rot(-z, RotAxis::z, UAngle::degree) * rot(theta, RotAxis::y, UAngle::degree)
           * rot(-zeta, RotAxis::z, UAngle::degree);
}

inline mat3d rot_earth_polar_motion(
    const JulianDate& jd,
    // TODO: add eop parameters here
    TimeScale scale_in = TimeScale::utc,
    TimeOffsets = TimeOffsets{},
    bool approx = true
) {
    // polar motion technically uses ut1 but eop data uses utc
    JulianDate jd_tai = jd;
    if (scale_in != TimeScale::utc) {
        jd_tai = jd_scale_convert(jd, scale_in, TimeScale::utc);
    }
    return mat3d1;
}

inline mat3d rot_earth_nutation(
    const JulianDate& jd,
    // TODO: add eop parameters here or in scope
    EarthNutationParams params,
    EarthNutationModel model = EarthNutationModel::IAU1980,
    TimeScale scale_in = TimeScale::tt,
    TimeOffsets offsets = TimeOffsets{},
    bool approx = true,
    i32 precision = 106
) {
    JulianDate jd_TT = jd;
    if (scale_in != TimeScale::tt) {
        JulianDate jd_TT = jd_scale_convert(jd, scale_in, TimeScale::tt);
    }
    f64 jd_TT_scalar = jd_to_scalar(jd_TT);

    f64 T = (jd_TT_scalar - 2451545.0) / 36525; // Centuries from J2000
    f64 T2 = T * T;
    f64 T3 = T * T * T;

    f64 epsilon = 23.43929111 * 3600 - 46.8150 * T - 0.00059 * T2 + 0.001813 * T3;
    f64 delta_Psi = 0.0;
    f64 delta_epsilon = 0.0;

    mat3d R = mat3d1;
    if (model == EarthNutationModel::IAU1980) {
        if (approx) {
            f64 l = 357.525 + 35999 * T;        // mean anomaly of the sun [deg]
            f64 F = 93.273 + 483202.019 * T;    // mean distance b/w nodes of moon [deg]
            f64 D = 297.850 + 445267.111 * T;   // mean distance b/w sun and moon [deg]
            f64 Omega = 125.045 - 1934.136 * T; // mean longitude of moon [deg]

            // Convert to radians
            l *= deg_to_rad;
            F *= deg_to_rad;
            D *= deg_to_rad;
            Omega *= deg_to_rad;

            delta_Psi = -17.200 * std::sin(Omega) + 0.202 * std::sin(2.0 * Omega)
                        - 1.319 * std::sin(2.0 * (F - D + Omega)) + 0.143 * std::sin(l)
                        - 0.227 * std::sin(2.0 * (F + Omega));
            delta_epsilon = 9.203 * std::cos(Omega) - 0.090 * std::cos(2.0 * Omega)
                            - 0.547 * std::cos(2.0 * (F - D + Omega))
                            + 0.098 * std::cos(2.0 * (F + Omega));

        } else {
            f64 l = DMS_to_deg(134, 57, 46.733) + DMS_to_deg(477198, 52, 2.633) * T
                    + DMS_to_deg(0, 0, 31.310) * T2 + DMS_to_deg(0, 0, 0.064) * T3;
            f64 l_prime = DMS_to_deg(357, 31, 39.804) + DMS_to_deg(35999, 3, 1.244) * T
                          - DMS_to_deg(0, 0, 0.577) * T2 - DMS_to_deg(0, 0, 0.012) * T3;
            f64 F = DMS_to_deg(93, 16, 18.877) + DMS_to_deg(483202, 1, 3.137) * T
                    - DMS_to_deg(0, 0, 13.257) * T2 + DMS_to_deg(0, 0, 0.011) * T3;
            f64 D = DMS_to_deg(297, 51, 1.307) + DMS_to_deg(445267, 6, 41.328) * T
                    - DMS_to_deg(0, 0, 6.891) * T2 + DMS_to_deg(0, 0, 0.019) * T3;
            f64 Omega = DMS_to_deg(125, 2, 40.280) - DMS_to_deg(1934, 8, 10.539) * T
                        + DMS_to_deg(0, 0, 7.455) * T2 + DMS_to_deg(0, 0, 0.008) * T3;

            for (i32 i = 0; i < precision; ++i) {
                f64 phi = params.lm(i) * l + params.ls(i) * l_prime + params.F(i) * F
                          + params.D(i) * D + params.Om(i) * Omega;
                delta_Psi
                    += (params.Psisin(i) + params.t_sin(i) * T) + sind(phi) / 10000.0;
                delta_epsilon
                    += (params.epscos(i) + params.t_cos(i) * T) * cosd(phi) / 10000.0;
            }
        }

        epsilon /= 3600.0;
        delta_Psi /= 3600.0;
        delta_epsilon /= 3600.0;

        R = rot(-epsilon - delta_epsilon, RotAxis::x, UAngle::degree)
            * rot(-delta_Psi, RotAxis::z, UAngle::degree)
            * rot(epsilon, RotAxis::x, UAngle::degree);

    } else if (model == EarthNutationModel::IAU1996) {
        R = mat3d1;
    }

    return R;
}

inline f64 earth_rot_angle_from_jd(
    const JulianDate& jd,
    TimeScale scale_in = TimeScale::ut1,
    TimeOffsets offsets = TimeOffsets{},
    UAngle angle_out = UAngle::degree
) {
    // Just a wrapper for GMST
    return gmst_from_jd(jd, scale_in, offsets, angle_out);
}

enum struct EarthFrame {
    ICRF = 1,
    GCRF = 1,
    J2000 = 1,
    EME2000 = 1,
    MOD = 2,
    TOD = 3,
    GTOD = 4,
    PEF = 5,
    TEME = 6,
    ITRF = 7
};

inline mat3d rot_eci_to_ecef(
    const JulianDate& jd,
    EarthFrame frame_in,
    EarthFrame frame_out,
    TimeScale scale_in = TimeScale::ut1,
    TimeOffsets offsets = TimeOffsets{},
    bool approx = true
) {
    mat3d R = mat3d1;
    f64 theta
        = earth_rot_angle_from_jd(jd, scale_in, offsets, UAngle::degree); // in degrees

    // Precession
    R = rot_earth_precession(jd, scale_in, offsets);

    // Nutation

    // Sidereal Rotation
    R = rot(theta, RotAxis::z, UAngle::degree) * R;

    // Polar Motion

    return R;
}

inline mat3d rot_ecef_to_eci(
    const JulianDate& jd,
    EarthFrame frame_in,
    EarthFrame frame_out,
    TimeScale scale_in = TimeScale::ut1,
    TimeOffsets offsets = TimeOffsets{},
    bool approx = true
) {
    return rot_eci_to_ecef(jd, frame_in, frame_out, scale_in, offsets, approx)
        .transpose();
}
