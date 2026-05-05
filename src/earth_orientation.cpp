#include "core/earth_orientation.hpp"

f64 gmst_from_jd(
    const JulianDate& jd,
    TimeOffsets offsets,
    TimeScale scale_in,
    UAngle angle_out) {
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

f64 lmst_from_jd(
    const JulianDate& jd,
    f64 longitude,
    TimeOffsets offsets,
    TimeScale scale_in,
    UAngle angle_in,
    UAngle angle_out) {
    if (angle_in != UAngle::degree) {
        longitude = convert_angle(longitude, angle_in, UAngle::degree);
    }

    f64 GMST
        = gmst_from_jd(jd, offsets, scale_in); // return degrees, convert in lmst scope
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

f64 earth_rot_angle_from_jd(
    const JulianDate& jd,
    TimeOffsets offsets,
    TimeScale scale_in,
    UAngle angle_out) {
    // Just a wrapper for GMST
    return gmst_from_jd(jd, offsets, scale_in, angle_out);
}

bool get_time_offsets(
    const JulianDate& jd,
    TimeOffsets& offsets,
    const LeapSecondParams& lsp,
    const EarthPolarMotionParams& pmp) {
    if (!lsp.loaded || !pmp.loaded) {
        return false;
    }

    i32 idx_ls = get_jd_index(jd, lsp.jd);
    offsets.tai_utc = lsp.leap_seconds(idx_ls);

    i32 idx_ut1_utc = get_jd_index(jd, pmp.jd);
    offsets.ut1_utc = pmp.ut1_utc(idx_ut1_utc);

    return true;
}

bool get_time_offsets(const JulianDate& jd, EarthOrientationParams& eop) {
    return get_time_offsets(jd, eop.offsets, eop.leap_seconds, eop.polar_motion);
}

bool load_nutation_model(
    std::string filename,
    EarthNutationParams& params,
    EarthNutationModel model,
    i32 precision,
    i32 lineskips,
    bool approx) {
    params.filename = filename;
    params.model = model;
    params.lineskips = lineskips;
    params.precision = precision;
    params.approx = approx;

    std::ifstream file(params.filename);
    if (!file) return false;
    std::string line;

    for (i32 i = 0; i < lineskips && std::getline(file, line); ++i) {}

    switch (model) {
    case EarthNutationModel::IAU1980: {
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
        params.lm.conservativeResize(i);
        params.ls.conservativeResize(i);
        params.F.conservativeResize(i);
        params.D.conservativeResize(i);
        params.Om.conservativeResize(i);
        params.period.conservativeResize(i);
        params.Psisin.conservativeResize(i);
        params.t_sin.conservativeResize(i);
        params.epscos.conservativeResize(i);
        params.t_cos.conservativeResize(i);
        params.precision = i;
    } break;
    case EarthNutationModel::IAU1996: {
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
        params.lm.conservativeResize(i);
        params.ls.conservativeResize(i);
        params.F.conservativeResize(i);
        params.D.conservativeResize(i);
        params.Om.conservativeResize(i);
        params.lme.conservativeResize(i);
        params.lv.conservativeResize(i);
        params.le.conservativeResize(i);
        params.lma.conservativeResize(i);
        params.lj.conservativeResize(i);
        params.lS.conservativeResize(i);
        params.lu.conservativeResize(i);
        params.ln.conservativeResize(i);
        params.pa.conservativeResize(i);
        params.period.conservativeResize(i);
        params.Psisin.conservativeResize(i);
        params.Psicos.conservativeResize(i);
        params.epssin.conservativeResize(i);
        params.epscos.conservativeResize(i);
        params.precision = i;
    } break;
    }
    params.loaded = true;
    return true;
}

bool load_nutation_model(EarthNutationParams& params) {
    return load_nutation_model(
        params.filename,
        params,
        params.model,
        params.precision,
        params.lineskips,
        params.approx
    );
}

bool load_leap_seconds(
    std::string filename,
    vecX<JulianDate>& jd,
    vecXd& leap_seconds,
    i32 lineskips) {
    i32 numlines = 50; // arbitrary
    jd.resize(numlines);
    leap_seconds.resize(numlines);

    // https://data.iana.org/time-zones/data/leap-seconds.list
    std::ifstream file(filename);
    if (!file) return false;
    std::string line;

    for (i32 i = 0; i < lineskips && std::getline(file, line); ++i) {}

    i32 i = 0;
    while (getline(file, line)) {
        if (i >= numlines) {
            numlines = i * 2;
            jd.conservativeResize(i * 2);
            leap_seconds.conservativeResize(i * 2);
        }

        std::istringstream iss(line);
        // NTP time in UTC
        f64 ntp_time;
        f64 tai_utc;

        bool failed = !(iss >> ntp_time >> tai_utc);
        if (failed) {
            continue;
        }

        jd(i) = mjd_to_jd(mjd_from_scalar(ntp_time / 86400 + 15020.0));
        leap_seconds(i) = tai_utc;

        ++i;
    }
    jd.conservativeResize(i);
    leap_seconds.conservativeResize(i);

    return true;
}

bool load_leap_seconds(
    std::string filename,
    LeapSecondParams& params,
    i32 lineskips) {
    params.loaded
        = load_leap_seconds(filename, params.jd, params.leap_seconds, lineskips);
    return params.loaded;
}

bool load_polar_motion_model_all(
    std::string filename,
    EarthPolarMotionParams& params,
    EarthPolarMotionModel model,
    const LeapSecondParams& lsp,
    i32 lineskips) {
    params.filename = filename;
    params.lineskips = lineskips;
    params.model = model;

    std::ifstream file(filename);
    if (!file) return false;
    std::string line;
    i32 numlines = 5000; // arbitrary

    params.jd.resize(numlines);
    params.xp.resize(numlines);
    params.yp.resize(numlines);
    params.ut1_utc.resize(numlines);

    for (i32 i = 0; i < lineskips && std::getline(file, line); ++i) {}

    switch (model) {
    case EarthPolarMotionModel::JPLEOP2: {
        i32 i = 0;
        while (getline(file, line)) {
            if (i >= numlines) {
                numlines = i * 2;
                params.jd.conservativeResize(i * 2);
                params.xp.conservativeResize(i * 2);
                params.yp.conservativeResize(i * 2);
                params.ut1_utc.conservativeResize(i * 2);
            }
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream iss(line);

            f64 mjd_temp;
            f64 tai_ut1;
            f64 xp, yp;

            bool failed = !(iss >> mjd_temp >> xp >> yp >> tai_ut1);
            if (failed) {
                continue;
            }

            // JPL EOP2 gives tai-ut1, need ut1-utc and tai-utc (leap seconds)
            params.xp(i) = xp;
            params.yp(i) = yp;
            params.jd(i) = mjd_to_jd(mjd_from_scalar(mjd_temp));
            i32 idx = get_jd_index(params.jd(i), lsp.jd);
            f64 leap_seconds = lsp.leap_seconds(idx);
            params.ut1_utc(i) = leap_seconds - tai_ut1; // tai-utc - (tai-ut1)

            ++i;
        }
        params.jd.conservativeResize(i);
        params.xp.conservativeResize(i);
        params.yp.conservativeResize(i);
        params.ut1_utc.conservativeResize(i);
    } break;
    case EarthPolarMotionModel::IAU1980: // same as below
    case EarthPolarMotionModel::IAU2000: {
        i32 i = 0;
        while (getline(file, line)) {
            if (i >= numlines) {
                numlines = i * 2;
                params.jd.conservativeResize(i * 2);
                params.xp.conservativeResize(i * 2);
                params.yp.conservativeResize(i * 2);
                params.ut1_utc.conservativeResize(i * 2);
            }
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream iss(line);
            f64 mjd_temp;
            f64 ut1_tai;
            f64 xp, yp;

            bool failed = !(iss >> mjd_temp >> xp >> yp >> ut1_tai);
            if (failed) {
                continue;
            }

            // JPL EOP2 gives tai-ut1, need ut1-utc and tai-utc (leap seconds)
            // TODO: complete this
            params.xp(i) = xp;
            params.yp(i) = yp;
            params.jd(i) = mjd_to_jd(mjd_from_scalar(mjd_temp));
            i32 idx = get_jd_index(params.jd(i), lsp.jd);
            f64 leap_seconds = lsp.leap_seconds(idx);
            params.ut1_utc(i) = leap_seconds - ut1_tai; // tai-utc + (ut1-tai)

            ++i;
        }
        params.jd.conservativeResize(i);
        params.xp.conservativeResize(i);
        params.yp.conservativeResize(i);
        params.ut1_utc.conservativeResize(i);
    } break;
    case EarthPolarMotionModel::IAU2000A: {
        i32 i = 0;
        while (getline(file, line)) {
            if (i >= numlines) {
                numlines = i * 2;
                params.jd.conservativeResize(i * 2);
                params.xp.conservativeResize(i * 2);
                params.yp.conservativeResize(i * 2);
                params.ut1_utc.conservativeResize(i * 2);
            }
            std::istringstream iss(line);
            f64 trash;
            f64 mjd_temp;
            f64 ut1_utc;
            f64 xp, yp;

            bool failed
                = !(iss >> trash >> trash >> trash >> trash >> mjd_temp >> xp >> yp
                    >> ut1_utc);
            if (failed) {
                continue;
            }

            params.xp(i) = xp / 1000.0; // convert to milliarcseconds
            params.yp(i) = yp / 1000.0;
            // JPL EOP2 gives tai-ut1, need ut1-utc and tai-utc (leap seconds)
            // TODO: complete this
            params.jd(i) = mjd_to_jd(mjd_from_scalar(mjd_temp));
            params.ut1_utc(i) = ut1_utc;

            ++i;
        }
        params.jd.conservativeResize(i);
        params.xp.conservativeResize(i);
        params.yp.conservativeResize(i);
        params.ut1_utc.conservativeResize(i);
    } break;
    }

    params.loaded = true;
    return true;
}

bool load_polar_motion_model_all(
    EarthPolarMotionParams& params,
    const LeapSecondParams& lsp) {
    return load_polar_motion_model_all(
        params.filename,
        params,
        params.model,
        lsp,
        params.lineskips
    );
}

i32 get_polar_motion_index(
    const JulianDate& jd,
    const EarthPolarMotionParams& params) {
    return get_jd_index(jd, params.jd);
}

bool load_all_eop(
    EarthOrientationParams& params,
    std::string leap_second_filename,
    i32 leap_second_lineskips,
    std::string polar_motion_filename,
    EarthPolarMotionModel polar_motion_model,
    i32 polar_motion_lineskips,
    std::string nutation_filename,
    EarthNutationModel nutation_model,
    i32 nutation_precision,
    bool nutation_approx,
    i32 nutation_lineskips) {
    LeapSecondParams leap_second_params;
    bool leap_seconds_ok = load_leap_seconds(
        leap_second_filename,
        leap_second_params,
        leap_second_lineskips
    );
    params.leap_seconds = leap_second_params;
    if (!leap_seconds_ok) return false;

    bool polar_motion_ok = false;
    if (leap_seconds_ok) {
        EarthPolarMotionParams polar_motion_params;
        polar_motion_ok = load_polar_motion_model_all(
            polar_motion_filename,
            polar_motion_params,
            polar_motion_model,
            leap_second_params,
            polar_motion_lineskips
        );
        params.polar_motion = polar_motion_params;
    }
    if (!polar_motion_ok) return false;

    EarthNutationParams nutation_params;
    bool nutation_ok = load_nutation_model(
        nutation_filename,
        nutation_params,
        nutation_model,
        nutation_precision,
        nutation_lineskips,
        nutation_approx
    );
    params.nutation = nutation_params;
    if (!nutation_ok) return false;

    params.loaded = leap_seconds_ok && polar_motion_ok && nutation_ok;
    return params.loaded;
}

bool load_all_eop(EarthOrientationParams& params) {
    bool leap_seconds_ok = load_leap_seconds(
        params.leap_seconds.filename,
        params.leap_seconds,
        params.leap_seconds.lineskips
    );
    if (!leap_seconds_ok) return false;

    bool polar_motion_ok = false;
    if (leap_seconds_ok) {
        polar_motion_ok = load_polar_motion_model_all(
            params.polar_motion.filename,
            params.polar_motion,
            params.polar_motion.model,
            params.leap_seconds,
            params.polar_motion.lineskips
        );
    }
    if (!polar_motion_ok) return false;

    bool nutation_ok = load_nutation_model(
        params.nutation.filename,
        params.nutation,
        params.nutation.model,
        params.nutation.precision,
        params.nutation.lineskips,
        params.nutation.approx
    );
    if (!nutation_ok) return false;

    params.loaded = leap_seconds_ok && polar_motion_ok && nutation_ok;
    return params.loaded;
}

mat3d rot_earth_precession(
    const JulianDate& jd,
    TimeOffsets offsets,
    TimeScale scale_in) {
    JulianDate jd_TT = jd;
    if (scale_in != TimeScale::tt) {
        jd_TT = jd_scale_convert(jd, scale_in, TimeScale::tt, offsets);
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

mat3d rot_earth_polar_motion(
    const JulianDate& jd,
    const EarthPolarMotionParams& params,
    TimeOffsets offsets,
    TimeScale scale_in) {
    // polar motion technically uses ut1 but eop data uses utc
    JulianDate jd_utc = jd;
    if (scale_in != TimeScale::utc) {
        jd_utc = jd_scale_convert(jd, scale_in, TimeScale::utc, offsets);
    }

    mat3d R = mat3d1;
    switch (params.model) {
    case EarthPolarMotionModel::JPLEOP2: {
        // NOTE: not sure if all four models use the same computation with different
        // values, if yes then get rid of switch
        i32 idx = get_polar_motion_index(jd_utc, params);
        f64 xp = params.xp(idx) / 3600.0 / 1000.0; // degrees
        f64 yp = params.yp(idx) / 3600.0 / 1000.0; // degrees

        if (params.approx) {
            xp *= deg_to_rad;
            yp *= deg_to_rad;
            R << 1.0, 0.0, xp, 0.0, 1.0, -yp, -xp, yp, 1.0;
        } else {
            R = rot(-xp, RotAxis::y, UAngle::degree)
                * rot(-yp, RotAxis::x, UAngle::degree);
        }

    } break;
    case EarthPolarMotionModel::IAU1980: {
    } break;
    case EarthPolarMotionModel::IAU2000: {
    } break;
    case EarthPolarMotionModel::IAU2000A: {
    } break;
    }

    return R;
}

mat3d rot_earth_nutation(
    const JulianDate& jd,
    // TODO: add eop parameters here or in scope
    const EarthNutationParams& params,
    TimeOffsets offsets,
    TimeScale scale_in) {
    JulianDate jd_TT = jd;
    if (scale_in != TimeScale::tt) {
        jd_TT = jd_scale_convert(jd, scale_in, TimeScale::tt, offsets);
    }
    f64 jd_TT_scalar = jd_to_scalar(jd_TT);

    f64 T = (jd_TT_scalar - 2451545.0) / 36525; // Centuries from J2000
    f64 T2 = T * T;
    f64 T3 = T * T * T;

    f64 epsilon = 23.43929111 * 3600 - 46.8150 * T - 0.00059 * T2 + 0.001813 * T3;
    f64 delta_Psi = 0.0;
    f64 delta_epsilon = 0.0;

    mat3d R = mat3d1;
    if (params.model == EarthNutationModel::IAU1980) {
        if (params.approx) {
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

            for (i32 i = 0; i < params.precision; ++i) {
                f64 phi = params.lm(i) * l + params.ls(i) * l_prime + params.F(i) * F
                          + params.D(i) * D + params.Om(i) * Omega;
                delta_Psi
                    += (params.Psisin(i) + params.t_sin(i) * T) * sind(phi) / 10000.0;
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

    } else if (params.model == EarthNutationModel::IAU1996) {
        // TODO: research 1996 model formula
        R = mat3d1;
    }

    return R;
}

EarthFrame frame_resolver(EarthFrame frame) {
    if (i32(frame) == 1) return EarthFrame::ICRF;
    if (i32(frame) == 4) return EarthFrame::GTOD;
    return frame;
}

i32 earth_frame_order(EarthFrame frame) {
    frame = frame_resolver(frame);

    switch (frame) {
    case EarthFrame::ICRF: return 1;
    case EarthFrame::MOD: return 2;
    case EarthFrame::TOD: return 3;
    case EarthFrame::GTOD: return 4;
    case EarthFrame::ITRS: return 6;
    }

    return -1;
}

bool is_same_earth_frame(EarthFrame a, EarthFrame b) {
    return earth_frame_order(a) == earth_frame_order(b);
}

bool is_forward_earth_frame_transform(EarthFrame source, EarthFrame target) {
    return earth_frame_order(source) < earth_frame_order(target);
}

bool is_reverse_earth_frame_transform(EarthFrame source, EarthFrame target) {
    return earth_frame_order(source) > earth_frame_order(target);
}

mat3d rot_earth_frame(
    const JulianDate& jd,
    EarthFrame frame_source,
    EarthFrame frame_target,
    const EarthOrientationParams& params,
    TimeOffsets offsets,
    TimeScale scale_in) {
    frame_source = frame_resolver(frame_source);
    frame_target = frame_resolver(frame_target);

    if (is_forward_earth_frame_transform(frame_source, frame_target)) {
        return earth_rot_helper::rot_earth(
            jd,
            frame_source,
            frame_target,
            params,
            offsets,
            scale_in
        );
    } else if (is_reverse_earth_frame_transform(frame_source, frame_target)) {
        return earth_rot_helper::rot_earth(
                   jd,
                   frame_target,
                   frame_source,
                   params,
                   offsets,
                   scale_in
        )
            .transpose();
    } else {
        return mat3d1;
    }
}

