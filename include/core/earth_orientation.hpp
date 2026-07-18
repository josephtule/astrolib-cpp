// Copyright 2025-2026 Joseph Tu Le
// SPDX-License-Identifier: Apache-2.0

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

// TODO: Store the current transform(s) somewhere

f64 gmst_from_jd(
    const JulianDate& jd,
    TimeOffsets offsets,
    TimeScale scale_in = TimeScale::ut1,
    UAngle angle_out = UAngle::degree
);

f64 lmst_from_jd(
    const JulianDate& jd,
    f64 longitude,
    TimeOffsets offsets,
    TimeScale scale_in = TimeScale::ut1,
    UAngle angle_in = UAngle::degree,
    UAngle angle_out = UAngle::degree
);

f64 earth_rot_angle_from_jd(
    const JulianDate& jd,
    TimeOffsets offsets,
    TimeScale scale_in = TimeScale::ut1,
    UAngle angle_out = UAngle::degree
);

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
enum struct EarthPolarMotionModel { JPLEOP2, IAU1980, IAU2000, IAU2000A };

struct LeapSecondParams {
    std::string filename;
    i32 lineskips = -1;
    vecX<JulianDate> jd;
    vecXd leap_seconds;
    bool loaded = false;
};
struct EarthPolarMotionParams {
    // Common
    std::string filename;
    i32 lineskips = -1;
    EarthPolarMotionModel model;
    bool loaded = false;

    vecX<JulianDate> jd;
    vecXd xp, yp;
    vecXd ut1_utc;
    bool approx = false;
};
struct EarthNutationParams {
    // Common
    std::string filename;
    i32 lineskips = -1;
    EarthNutationModel model = EarthNutationModel::IAU1980;
    bool approx = false;
    i32 precision = 106;
    bool loaded = false;
    vecXd lm, ls, F, D, Om, period;
    // IAU1980
    vecXd Psisin, t_sin, epscos, t_cos;
    // IAU1996
    // https://iers-conventions.obspm.fr/content/tn36.pdf
    vecXd lme, lv, le, lma, lj, lS, lu, ln, pa, Psicos, epssin;
};
struct EarthOrientationParams {
    JulianDate jd;
    TimeOffsets offsets{};
    LeapSecondParams leap_seconds{};
    bool loaded = false;

    // Nutation
    EarthNutationParams nutation{};

    // Polar Motion
    EarthPolarMotionParams polar_motion;

    // Precision
};

bool get_time_offsets(
    const JulianDate& jd,
    TimeOffsets& offsets,
    const LeapSecondParams& lsp,
    const EarthPolarMotionParams& pmp
);
bool get_time_offsets(const JulianDate& jd, EarthOrientationParams& eop);

bool load_nutation_model(
    std::string filename,
    EarthNutationParams& params,
    EarthNutationModel model,
    i32 precision,
    i32 lineskips,
    bool approx = false
);
bool load_nutation_model(EarthNutationParams& params);

bool load_leap_seconds(
    std::string filename,
    vecX<JulianDate>& jd,
    vecXd& leap_seconds,
    i32 lineskips
);
bool load_leap_seconds(
    std::string filename,
    LeapSecondParams& params,
    i32 lineskips
);

bool load_polar_motion_model_all(
    std::string filename,
    EarthPolarMotionParams& params,
    EarthPolarMotionModel model,
    const LeapSecondParams& lsp,
    i32 lineskips
);
bool load_polar_motion_model_all(
    EarthPolarMotionParams& params,
    const LeapSecondParams& lsp
);

i32 get_polar_motion_index(
    const JulianDate& jd,
    const EarthPolarMotionParams& params
);

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
    i32 nutation_lineskips
);

bool load_all_eop(EarthOrientationParams& params);

mat3d rot_earth_precession(
    const JulianDate& jd,
    TimeOffsets offsets,
    TimeScale scale_in = TimeScale::tt
);

mat3d rot_earth_polar_motion(
    const JulianDate& jd,
    const EarthPolarMotionParams& params,
    TimeOffsets offsets,
    TimeScale scale_in = TimeScale::utc
);

mat3d rot_earth_nutation(
    const JulianDate& jd,
    // TODO: add eop parameters here or in scope
    const EarthNutationParams& params,
    TimeOffsets offsets,
    TimeScale scale_in = TimeScale::tt
);

enum struct EarthFrame : i32 {
    // Ordered by frame hierarchy:
    // GCRF/ICRF/J2000/EME2000 -> MOD -> TOD -> GTOD/PEF -> ITRS.

    // ECI frames
    ICRF = 1,
    GCRF = 1,
    J2000 = 1,
    EME2000 = 1, // ICRF/GCRF \approx J2000/EME2000
    MOD = 2,
    TOD = 3,

    // ECEF frames
    GTOD = 4,
    PEF = 4,
    // TEME = 5,
    ITRS = 6,
};

// NOTE:
// These rotations are passive rotations
// Unless noted with a "to", all rotations are R_A_B : B -> A
EarthFrame frame_resolver(EarthFrame frame);

i32 earth_frame_order(EarthFrame frame);

bool is_same_earth_frame(EarthFrame a, EarthFrame b);

bool is_forward_earth_frame_transform(EarthFrame source, EarthFrame target);

bool is_reverse_earth_frame_transform(EarthFrame source, EarthFrame target);

namespace earth_rot_helper {
inline mat3d rot_earth(
    const JulianDate& jd,
    EarthFrame frame_source,
    EarthFrame frame_target,
    const EarthOrientationParams& params,
    const TimeOffsets& offsets,
    TimeScale scale_in = TimeScale::ut1
) {
    mat3d R = mat3d1;

    if (is_reverse_earth_frame_transform(frame_source, frame_target)) return R;

    EarthFrame frame = frame_source;

    // ECI ---------------------------------------------------------------------
    if (frame == EarthFrame::ICRF || frame == EarthFrame::GCRF
        || frame == EarthFrame::J2000 || frame == EarthFrame::EME2000) {
        // ICRF/GCRF to MOD
        // Precession
        R = rot_earth_precession(jd, offsets, scale_in) * R;
        frame = EarthFrame::MOD;
        if (frame == frame_target) return R;
    }

    if (frame == EarthFrame::MOD) {
        // MOD to TOD
        // Nutation
        R = rot_earth_nutation(jd, params.nutation, offsets, scale_in) * R;
        frame = EarthFrame::TOD;
        if (frame == frame_target) return R;
    }

    // ECEF --------------------------------------------------------------------
    if (frame == EarthFrame::TOD) {
        // TOD to GTOD/PEF
        // TODO: add TOD to TEME and TEME to GTOD/PEF
        // Sidereal Rotation
        f64 theta = earth_rot_angle_from_jd(
            jd,
            offsets,
            scale_in,
            UAngle::degree
        ); // in degrees
        R = rot(theta, RotAxis::z, UAngle::degree) * R;
        frame = EarthFrame::GTOD;
        if (frame == frame_target) return R;
    }

    if (frame == EarthFrame::GTOD || frame == EarthFrame::PEF) {
        // GTOD/PEF to ITRF
        // Polar Motion
        R = rot_earth_polar_motion(jd, params.polar_motion, offsets, scale_in) * R;
        frame = EarthFrame::ITRS;
        if (frame == frame_target) return R;
    }

    return R;
}

} // namespace earth_rot_helper


mat3d rot_earth_frame(
    const JulianDate& jd,
    EarthFrame frame_source,
    EarthFrame frame_target,
    const EarthOrientationParams& params,
    TimeOffsets offsets,
    TimeScale scale_in = TimeScale::ut1
);
