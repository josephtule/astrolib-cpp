#pragma once

#include "core/earth_orientation.hpp"
#include "core/state.hpp"
#include "core/time.hpp"
#include "core/transform.hpp"
#include "util/constants.hpp"
#include "util/math.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

// TODO: some functions will have earth based coordinates and others generic planet/body
// based, need to change to planet/body later on, but these should be general

inline mat3d rot_enu_from_bcbf(
    f64 latitude, // planetodetic latitude
    f64 longitude,
    UAngle angle_in = UAngle::degree
) {
    // Passive rotation from BCBF to ENU of a station
    if (angle_in != UAngle::radian) {
        latitude = convert_angle(latitude, angle_in, UAngle::radian);
        longitude = convert_angle(longitude, angle_in, UAngle::radian);
    }

    mat3d R_ENU_BCBF = mat3d1;

    R_ENU_BCBF.row(0) = vec3d{-sin(longitude), cos(longitude), 0.0};
    R_ENU_BCBF.row(1) = vec3d{
        -sin(latitude) * cos(longitude),
        -sin(latitude) * sin(longitude),
        cos(latitude)
    };
    R_ENU_BCBF.row(2) = vec3d{
        cos(latitude) * cos(longitude),
        cos(latitude) * sin(longitude),
        sin(latitude)
    };

    return R_ENU_BCBF;
}

inline vec3d bcbf_rel_to_enu(
    const vec3d& r_rel_bcbf,
    f64 latitude, // planetodetic
    f64 longitude,
    UAngle angle_in = UAngle::degree
) {
    // Takes a relative vector (in BCBF frame) and gives local ENU coordinates depending
    // on planetodetic latitude and longitude
    vec3d rho_enu = rot_enu_from_bcbf(latitude, longitude, angle_in) * r_rel_bcbf;

    return rho_enu;
}

// Azimuth is measured clockwise from the north
inline vec3d azel_from_enu(
    const vec3d& r_rel_enu, // must be ENU
    UAngle angle_out = UAngle::degree,
    f64 tol = tol12
) {
    vec3d azel = vec3d0;

    f64 range = r_rel_enu.norm(); // range
    if (range <= tol) {
        // Azimuth/elevation undefined for zero range.
        return vec3d0;
    }
    f64 elevation = std::asin(r_rel_enu(2) / range); // elevation
    f64 azimuth = wrap_angle(
        std::atan2(r_rel_enu(0), r_rel_enu(1)),
        0.0,
        2.0 * pi,
        UAngle::radian,
        UAngle::radian
    );

    elevation = convert_angle(elevation, UAngle::radian, angle_out);
    azimuth = convert_angle(azimuth, UAngle::radian, angle_out);

    azel = vec3d{azimuth, elevation, range};

    return azel;
}

inline vec3d enu_to_sez(const vec3d& r_enu) {
    // ENU = [E, N, U]
    // SEZ = [S, E, Z] = [-N, E, U]
    return vec3d{-r_enu(1), r_enu(0), r_enu(2)};
}
inline vec3d sez_to_enu(const vec3d& r_sez) {
    return vec3d{r_sez(1), -r_sez(0), r_sez(2)};
}

inline vec3d station_rel_bcbf(const vec3d& r_target_bcbf, const vec3d& r_station_bcbf) {
    return r_target_bcbf - r_station_bcbf;
}

inline vec3d azel_from_bcbf(
    const vec3d& r_target_bcbf,
    const vec3d& r_station_bcbf,
    f64 station_lat,
    f64 station_lon,
    UAngle angle_in = UAngle::degree,
    UAngle angle_out = UAngle::degree,
    f64 tol = tol12
) {
    vec3d r_rel_bcbf = station_rel_bcbf(r_target_bcbf, r_station_bcbf);
    vec3d r_rel_enu = bcbf_rel_to_enu(r_rel_bcbf, station_lat, station_lon, angle_in);
    return azel_from_enu(r_rel_enu, angle_out, tol);
}

inline vec3d azel_from_source_frame(
    const vec3d& r_target_source,
    const vec3d& r_station_source,
    const mat3d& R_BCBF_source, // R: source frame -> bcbf
    f64 station_lat,
    f64 station_lon,
    UAngle angle_in = UAngle::degree,
    UAngle angle_out = UAngle::degree,
    f64 tol = tol12
) {
    vec3d r_target_bcbf = R_BCBF_source * r_target_source;
    vec3d r_station_bcbf = R_BCBF_source * r_station_source;

    return azel_from_bcbf(
        r_target_bcbf,
        r_station_bcbf,
        station_lat,
        station_lon,
        angle_in,
        angle_out,
        tol
    );
}

inline vec3d azel_from_source_frame(
    const vec3d& r_target_source,
    const vec3d& r_station_source,
    const vec4d& q_BCBF_source,
    f64 station_lat,
    f64 station_lon,
    UAngle angle_in = UAngle::degree,
    UAngle angle_out = UAngle::degree,
    f64 tol = tol12
) {
    vec3d r_target_bcbf = ep_rotate_fast_passive(q_BCBF_source, r_target_source);
    vec3d r_station_bcbf = ep_rotate_fast_passive(q_BCBF_source, r_station_source);

    return azel_from_bcbf(
        r_target_bcbf,
        r_station_bcbf,
        station_lat,
        station_lon,
        angle_in,
        angle_out,
        tol
    );
}

inline vec3d azel_from_earth_frame(
    const vec3d& r_target_source,
    const vec3d& r_station_source,
    f64 station_lat,
    f64 station_lon,
    const JulianDate& jd,
    EarthFrame frame_source,
    const EarthOrientationParams& eop,
    TimeOffsets offsets,
    TimeScale scale_in = TimeScale::utc,
    UAngle angle_in = UAngle::degree,
    UAngle angle_out = UAngle::degree,
    f64 tol = tol12
) {
    mat3d R_ITRS_source
        = rot_earth_frame(jd, frame_source, EarthFrame::ITRS, eop, offsets, scale_in);

    vec3d r_target_bcbf = R_ITRS_source * r_target_source;
    vec3d r_station_bcbf = R_ITRS_source * r_station_source;

    return azel_from_bcbf(
        r_target_bcbf,
        r_station_bcbf,
        station_lat,
        station_lon,
        angle_in,
        angle_out,
        tol
    );
}

inline vec3d azel_rates_from_enu(
    const vec3d& r_rel_enu,
    const vec3d& v_rel_enu,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
) {
    vec3d azel_dot = vec3d0;

    f64 E = r_rel_enu(0);
    f64 N = r_rel_enu(1);
    f64 U = r_rel_enu(2);
    f64 E_dot = v_rel_enu(0);
    f64 N_dot = v_rel_enu(1);
    f64 U_dot = v_rel_enu(2);

    f64 rho2 = r_rel_enu.squaredNorm();
    f64 rho = std::sqrt(rho2);
    f64 rhoEN2 = (r_rel_enu.segment<2>(0)).squaredNorm();
    f64 rhoEN = std::sqrt(rhoEN2);
    if (rho <= tol || rhoEN <= tol) return azel_dot;

    f64 rho_dot = r_rel_enu.dot(v_rel_enu) / rho;
    f64 az_dot = (N * E_dot - E * N_dot) / (rhoEN2);
    f64 el_dot
        = (U_dot * rho - U * rho_dot) / (rho2 * std::sqrt(1 - std::pow(U / rho, 2)));

    if (angle_out != UAngle::radian) {
        az_dot = convert_angle(az_dot, UAngle::radian, angle_out);
        el_dot = convert_angle(el_dot, UAngle::radian, angle_out);
    }
    azel_dot = vec3d{az_dot, el_dot, rho_dot};
    return azel_dot;
}

inline vec3d azel_rates_from_bcbf(
    const vec3d& r_target_bcbf,
    const vec3d& v_target_bcbf,
    const vec3d& r_station_bcbf,
    const vec3d& v_station_bcbf,
    f64 station_lat,
    f64 station_lon,
    UAngle angle_in = UAngle::degree,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
) {
    // NOTE: both r_rel and v_rel need to be in the same frame, BCI -> BCBF requires
    // BKE/transport theorem addition before rotating into the BCBF frame if BCBF is
    // rotating within BCI
    // General form:
    // B_v = A_v + B_w_A x r
    // where B_w_A is the angular velocity of the A-frame in the B-frame, all left-scripts
    // denote derivatives taken in the scripted frame (i.e. B_v is the derivative of r in
    // the B frame), r is an arbitrary vector.

    // The components of these vectors are arbitrary if we want to express them in the B
    // frame but A_v and r are in the A frame we have to rotate them into the correct
    // frame (B_w_A should be too for ease of computation, but is not necessary as the
    // cross-product takes into also "crosses" the frame unit-vectors as well): B_v (in
    // B-coordinates) = R_B_A * (A_v + B_w_A x r) (in A-coordinates) where R_B_A is the
    // rotation matrix from A -> B (R_B_A: A -> B)

    // For our case:
    // v_bcbf = R_BCBF_BCI * (v_bci + BCBF_w_BCI_bci x r_bci)
    // where BCBF_w_BCI_bci is angular velocity of BCI wrt BCBF, expressed in BCI

    vec3d r_rel_bcbf = r_target_bcbf - r_station_bcbf;
    vec3d v_rel_bcbf = v_target_bcbf - v_station_bcbf;

    mat3d R_ENU_BCBF = rot_enu_from_bcbf(station_lat, station_lon, angle_in);

    vec3d r_rel_enu = R_ENU_BCBF * r_rel_bcbf;
    vec3d v_rel_enu = R_ENU_BCBF * v_rel_bcbf;

    return azel_rates_from_enu(r_rel_enu, v_rel_enu, angle_out, tol);
}

inline vec3d stat_r_bcbf_from_detic(
    const vec3d& llh,
    const Celestial& body,
    UAngle angle_in = UAngle::degree
) {
    if (body.semimajor_axis <= 0.0 || body.semiminor_axis <= 0.0) return vec3d0;
    return detic_to_bcbf(llh, body, angle_in);
}

inline vec3d stat_target_rel_enu(
    const vec3d& r_station_body,
    const vec3d& r_target_body,
    f64 lat, // planetodetic
    f64 lon,
    UAngle angle_in = UAngle::degree
) {
    vec3d r_rel_body = r_target_body - r_station_body;
    vec3d r_rel_enu = bcbf_rel_to_enu(r_rel_body, lat, lon, angle_in);
    return r_rel_enu;
}

inline mat3d stat_rot_enu_from_detic(const vec3d& llh, UAngle angle_in = UAngle::radian) {
    f64 lat = llh(0), lon = llh(1);

    return rot_enu_from_bcbf(lat, lon, angle_in);
}

inline vec3d stat_target_rel_enu(
    const vec3d& r_station_body,
    const vec3d& r_target_body,
    const vec3d& stat_llh_body,
    UAngle angle_in = UAngle::radian
) {
    return stat_target_rel_enu(
        r_station_body,
        r_target_body,
        stat_llh_body(0),
        stat_llh_body(1),
        angle_in
    );
}

inline vec4d stat_att_enu_from_detic(
    const StateAtt& x_anchor_att,
    const vec3d& llh,
    UAngle angle_in = UAngle::radian
) {
    // rotation to position on anchor, from BCBF (in lat/lon) to ENU
    mat3d R_ENU_BCBF = stat_rot_enu_from_detic(llh, angle_in);
    // rotation of anchor orientation, from sim inertial to BCBF
    mat3d R_BCBF_I = ep_to_dcm(x_anchor_att.q);
    mat3d R_ENU_I = R_ENU_BCBF * R_BCBF_I; // compose rotations

    return dcm_to_ep(R_ENU_I);
}

// inline vec3d azel_from_radec(){
//
// }

// TODO: add quaternion versions (faster with StateAtt.q vs getting DCM)