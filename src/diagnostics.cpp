// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/diagnostics.hpp"
#include "core/body.hpp"
#include "core/dynamics_rotational.hpp"
#include "core/entity.hpp"
#include "core/estimation_batch.hpp"
#include "core/estimation_common.hpp"
#include "core/estimation_recursive.hpp"
#include "core/estimation_world.hpp"
#include "core/history_io.hpp"
#include "core/ingest.hpp"
#include "core/integrator_adaptive.hpp"
#include "core/integrator_fixed.hpp"
#include "core/interpolation.hpp"
#include "core/measurement.hpp"
#include "core/measurement_world.hpp"
#include "core/observation_type.hpp"
#include "core/observations.hpp"
#include "core/od_dynamics.hpp"
#include "core/orbit_determination.hpp"
#include "core/orbital_elements.hpp"
#include "core/planets.hpp"
#include "core/scenario_io.hpp"
#include "core/state.hpp"
#include "core/station_geometry.hpp"
#include "core/time.hpp"
#include "core/transform.hpp"
#include "core/world.hpp"
#include "core/world_history.hpp"
#include "core/world_stepper.hpp"

#include "graphics/render_loop.hpp"
#include "util/constants.hpp"
#include "util/printing.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include "examples/scenarios.hpp"

#include "graphics/raygen.hpp"
#include "graphics/rdraw.hpp"
#include "graphics/renderer.hpp"
#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"

#include <chrono>
#include <filesystem>
#include <limits>
#include <memory>
#include <random>

void print_diag_title(const std::string& title = "") {
    std::string line = "-----------------------------------------------------------";
    size_t width = line.size();
    std::string spacer = " ";
    if (title.empty()) {
        spacer = "-";
    }
    std::string withSpace = title + spacer;
    std::string trimmed = withSpace.substr(0, width);
    line.replace(0, trimmed.size(), trimmed);
    std::println("{}", line);
}

void run_gravity_diag(
    World& world,
    EntityId earth_id,
    EntityId urath_id,
    EntityId sat_id,
    EntityId stat_id
) {
    Celestial* earth = world.celestial(earth_id);
    Celestial* urath = world.celestial(urath_id);
    Satellite* sat = world.satellite(sat_id);
    Station* stat = world.station(stat_id);

    vec3d a_from_earth = world.gravity_accel_from(sat_id, earth_id);
    vec3d a_from_urath = world.gravity_accel_from(sat_id, urath_id);
    vec3d a_total = world.gravity_accel_on(sat_id);

    std::println("Satellite Gravity Diagnostic -------------------------------");
    std::println("a_earth = {}, mag = {}", a_from_earth, a_from_earth.norm());
    std::println("a_urath = {}, mag = {}", a_from_urath, a_from_urath.norm());
    std::println("a_total = {}, mag = {}", a_total, a_total.norm());
}

void run_station_geo_diag(
    World& world,
    EntityId earth_id,
    EntityId urath_id,
    EntityId sat_id,
    EntityId stat_id
) {
    Celestial* earth = world.celestial(earth_id);
    Celestial* urath = world.celestial(urath_id);
    Satellite* sat = world.satellite(sat_id);
    Station* stat = world.station(stat_id);

    std::println("Station Geometry Diagnostic --------------------------------");
    std::println("stat_r_inertial: {}", world.stat_r_inertial(stat_id));
    std::println("stat_v_inertial: {}", world.stat_v_inertial(stat_id));
    std::println("z inertial earth: {}", world.body_z_inertial(earth_id));
    std::println("w inertial earth: {}", world.body_w_inertial(earth_id));
    std::println("w body earth: {}", earth->x_att.w);

    // std::println("stat_r_inertial: {}", world.stat_x_tr_inertial(stat_id).r);
    // std::println("stat_v_inertial: {}", world.stat_x_tr_inertial(stat_id).v);
}

void run_zonal_orientation_diag(
    World& world,
    EntityId earth_id,
    EntityId urath_id,
    EntityId sat_id,
    EntityId stat_id
) {
    Celestial* earth = world.celestial(earth_id);
    Celestial* urath = world.celestial(urath_id);
    Satellite* sat = world.satellite(sat_id);
    Station* stat = world.station(stat_id);

    sat->x_tr.r = 7000.0 * vec3d{std::cos(pio2), 0.0, std::sin(pio2)};

    vec4d earth_q0 = earth->x_att.q;

    // pointmass (no orientation)
    earth->x_att.q = q_identity;
    earth->gravity_model = GravityModel::pointmass;
    vec3d a_point = world.gravity_accel_from(sat_id, earth_id);

    // pointmass (tilted)
    earth->x_att.q = earth_q0;
    vec3d a_point_tilt = world.gravity_accel_from(sat_id, earth_id);

    // zonal (no orientation)
    earth->gravity_model = GravityModel::zonal;
    earth->x_att.q = q_identity;
    vec3d a_zonal = world.gravity_accel_from(sat_id, earth_id);

    // zonal (tilted)
    earth->x_att.q = earth_q0;
    vec3d a_zonal_tilt = world.gravity_accel_from(sat_id, earth_id);

    std::println("Zonal Orientation Diagnostic -------------------------------");
    std::println("a_point = {} (mag = {})", a_point, a_point.norm());
    std::println("a_point_tilt = {} (mag = {})", a_point_tilt, a_point_tilt.norm());
    std::println("a_zonal = {} (mag = {})", a_zonal, a_zonal.norm());
    std::println("a_zonal_tilt = {} (mag = {})", a_zonal_tilt, a_zonal_tilt.norm());
    std::println(
        "a_point tilt error = {} (mag = {})",
        a_point - a_point_tilt,
        (a_point - a_point_tilt).norm()
    );
    std::println(
        "a_zonal tilt error = {} (mag = {})",
        a_zonal - a_zonal_tilt,
        (a_zonal - a_zonal_tilt).norm()
    );
}

void run_spherical_harmonics_diag(
    World& world,
    EntityId earth_id,
    EntityId urath_id,
    EntityId sat_id,
    EntityId stat_id
) {
    Celestial* earth = world.celestial(earth_id);
    Celestial* urath = world.celestial(urath_id);
    Satellite* sat = world.satellite(sat_id);
    Station* stat = world.station(stat_id);

    i32 degree = 6;

    sat->x_tr.r = 7000.0 * vec3d{std::cos(pio2), 0.0, std::sin(pio2)};

    mat3d earth_tilt_dcm = rotX(23.44 * deg_to_rad);
    earth->x_att.q = dcm_to_ep(earth_tilt_dcm);

    // Pointmass
    earth->gravity_model = GravityModel::pointmass;
    vec3d a_point = world.gravity_accel_from(sat_id, earth_id);

    // Zonal, J_degree
    earth->gravity_model = GravityModel::zonal;
    earth->degree = degree;
    vec3d a_zonal = world.gravity_accel_from(sat_id, earth_id);

    // Spherical Harmonics, n = degree, m = 0
    earth->gravity_model = GravityModel::spherical_harmonics;
    earth->degree = degree;
    earth->order = 0;
    read_sphh_coefs(
        std::string(PROJECT_ROOT) + "/assets/egm2008_120.txt",
        earth->C,
        earth->S,
        earth->degree,
        earth->order
    );
    vec3d a_sphh = world.gravity_accel_from(sat_id, earth_id);

    std::println("Spherical Harmonics Diagnostic -----------------------------");
    std::println("a_point = {} (mag = {})", a_point, a_point.norm());
    std::println("a_zonal = {} (mag = {})", a_zonal, a_zonal.norm());
    std::println("a_sphh = {} (mag = {})", a_sphh, a_sphh.norm());
    vec3d error = a_zonal - a_sphh;
    std::println("error= {} (mag = {})", error, error.norm());
}

void run_sphh_longitude_diag(
    World& world,
    EntityId earth_id,
    EntityId urath_id,
    EntityId sat_id,
    EntityId stat_id
) {
    Celestial* earth = world.celestial(earth_id);
    Celestial* urath = world.celestial(urath_id);
    Satellite* sat = world.satellite(sat_id);
    Station* stat = world.station(stat_id);

    mat3d earth_tilt_dcm = rotX(23.44 * deg_to_rad);
    earth->x_att.q = dcm_to_ep(earth_tilt_dcm);

    i32 degree = 6;
    i32 order = 6;
    earth->gravity_model = GravityModel::spherical_harmonics;
    earth->degree = degree;
    earth->order = order;
    read_sphh_coefs(
        std::string(PROJECT_ROOT) + "/assets/egm2008_120.txt",
        earth->C,
        earth->S,
        earth->degree,
        earth->order
    );

    // Spherical Harmonics, zonal only
    vec3d llh = vec3d{-45.0, 0, 7000.0};       // in body fixed frame
    vec3d r_body = detic_to_bcbf(llh, *earth); // in body fixed frame
    sat->x_tr.r
        = ep_rotate_fast_passive(ep_conj(earth->x_att.q), r_body); // back to inertial
    earth->degree = degree;
    earth->order = 0;
    vec3d a_sphh_zonal_long1 = world.gravity_accel_from(sat_id, earth_id);

    // Spherical Harmonics, zonal only, different longitude
    llh = vec3d{-45.0, 57.0, 7000.0};
    r_body = detic_to_bcbf(llh, *earth);
    sat->x_tr.r = ep_rotate_fast_passive(ep_conj(earth->x_att.q), r_body);
    vec3d a_sphh_zonal_long2 = world.gravity_accel_from(sat_id, earth_id);

    // Spherical Harmonics
    llh = vec3d{-45.0, 0.0, 7000.0};
    r_body = detic_to_bcbf(llh, *earth);
    sat->x_tr.r = ep_rotate_fast_passive(ep_conj(earth->x_att.q), r_body);
    earth->order = order;
    vec3d a_sphh_long1 = world.gravity_accel_from(sat_id, earth_id);

    // Spherical Harmonics, different longitude
    llh = vec3d{-45.0, 57.0, 7000.0};
    r_body = detic_to_bcbf(llh, *earth);
    sat->x_tr.r = ep_rotate_fast_passive(ep_conj(earth->x_att.q), r_body);
    vec3d a_sphh_long2 = world.gravity_accel_from(sat_id, earth_id);

    std::println("Spherical Harmonics Diagnostic -----------------------------");
    std::println(
        "a_sphh_zonal_long1 = {} (mag = {})",
        a_sphh_zonal_long1,
        a_sphh_zonal_long1.norm()
    );
    std::println(
        "a_sphh_zonal_long2 = {} (mag = {})",
        a_sphh_zonal_long2,
        a_sphh_zonal_long2.norm()
    );
    f64 zonal_long_mag_error = a_sphh_zonal_long2.norm() - a_sphh_zonal_long1.norm();
    std::println("zonal_long_error = {}", zonal_long_mag_error);

    std::println("a_sphh_long1 = {} (mag = {})", a_sphh_long1, a_sphh_long1.norm());
    std::println("a_sphh_long2 = {} (mag = {})", a_sphh_long2, a_sphh_long2.norm());
    vec3d sphh_long_error = a_sphh_long2 - a_sphh_long1;
    std::println(
        "sphh_long_error = {} (mag = {})",
        sphh_long_error,
        sphh_long_error.norm()
    );
}

void run_epkde(
    World& world,
    EntityId earth_id,
    EntityId urath_id,
    EntityId sat_id,
    EntityId stat_id
) {
    Celestial* earth = world.celestial(earth_id);

    vec4d q = q_identity;
    vec3d w = 7.292115000000000e-05 * axis_z;

    vec4d q_dot = k_eulerparams(q, w);
    std::println("q = {}", q);
    std::println("w = {}", w);
    std::println("q_dot = {}", q_dot);
}

void run_earth_rot_diag(const JulianDate& jd, const EarthOrientationParams& eop) {
    mat3d R_ITRS_J2000 = rot_earth_frame(
        jd,
        EarthFrame::J2000,
        EarthFrame::ITRS,
        eop,
        eop.offsets,
        TimeScale::utc
    );

    mat3d R_J2000_ITRS = rot_earth_frame(
        jd,
        EarthFrame::ITRS,
        EarthFrame::J2000,
        eop,
        eop.offsets,
        TimeScale::utc
    );

    mat3d inv_err = R_J2000_ITRS * R_ITRS_J2000 - mat3d1;
    std::println("inverse error = {}", inv_err.norm());

    mat3d orth_err = R_ITRS_J2000.transpose() * R_ITRS_J2000 - mat3d1;
    std::println("orthogonality error = {}", orth_err.norm());

    f64 det = R_ITRS_J2000.determinant();
    std::println("R_ITRS_J2000 determinant = {}", det);

    mat3d R_TOD_J2000 = rot_earth_frame(
        jd,
        EarthFrame::J2000,
        EarthFrame::TOD,
        eop,
        eop.offsets,
        TimeScale::utc
    );
    mat3d R_ITRS_TOD = rot_earth_frame(
        jd,
        EarthFrame::TOD,
        EarthFrame::ITRS,
        eop,
        eop.offsets,
        TimeScale::utc
    );
    f64 inter_err = (R_ITRS_TOD * R_TOD_J2000 - R_ITRS_J2000).norm();
    std::println("Intermediate error = {}", inter_err);

    f64 theta = earth_rot_angle_from_jd(jd, eop.offsets, TimeScale::utc, UAngle::degree);
    mat3d R_GTOD_TOD = rot(theta, RotAxis::z, UAngle::degree);
    vec3d x_GTOD = R_GTOD_TOD * axis_x;
    mat3d R_GTOD_TOD_rotZ = rotZ(theta, UAngle::degree);
    vec3d x_GTOD_rotZ = R_GTOD_TOD_rotZ * axis_x;
    f64 x_err = (x_GTOD - x_GTOD_rotZ).norm();
    std::println("x_err = {}", x_err);
}

void run_station_obs_geom_diag(
    Celestial& earth,
    const JulianDate& jd,
    const EarthOrientationParams& eop
) {
    f64 lat = 0.0;
    f64 lon = 0.0;

    // equator / prime meridian
    // BCBF +y -> east, +z -> north, +x -> up
    vec3d east_enu = bcbf_rel_to_enu(axis_y, lat, lon);
    vec3d north_enu = bcbf_rel_to_enu(axis_z, lat, lon);
    vec3d up_enu = bcbf_rel_to_enu(axis_x, lat, lon);

    std::println("east error = {}", (east_enu - axis_x).norm());
    std::println("north error = {}", (north_enu - axis_y).norm());
    std::println("up error = {}", (up_enu - axis_z).norm());

    // body-fixed station observation
    vec3d r_station_bcbf = vec3d{earth.semimajor_axis, 0.0, 0.0};

    vec3d r_target_bcbf = r_station_bcbf + axis_y; // due east
    vec3d azel = azel_from_bcbf(r_target_bcbf, r_station_bcbf, lat, lon);
    std::println("bcbf due east azel = {}", azel);

    r_target_bcbf = r_station_bcbf + axis_z; // due north
    azel = azel_from_bcbf(r_target_bcbf, r_station_bcbf, lat, lon);
    std::println("bcbf due north azel = {}", azel);

    r_target_bcbf = r_station_bcbf + axis_x; // overhead
    azel = azel_from_bcbf(r_target_bcbf, r_station_bcbf, lat, lon);
    std::println("bcbf overhead azel = {}", azel);

    // due east
    r_target_bcbf = r_station_bcbf + axis_y;
    azel = azel_from_bcbf(
        r_target_bcbf,
        r_station_bcbf,
        lat,
        lon,
        UAngle::degree,
        UAngle::radian
    );
    std::println("bcbf due east azel rad = {}", azel);

    // nonzero longitude
    lat = 0.0;
    lon = 90.0;
    r_station_bcbf = vec3d{0.0, earth.semimajor_axis, 0.0};
    r_target_bcbf = r_station_bcbf - axis_x;
    azel = azel_from_bcbf(r_target_bcbf, r_station_bcbf, lat, lon);
    std::println("bcbf lon 90 due east azel = {}", azel);

    // at north pole
    lat = 90.0;
    lon = 0.0;
    r_station_bcbf = vec3d{0.0, 0.0, earth.semiminor_axis};
    r_target_bcbf = r_station_bcbf + axis_z;
    azel = azel_from_bcbf(r_target_bcbf, r_station_bcbf, lat, lon);
    std::println("bcbf north pole overhead azel = {}", azel);

    // source frame identity
    lat = 0.0;
    lon = 0.0;
    r_station_bcbf = vec3d{earth.semimajor_axis, 0.0, 0.0};
    r_target_bcbf = r_station_bcbf + axis_y;
    azel = azel_from_source_frame(r_target_bcbf, r_station_bcbf, mat3d1, lat, lon);
    std::println("source identity due east azel = {}", azel);

    // source frame rotated
    mat3d R_source_bcbf = rotZ(30.0, UAngle::degree);
    mat3d R_bcbf_source = R_source_bcbf.transpose();
    vec3d r_station_source = R_source_bcbf * r_station_bcbf;
    vec3d r_target_source = R_source_bcbf * r_target_bcbf;
    azel = azel_from_source_frame(
        r_target_source,
        r_station_source,
        R_bcbf_source,
        lat,
        lon
    );
    std::println("source rotated due east azel = {}", azel);

    // earth frame identity
    azel = azel_from_earth_frame(
        r_target_bcbf,
        r_station_bcbf,
        lat,
        lon,
        jd,
        EarthFrame::ITRS,
        eop,
        eop.offsets,
        TimeScale::utc
    );
    std::println("earth ITRS due east azel = {}", azel);

    // earth frame J2000
    mat3d R_J2000_ITRS = rot_earth_frame(
        jd,
        EarthFrame::ITRS,
        EarthFrame::J2000,
        eop,
        eop.offsets,
        TimeScale::utc
    );
    vec3d r_station_j2000 = R_J2000_ITRS * r_station_bcbf;
    vec3d r_target_j2000 = R_J2000_ITRS * r_target_bcbf;
    azel = azel_from_earth_frame(
        r_target_j2000,
        r_station_j2000,
        lat,
        lon,
        jd,
        EarthFrame::J2000,
        eop,
        eop.offsets,
        TimeScale::utc
    );
    std::println("earth J2000 due east azel = {}", azel);

    // ENU angle checks
    azel = azel_from_enu(axis_x);
    std::println("enu east azel = {}", azel);

    azel = azel_from_enu(axis_y);
    std::println("enu north azel = {}", azel);

    azel = azel_from_enu(axis_z);
    std::println("enu up azel = {}", azel);

    azel = azel_from_enu(vec3d0);
    std::println("enu zero azel = {}", azel);

    // ENU rate checks
    vec3d azel_dot = azel_rates_from_enu(axis_y, axis_x);
    std::println("enu rate north/east azel_dot = {}", azel_dot);

    // BCBF rate checks
    lat = 0.0;
    lon = 0.0;
    r_station_bcbf = vec3d{earth.semimajor_axis, 0.0, 0.0};
    vec3d v_station_bcbf = vec3d0;
    r_target_bcbf = r_station_bcbf + axis_z;
    vec3d v_target_bcbf = axis_y;
    azel_dot = azel_rates_from_bcbf(
        r_target_bcbf,
        v_target_bcbf,
        r_station_bcbf,
        v_station_bcbf,
        lat,
        lon
    );
    std::println("bcbf rate north/east azel_dot = {}", azel_dot);
}

void run_rv_coe_diag(const Celestial& body) {
    f64 R_body = body.mean_radius;
    f64 mu = body.mu;

    // non-circular non-equatorial
    f64 theta = pio4;
    f64 st = std::sin(theta);
    f64 ct = std::cos(theta);
    vec3d r = vec3d{1.5, 0.0, 0.0} * R_body;
    vec3d v = vec3d{0.0, 1.5 * ct, st} * std::sqrt(mu / r.norm());
    OEClassical coe = rv_to_classical(r, v, mu);
    StateTr x = classical_to_rv(coe, mu);
    f64 r_err = (r - x.r).norm();
    f64 v_err = (v - x.v).norm();
    std::println("Non-circular, Non-equatorial");
    std::println("Position error: {}", r_err);
    std::println("Velocity error: {}", v_err);

    // equatorial circular
    theta = 0.0;
    st = std::sin(theta);
    ct = std::cos(theta);
    r = vec3d{1.5, 0.0, 0.0} * R_body;
    v = vec3d{0.0, ct, st} * std::sqrt(mu / r.norm());
    coe = rv_to_classical(r, v, mu);
    x = classical_to_rv(coe, mu);
    r_err = (r - x.r).norm();
    v_err = (v - x.v).norm();
    std::println("Circular, equatorial");
    std::println("Position error: {}", r_err);
    std::println("Velocity error: {}", v_err);

    // equatorial elliptical
    theta = 0.0;
    st = std::sin(theta);
    ct = std::cos(theta);
    r = vec3d{1.5, 0.0, 0.0} * R_body;
    v = vec3d{0.0, 3.0 * ct, st} * std::sqrt(mu / r.norm());
    coe = rv_to_classical(r, v, mu);
    x = classical_to_rv(coe, mu);
    r_err = (r - x.r).norm();
    v_err = (v - x.v).norm();
    std::println("Non-Circular, equatorial");
    std::println("Position error: {}", r_err);
    std::println("Velocity error: {}", v_err);

    // circular non-equatorial
    theta = pi / 3.0;
    st = std::sin(theta);
    ct = std::cos(theta);
    r = vec3d{1.5, 0.0, 0.0} * R_body;
    v = vec3d{0.0, ct, st} * std::sqrt(mu / r.norm());
    coe = rv_to_classical(r, v, mu);
    x = classical_to_rv(coe, mu);
    r_err = (r - x.r).norm();
    v_err = (v - x.v).norm();
    std::println("Circular, Non-equatorial");
    std::println("Position error: {}", r_err);
    std::println("Velocity error: {}", v_err);

    f64 a = 1.5 * R_body;

    auto run_coe_base = [mu](const char* label, OEClassical coe0) {
        StateTr x0 = classical_to_rv(coe0, mu);
        OEClassical coe1 = rv_to_classical(x0, mu);
        StateTr x1 = classical_to_rv(coe1, mu);

        std::println("{}", label);
        std::println("Position error: {}", (x0.r - x1.r).norm());
        std::println("Velocity error: {}", (x0.v - x1.v).norm());
        std::println(
            "coe = [{}, {}, {}, {}, {}, {}]",
            coe1.sma,
            coe1.ecc,
            coe1.inc,
            coe1.raan,
            coe1.aop,
            coe1.ta
        );
    };

    // coe seed, non-circular non-equatorial
    run_coe_base(
        "COE seed: Non-circular, Non-equatorial",
        OEClassical{a, 0.2, 0.7, 0.4, 0.3, 0.8}
    );

    // coe seed, circular equatorial
    run_coe_base(
        "COE seed: Circular, equatorial",
        OEClassical{a, 0.0, 0.0, 0.0, 0.0, 0.8}
    );

    // coe seed, non-circular equatorial
    run_coe_base(
        "COE seed: Non-circular, equatorial",
        OEClassical{a, 0.2, 0.0, 0.0, 0.3, 0.8}
    );

    // coe seed, circular non-equatorial
    run_coe_base(
        "COE seed: Circular, Non-equatorial",
        OEClassical{a, 0.0, 0.7, 0.4, 0.0, 0.8}
    );
}

void run_radec_diag() {
    vec3d radec;
    vec3d radec_dot;
    vec3d r_rel;
    vec3d v_rel;

    r_rel = axis_x;
    radec = radec_from_rel(r_rel, UAngle::degree);
    std::println("r_rel = {}\nRA = {}\nDec = {}", r_rel, radec(0), radec(1));

    r_rel = axis_y;
    radec = radec_from_rel(r_rel, UAngle::degree);
    std::println("r_rel = {}\nRA = {}\nDec = {}", r_rel, radec(0), radec(1));

    r_rel = axis_z;
    radec = radec_from_rel(r_rel, UAngle::degree);
    std::println("r_rel = {}\nRA = {}\nDec = {}", r_rel, radec(0), radec(1));

    r_rel = axis_x;
    v_rel = axis_y;
    radec_dot = radec_rates_from_rel(r_rel, v_rel);
    std::println(
        "r_rel = {}\nv_rel = {}\nRA_dot = {}\nDec_dot = {}",
        r_rel,
        v_rel,
        radec_dot(0),
        radec_dot(1)
    );
}

f64 tof_elliptic_ta(
    f64 sma,
    f64 ecc,
    f64 ta1,
    f64 ta2,
    f64 mu,
    UAngle angle_in = UAngle::radian,
    i32 n_rev = 0,
    f64 tol = tol12
) {
    // temp, for diags only
    if (mu <= 0.0 || sma <= 0.0 || ecc < 0.0 || ecc >= 1.0 || n_rev < 0) {
        return std::numeric_limits<f64>::quiet_NaN();
    }

    if (angle_in != UAngle::radian) {
        ta1 = convert_angle(ta1, angle_in, UAngle::radian);
        ta2 = convert_angle(ta2, angle_in, UAngle::radian);
    }

    auto mean_anomaly_from_ta = [ecc](f64 ta) {
        f64 half_ta = 0.5 * ta;
        f64 E = 2.0
                * std::atan2(
                    std::sqrt(1.0 - ecc) * std::sin(half_ta),
                    std::sqrt(1.0 + ecc) * std::cos(half_ta)
                );
        E = wrap_angle(E, 0.0, twopi);
        return E - ecc * std::sin(E);
    };

    f64 M1 = mean_anomaly_from_ta(ta1);
    f64 M2 = mean_anomaly_from_ta(ta2);
    f64 dM = wrap_angle(M2 - M1, 0.0, twopi) + static_cast<f64>(n_rev) * twopi;
    f64 n = std::sqrt(mu / (sma * sma * sma));
    if (n <= tol) return std::numeric_limits<f64>::quiet_NaN();

    return dM / n;
}

void run_iod_diag(const Celestial& body) {
    OEClassical coes{
        .sma = body.mean_radius * 2.0,
        .ecc = 0.1,
        .inc = pio4,
        .raan = 0.0,
        .aop = 0.0
    };

    f64 ta1 = 0.0;
    f64 ta2 = 0.2;
    f64 ta3 = 0.4;

    coes.ta = ta1;
    StateTr x1 = classical_to_rv(coes, body.mu, UAngle::radian);
    coes.ta = ta2;
    StateTr x2 = classical_to_rv(coes, body.mu, UAngle::radian);
    coes.ta = ta3;
    StateTr x3 = classical_to_rv(coes, body.mu, UAngle::radian);

    f64 t1 = 0.0;
    f64 t2 = tof_elliptic_ta(coes.sma, coes.ecc, ta1, ta2, body.mu);
    f64 t3 = tof_elliptic_ta(coes.sma, coes.ecc, ta1, ta3, body.mu);

    vec3d R1{0.0, 0.0, body.mean_radius};
    vec3d R2 = R1;
    vec3d R3 = R1;

    // relative positions
    vec3d r_rel1 = x1.r - R1;
    vec3d r_rel2 = x2.r - R2;
    vec3d r_rel3 = x3.r - R3;

    vec3d radec1 = radec_from_rel(r_rel1);
    vec3d radec2 = radec_from_rel(r_rel2);
    vec3d radec3 = radec_from_rel(r_rel3);

    vec3d t{t1, t2, t3};

    vec3d ra{radec1(0), radec2(0), radec3(0)};
    vec3d dec{radec1(1), radec2(1), radec3(1)};

    mat3d R;
    R.col(0) = R1;
    R.col(1) = R2;
    R.col(2) = R3;

    IODAnglesObs3 arc = iod_angles3_from_radec(t, ra, dec, R);

    // Gauss
    IODResult result_gauss = iod_gauss(arc, body.mu);
    std::println("Gauss success = {}", result_gauss.success);
    if (result_gauss.success) {
        std::println("Gauss status = {}", i32(result_gauss.status));
        std::println("Gauss iterations = {}", result_gauss.iterations);
        std::println("r2 error = {}", (result_gauss.x.r - x2.r).norm());
        std::println("v2 error = {}", (result_gauss.x.v - x2.v).norm());
    }

    // Gibbs
    IODResult result_gibbs = iod_gibbs(x1.r, x2.r, x3.r, body.mu);
    if (result_gibbs.success) {
        std::println("Gibbs success = {}", result_gibbs.success);
        std::println("Gibbs status = {}", i32(result_gibbs.status));
        std::println("r2 error = {}", (result_gibbs.x.r - x2.r).norm());
        std::println("v2 error = {}", (result_gibbs.x.v - x2.v).norm());
    }

    // Herrick-Gibbs
    IODResult result_herrickgibbs
        = iod_herrickgibbs(t1, t2, t3, x1.r, x2.r, x3.r, body.mu);
    if (result_herrickgibbs.success) {
        std::println("Herrick-Gibbs success = {}", result_herrickgibbs.success);
        std::println("Herrick-Gibbs status = {}", i32(result_herrickgibbs.status));
        std::println("r2 error = {}", (result_herrickgibbs.x.r - x2.r).norm());
        std::println("v2 error = {}", (result_herrickgibbs.x.v - x2.v).norm());
    }
}

void run_od_prop_diag(const Celestial& body) {
    OEClassical coes{
        .sma = body.mean_radius * 2.0,
        .ecc = 0.1,
        .inc = pio4,
        .raan = 0.0,
        .aop = 0.0
    };

    f64 ta1 = 0.0;
    f64 ta2 = 0.2;

    coes.ta = ta1;
    StateTr x1 = classical_to_rv(coes, body.mu, UAngle::radian);
    coes.ta = ta2;
    StateTr x2_truth = classical_to_rv(coes, body.mu, UAngle::radian);
    f64 tof = tof_elliptic_ta(coes.sma, coes.ecc, ta1, ta2, body.mu);

    i32 n_steps = 100;
    f64 dt = tof / n_steps;

    ODDynamicsConfig cfg{.tr_model = ODTrDynamicsModel::two_body, .mu = body.mu};
    StateTr x2_prop = propagate_tr_od(0.0, x1, dt, n_steps, cfg);

    std::println("OD Propagation Error");
    StateTr x2_err = x2_truth - x2_prop;
    std::println("Position error = {}", x2_err.r.norm());
    std::println("Velocity error = {}", x2_err.v.norm());

    // STM propagation

    std::println("STM Propagation Error");
    StateTr x0 = x1;
    VarStateTr y0;
    y0.x = x0;
    y0.Phi = mat6d1;

    VarStateTr yf = propagate_var_tr_od(0.0, y0, dt, n_steps, cfg);
    StateTr xf = propagate_tr_od(0.0, x0, dt, n_steps, cfg);

    std::println("State/var position error = {}", (yf.x.r - xf.r).norm());
    std::println("State/var velocity error = {}", (yf.x.v - xf.v).norm());

    vec6d x0_vec = statetr_to_vec6d(x0);
    vec6d xf_vec = statetr_to_vec6d(xf);

    f64 eps_pos = 1e-3; // km
    f64 eps_vel = 1e-6; // km/s
    std::println("FD eps position = {}", eps_pos);
    std::println("FD eps velocity = {}", eps_vel);

    mat6d Phi_fd;
    // Build each STM column by perturbing one initial state component,
    // propagating it, then differencing against the nominal final state
    // forward finite difference:
    // https://www.dam.brown.edu/people/alcyew/handouts/numdiff.pdf
    // TODO: maybe add jacobian-free (use finite diff) version for od?
    for (i32 i = 0; i < 6; ++i) {
        f64 eps_i = (i < 3) ? eps_pos : eps_vel;

        vec6d x0_pert_vec = x0_vec;
        x0_pert_vec(i) += eps_i;
        StateTr x0_pert = vec6d_to_statetr(x0_pert_vec);
        StateTr xf_pert = propagate_tr_od(0.0, x0_pert, dt, n_steps, cfg);
        vec6d xf_pert_vec = statetr_to_vec6d(xf_pert);
        Phi_fd.col(i) = (xf_pert_vec - xf_vec) / eps_i;
    }

    mat6d Phi_err = yf.Phi - Phi_fd;
    f64 Phi_err_rel = Phi_err.norm() / yf.Phi.norm();
    std::println("STM fd error norm = {}", Phi_err.norm());
    std::println("STM fd max error = {}", Phi_err.cwiseAbs().maxCoeff());
    std::println("STM fd relative error = {}", Phi_err_rel);

    std::println("Phi norm = {}", yf.Phi.norm());
    std::println("Phi fd norm = {}", Phi_fd.norm());
}

void run_od_prop_adaptive_diag(const Celestial& body) {
    print_diag_title("OD Propagation Error (Adaptive)");
    OEClassical coes{
        .sma = body.mean_radius * 2.0,
        .ecc = 0.1,
        .inc = pio4,
        .raan = 0.0,
        .aop = 0.0
    };

    f64 ta1 = 0.0;
    f64 ta2 = 0.2;

    coes.ta = ta1;
    StateTr x1 = classical_to_rv(coes, body.mu, UAngle::radian);
    coes.ta = ta2;
    StateTr x2_truth = classical_to_rv(coes, body.mu, UAngle::radian);
    f64 tof = tof_elliptic_ta(coes.sma, coes.ecc, ta1, ta2, body.mu);

    i32 n_steps = 24;
    f64 dt = tof / n_steps;

    f64 t0 = 0.0;

    ODDynamicsConfig cfg{.tr_model = ODTrDynamicsModel::two_body, .mu = body.mu};
    f64 t_interval = tof;
    StateTr x2_prop = propagate_tr_od(t0, x1, t_interval, n_steps, cfg);

    AdaptiveIntegratorConfig adaptive_cfg{};
    f64 tf = t0 + t_interval;

    auto run_case = [&](const AdaptiveIntegratorConfig& integrator_cfg,
                        const StateTr& x0,
                        f64 interval,
                        const ODDynamicsConfig* dyn_cfg = nullptr) {
        const ODDynamicsConfig& selected_dyn_cfg = dyn_cfg ? *dyn_cfg : cfg;
        auto time_start = std::chrono::steady_clock::now();
        auto result = propagate_tr_od_adaptive(
            t0,
            x0,
            interval,
            selected_dyn_cfg,
            integrator_cfg
        );
        auto time_end = std::chrono::steady_clock::now();
        f64 runtime_us
            = std::chrono::duration<f64, std::micro>(time_end - time_start).count();
        return std::pair{result, runtime_us};
    };

    auto print_result = [&](const std::string& label,
                            const AdaptivePropagationResult<StateTr>& result,
                            f64 runtime_us,
                            f64 requested_tf,
                            const StateTr* truth = nullptr) {
        std::println("{}", label);
        std::println("Status: {}", status_string(result.status));
        std::println("Requested/returned time: {} / {}", requested_tf, result.t);
        std::println("Endpoint error: {}", std::abs(result.t - requested_tf));
        if (truth && result.status == StatusCode::ok) {
            StateTr error = *truth - result.x;
            std::println("Position error: {}", error.r.norm());
            std::println("Velocity error: {}", error.v.norm());
        }
        std::println("Final error norm: {}", result.final_error_norm);
        std::println(
            "Attempts/accepted/rejected: {} / {} / {}",
            result.stats.attempted_steps,
            result.stats.accepted_steps,
            result.stats.rejected_steps
        );
        std::println("Derivative evaluations: {}", result.stats.deriv_evals);
        std::println(
            "Min/max/final dt: {} / {} / {}",
            result.stats.min_accepted_dt,
            result.stats.max_accepted_dt,
            result.stats.final_accepted_dt
        );
        std::println("Runtime (us): {}", runtime_us);
    };

    auto print_status_check = [](StatusCode actual, StatusCode expected) {
        std::println(
            "Expected status: {}, passed: {}",
            status_string(expected),
            actual == expected
        );
    };

    auto [adaptive_res, adaptive_runtime_us] = run_case(adaptive_cfg, x1, t_interval);

    StateTr x2_err = x2_truth - x2_prop;

    std::println("\nBaseline ---------------------------------------------------");
    std::println("Position error (fixed) = {}", x2_err.r.norm());
    std::println("Velocity error (fixed) = {}", x2_err.v.norm());
    std::println(
        "Fixed Stepper Stats:\nFixed dt: {},\nDerivative Evals: {}",
        dt,
        n_steps * 4
    );
    print_result("Adaptive result:", adaptive_res, adaptive_runtime_us, tf, &x2_truth);
    std::println(
        "Step growth observed: {}",
        adaptive_res.stats.max_accepted_dt > adaptive_cfg.dt_initial
    );
    std::println(
        "Stats consistent: {}",
        adaptive_res.stats.attempted_steps
                == adaptive_res.stats.accepted_steps + adaptive_res.stats.rejected_steps
            && adaptive_res.stats.deriv_evals == adaptive_res.stats.attempted_steps * 7
    );

    std::println("\nOversized Initial Step -------------------------------------");
    AdaptiveIntegratorConfig oversized_cfg = adaptive_cfg;
    oversized_cfg.rel_tol = 1e-12;
    oversized_cfg.abs_tol_r = 1e-12;
    oversized_cfg.abs_tol_v = 1e-15;
    oversized_cfg.dt_initial = tof;
    oversized_cfg.dt_max = tof;
    auto [oversized_res, oversized_runtime_us] = run_case(oversized_cfg, x1, t_interval);
    print_result(
        "Oversized-step result:",
        oversized_res,
        oversized_runtime_us,
        tf,
        &x2_truth
    );
    std::println(
        "Rejected then recovered: {}",
        oversized_res.status == StatusCode::ok && oversized_res.stats.rejected_steps > 0
    );

    std::println("\nTolerance Sweep --------------------------------------------");
    struct ToleranceCase {
        std::string label;
        f64 rel_tol;
        f64 abs_tol_r;
        f64 abs_tol_v;
    };
    array<ToleranceCase, 3> tolerance_cases{{
        {.label = "Loose", .rel_tol = 1e-6, .abs_tol_r = 1e-6, .abs_tol_v = 1e-9},
        {.label = "Medium", .rel_tol = 1e-9, .abs_tol_r = 1e-9, .abs_tol_v = 1e-12},
        {.label = "Tight", .rel_tol = 1e-12, .abs_tol_r = 1e-12, .abs_tol_v = 1e-15},
    }};

    array<f64, 3> tolerance_position_errors{};
    for (size_t i = 0; i < tolerance_cases.size(); ++i) {
        AdaptiveIntegratorConfig tolerance_cfg = adaptive_cfg;
        tolerance_cfg.rel_tol = tolerance_cases[i].rel_tol;
        tolerance_cfg.abs_tol_r = tolerance_cases[i].abs_tol_r;
        tolerance_cfg.abs_tol_v = tolerance_cases[i].abs_tol_v;
        auto [tolerance_res, tolerance_runtime_us]
            = run_case(tolerance_cfg, x1, t_interval);
        if (tolerance_res.status == StatusCode::ok) {
            tolerance_position_errors[i] = (x2_truth - tolerance_res.x).r.norm();
        } else {
            tolerance_position_errors[i] = inf<f64>;
        }
        print_result(
            tolerance_cases[i].label + " tolerance:",
            tolerance_res,
            tolerance_runtime_us,
            tf,
            &x2_truth
        );
    }
    std::println(
        "Tighter tolerances reduce position error: {}",
        tolerance_position_errors[1] < tolerance_position_errors[0]
            && tolerance_position_errors[2] < tolerance_position_errors[1]
    );

    std::println("\nBackward Propagation ---------------------------------------");
    AdaptivePropagationResult<StateTr> backward_res{};
    f64 backward_runtime_us = 0.0;
    if (adaptive_res.status == StatusCode::ok) {
        auto time_start = std::chrono::steady_clock::now();
        backward_res = propagate_tr_od_adaptive(
            tf,
            adaptive_res.x,
            -t_interval,
            cfg,
            adaptive_cfg
        );
        auto time_end = std::chrono::steady_clock::now();
        backward_runtime_us
            = std::chrono::duration<f64, std::micro>(time_end - time_start).count();
    }
    print_result("Forward/backward closure:", backward_res, backward_runtime_us, t0, &x1);

    std::println("\nInvalid Configuration --------------------------------------");
    auto run_invalid_config
        = [&](const std::string& label, const AdaptiveIntegratorConfig& invalid_cfg) {
              auto [result, runtime_us] = run_case(invalid_cfg, x1, t_interval);
              print_result(label, result, runtime_us, tf);
              print_status_check(result.status, StatusCode::invalid_input);
          };

    AdaptiveIntegratorConfig zero_tolerance_cfg = adaptive_cfg;
    zero_tolerance_cfg.rel_tol = 0.0;
    run_invalid_config("Zero relative tolerance:", zero_tolerance_cfg);

    AdaptiveIntegratorConfig invalid_bounds_cfg = adaptive_cfg;
    invalid_bounds_cfg.dt_min = 20.0;
    invalid_bounds_cfg.dt_max = 10.0;
    run_invalid_config("Invalid step bounds:", invalid_bounds_cfg);

    AdaptiveIntegratorConfig zero_initial_step_cfg = adaptive_cfg;
    zero_initial_step_cfg.dt_initial = 0.0;
    run_invalid_config("Zero initial step:", zero_initial_step_cfg);

    std::println("\nFailure Statuses -------------------------------------------");
    StateTr non_finite_state = x1;
    non_finite_state.r(0) = qNaN<f64>;
    auto [non_finite_state_res, non_finite_state_runtime_us]
        = run_case(adaptive_cfg, non_finite_state, t_interval);
    print_result(
        "Non-finite initial state:",
        non_finite_state_res,
        non_finite_state_runtime_us,
        tf
    );
    print_status_check(non_finite_state_res.status, StatusCode::invalid_state);

    ODDynamicsConfig non_finite_dyn_cfg = cfg;
    non_finite_dyn_cfg.mu = qNaN<f64>;
    auto [non_finite_result_res, non_finite_result_runtime_us]
        = run_case(adaptive_cfg, x1, t_interval, &non_finite_dyn_cfg);
    print_result(
        "Non-finite dynamics result:",
        non_finite_result_res,
        non_finite_result_runtime_us,
        tf
    );
    print_status_check(non_finite_result_res.status, StatusCode::non_finite_result);

    AdaptiveIntegratorConfig rejection_limit_cfg = oversized_cfg;
    rejection_limit_cfg.max_rejections = 1;
    auto [rejection_limit_res, rejection_limit_runtime_us]
        = run_case(rejection_limit_cfg, x1, t_interval);
    print_result("Rejection limit:", rejection_limit_res, rejection_limit_runtime_us, tf);
    print_status_check(rejection_limit_res.status, StatusCode::max_rejections_reached);

    AdaptiveIntegratorConfig attempt_limit_cfg = oversized_cfg;
    attempt_limit_cfg.max_attempts = 1;
    attempt_limit_cfg.max_rejections = 1000;
    auto [attempt_limit_res, attempt_limit_runtime_us]
        = run_case(attempt_limit_cfg, x1, t_interval);
    print_result("Attempt limit:", attempt_limit_res, attempt_limit_runtime_us, tf);
    print_status_check(attempt_limit_res.status, StatusCode::max_steps_reached);

    AdaptiveIntegratorConfig underflow_cfg = oversized_cfg;
    underflow_cfg.dt_min = tof * 0.5;
    auto [underflow_res, underflow_runtime_us] = run_case(underflow_cfg, x1, t_interval);
    print_result("Step underflow:", underflow_res, underflow_runtime_us, tf);
    print_status_check(underflow_res.status, StatusCode::step_size_underflow);

    std::println("\nZero Duration ----------------------------------------------");
    auto [zero_duration_res, zero_duration_runtime_us] = run_case(adaptive_cfg, x1, 0.0);
    print_result(
        "Zero-duration result:",
        zero_duration_res,
        zero_duration_runtime_us,
        t0,
        &x1
    );
    std::println(
        "Zero-duration counts unchanged: {}",
        zero_duration_res.stats.attempted_steps == 0
            && zero_duration_res.stats.accepted_steps == 0
            && zero_duration_res.stats.rejected_steps == 0
            && zero_duration_res.stats.deriv_evals == 0
    );
}

void run_measurement_jacobian_diag() {
    MeasurementContext ctx;
    // zero observer state
    ctx.x_tr_observer.r = vec3d{0.0, 0.0, 0.0};
    ctx.x_tr_observer.v = vec3d{0.0, 0.0, 0.0};

    // arbitrary target state
    ctx.x_tr_target.r = vec3d{7000.0, 1000.0, 1300.0};
    ctx.x_tr_target.v = vec3d{-0.5, 7.2, 1.0};

    matXd H_fd = jacobian_fd_measurement(ObservationType::radec, ctx);
    matXd H_an = jacobian_radec(ctx);

    matXd H_err = H_fd - H_an;

    std::println("RA/Dec H fd/an error norm: {}", H_err.norm());
    std::println("RA/Dec H fd/an max error: {}", H_err.cwiseAbs().maxCoeff());

    // should be zero
    std::println("RA/Dec H fd velocity cols norm: {}", H_fd.block(0, 3, 2, 3).norm());

    // degree output check
    matXd H_fd_deg = jacobian_fd_measurement(
        ObservationType::radec,
        ctx,
        UAngle::radian,
        UAngle::degree
    );

    matXd H_an_deg = jacobian_radec(ctx, UAngle::degree);
    matXd H_err_deg = H_fd_deg - H_an_deg;

    std::println("RA/Dec H deg fd/an error norm: {}", H_err_deg.norm());
    std::println("RA/Dec H deg fd/an max error: {}", H_err_deg.cwiseAbs().maxCoeff());
}

void run_batch_od_diag() {

    vec3d llh = vec3d{0.0, 0.0, 0.0}; // [lat, lon, h] = [deg, deg, sim units]
    EarthStationSatScenario scenario = make_earth_station_sat_scenario(llh);
    World& world = scenario.world;
    EntityId earth_id = scenario.earth_id;
    EntityId stat_id = scenario.stat_id;
    EntityId sat_id = scenario.sat_id;
    Celestial* earth = world.celestial(earth_id);
    Satellite* sat = world.satellite(sat_id);

    StateTr x0_truth;
    x0_truth.r = vec3d{7000.0, 1000.0, 1300.0};
    x0_truth.v = vec3d{-0.5, 7.2, 1.0};
    sat->x_tr = x0_truth;

    f64 t_meas0 = 0.0;
    f64 t_measf = 600.0;
    f64 N_meas = 50;
    vecXd t_meas = vecXd::LinSpaced(N_meas, t_meas0, t_measf);
    f64 sigma_rad = 1e-6;
    f64 sigma_range = 1e-3;
    f64 sigma_range_rate = 1e-6;
    f64 sigma_r = 1e-3;
    f64 sigma_v = 1e-6;

    struct BatchDiagCase {
        const char* name = "";
        bool use_radec = false;
        bool use_range = false;
        bool use_range_rate = false;
        bool use_pos = false;
        bool use_pos_vel = false;
    };

    auto run_case = [&](const BatchDiagCase& diag_case) {
        ODBatchInput input;
        input.x0_guess = x0_truth;
        input.dyn_config.tr_model = ODTrDynamicsModel::two_body;
        input.dyn_config.mu = earth->mu;
        input.t0 = 0.0;
        input.max_iters = 10;
        input.prop_steps = 200;
        input.tol_dx = 1e-6;
        input.tol_residual = 1e-8;

        // add perturbations
        input.x0_guess.r += vec3d{1.0, -1.0, 0.5};
        input.x0_guess.v += vec3d{1e-3, -1e-3, 5e-4};

        std::mt19937_64 rng(12345);
        std::normal_distribution<f64> noise_unit(0.0, 1.0);

        for (i32 i = 0; i < N_meas; ++i) {
            // propagate to get measurements
            sat->x_tr = propagate_tr_od(
                0.0,
                x0_truth,
                t_meas(i),
                input.prop_steps,
                input.dyn_config
            );
            world.set_stat_anchor_detic(stat_id, earth_id, llh);
            StateTr x_tr_obsv = world.stat_x_tr_inertial(stat_id);

            Measurement meas;
            meas.target_id = sat_id;
            meas.observer_id = stat_id;

            if (diag_case.use_radec) {
                // get measurement (RADec)
                meas.t = t_meas(i);
                meas.type = ObservationType::radec;
                meas.z = world_predict_measurement(
                    world,
                    meas.type,
                    stat_id,
                    sat_id,
                    UAngle::radian,
                    UAngle::radian
                );
                meas.z(0) += noise_unit(rng) * sigma_rad;
                meas.z(1) += noise_unit(rng) * sigma_rad;
                meas.R = matXd1<2> * sigma_rad * sigma_rad;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_tr_obsv);
            }
            if (diag_case.use_range) {
                // get measurement (range)
                meas.t = t_meas(i);
                meas.type = ObservationType::range;
                meas.z = world_predict_measurement(
                    world,
                    meas.type,
                    stat_id,
                    sat_id,
                    UAngle::radian,
                    UAngle::radian
                );
                meas.z(0) += noise_unit(rng) * sigma_range;
                meas.R = matXd1<1> * sigma_range * sigma_range;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_tr_obsv);
            }
            if (diag_case.use_range_rate) {
                // get measurement (range_rate)
                meas.t = t_meas(i);
                meas.type = ObservationType::range_rate;
                meas.z = world_predict_measurement(
                    world,
                    meas.type,
                    stat_id,
                    sat_id,
                    UAngle::radian,
                    UAngle::radian
                );
                meas.z(0) += noise_unit(rng) * sigma_range_rate;
                meas.R = matXd1<1> * sigma_range_rate * sigma_range_rate;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_tr_obsv);
            }
            if (diag_case.use_pos) {
                // get measurement (pos)
                meas.t = t_meas(i);
                meas.type = ObservationType::pos;
                meas.z = world_predict_measurement(
                    world,
                    meas.type,
                    stat_id,
                    sat_id,
                    UAngle::radian,
                    UAngle::radian
                );
                meas.z(0) += noise_unit(rng) * sigma_r;
                meas.z(1) += noise_unit(rng) * sigma_r;
                meas.z(2) += noise_unit(rng) * sigma_r;
                meas.R = matXd::Identity(3, 3) * sigma_r * sigma_r;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_tr_obsv);
            }
            if (diag_case.use_pos_vel) {
                // get measurement (posvel)
                meas.t = t_meas(i);
                meas.type = ObservationType::pos_vel;
                meas.z = world_predict_measurement(
                    world,
                    meas.type,
                    stat_id,
                    sat_id,
                    UAngle::radian,
                    UAngle::radian
                );
                meas.z(0) += noise_unit(rng) * sigma_r;
                meas.z(1) += noise_unit(rng) * sigma_r;
                meas.z(2) += noise_unit(rng) * sigma_r;
                meas.z(3) += noise_unit(rng) * sigma_v;
                meas.z(4) += noise_unit(rng) * sigma_v;
                meas.z(5) += noise_unit(rng) * sigma_v;
                meas.R = matXd::Zero(6, 6);
                meas.R.block(0, 0, 3, 3) = mat3d1 * sigma_r * sigma_r;
                meas.R.block(3, 3, 3, 3) = mat3d1 * sigma_v * sigma_v;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_tr_obsv);
            }
        }

        // solve
        ODBatchResult result = od_batch_lumve(input);
        f64 initial_err = (statetr_to_vec6d(input.x0_guess - x0_truth)).norm();
        f64 final_err = (statetr_to_vec6d(result.x0_est - x0_truth)).norm();
        f64 final_r_err = (result.x0_est.r - x0_truth.r).norm();
        f64 final_v_err = (result.x0_est.v - x0_truth.v).norm();
        i32 meas_per_epoch = static_cast<i32>(diag_case.use_radec)
                             + static_cast<i32>(diag_case.use_range)
                             + static_cast<i32>(diag_case.use_range_rate)
                             + static_cast<i32>(diag_case.use_pos)
                             + static_cast<i32>(diag_case.use_pos_vel);

        std::println("-----------------------------------------------------------");
        std::println("LUMVE Case: {}", diag_case.name);
        std::println(
            "RADec: {}, Range: {}, Range Rate: {}, Position: {}, State: {}",
            diag_case.use_radec,
            diag_case.use_range,
            diag_case.use_range_rate,
            diag_case.use_pos,
            diag_case.use_pos_vel
        );
        std::println("LUMVE Epochs = {}", N_meas);
        std::println("LUMVE Measurements Per Epoch = {}", meas_per_epoch);
        std::println("LUMVE Initial Error = {}", initial_err);
        std::println("LUMVE Final Error = {}", final_err);
        std::println("LUMVE Final Position Error = {}", final_r_err);
        std::println("LUMVE Final Velocity Error = {}", final_v_err);
        std::println("LUMVE Success = {}", od_status_success(result.status));
        std::println("LUMVE Status: {}", status_string(result.status));
        std::println("LUMVE Measurements = {}", input.measurements.size());
        std::println("LUMVE Iterations = {}", result.iterations);
        std::println("LUMVE Residual Norm = {}", result.residual_norm);
        std::println("LUMVE Raw Residual Norm = {}", result.raw_residual_norm);
        std::println("LUMVE Delta x Norm = {}", result.dx_norm);
        std::println("LUMVE Covariance Norm = {}", result.covariance.norm());
    };

    svec<BatchDiagCase> cases{
        {.name = "RADec Only", .use_radec = true},
        {.name = "RADec + Range", .use_radec = true, .use_range = true},
        {.name = "RADec + Range + Range Rate",
         .use_radec = true,
         .use_range = true,
         .use_range_rate = true},
        {.name = "Position Only", .use_pos = true},
        {.name = "State Only", .use_pos_vel = true},
        {.name = "All Measurements",
         .use_radec = true,
         .use_range = true,
         .use_range_rate = true,
         .use_pos = true,
         .use_pos_vel = true}
    };

    for (const BatchDiagCase& diag_case : cases) {
        run_case(diag_case);
    }
}

void run_checkpoint_diag() {
    World world;
    EntityId cel_id = world.spawn_celestial();
    EntityId sat_id = world.spawn_satellite();
    EntityId stat_id = world.spawn_station();
    world.reset_time(0.0);
    Satellite* sat = world.satellite(sat_id);
    sat->x_tr.r = vec3d{1.0, 1.0, 1.0} * 7000;
    sat->x_tr.v = vec3d{0.0, 7.5, 0.0};
    std::println("Initial World:");
    std::println("Sim Time: {}", world.t_sim());
    std::println("Sat Pos = {}", sat->x_tr.r);
    std::println("Sat Vel = {}", sat->x_tr.v);
    WorldStateSnapshot world_state_snapshot = world.capture_checkpoint();

    world.reset_time(100.0);
    sat->x_tr.r = vec3d{1.0, 1.0, 1.0} * 1000;
    sat->x_tr.v = vec3d{0.0, 0.0, 0.0};
    std::println("Mutated World:");
    std::println("Sim Time: {}", world.t_sim());
    std::println("Sat Pos = {}", sat->x_tr.r);
    std::println("Sat Vel = {}", sat->x_tr.v);

    bool restored = world.restore_checkpoint_state(world_state_snapshot);
    std::println("Restored World: {}", restored);
    std::println("Sim Time: {}", world.t_sim());
    std::println("Sat Pos = {}", sat->x_tr.r);
    std::println("Sat Vel = {}", sat->x_tr.v);
}

void run_station_anchor_diag() {
    vec3d llh = vec3d{0.0, 0.0, 0.0}; // [lat, lon, h] = [deg, deg, sim units]
    EarthStationSatScenario scenario = make_earth_station_sat_scenario(llh);
    World& world = scenario.world;
    EntityId earth_id = scenario.earth_id;
    EntityId stat_id = scenario.stat_id;
    EntityId sat_id = scenario.sat_id;
    Celestial* earth = world.celestial(earth_id);
    Satellite* sat = world.satellite(sat_id);
    Station* stat = world.station(stat_id);

    std::println("Station Anchor Diagnostic -------------------------------");
    std::println("station anchor id = {}", stat->anchor_id);
    std::println("station r_body_BCBF = {}", stat->r_body_BCBF);
    std::println("station llh_BCBF = {}", stat->llh_BCBF);

    vec3d r_stat_expected_BCBF = vec3d{earth->semimajor_axis, 0.0, 0.0};
    vec3d r_stat_BCI = world.stat_r_inertial(stat_id);
    std::println("station inertial pos = {}", r_stat_BCI);
    std::println(
        "station inertial pos error = {}",
        (r_stat_BCI - r_stat_expected_BCBF).norm()
    );

    mat3d R_ENU_BCBF = world.stat_rot_enu_from_body(stat_id);
    vec3d y_BCBF_to_ENU = R_ENU_BCBF * axis_y;
    vec3d z_BCBF_to_ENU = R_ENU_BCBF * axis_z;
    vec3d x_BCBF_to_ENU = R_ENU_BCBF * axis_x;
    std::println(
        "+Y_BCBF -> ENU = {}, error = {}",
        y_BCBF_to_ENU,
        (y_BCBF_to_ENU - axis_x).norm()
    );
    std::println(
        "+Z_BCBF -> ENU = {}, error = {}",
        z_BCBF_to_ENU,
        (z_BCBF_to_ENU - axis_y).norm()
    );
    std::println(
        "+X_BCBF -> ENU = {}, error = {}",
        x_BCBF_to_ENU,
        (x_BCBF_to_ENU - axis_z).norm()
    );

    sat->x_tr.r = earth->x_tr.r + vec3d{earth->semimajor_axis + 1000.0, 0.0, 0.0};
    vec3d r_rel_up_ENU = world.stat_rel_enu(stat_id, sat_id);
    std::println("target overhead ENU = {}", r_rel_up_ENU);
    std::println(
        "target overhead ENU error = {}",
        (r_rel_up_ENU - 1000.0 * axis_z).norm()
    );

    sat->x_tr.r = earth->x_tr.r + vec3d{earth->semimajor_axis, 1000.0, 0.0};
    vec3d r_rel_east_ENU = world.stat_rel_enu(stat_id, sat_id);
    std::println("target east ENU = {}", r_rel_east_ENU);
    std::println("target east ENU error = {}", (r_rel_east_ENU - 1000.0 * axis_x).norm());

    // different cases for llh
    const svec<vec3d> llh_cases = {
        vec3d{0.0, 90.0, 0.0}, // equator
        vec3d{45.0, 0.0, 0.0}, // mid lat
        vec3d{90.0, 0.0, 0.0}, // north pole
    };
    for (const vec3d& llh_case : llh_cases) {
        bool set_case = world.set_stat_anchor_detic(stat_id, earth_id, llh_case);
        std::println("LLH case [deg, deg, sim] = {}, set = {}", llh_case, set_case);
        std::println("station r_body_BCBF = {}", stat->r_body_BCBF);
        std::println("station llh_BCBF = {}", stat->llh_BCBF);
        mat3d R_case_ENU_BCBF = world.stat_rot_enu_from_body(stat_id);
        std::println("station R_ENU_BCBF row 0 = {}", vec3d{R_case_ENU_BCBF.row(0)});
        std::println("station R_ENU_BCBF row 1 = {}", vec3d{R_case_ENU_BCBF.row(1)});
        std::println("station R_ENU_BCBF row 2 = {}", vec3d{R_case_ENU_BCBF.row(2)});
    }

    // station on rotating body, expect v_BCI = w x r_BCI
    world.set_stat_anchor_detic(stat_id, earth_id, llh);
    earth->x_att.q = vec4d{0.0, 0.0, 0.0, 1.0};
    earth->x_att.w = vec3d{0.0, 0.0, 7.292115000000000e-05};
    vec3d v_stat_BCI_expected
        = earth->x_att.w.cross(vec3d{earth->semimajor_axis, 0.0, 0.0});
    vec3d v_stat_BCI = world.stat_v_inertial(stat_id);
    std::println("station spin velocity BCI = {}", v_stat_BCI);
    std::println(
        "station spin velocity error = {}",
        (v_stat_BCI - v_stat_BCI_expected).norm()
    );

    // check az/el
    // put sat overhead station, el = 90deg expected
    sat->x_tr.r = earth->x_tr.r + vec3d{earth->semimajor_axis + 1000.0, 0.0, 0.0};
    vec3d azel_overhead = azel_from_enu(world.stat_rel_enu(stat_id, sat_id));
    std::println("target overhead azelrho = {}", azel_overhead);
    std::println(
        "target overhead elevation error = {}",
        std::abs(azel_overhead(1) - 90.0)
    );
    // put sat directly east, el = 0deg and az = 90deg expected
    sat->x_tr.r = earth->x_tr.r + vec3d{earth->semimajor_axis, 1000.0, 0.0};
    vec3d azel_east = azel_from_enu(world.stat_rel_enu(stat_id, sat_id));
    std::println("target east azelrho = {}", azel_east);
    std::println("target east azimuth error = {}", std::abs(azel_east(0) - 90.0));
    std::println("target east elevation error = {}", std::abs(azel_east(1)));
    std::println("-----------------------------------------------------------");
}

void run_world_measurement_diag() {
    vec3d llh = vec3d{0.0, 0.0, 0.0}; // [lat, lon, h] = [deg, deg, sim units]
    EarthStationSatScenario scenario = make_earth_station_sat_scenario(llh);
    World& world = scenario.world;
    EntityId earth_id = scenario.earth_id;
    EntityId stat_id = scenario.stat_id;
    EntityId sat_id = scenario.sat_id;
    Celestial* earth = world.celestial(earth_id);
    Satellite* sat = world.satellite(sat_id);

    // satellite + measurements
    sat->x_tr.v = earth->x_tr.v + vec3d{0.1, 0.2, 0.3};

    // overhead
    sat->x_tr.r = earth->x_tr.r + vec3d{earth->semimajor_axis + 1000.0, 0.0, 0.0};
    vecXd z_azel_over = world_predict_measurement(
        world,
        ObservationType::azel,
        stat_id,
        sat_id,
        UAngle::radian,
        UAngle::degree
    );
    vecXd z_range_over = world_predict_measurement(
        world,
        ObservationType::range,
        stat_id,
        sat_id,
        UAngle::radian,
        UAngle::degree
    );
    vecXd z_range_rate_over = world_predict_measurement(
        world,
        ObservationType::range_rate,
        stat_id,
        sat_id,
        UAngle::radian,
        UAngle::degree
    );
    StateTr x_stat_over = world.stat_x_tr_inertial(stat_id);
    StateTr x_rel_over = sat->x_tr - x_stat_over;
    f64 rho_dot_over_expected = x_rel_over.r.dot(x_rel_over.v) / x_rel_over.r.norm();

    std::println("World Measurement Diagnostic ------------------------------");
    std::println("overhead azel = {}", z_azel_over);
    std::println("overhead elevation error = {}", std::abs(z_azel_over(1) - 90.0));
    std::println("overhead range = {}", z_range_over);
    std::println("overhead range error = {}", std::abs(z_range_over(0) - 1000.0));
    std::println("overhead range-rate = {}", z_range_rate_over);
    std::println(
        "overhead range-rate error = {}",
        std::abs(z_range_rate_over(0) - rho_dot_over_expected)
    );

    sat->x_tr.r = earth->x_tr.r + vec3d{earth->semimajor_axis, 1000.0, 0.0};
    vecXd z_azel_east = world_predict_measurement(
        world,
        ObservationType::azel,
        stat_id,
        sat_id,
        UAngle::radian,
        UAngle::degree
    );
    vecXd z_range_east = world_predict_measurement(
        world,
        ObservationType::range,
        stat_id,
        sat_id,
        UAngle::radian,
        UAngle::degree
    );
    vecXd z_range_rate_east = world_predict_measurement(
        world,
        ObservationType::range_rate,
        stat_id,
        sat_id,
        UAngle::radian,
        UAngle::degree
    );
    StateTr x_stat_east = world.stat_x_tr_inertial(stat_id);
    StateTr x_rel_east = sat->x_tr - x_stat_east;
    f64 rho_dot_east_expected = x_rel_east.r.dot(x_rel_east.v) / x_rel_east.r.norm();

    std::println("east azel = {}", z_azel_east);
    std::println("east azimuth error = {}", std::abs(z_azel_east(0) - 90.0));
    std::println("east elevation error = {}", std::abs(z_azel_east(1)));
    std::println("east range = {}", z_range_east);
    std::println("east range error = {}", std::abs(z_range_east(0) - 1000.0));
    std::println("east range-rate = {}", z_range_rate_east);
    std::println(
        "east range-rate error = {}",
        std::abs(z_range_rate_east(0) - rho_dot_east_expected)
    );

    vecXd z_radec_world = world_predict_measurement(
        world,
        ObservationType::radec,
        stat_id,
        sat_id,
        UAngle::radian,
        UAngle::degree
    );
    vecXd z_radec_explicit = predict_measurement(
        ObservationType::radec,
        sat->x_tr,
        world.stat_x_tr_inertial(stat_id),
        UAngle::degree
    );
    std::println("east radec world = {}", z_radec_world);
    std::println("east radec explicit = {}", z_radec_explicit);
    std::println("east radec error = {}", (z_radec_world - z_radec_explicit).norm());
    std::println("-----------------------------------------------------------");
}

void run_ekf_world_diag() {
    vec3d llh = vec3d{0.0, 0.0, 0.0}; // [lat, lon, h] = [deg, deg, sim units]
    EarthStationSatScenario scenario = make_earth_station_sat_scenario(llh);
    World& world = scenario.world;
    EntityId earth_id = scenario.earth_id;
    EntityId stat_id = scenario.stat_id;
    EntityId sat_id = scenario.sat_id;
    Celestial* earth = world.celestial(earth_id);
    Satellite* sat = world.satellite(sat_id);

    StateTr x0_truth;
    x0_truth.r = vec3d{7000.0, 1000.0, 1300.0};
    x0_truth.v = vec3d{-0.5, 7.2, 1.0};
    sat->x_tr = x0_truth;

    f64 t_meas0 = 0.0;
    f64 t_measf = 600.0;
    f64 N_meas = 50;
    vecXd t_meas = vecXd::LinSpaced(N_meas, t_meas0, t_measf);
    f64 sigma_rad = 1e-6;
    f64 sigma_range = 1e-3;
    f64 sigma_range_rate = 1e-6;
    f64 sigma_r = 1e-3;
    f64 sigma_v = 1e-6;

    struct EKFDiagCase {
        const char* name = "";
        bool use_radec = false;
        bool use_range = false;
        bool use_range_rate = false;
        bool use_pos = false;
        bool use_pos_vel = false;
    };

    auto run_case = [&](const EKFDiagCase& diag_case) {
        ODEKFOfflineInput input;
        input.initial_filter.x = x0_truth;
        input.initial_filter.t = 0.0;
        input.initial_filter.P = mat6d1;
        input.dyn_config.tr_model = ODTrDynamicsModel::two_body;
        input.dyn_config.mu = earth->mu;
        input.prop_steps = 200;
        input.Q = mat6d1 * 1e-4;

        // add perturbations
        input.initial_filter.x.r += vec3d{1.0, -1.0, 0.5};
        input.initial_filter.x.v += vec3d{1e-3, -1e-3, 5e-4};

        std::mt19937_64 rng(12345);
        std::normal_distribution<f64> noise_unit(0.0, 1.0);

        for (i32 i = 0; i < N_meas; ++i) {
            // propagate to get measurements
            sat->x_tr = propagate_tr_od(
                0.0,
                x0_truth,
                t_meas(i),
                input.prop_steps,
                input.dyn_config
            );

            // use stationary observer for now
            world.set_stat_anchor_detic(stat_id, earth_id, llh);
            StateTr x_tr_obsv = world.stat_x_tr_inertial(stat_id);

            Measurement meas;
            meas.target_id = sat_id;
            meas.observer_id = stat_id;

            if (diag_case.use_radec) {
                // get measurement (RADec)
                meas.t = t_meas(i);
                meas.type = ObservationType::radec;
                meas.z = world_predict_measurement(
                    world,
                    meas.type,
                    stat_id,
                    sat_id,
                    UAngle::radian,
                    UAngle::radian
                );
                meas.z(0) += noise_unit(rng) * sigma_rad;
                meas.z(1) += noise_unit(rng) * sigma_rad;
                meas.R = matXd1<2> * sigma_rad * sigma_rad;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_tr_obsv);
            }
            if (diag_case.use_range) {
                // get measurement (range)
                meas.t = t_meas(i);
                meas.type = ObservationType::range;
                meas.z = world_predict_measurement(
                    world,
                    meas.type,
                    stat_id,
                    sat_id,
                    UAngle::radian,
                    UAngle::radian
                );
                meas.z(0) += noise_unit(rng) * sigma_range;
                meas.R = matXd1<1> * sigma_range * sigma_range;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_tr_obsv);
            }
            if (diag_case.use_range_rate) {
                // get measurement (range_rate)
                meas.t = t_meas(i);
                meas.type = ObservationType::range_rate;
                meas.z = world_predict_measurement(
                    world,
                    meas.type,
                    stat_id,
                    sat_id,
                    UAngle::radian,
                    UAngle::radian
                );
                meas.z(0) += noise_unit(rng) * sigma_range_rate;
                meas.R = matXd1<1> * sigma_range_rate * sigma_range_rate;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_tr_obsv);
            }
            if (diag_case.use_pos) {
                // get measurement (pos)
                meas.t = t_meas(i);
                meas.type = ObservationType::pos;
                meas.z = world_predict_measurement(
                    world,
                    meas.type,
                    stat_id,
                    sat_id,
                    UAngle::radian,
                    UAngle::radian
                );
                meas.z(0) += noise_unit(rng) * sigma_r;
                meas.z(1) += noise_unit(rng) * sigma_r;
                meas.z(2) += noise_unit(rng) * sigma_r;
                meas.R = matXd::Identity(3, 3) * sigma_r * sigma_r;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_tr_obsv);
            }
            if (diag_case.use_pos_vel) {
                // get measurement (posvel)
                meas.t = t_meas(i);
                meas.type = ObservationType::pos_vel;
                meas.z = world_predict_measurement(
                    world,
                    meas.type,
                    stat_id,
                    sat_id,
                    UAngle::radian,
                    UAngle::radian
                );
                meas.z(0) += noise_unit(rng) * sigma_r;
                meas.z(1) += noise_unit(rng) * sigma_r;
                meas.z(2) += noise_unit(rng) * sigma_r;
                meas.z(3) += noise_unit(rng) * sigma_v;
                meas.z(4) += noise_unit(rng) * sigma_v;
                meas.z(5) += noise_unit(rng) * sigma_v;
                meas.R = matXd::Zero(6, 6);
                meas.R.block(0, 0, 3, 3) = mat3d1 * sigma_r * sigma_r;
                meas.R.block(3, 3, 3, 3) = mat3d1 * sigma_v * sigma_v;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_tr_obsv);
            }
        }

        // run ekf
        ODEKFResult result = od_ekf_offline(input);

        f64 initial_err = (statetr_to_vec6d(input.initial_filter.x - x0_truth)).norm();
        f64 final_err = (statetr_to_vec6d(result.filter.x - sat->x_tr)).norm();
        f64 final_r_err = (result.filter.x.r - sat->x_tr.r).norm();
        f64 final_v_err = (result.filter.x.v - sat->x_tr.v).norm();
        i32 meas_per_epoch = static_cast<i32>(diag_case.use_radec)
                             + static_cast<i32>(diag_case.use_range)
                             + static_cast<i32>(diag_case.use_range_rate)
                             + static_cast<i32>(diag_case.use_pos)
                             + static_cast<i32>(diag_case.use_pos_vel);

        std::println("-----------------------------------------------------------");
        std::println("MIXED EKF Case: {}", diag_case.name);
        std::println(
            "RADec: {}, Range: {}, Range Rate: {}, Position: {}, State: {}",
            diag_case.use_radec,
            diag_case.use_range,
            diag_case.use_range_rate,
            diag_case.use_pos,
            diag_case.use_pos_vel
        );
        std::println("MIXED EKF Epochs = {}", N_meas);
        std::println("MIXED EKF Measurements Per Epoch = {}", meas_per_epoch);
        std::println("MIXED EKF Initial Error = {}", initial_err);
        std::println("MIXED EKF Final Error = {}", final_err);
        std::println("MIXED EKF Final Position Error = {}", final_r_err);
        std::println("MIXED EKF Final Velocity Error = {}", final_v_err);
        std::println("MIXED EKF Final Time = {}", result.filter.t);
        std::println("MIXED EKF Success = {}", od_status_success(result.status));
        std::println("MIXED EKF Status: {}", status_string(result.status));
        std::println(
            "MIXED EKF Processed Measurements = {}",
            result.processed_measurements
        );
        std::println("MIXED EKF Total Measurements = {}", input.measurements.size());
        std::println("MIXED EKF Residual Norm = {}", result.residual_norm);
        std::println("MIXED EKF Raw Residual Norm = {}", result.raw_residual_norm);
        std::println("MIXED EKF Final Covariance Norm = {}", result.filter.P.norm());
    };

    svec<EKFDiagCase> cases{
        {.name = "RADec Only", .use_radec = true},
        {.name = "RADec + Range", .use_radec = true, .use_range = true},
        {.name = "RADec + Range + Range Rate",
         .use_radec = true,
         .use_range = true,
         .use_range_rate = true},
        {.name = "Position Only", .use_pos = true},
        {.name = "State Only", .use_pos_vel = true},
        {.name = "All Measurements",
         .use_radec = true,
         .use_range = true,
         .use_range_rate = true,
         .use_pos = true,
         .use_pos_vel = true}
    };

    for (const EKFDiagCase& diag_case : cases) {
        run_case(diag_case);
    }
}

void run_world_stepper_diag() {
    print_diag_title("World Stepper");

    auto run_case = [&](IntegratorTypeFixed integrator_tr) {
        vec3d llh = vec3d{0.0, 0.0, 0.0}; // [lat, lon, h] = [deg, deg, sim units]
        EarthStationSatScenario scenario = make_earth_station_sat_scenario(llh);
        World& world = scenario.world;
        EntityId earth_id = scenario.earth_id;
        EntityId sat_id = scenario.sat_id;
        EntityId stat_id = scenario.stat_id;
        Celestial* earth = world.celestial(earth_id);
        Satellite* sat = world.satellite(sat_id);
        Station* stat = world.station(stat_id);

        // earth
        earth->propagate_att = true;
        earth->set_spin_rate(0.1);

        // station
        // shouldn't propagate attitude because anchored
        stat->anchored = true;               // for testing
        stat->mass_properties.active = true; // for testing
        stat->propagate_att = true;          // for testing
        // need false, true, true to propagate station attitude
        stat->mass_properties.mass = 500.0;
        stat->mass_properties.principal_axes = true;
        stat->mass_properties.I.diagonal() = vec3d{100.0, 200.0, 300.0};
        stat->x_att.q = q_identity;
        stat->x_att.w = {0, 0, 0.01};
        stat->propagate_tr = true;
        stat->propagate_att = true;

        // satellite
        StateTr x0;
        x0.r = vec3d{7000.0, 0.0, 0.0};
        x0.v = vec3d{0.0, 7.546053290107541, 0.0};

        sat->mass_properties.active = true;
        sat->mass_properties.mass = 500.0;
        sat->mass_properties.principal_axes = true;
        sat->mass_properties.I.diagonal() = vec3d{100.0, 200.0, 300.0};
        sat->x_att.q = q_identity;
        sat->x_att.w = {0, 0, 0.01};
        sat->x_tr = x0;
        sat->propagate_tr = true;
        sat->propagate_att = true;

        WorldStepperConfig cfg;
        cfg.step_tr = true;
        cfg.step_att = true;
        cfg.substeps = 1;
        cfg.ticks = 10;
        cfg.dt_scale = 1.0 / cfg.ticks;
        cfg.integrator_tr = integrator_tr;
        cfg.integrator_att = IntegratorTypeFixed::rk4;

        f64 t_span = 100.0;
        f64 t0 = 0.0;
        i32 n_steps = 500;
        f64 dt = t_span / n_steps;
        world.reset_time(t0);

        ODDynamicsConfig ODcfg;
        ODcfg.R_cb_ref = earth->mean_radius;
        ODcfg.integrator = integrator_tr;
        ODcfg.mu = earth->mu;
        ODcfg.tr_model = ODTrDynamicsModel::two_body;

        // "truth"
        i32 n_ref_steps = n_steps * cfg.ticks * cfg.substeps;
        StateTr x_OD = propagate_tr_od(t0, x0, t_span, n_ref_steps, ODcfg);

        // world stepper
        WorldStepperStats stats;
        StatusCode step_status = StatusCode::ok;
        for (i32 i = 0; i < n_steps; ++i) {
            WorldStepResult step_result = step_world(world, dt, cfg);
            stats += step_result.stats;
            step_status = step_result.status;
            if (step_status != StatusCode::ok) break;
        }

        StateTr x_err = x_OD - sat->x_tr;
        f64 pos_err = x_err.r.norm();
        f64 vel_err = x_err.v.norm();
        f64 pos_rel_err = pos_err / x_OD.r.norm();
        f64 vel_rel_err = vel_err / x_OD.v.norm();

        std::println(
            "{} Success = {}",
            integrator_name(integrator_tr),
            step_status == StatusCode::ok
        );
        std::println("{} Position Error = {}", integrator_name(integrator_tr), pos_err);
        std::println("{} Velocity Error = {}", integrator_name(integrator_tr), vel_err);
        std::println(
            "{} Relative Position Error = {}",
            integrator_name(integrator_tr),
            pos_rel_err
        );
        std::println(
            "{} Relative Velocity Error = {}",
            integrator_name(integrator_tr),
            vel_rel_err
        );
        std::println("{} Final Time = {}", integrator_name(integrator_tr), world.t_sim());
        std::println();
    };

    run_case(IntegratorTypeFixed::rk1);
    run_case(IntegratorTypeFixed::rk2);
    run_case(IntegratorTypeFixed::rk2_heun);
    run_case(IntegratorTypeFixed::rk2_ralston);
    run_case(IntegratorTypeFixed::rk3);
    run_case(IntegratorTypeFixed::rk4);
}

void run_body_fixed_gravity_timing_diag() {
    print_diag_title("Body Fixed Gravity Timing");

    // celestial simple spin
    {
        World world;
        EntityId earth_id = wgs84(world);
        Celestial* earth = world.celestial(earth_id);

        f64 spin_rate = 0.1;
        f64 t_span = 10.0;
        f64 dt = 0.01;
        i32 n_steps = static_cast<i32>(t_span / dt);

        earth->x_att.q = q_identity;
        earth->x_att.w = {0.0, 0.0, spin_rate};
        earth->attitude_model = CelestialAttitudeModel::simple_spin;
        earth->propagate_att = true;

        WorldStepperConfig cfg;
        cfg.step_tr = false;
        cfg.step_att = true;
        cfg.integrator_tr = IntegratorTypeFixed::rk1;
        cfg.integrator_att = IntegratorTypeFixed::rk4;
        cfg.substeps = 1;
        cfg.ticks = 1;
        cfg.dt_scale = 1.0;

        WorldStepperStats stats;
        StatusCode step_status = StatusCode::ok;
        for (i32 i = 0; i < n_steps; ++i) {
            WorldStepResult step_result = step_world(world, dt, cfg);
            stats += step_result.stats;
            step_status = step_result.status;
            if (step_status != StatusCode::ok) break;
        }

        f64 theta = spin_rate * world.t_sim();
        vec4d q_expected{0.0, 0.0, std::sin(theta / 2.0), std::cos(theta / 2.0)};

        std::println("Simple Spin Success = {}", step_status == StatusCode::ok);
        std::println("Simple Spin Final Time = {}", world.t_sim());
        std::println("Simple Spin q = {}", earth->x_att.q);
        std::println("Simple Spin Expected q = {}", q_expected);
        std::println("Simple Spin q Error = {}", (earth->x_att.q - q_expected).norm());
        std::println("Simple Spin q Norm = {}", earth->x_att.q.norm());
        std::println();
    }

    // body-fixed source orientation diff
    {
        World world;
        EntityId earth_id = wgs84(world);
        EntityId sat_id = world.spawn_satellite();

        Celestial* earth = world.celestial(earth_id);
        Satellite* sat = world.satellite(sat_id);

        // set earth/satellite state
        earth->propagate_tr = false;
        earth->propagate_att = true;
        earth->x_tr.r = vec3d0;
        earth->x_att.q = q_identity;
        earth->x_att.w = {0.0, 0.0, 0.1};
        sat->x_tr.r = {7000.0, 500.0, 1000.0};

        // different orientation
        f64 theta = 1.0;
        vec4d q_spin{0.0, 0.0, std::sin(theta / 2.0), std::cos(theta / 2.0)};

        // integration options
        WorldStepperConfig cfg;
        cfg.step_tr = true;
        cfg.step_att = true;
        cfg.substeps = 1;
        cfg.ticks = 10;
        cfg.dt_scale = 1.0 / cfg.ticks;
        cfg.integrator_tr = IntegratorTypeFixed::rk1;
        cfg.integrator_att = IntegratorTypeFixed::rk4;

        f64 t_span = 100.0;
        f64 t0 = 0.0;
        i32 n_steps = 1000;
        f64 dt = t_span / n_steps;
        world.reset_time(t0);

        WorldStateSnapshot world_snapshot = world.capture_checkpoint();

        WorldStepperStats stats;
        StatusCode step_status = StatusCode::ok;
        auto prop = [&]() {
            stats = WorldStepperStats{};
            step_status = StatusCode::ok;
            for (i32 i = 0; i < n_steps; ++i) {
                WorldStepResult step_result = step_world(world, dt, cfg);
                stats += step_result.stats;
                step_status = step_result.status;
                if (step_status != StatusCode::ok) break;
            }
        };

        // point-mass
        earth->gravity_model = GravityModel::pointmass;
        earth->x_att.q = q_identity;
        vec3d a_pm_0 = world.gravity_accel_on(sat_id);
        vec3d a_pm_stage = world.gravity_accel_on(sat_id, sat->x_tr);
        earth->x_att.q = q_spin;
        vec3d a_pm_1_spin = world.gravity_accel_on(sat_id);
        vec3d r_pm_1_spin = sat->x_tr.r;
        world.restore_checkpoint_state(world_snapshot);
        earth->attitude_model = CelestialAttitudeModel::simple_spin;
        prop();
        vec3d r_pm_1_prop = sat->x_tr.r;
        vec3d v_pm_1_prop = sat->x_tr.v;
        world.restore_checkpoint_state(world_snapshot);
        earth->attitude_model = CelestialAttitudeModel::fixed;
        prop();
        vec3d r_pm_1_prop_fixed = sat->x_tr.r;
        vec3d v_pm_1_prop_fixed = sat->x_tr.v;

        // zonal (j2)
        world.restore_checkpoint_state(world_snapshot);
        earth->gravity_model = GravityModel::zonal;
        earth->degree = 2;
        earth->x_att.q = q_identity;
        vec3d a_zonal_0 = world.gravity_accel_on(sat_id);
        vec3d a_zonal_stage = world.gravity_accel_on(sat_id, sat->x_tr);
        earth->x_att.q = q_spin;
        vec3d a_zonal_1_spin = world.gravity_accel_on(sat_id);
        world.restore_checkpoint_state(world_snapshot);
        earth->attitude_model = CelestialAttitudeModel::simple_spin;
        prop();
        vec3d r_zonal_1_prop = sat->x_tr.r;
        vec3d v_zonal_1_prop = sat->x_tr.v;
        world.restore_checkpoint_state(world_snapshot);
        earth->attitude_model = CelestialAttitudeModel::fixed;
        prop();
        vec3d r_zonal_1_prop_fixed = sat->x_tr.r;
        vec3d v_zonal_1_prop_fixed = sat->x_tr.v;

        // spherical harmonics (2x2)
        world.restore_checkpoint_state(world_snapshot);
        earth->gravity_model = GravityModel::spherical_harmonics;
        earth->degree = 2;
        earth->order = 2;
        earth->C = matXd::Zero(earth->degree + 1, earth->order + 1);
        earth->S = matXd::Zero(earth->degree + 1, earth->order + 1);
        earth->C(2, 2) = 1.0e-6;
        earth->S(2, 1) = -5.0e-7;
        earth->x_att.q = q_identity;
        vec3d a_sphh_0 = world.gravity_accel_on(sat_id);
        vec3d a_sphh_stage = world.gravity_accel_on(sat_id, sat->x_tr);
        earth->x_att.q = q_spin;
        vec3d a_sphh_1_spin = world.gravity_accel_on(sat_id);
        world.restore_checkpoint_state(world_snapshot);
        earth->attitude_model = CelestialAttitudeModel::simple_spin;
        prop();
        vec3d r_sphh_1_prop = sat->x_tr.r;
        vec3d v_sphh_1_prop = sat->x_tr.v;
        world.restore_checkpoint_state(world_snapshot);
        earth->attitude_model = CelestialAttitudeModel::fixed;
        prop();
        vec3d r_sphh_1_prop_fixed = sat->x_tr.r;
        vec3d v_sphh_1_prop_fixed = sat->x_tr.v;

        // difference in acceleration due to spin
        std::println("Pointmass Spin Diff = {}", (a_pm_1_spin - a_pm_0).norm());
        std::println("Zonal Spin Diff = {}", (a_zonal_1_spin - a_zonal_0).norm());
        std::println(
            "Spherical Harmonics Spin Diff = {}",
            (a_sphh_1_spin - a_sphh_0).norm()
        );
        std::println("Pointmass Accel = {}", a_pm_0);
        std::println("Zonal Accel = {}", a_zonal_0);
        std::println("Spherical Harmonics Accel = {}", a_sphh_0);
        std::println("Spherical Harmonics Rotated Accel = {}", a_sphh_1_spin);
        std::println(
            "Pointmass Same-State Staging Diff = {}",
            (a_pm_stage - a_pm_0).norm()
        );
        std::println(
            "Zonal Same-State Staging Diff = {}",
            (a_zonal_stage - a_zonal_0).norm()
        );
        std::println(
            "Spherical Harmonics Same-State Staging Diff = {}",
            (a_sphh_stage - a_sphh_0).norm()
        );

        std::println();

        // difference in state simple spin vs fixed
        std::println(
            "Pointmass Position Error = {}",
            (r_pm_1_prop_fixed - r_pm_1_prop).norm()
        );
        std::println(
            "Pointmass Velocity Error = {}",
            (v_pm_1_prop_fixed - v_pm_1_prop).norm()
        );
        std::println(
            "Zonal Position Error = {}",
            (r_zonal_1_prop_fixed - r_zonal_1_prop).norm()
        );
        std::println(
            "Zonal Velocity Error = {}",
            (v_zonal_1_prop_fixed - v_zonal_1_prop).norm()
        );
        std::println(
            "Spherical Harmonics Position Error = {}",
            (r_sphh_1_prop_fixed - r_sphh_1_prop).norm()
        );
        std::println(
            "Spherical Harmonics Velocity Error = {}",
            (v_sphh_1_prop_fixed - v_sphh_1_prop).norm()
        );
    }
}

void run_moving_source_world_diag() {
    World world;

    // Earth
    EntityId earth_id = wgs84(world);
    Celestial* earth = world.celestial(earth_id);
    earth->gravity_model = GravityModel::zonal;
    earth->name = "Earth";
    earth->attitude_model = CelestialAttitudeModel::fixed;
    earth->propagate_att = false;
    earth->degree = 4;

    // Urath
    EntityId urath_id = wgs84(world);
    Celestial* urath = world.celestial(urath_id);
    urath->gravity_model = GravityModel::pointmass;
    urath->name = "Urath";
    urath->attitude_model = CelestialAttitudeModel::fixed;
    urath->propagate_att = false;
    urath->degree = 4;
    // urath->mu /= 10.0;

    // satellite
    EntityId sat_id = world.spawn_satellite();
    Satellite* sat = world.satellite(sat_id);
    sat->propagate_att = false;

    // binary body orbit initial conditions
    f64 separation = 100'000.0;
    f64 r_mag = separation / 2.0;
    f64 v_mag = std::sqrt(earth->mu / (2.0 * separation));
    f64 period = (2.0 * pi * r_mag) / v_mag;
    earth->x_tr.r = vec3d{r_mag, 0.0, 0.0};
    earth->x_tr.v = vec3d{0.0, -v_mag, 0.0};
    urath->x_tr.r = vec3d{-r_mag, 0.0, 0.0};
    urath->x_tr.v = vec3d{0.0, v_mag, 0.0};
    // sat orbit around urath
    f64 r_mag_sat = urath->semimajor_axis + 1000.0;
    f64 v_mag_sat = std::sqrt(urath->mu / r_mag_sat);
    sat->x_tr.r = urath->x_tr.r + vec3d{r_mag_sat, 0.0, 0.0};
    sat->x_tr.v = urath->x_tr.v + vec3d{0.0, v_mag_sat, 0.0};

    // integration options
    WorldStepperConfig cfg;
    cfg.step_tr = true;
    cfg.step_att = true;
    cfg.substeps = 10;
    cfg.ticks = 1;
    cfg.dt_scale = 1.0;
    cfg.integrator_tr = IntegratorTypeFixed::rk4;
    cfg.integrator_att = IntegratorTypeFixed::rk4;

    f64 t_span = period;
    f64 t0 = 0.0;
    i32 n_steps = 100000;
    f64 dt = t_span / n_steps;
    world.reset_time(t0);

    WorldStateSnapshot world_snapshot = world.capture_checkpoint();

    bool output_csv = false;
    std::ofstream file(std::string(PROJECT_ROOT) + "/assets/world_traj.csv");

    if (output_csv) {
        file << "i,"
             << "earth_x,earth_y,earth_z,"
             << "urath_x,urath_y,urath_z,"
             << "sat_x,sat_y,sat_z\n";
    }

    StateTr x_earth0 = earth->x_tr;
    StateTr x_urath0 = urath->x_tr;
    StateTr x_sat0 = sat->x_tr;

    i32 print_i = 25000;
    i32 csv_i = 10;
    WorldStepperStats stats;
    StatusCode step_status = StatusCode::ok;
    auto prop = [&](bool write_csv, bool print_progress) {
        stats = WorldStepperStats{};
        step_status = StatusCode::ok;
        for (i32 i = 0; i < n_steps; ++i) {
            if (write_csv && output_csv && i % csv_i == 0) {
                file << i << "," << earth->x_tr.r(0) << "," << earth->x_tr.r(1) << ","
                     << earth->x_tr.r(2) << "," << urath->x_tr.r(0) << ","
                     << urath->x_tr.r(1) << "," << urath->x_tr.r(2) << ","
                     << sat->x_tr.r(0) << "," << sat->x_tr.r(1) << "," << sat->x_tr.r(2)
                     << "\n";
            }
            if (print_progress && i % print_i == 0) {
                std::println("i: {}", i);
                std::println("Earth Position = {}", earth->x_tr.r);
                std::println("Urath Position = {}", urath->x_tr.r);
                std::println("Satellite Position = {}", sat->x_tr.r);
            }
            WorldStepResult step_result = step_world(world, dt, cfg);
            stats += step_result.stats;
            step_status = step_result.status;
            if (step_status != StatusCode::ok) break;
        }
    };

    // moving earth and urath
    earth->propagate_tr = true;
    urath->propagate_tr = true;
    sat->propagate_tr = true;
    prop(false, false);
    StateTr x_earth_moving = earth->x_tr;
    StateTr x_urath_moving = urath->x_tr;
    StateTr x_sat_moving = sat->x_tr;
    f64 t_moving = world.t_sim();
    WorldStepperStats stats_moving = stats;
    StatusCode status_moving = step_status;

    // fixed earth and urath
    world.restore_checkpoint_state(world_snapshot);
    earth->propagate_tr = false;
    urath->propagate_tr = false;
    sat->propagate_tr = true;
    prop(false, false);
    StateTr x_earth_fixed = earth->x_tr;
    StateTr x_urath_fixed = urath->x_tr;
    StateTr x_sat_fixed = sat->x_tr;
    f64 t_fixed = world.t_sim();
    WorldStepperStats stats_fixed = stats;
    StatusCode status_fixed = step_status;

    std::println("Moving Source Success = {}", status_moving == StatusCode::ok);
    std::println("Fixed Source Success = {}", status_fixed == StatusCode::ok);
    std::println("Moving Source Final Time = {}", t_moving);
    std::println("Fixed Source Final Time = {}", t_fixed);
    std::println("Earth Displacement = {}", (x_earth_moving.r - x_earth0.r).norm());
    std::println("Urath Displacement = {}", (x_urath_moving.r - x_urath0.r).norm());
    std::println("Satellite Displacement = {}", (x_sat_moving.r - x_sat0.r).norm());
    std::println(
        "Fixed Source Earth Displacement = {}",
        (x_earth_fixed.r - x_earth0.r).norm()
    );
    std::println(
        "Fixed Source Urath Displacement = {}",
        (x_urath_fixed.r - x_urath0.r).norm()
    );
    std::println(
        "Moving Source Satellite Position Diff = {}",
        (x_sat_moving.r - x_sat_fixed.r).norm()
    );
    std::println(
        "Moving Source Satellite Velocity Diff = {}",
        (x_sat_moving.v - x_sat_fixed.v).norm()
    );
    std::println("Moving Source Final Earth Position = {}", x_earth_moving.r);
    std::println("Moving Source Final Urath Position = {}", x_urath_moving.r);
    std::println("Moving Source Final Satellite Position = {}", x_sat_moving.r);
    std::println("Fixed Source Final Earth Position = {}", x_earth_fixed.r);
    std::println("Fixed Source Final Urath Position = {}", x_urath_fixed.r);
    std::println("Fixed Source Final Satellite Position = {}", x_sat_fixed.r);
}

void print_tle_data_summary(const std::string& label, const TLEData& tle) {
    std::println(
        "{}: sat_id={}, sat_num={}, name='{}', converted={}, units={}, sma={}, "
        "n_rad_s={}, jd={}",
        label,
        tle.sat_id,
        tle.sat_num,
        tle.name,
        tle.converted,
        i32(tle.units_angle),
        tle.sma,
        tle.n_rad_s,
        jd_to_scalar(tle.jd_utc)
    );
}

void run_tle_status_reader_diag() {
    print_diag_title("TLE Status Reader");

    Celestial earth = wgs84();
    TLEReadOptions raw_opts{.millennium = 2000};
    TLEReadOptions conv_opts{
        .millennium = 2000,
        .convert = true,
        .mu = earth.mu,
        .angle_out = UAngle::radian
    };

    std::string tle_single_file = pwd + "/assets/tle_iss.txt";
    TLEData tle_single;
    TLEStatus single_status
        = read_tle_data_single(tle_single_file, tle_single, conv_opts);
    std::println("Single Status: {}", tle_status_string(single_status));
    if (single_status == TLEStatus::ok) {
        print_tle_data_summary("Single", tle_single);
    }

    std::println();

    std::string tle_count_file = pwd + "/assets/tle_all.txt";
    svec<TLEData> tles_count;
    i32 count = 3;
    TLEStatus count_status
        = read_tle_data_count(tle_count_file, tles_count, count, conv_opts);
    std::println(
        "Count Status: {} --- count = {}",
        tle_status_string(count_status),
        tles_count.size()
    );
    if (count_status == TLEStatus::ok) {
        for (i32 i = 0; i < tles_count.size(); ++i) {
            print_tle_data_summary("Count[" + std::to_string(i) + "]", tles_count[i]);
        }
    }

    std::println();

    std::string tle_all_file = pwd + "/assets/tle_trunc.txt";
    svec<TLEData> tles_all;
    TLEStatus all_status = read_tle_data_all(tle_all_file, tles_all, raw_opts);
    std::println(
        "All Small Status: {} --- count = {}",
        tle_status_string(all_status),
        tles_all.size()
    );
    if (all_status == TLEStatus::ok && !tles_all.empty()) {
        print_tle_data_summary("All[0]", tles_all.front());
        print_tle_data_summary("All[last]", tles_all.back());
    }

    std::println();

    std::string tle_idx_file = pwd + "/assets/tle_all.txt";
    svec<TLEData> tles_idx;
    svec<i32> idx = {0, 2, 4};
    TLEStatus idx_status = read_tle_data_index(tle_idx_file, tles_idx, idx, conv_opts);
    std::println(
        "Index Status: {} --- count = {}",
        tle_status_string(idx_status),
        tles_idx.size()
    );
    if (idx_status == TLEStatus::ok) {
        for (i32 i = 0; i < tles_idx.size(); ++i) {
            print_tle_data_summary("Index[" + std::to_string(i) + "]", tles_idx[i]);
        }
    }

    std::println();

    std::string tle_satnum_file = pwd + "/assets/tle_all.txt";
    TLEData tle_satnum;
    i32 satnum = 25544;
    TLEStatus satnum_status
        = read_tle_data_single_satnum(tle_satnum_file, tle_satnum, satnum, conv_opts);
    std::println("Single Satnum Status: {}", tle_status_string(satnum_status));
    if (satnum_status == TLEStatus::ok) {
        print_tle_data_summary("Satnum", tle_satnum);
    }

    svec<TLEData> tles_satnums;
    svec<i32> satnums = {5, 11, 25544};
    TLEStatus satnums_status
        = read_tle_data_satnums(tle_satnum_file, tles_satnums, satnums, conv_opts);
    std::println(
        "Satnums Status: {} --- count = {}",
        tle_status_string(satnums_status),
        tles_satnums.size()
    );
    if (satnums_status == TLEStatus::ok) {
        for (i32 i = 0; i < tles_satnums.size(); ++i) {
            print_tle_data_summary("Satnums[" + std::to_string(i) + "]", tles_satnums[i]);
        }
    }

    std::println();

    TLEData tle_satid;
    std::string satid = "T0058";
    TLEStatus satid_status
        = read_tle_data_single_satid(tle_satnum_file, tle_satid, satid, conv_opts);
    std::println("Single Satid Status: {}", tle_status_string(satid_status));
    if (satid_status == TLEStatus::ok) {
        print_tle_data_summary("Satid", tle_satid);
    }

    svec<TLEData> tles_satids;
    svec<std::string> satids = {"T0058", "T0059"};
    TLEStatus satids_status
        = read_tle_data_satids(tle_satnum_file, tles_satids, satids, conv_opts);
    std::println(
        "Satids Status: {} --- count = {}",
        tle_status_string(satids_status),
        tles_satids.size()
    );
    if (satids_status == TLEStatus::ok) {
        for (i32 i = 0; i < tles_satids.size(); ++i) {
            print_tle_data_summary("Satids[" + std::to_string(i) + "]", tles_satids[i]);
        }
    }

    std::println();

    TLEData tle_missing;
    i32 satnum_missing = 99999999;
    TLEStatus missing_status = read_tle_data_single_satnum(
        tle_satnum_file,
        tle_missing,
        satnum_missing,
        raw_opts
    );
    std::println(
        "Missing Satnum Status: {} (should fail)",
        tle_status_string(missing_status)
    );

    std::println();

    std::string tle_fail_file = pwd + "/assets/tle_iss_fail.txt";
    TLEData tle_bad;
    TLEStatus bad_status
        = read_tle_data_single(pwd + "/assets/tle_iss_fail.txt", tle_bad, raw_opts);
    std::println("Bad Checksum Status: {} (should fail)", tle_status_string(bad_status));
}

void run_staged_attitude_gravity_diag() {
    print_diag_title("Staged Attitude Diagnostic");

    World world;

    // Earth
    EntityId earth_id = wgs84(world);
    Celestial* earth = world.celestial(earth_id);
    earth->name = "Earth";
    earth->gravity_model = GravityModel::spherical_harmonics;
    earth->degree = 4;
    earth->order = 4;
    bool gfc_ok = read_gfc(
        pwd + "/assets/earth/EGM2008.gfc.txt",
        earth->C,
        earth->S,
        earth->degree,
        earth->order
    );
    earth->propagate_tr = false;
    earth->propagate_att = false;
    earth->attitude_model = CelestialAttitudeModel::simple_spin;
    earth->x_att.q = dcm_to_ep(rotX(23.44, UAngle::degree));
    earth->set_spin_rate(earth->spin_rate() * 4);

    // satellite
    EntityId sat_id = world.spawn_satellite();
    Satellite* sat = world.satellite(sat_id);
    sat->propagate_tr = true;
    sat->propagate_att = false;
    sat->x_tr
        = classical_to_rv(8000, 0.15, 10.0, 10.0, 10.0, 0, earth->mu, UAngle::degree);
    sat->set_I(vec3d{100.0, 200.0, 300.0}.asDiagonal());
    sat->x_att.w = vec3d{0.000001, 0.0025, 0.000001}; // intermediate axis rotation

    // integration options
    WorldStepperConfig cfg;
    cfg.step_tr = true;
    cfg.step_att = true;
    cfg.substeps = 1;
    cfg.ticks = 1;
    cfg.dt_scale = 1.0;
    cfg.integrator_tr = IntegratorTypeFixed::rk4;
    cfg.integrator_att = IntegratorTypeFixed::rk4;

    f64 t_span = 10000;
    f64 t0 = 0.0;
    i32 n_steps = 10000;
    f64 dt = t_span / n_steps;
    world.reset_time(t0);

    WorldStateSnapshot world_snapshot = world.capture_checkpoint();

    WorldStepperStats stats;
    StatusCode step_status = StatusCode::ok;
    auto prop = [&]() {
        stats = WorldStepperStats{};
        step_status = StatusCode::ok;
        for (i32 i = 0; i < n_steps; ++i) {
            WorldStepResult step_result = step_world(world, dt, cfg);
            stats += step_result.stats;
            step_status = step_result.status;
            if (step_status != StatusCode::ok) break;
        }
    };

    auto run_case = [&](GravityModel gravity_model, std::string name) {
        world.restore_checkpoint_state(world_snapshot);
        earth->gravity_model = gravity_model;
        earth->attitude_model = CelestialAttitudeModel::fixed;
        earth->propagate_att = false;
        prop();
        StateTr xf_fixed = sat->x_tr;
        WorldStepperStats stats_fixed = stats;
        StatusCode status_fixed = step_status;

        world.restore_checkpoint_state(world_snapshot);
        earth->gravity_model = gravity_model;
        earth->attitude_model = CelestialAttitudeModel::simple_spin;
        earth->propagate_att = true;
        prop();
        StateTr xf_rot = sat->x_tr;
        WorldStepperStats stats_rot = stats;
        StatusCode status_rot = step_status;

        StateTr x_err = xf_fixed - xf_rot;

        std::println("{} Fixed Success = {}", name, status_fixed == StatusCode::ok);
        std::println("{} Rotating Success = {}", name, status_rot == StatusCode::ok);
        std::println("{} Position Difference = {}", name, x_err.r.norm());
        std::println("{} Velocity Difference = {}", name, x_err.v.norm());
        std::println();
    };

    run_case(GravityModel::pointmass, "Pointmass");
    run_case(GravityModel::spherical_harmonics, "Spherical Harmonics");
    run_case(GravityModel::zonal, "Zonal");
}

void run_render_pipeline_diag() {
    print_diag_title("Render Pipeline Diag");
    World world;

    // Earth
    EntityId earth_id = wgs84(world);
    Celestial* earth = world.celestial(earth_id);
    earth->name = "Earth";
    earth->gravity_model = GravityModel::zonal;
    earth->degree = 4;
    earth->order = 4;
    bool gfc_ok = read_gfc(
        pwd + "/assets/earth/EGM2008.gfc.txt",
        earth->C,
        earth->S,
        earth->degree,
        earth->order
    );
    if (!gfc_ok) {
        std::println("GFC Read Failed");
        return;
    }
    earth->propagate_tr = false;
    earth->propagate_att = true;
    earth->attitude_model = CelestialAttitudeModel::simple_spin;
    earth->x_att.q = dcm_to_ep(rotX(23.44, UAngle::degree));
    earth->set_spin_rate(earth->spin_rate() * 4);
    BuiltinRenderAssets assets;

    // satellite
    EntityId sat_id = world.spawn_satellite();
    Satellite* sat = world.satellite(sat_id);
    sat->propagate_tr = true;
    sat->propagate_att = true;
    sat->x_tr
        = classical_to_rv(8000, 0.15, 10.0, 10.0, 10.0, 0, earth->mu, UAngle::degree);
    sat->set_I(vec3d{100.0, 200.0, 300.0}.asDiagonal());
    sat->x_att.w = vec3d{0.000001, 0.0025, 0.000001}; // intermediate axis rotation

    // tle satellites
    bool load_tle_sats = true;
    if (load_tle_sats) {
        TLEReadOptions tle_opts{
            .convert = true,
            .millennium = 2000,
            .angle_out = UAngle::radian,
            .mu = earth->mu
        };
        svec<TLEData> tles;
        svec<uptr<Satellite>> tle_sats;
        // svec<i32> sat_nums = {25544, 11, 5};
        // TLEStatus tle_status = read_tle_data_satnums(
        //     pwd + "/assets/tle_all.txt",
        //     tles,
        //     sat_nums,
        //     tle_opts
        // );
        i32 num_sats = 10;
        TLEStatus tle_status
            = read_tle_data_count(pwd + "/assets/tle_all.txt", tles, num_sats, tle_opts);
        if (tle_status != TLEStatus::ok) {
            std::println("TLE Read Failure: {}", tle_status_string(tle_status));
            return;
        }
        tle_status = sats_from_tle_data(tle_sats, tles, tle_opts);
        if (tle_status != TLEStatus::ok) {
            std::println("TLE Read Failure: {}", tle_status_string(tle_status));
            return;
        }
        for (auto& tle_sat : tle_sats) {
            tle_sat->propagate_tr = true;
            tle_sat->propagate_att = false;
        }
        world.insert_satellites(std::move(tle_sats));
        // NOTE: these are in TEME and not rotated into sim inertial, just for visual
    }

    // station
    EntityId stat1_id = world.spawn_station();
    Station* stat1 = world.station(stat1_id);
    world.set_stat_anchor_detic(stat1_id, earth_id, vec3d{0.0, 0.0, 0.0});

    EntityId stat2_id = world.spawn_station();
    Station* stat2 = world.station(stat2_id);
    world.set_stat_anchor_detic(stat2_id, earth_id, vec3d{45.0, 45.0, 0.0});

    EntityId stat3_id = world.spawn_station();
    Station* stat3 = world.station(stat3_id);
    world.set_stat_anchor_detic(stat3_id, earth_id, vec3d{90.0, 90.0, 0.0});

    // integrator
    WorldStepperConfig cfg;
    cfg.step_tr = true;
    cfg.step_att = true;
    cfg.substeps = 1;
    cfg.ticks = 100;
    cfg.dt_scale = 2.0;
    cfg.integrator_tr = IntegratorTypeFixed::rk4;
    cfg.integrator_att = IntegratorTypeFixed::rk4;

    f64 t0 = 0.0;
    f64 dt = 1.0 / cfg.ticks * cfg.dt_scale;
    world.reset_time(t0);
    WorldStepperStats stats;
    StatusCode step_status = StatusCode::ok;
    WorldStepperWorkspace wksp;

    // windowing and graphics
    int screenWidth = 960;
    int screenHeight = 800;

    SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_MSAA_4X_HINT | FLAG_WINDOW_HIGHDPI);

    InitWindow(screenWidth, screenHeight, "Basic Window");

    Camera3D camera;
    vec3f cam_pos = vec3f{1.0, 1.0, 1.0} * 25000;
    camera.position = eig_to_rl(cam_pos);
    camera.fovy = 45;
    camera.projection = CAMERA_PERSPECTIVE;
    camera.up = eig_to_rl(axis_z);
    camera.target = eig_to_rl(origin);
    rlSetClipPlanes(1.0e3, 1.0e6);

    i32 sphere_i = 32;
    Image checker_pattern = GenImageChecked(sphere_i, sphere_i, 1, 1, RAYWHITE, GRAY);
    Texture2D texture = LoadTextureFromImage(checker_pattern);

    Mesh sphere_mesh = GenMeshSphere(1., sphere_i, sphere_i);
    Model sphere_model = LoadModelFromMesh(sphere_mesh);
    sphere_model.transform = MatrixIdentity();
    sphere_model.materials[0] = LoadMaterialDefault();
    sphere_model.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = texture;

    checker_pattern = GenImageChecked(2, 2, 1, 1, RAYWHITE, GRAY);
    texture = LoadTextureFromImage(checker_pattern);
    Mesh cube_mesh = GenMeshCube(1.0f, 1.0f, 1.0f);
    Model cube_model = LoadModelFromMesh(cube_mesh);
    cube_model.transform = MatrixIdentity();
    cube_model.materials[0] = LoadMaterialDefault();
    cube_model.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = texture;

    Mesh cylinder_mesh = GenMeshCylinder(1.0f, 1.0f, sphere_i);
    Model cylinder_model = LoadModelFromMesh(cylinder_mesh);
    cylinder_model.transform = MatrixIdentity();
    cylinder_model.materials[0] = LoadMaterialDefault();
    cylinder_model.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = texture;

    UnloadImage(checker_pattern);

    SetTargetFPS(60);
    bool paused = false;
    while (!WindowShouldClose()) {
        // std::cout << 1.0f/GetFrameTime() << '\n';
        BeginDrawing();
        BeginMode3D(camera);
        ClearBackground(Color({30, 30, 30, 255}));

        RenderSceneSnapshot scene = build_render_scene_snapshot(world, assets);
        render_scene_snapshot(scene, sphere_model, cube_model, cylinder_model);

        if (IsKeyPressed(KEY_SPACE)) {
            paused = !paused;
        }
        if (!paused) {
            WorldStepResult step_result = step_world(world, dt, cfg, wksp);
            stats += step_result.stats;
            step_status = step_result.status;
            if (step_status != StatusCode::ok) {
                std::println("Simulation Failure");
                break;
            }
        }
        svec<EntityId> stat_ids = world.active_station_ids();
        for (EntityId stat_id : stat_ids) {
            draw_axes(
                world.stat_x_tr_inertial(stat_id),
                world.stat_x_att_inertial(stat_id),
                1000.0f
            );
        }

        EndMode3D();
        EndDrawing();
    }
    print_diag_title();
}

void run_world_workspace_diag() {
    print_diag_title("World Stepper Workspace");

    World world;

    // Earth
    EntityId earth_id = wgs84(world);
    Celestial* earth = world.celestial(earth_id);
    earth->name = "Earth";
    earth->gravity_model = GravityModel::spherical_harmonics;
    earth->degree = 4;
    earth->order = 4;
    bool gfc_ok = read_gfc(
        pwd + "/assets/earth/EGM2008.gfc.txt",
        earth->C,
        earth->S,
        earth->degree,
        earth->order
    );
    if (!gfc_ok) {
        std::println("GFC Load Failed");
        return;
    }
    earth->propagate_tr = false;
    earth->propagate_att = true;
    earth->attitude_model = CelestialAttitudeModel::simple_spin;
    earth->x_att.q = dcm_to_ep(rotX(23.44, UAngle::degree));
    earth->set_spin_rate(earth->spin_rate() * 4);

    // satellite
    EntityId sat_id = world.spawn_satellite();
    Satellite* sat = world.satellite(sat_id);
    sat->propagate_tr = true;
    sat->propagate_att = true;
    sat->x_tr
        = classical_to_rv(8000, 0.15, 10.0, 10.0, 10.0, 0, earth->mu, UAngle::degree);
    sat->set_I(vec3d{100.0, 200.0, 300.0}.asDiagonal());
    sat->x_att.w = vec3d{0.000001, 0.0025, 0.000001}; // intermediate axis rotation

    // station
    EntityId stat1_id = world.spawn_station();
    world.set_stat_anchor_detic(stat1_id, earth_id, vec3d{0.0, 0.0, 0.0});

    // integrator
    WorldStepperConfig cfg;
    cfg.step_tr = true;
    cfg.step_att = true;
    cfg.substeps = 1;
    cfg.ticks = 1;
    cfg.dt_scale = 1.0;
    cfg.integrator_tr = IntegratorTypeFixed::rk4;
    cfg.integrator_att = IntegratorTypeFixed::rk4;

    f64 t_span = 10000;
    f64 t0 = 0.0;
    i32 n_steps = 10000;
    f64 dt = t_span / n_steps;
    world.reset_time(t0);

    WorldStepperStats stats;
    StatusCode step_status = StatusCode::ok;

    WorldStepperWorkspace wksp;

    auto print_workspace = [&](const std::string& label) {
        std::println("{} Propagated TR IDs = {}", label, wksp.propagated_tr_ids.size());
        std::println("{} Propagated ATT IDs = {}", label, wksp.propagated_att_ids.size());
        std::println("{} Celestial ATT IDs = {}", label, wksp.celestial_att_ids.size());
        std::println("{} Gravity Source IDs = {}", label, wksp.gravity_source_ids.size());
        std::println("{} Source ATT IDs = {}", label, wksp.source_att_ids.size());
        std::println("{} Dirty = {}", label, wksp.dirty);
    };

    rebuild_world_stepper_workspace(world, wksp);
    print_workspace("WORLD WKSP Initial");

    auto contains_id = [](const svec<EntityId>& ids, EntityId id) {
        for (EntityId current_id : ids) {
            if (current_id == id) return true;
        }
        return false;
    };

    i32 workspace_test_steps = 0;

    // adding a propagated body requires the cached ID sets to be rebuilt
    EntityId added_sat_id = world.spawn_satellite();
    Satellite* added_sat = world.satellite(added_sat_id);
    if (added_sat == nullptr) {
        std::println("Added Satellite Lookup Failed");
        return;
    }
    added_sat->propagate_tr = true;
    added_sat->propagate_att = false;
    added_sat->x_tr
        = classical_to_rv(9000, 0.05, 20.0, 0.0, 0.0, 0.0, earth->mu, UAngle::degree);

    invalidate_stepper_wksp(wksp);
    bool add_marked_dirty = wksp.dirty;
    WorldStepResult add_step_result = step_world(world, dt, cfg, wksp);
    stats += add_step_result.stats;
    step_status = add_step_result.status;
    ++workspace_test_steps;
    bool add_rebuilt = !wksp.dirty
        && contains_id(wksp.propagated_tr_ids, added_sat_id);

    if (step_status != StatusCode::ok) {
        std::println("Added Body Step Status = {}", status_string(step_status));
        return;
    }

    // deactivating the body removes it from the propagated ID cache
    if (!world.make_inactive(added_sat_id)) {
        std::println("Added Satellite Deactivation Failed");
        return;
    }
    invalidate_stepper_wksp(wksp);
    bool inactive_marked_dirty = wksp.dirty;
    WorldStepResult inactive_step_result = step_world(world, dt, cfg, wksp);
    stats += inactive_step_result.stats;
    step_status = inactive_step_result.status;
    ++workspace_test_steps;
    bool inactive_rebuilt = !wksp.dirty
        && !contains_id(wksp.propagated_tr_ids, added_sat_id);

    if (step_status != StatusCode::ok) {
        std::println("Inactive Body Step Status = {}", status_string(step_status));
        return;
    }

    // state values do not change workspace membership
    sat->x_tr.r(0) += 1.0;
    bool state_change_kept_clean = !wksp.dirty;

    std::println("WORLD WKSP Invalidation Checks");
    std::println("    Add Marked Dirty = {}", add_marked_dirty);
    std::println("    Added ID Present After Rebuild = {}", add_rebuilt);
    std::println("    Inactive Marked Dirty = {}", inactive_marked_dirty);
    std::println("    Inactive ID Removed After Rebuild = {}", inactive_rebuilt);
    std::println("    State Change Kept Cache Clean = {}", state_change_kept_clean);
    std::println(
        "    Passed = {}",
        add_marked_dirty && add_rebuilt && inactive_marked_dirty
            && inactive_rebuilt && state_change_kept_clean
    );

    for (i32 i = 0; i < n_steps; ++i) {
        WorldStepResult step_result = step_world(world, dt, cfg, wksp);
        stats += step_result.stats;
        step_status = step_result.status;
        if (step_status != StatusCode::ok) {
            std::println("Integration Failed");
            break;
        }
    }

    print_workspace("WORLD WKSP Final");
    std::println("WORLD WKSP Success = {}", step_status == StatusCode::ok);
    std::println("WORLD WKSP Ticks Completed = {}", stats.ticks_completed);
    std::println("WORLD WKSP Substeps Completed = {}", stats.substeps_completed);
    std::println("WORLD WKSP Time Advanced = {}", stats.dt_sim_advanced);
    std::println("WORLD WKSP Final Time = {}", world.t_sim());
    std::println("WORLD WKSP Expected Time = {}", t0 + t_span + workspace_test_steps * dt);
}

void run_world_measurement_context_diag() {
    World world;

    // Earth
    EntityId earth_id = wgs84(world);
    Celestial* earth = world.celestial(earth_id);
    earth->name = "Earth";

    // satellite
    EntityId sat_id = world.spawn_satellite();
    Satellite* sat = world.satellite(sat_id);
    sat->x_tr
        = classical_to_rv(8000, 0.15, 10.0, 10.0, 10.0, 0, earth->mu, UAngle::degree);

    // station
    EntityId stat_id = world.spawn_station();
    world.set_stat_anchor_detic(stat_id, earth_id, vec3d{0.0, 0.0, 0.0});

    print_diag_title("World Measurement Context Diagnostic");

    struct MeasurementContextDiagCase {
        const char* name = "";
        ObservationType type = ObservationType::radec;
        UAngle angle_out = UAngle::radian;
    };

    svec<MeasurementContextDiagCase> cases{
        {.name = "RA/Dec", .type = ObservationType::radec, .angle_out = UAngle::degree},
        {.name = "Range", .type = ObservationType::range, .angle_out = UAngle::radian},
        {.name = "Range Rate",
         .type = ObservationType::range_rate,
         .angle_out = UAngle::radian},
        {.name = "Az/El", .type = ObservationType::azel, .angle_out = UAngle::degree}
    };

    for (const auto& diag_case : cases) {
        vecXd z_world = world_predict_measurement(
            world,
            diag_case.type,
            stat_id,
            sat_id,
            UAngle::radian,
            diag_case.angle_out
        );

        MeasurementContext ctx;
        StatusCode ctx_status = make_world_station_measurement_context(
            world,
            ctx,
            stat_id,
            sat->x_tr,
            diag_case.type
        );

        vecXd z_ctx = predict_measurement(
            diag_case.type,
            ctx,
            UAngle::radian,
            diag_case.angle_out
        );

        std::println("{} Status: {}", diag_case.name, status_string(ctx_status));
        std::println("{} via World: {}", diag_case.name, z_world);
        std::println("{} via Context: {}", diag_case.name, z_ctx);
        std::println("{} Error: {}", diag_case.name, (z_world - z_ctx).norm());
    }
}
// -------------------------------------------------------------------------------------------------
void run_iod_lumve_ekf_init_diag() {
    World world;

    // Diagnostic parameters
    GravityModel truth_model = GravityModel::spherical_harmonics;
    i32 truth_deg_ord = 32;
    GravityModel est_model = GravityModel::zonal;
    i32 est_degree = 4;
    bool use_iod_initial_guess = true; // keep IOD optional for now

    // Earth
    EntityId earth_id = wgs84(world);
    Celestial* earth = world.celestial(earth_id);
    earth->name = "Earth";
    earth->gravity_model = truth_model;
    earth->degree = truth_deg_ord;
    earth->order = truth_deg_ord;
    bool gfc_ok = read_gfc(
        pwd + "/assets/earth/EGM2008.gfc.txt",
        earth->C,
        earth->S,
        earth->degree,
        earth->order
    );
    if (!gfc_ok) {
        std::println("GFC Load Failed");
        return;
    }
    earth->propagate_tr = true;
    earth->propagate_att = true;
    earth->attitude_model = CelestialAttitudeModel::simple_spin;
    earth->x_att.q = dcm_to_ep(rotX(23.44, UAngle::degree));
    earth->set_spin_rate(earth->spin_rate());

    // Stations
    vec3d llh1 = vec3d{0.0, 0.0, 0.0}; // [lat, lon, h] = [deg, deg, sim units]
    vec3d llh2 = vec3d{45.0, 30.0, 0.0};
    EntityId stat1_id = world.spawn_station();
    Station* stat1 = world.station(stat1_id);
    world.set_stat_anchor_detic(stat1_id, earth_id, llh1);

    EntityId stat2_id = world.spawn_station();
    Station* stat2 = world.station(stat2_id);
    world.set_stat_anchor_detic(stat2_id, earth_id, llh2);

    // Satellite
    EntityId sat_id = world.spawn_satellite();
    Satellite* sat = world.satellite(sat_id);
    StateTr x0_truth;
    x0_truth.r = vec3d{7000.0, 1000.0, 1300.0};
    x0_truth.v = vec3d{-0.5, 7.2, 1.0};
    sat->x_tr = x0_truth;

    // Orbit Determination Config
    ODDynamicsConfig dyn_config = make_od_cfg_from_celestial(*earth);
    dyn_config.tr_model = worldtrmodel_to_odtrmodel(
        est_model
    ); // have od model and "truth" model discrepancy
    dyn_config.zonal_degree = est_degree;
    dyn_config.update_body_attitude = true;

    f64 t_meas0 = 0.0;
    f64 t_measf = 600.0;
    i32 N_meas = 100 * 2; // half lumve, half ekf
    i32 split_idx = N_meas / 2;
    vecXd t_meas = vecXd::LinSpaced(N_meas, t_meas0, t_measf);
    f64 sigma_rad = 1e-5;
    f64 sigma_range = 1e-3;
    f64 sigma_range_rate = 1e-6;
    f64 sigma_r = 1e-3;
    f64 sigma_v = 1e-6;

    WorldStepperConfig cfg;
    cfg.step_tr = true;
    cfg.step_att = true;
    cfg.substeps = 1;
    cfg.ticks = 10;
    cfg.dt_scale = 1.0 / cfg.ticks;
    cfg.integrator_tr = IntegratorTypeFixed::rk4;
    cfg.integrator_att = IntegratorTypeFixed::rk4;
    // radec + range measurements only

    WorldStateSnapshot world_snapshot = world.capture_checkpoint();

    i32 n_steps = 100;
    f64 dt;

    WorldStepperStats stats;
    StatusCode step_status = StatusCode::ok;
    auto prop = [&]() {
        stats = WorldStepperStats{};
        step_status = StatusCode::ok;
        for (i32 i = 0; i < n_steps; ++i) {
            WorldStepResult step_result = step_world(world, dt, cfg);
            stats += step_result.stats;
            step_status = step_result.status;
            if (step_status != StatusCode::ok) break;
        }
    };

    // IOD
    std::array<f64, 3> ts;
    std::array<vec2d, 3> radecs;
    std::array<vec3d, 3> Rs;
    std::array<vec3d, 3> rs;

    f64 t0 = 0.0;
    ts[1] = t0;
    radecs[1]
        = world_predict_measurement(world, ObservationType::radec, stat1_id, sat_id);
    Rs[1] = world.stat_x_tr_inertial(stat1_id).r; // ignore checks for now
    rs[1] = sat->x_tr.r;

    f64 t_span = -30.0;
    dt = t_span / n_steps;
    prop();
    ts[0] = t0 + t_span;
    radecs[0]
        = world_predict_measurement(world, ObservationType::radec, stat1_id, sat_id);
    Rs[0] = world.stat_x_tr_inertial(stat1_id).r; // ignore checks for now
    rs[0] = sat->x_tr.r;

    world.restore_checkpoint_state(world_snapshot);

    t_span = 30.0;
    dt = t_span / n_steps;
    prop();
    ts[2] = t0 + t_span;
    radecs[2]
        = world_predict_measurement(world, ObservationType::radec, stat1_id, sat_id);
    Rs[2] = world.stat_x_tr_inertial(stat1_id).r; // ignore checks for now
    rs[2] = sat->x_tr.r;

    world.restore_checkpoint_state(world_snapshot);

    IODAnglesObs3 gauss_input = iod_angles3_from_radec(ts, radecs, Rs, UAngle::radian);
    IODResult iod_res = iod_gauss(gauss_input, earth->mu);
    // IODResult iod_res
    //     = iod_herrickgibbs(ts[0], ts[1], ts[2], rs[0], rs[1], rs[2], earth->mu);

    // LUMVE
    ODBatchInput lumve_input;
    bool iod_guess_ok = iod_res.success && statetr_to_vec6d(iod_res.x).allFinite();
    if (use_iod_initial_guess && iod_guess_ok) {
        lumve_input.x0_guess = iod_res.x;
    } else {
        lumve_input.x0_guess = x0_truth;
        lumve_input.x0_guess.r += vec3d{1.0, -1.0, 0.5};
        lumve_input.x0_guess.v += vec3d{1e-4, -1e-4, 2.5e-4};
    }

    lumve_input.dyn_config = dyn_config;
    lumve_input.t0 = 0.0;
    lumve_input.max_iters = 10;
    lumve_input.prop_steps = 100;
    lumve_input.tol_dx = 1e-6;
    lumve_input.tol_residual = 1e-8;

    std::mt19937_64 rng(12345);
    std::normal_distribution<f64> noise_unit(0.0, 1.0);

    // create measurements
    for (i32 i = 0; i < split_idx; ++i) {
        dt = t_meas(i) - world.t_sim();
        WorldStepResult step_result = step_world(world, dt, cfg);
        stats += step_result.stats;
        step_status = step_result.status;
        if (step_status != StatusCode::ok) break;

        // stat1 -> radec + range
        StateTr x_tr_obsv1 = world.stat_x_tr_inertial(stat1_id);

        // stat2 -> range-rate
        StateTr x_tr_obsv2 = world.stat_x_tr_inertial(stat2_id);

        Measurement meas;
        meas.target_id = sat_id;
        meas.observer_id = stat1_id;

        // get measurement (RADec)
        meas.t = t_meas(i);
        meas.type = ObservationType::radec;
        meas.z = world_predict_measurement(
            world,
            meas.type,
            stat1_id,
            sat_id,
            UAngle::radian,
            UAngle::radian
        );
        meas.z(0) += noise_unit(rng) * sigma_rad;
        meas.z(1) += noise_unit(rng) * sigma_rad;
        meas.R = matXd1<2> * sigma_rad * sigma_rad;
        lumve_input.measurements.push_back(meas);
        lumve_input.observer_states.push_back(x_tr_obsv1);

        // get measurement (range)
        meas.t = t_meas(i);
        meas.type = ObservationType::range;
        meas.z = world_predict_measurement(
            world,
            meas.type,
            stat1_id,
            sat_id,
            UAngle::radian,
            UAngle::radian
        );
        meas.z(0) += noise_unit(rng) * sigma_range;
        meas.R = matXd1<1> * sigma_range * sigma_range;
        lumve_input.measurements.push_back(meas);
        lumve_input.observer_states.push_back(x_tr_obsv1);

        // get measurement (range-rate)
        meas.t = t_meas(i);
        meas.type = ObservationType::range_rate;
        meas.z = world_predict_measurement(
            world,
            meas.type,
            stat2_id,
            sat_id,
            UAngle::radian,
            UAngle::radian
        );
        meas.z(0) += noise_unit(rng) * sigma_range_rate;
        meas.R = matXd1<1> * sigma_range_rate * sigma_range_rate;
        lumve_input.measurements.push_back(meas);
        lumve_input.observer_states.push_back(x_tr_obsv2);
    }

    // solve LUMVE
    ODBatchResult lumve_result = od_batch_lumve(lumve_input);

    // offline EKF
    ODEKFOfflineInput ekf_input;
    bool lumve_state_ok = od_status_success(lumve_result.status)
                          && statetr_to_vec6d(lumve_result.x0_est).allFinite();
    ekf_input.initial_filter.x
        = lumve_state_ok ? lumve_result.x0_est : lumve_input.x0_guess;
    ekf_input.initial_filter.t = lumve_input.t0;

    bool lumve_cov_ok = lumve_result.covariance.allFinite()
                        && lumve_result.covariance.diagonal().minCoeff() > 0.0;
    ekf_input.initial_filter.P = lumve_cov_ok ? lumve_result.covariance : mat6d1;

    ekf_input.dyn_config = dyn_config;
    ekf_input.prop_steps = 200;
    ekf_input.Q = mat6d1 * 1e-4;

    for (i32 i = split_idx; i < N_meas; ++i) {
        dt = t_meas(i) - world.t_sim();
        WorldStepResult step_result = step_world(world, dt, cfg);
        stats += step_result.stats;
        step_status = step_result.status;
        if (step_status != StatusCode::ok) break;
        // stat1 -> radec + range
        StateTr x_tr_obsv1 = world.stat_x_tr_inertial(stat1_id);

        // stat2 -> range-rate
        StateTr x_tr_obsv2 = world.stat_x_tr_inertial(stat2_id);

        Measurement meas;
        meas.target_id = sat_id;
        meas.observer_id = stat1_id;

        // get measurement (RADec)
        meas.t = t_meas(i);
        meas.type = ObservationType::radec;
        meas.z = world_predict_measurement(
            world,
            meas.type,
            stat1_id,
            sat_id,
            UAngle::radian,
            UAngle::radian
        );
        meas.z(0) += noise_unit(rng) * sigma_rad;
        meas.z(1) += noise_unit(rng) * sigma_rad;
        meas.R = matXd1<2> * sigma_rad * sigma_rad;
        ekf_input.measurements.push_back(meas);
        ekf_input.observer_states.push_back(x_tr_obsv1);

        // get measurement (range)
        meas.t = t_meas(i);
        meas.type = ObservationType::range;
        meas.z = world_predict_measurement(
            world,
            meas.type,
            stat1_id,
            sat_id,
            UAngle::radian,
            UAngle::radian
        );
        meas.z(0) += noise_unit(rng) * sigma_range;
        meas.R = matXd1<1> * sigma_range * sigma_range;
        ekf_input.measurements.push_back(meas);
        ekf_input.observer_states.push_back(x_tr_obsv1);

        // get measurement (range-rate)
        meas.t = t_meas(i);
        meas.type = ObservationType::range_rate;
        meas.z = world_predict_measurement(
            world,
            meas.type,
            stat2_id,
            sat_id,
            UAngle::radian,
            UAngle::radian
        );
        meas.z(0) += noise_unit(rng) * sigma_range_rate;
        meas.R = matXd1<1> * sigma_range_rate * sigma_range_rate;
        ekf_input.measurements.push_back(meas);
        ekf_input.observer_states.push_back(x_tr_obsv2);
    }

    // solve offline ekf
    ODEKFResult ekf_result = od_ekf_offline(ekf_input);

    StateTr x_truth_final = propagate_tr_od(
        0.0,
        x0_truth,
        t_meas(N_meas - 1),
        ekf_input.prop_steps,
        ekf_input.dyn_config
    );

    f64 iod_err = iod_guess_ok ? statetr_to_vec6d(iod_res.x - x0_truth).norm() : 0.0;
    f64 lumve_initial_err = statetr_to_vec6d(lumve_input.x0_guess - x0_truth).norm();
    f64 lumve_final_err = statetr_to_vec6d(lumve_result.x0_est - x0_truth).norm();
    f64 ekf_initial_err = statetr_to_vec6d(ekf_input.initial_filter.x - x0_truth).norm();
    f64 ekf_final_err = statetr_to_vec6d(ekf_result.filter.x - x_truth_final).norm();
    f64 ekf_final_r_err = (ekf_result.filter.x.r - x_truth_final.r).norm();
    f64 ekf_final_v_err = (ekf_result.filter.x.v - x_truth_final.v).norm();

    print_diag_title("IOD -> LUMVE -> EKF Pipeline Diagnostic");
    std::println("Truth Model = {}", gravity_model_str(truth_model));
    std::println("Truth Degree/Order = {}", truth_deg_ord);
    std::println("Estimator Model = {}", gravity_model_str(est_model));
    std::println("Estimator Zonal Degree = {}", est_degree);
    std::println();
    std::println("IOD Used As Initial Guess = {}", use_iod_initial_guess && iod_guess_ok);
    std::println("IOD Success = {}", iod_res.success);
    std::println("IOD Status: {}", istatus_string(iod_res.status));
    std::println("IOD Initial Error = {}", iod_err);

    std::println();

    std::println("LUMVE Success = {}", od_status_success(lumve_result.status));
    std::println("LUMVE Status: {}", status_string(lumve_result.status));
    std::println("LUMVE Initial Error = {}", lumve_initial_err);
    std::println("LUMVE Final Error = {}", lumve_final_err);
    std::println("LUMVE Measurements = {}", lumve_input.measurements.size());
    std::println("LUMVE Iterations = {}", lumve_result.iterations);
    std::println("LUMVE Residual Norm = {}", lumve_result.residual_norm);
    std::println("LUMVE Raw Residual Norm = {}", lumve_result.raw_residual_norm);
    std::println("LUMVE Delta x Norm = {}", lumve_result.dx_norm);
    std::println("LUMVE Covariance Used = {}", lumve_cov_ok);
    std::println("LUMVE Covariance Norm = {}", lumve_result.covariance.norm());

    std::println();

    std::println("EKF Success = {}", od_status_success(ekf_result.status));
    std::println("EKF Status: {}", status_string(ekf_result.status));
    std::println("EKF Initial Error = {}", ekf_initial_err);
    std::println("EKF Final Error = {}", ekf_final_err);
    std::println("EKF Final Position Error = {}", ekf_final_r_err);
    std::println("EKF Final Velocity Error = {}", ekf_final_v_err);
    std::println("EKF Final Time = {}", ekf_result.filter.t);
    std::println("EKF Expected Final Time = {}", t_meas(N_meas - 1));
    std::println("EKF Processed Measurements = {}", ekf_result.processed_measurements);
    std::println("EKF Total Measurements = {}", ekf_input.measurements.size());
    std::println("EKF Residual Norm = {}", ekf_result.residual_norm);
    std::println("EKF Raw Residual Norm = {}", ekf_result.raw_residual_norm);
    std::println("EKF Final Covariance Norm = {}", ekf_result.filter.P.norm());

    std::println();
    std::println("World Steps: {}", stats.substeps_completed * stats.ticks_completed);
}

void run_od_zonal_jacobian_diag() {
    World world;

    EntityId earth_id = wgs84(world);
    Celestial* earth = world.celestial(earth_id);
    if (!earth) {
        std::println("OD Zonal Jacobian Diagnostic: Earth Spawn Failed");
        return;
    }

    earth->name = "Earth";
    earth->attitude_model = CelestialAttitudeModel::simple_spin;
    earth->x_att.q = dcm_to_ep(rotX(23.44, UAngle::degree));
    earth->set_spin_rate(earth->spin_rate());

    StateTr x0;
    x0.r = vec3d{7000.0, 1000.0, 1300.0};
    x0.v = vec3d{-0.5, 7.2, 1.0};

    auto print_case = [&](const std::string& name, f64 t, const ODDynamicsConfig& cfg) {
        mat6d G = jacobian_tr_od(t, x0, cfg);
        mat6d G_fd = jacobian_fd_od_dynamics(t, x0, cfg);

        mat6d G_err = G - G_fd;
        mat3d A_err = G.block<3, 3>(3, 0) - G_fd.block<3, 3>(3, 0); // accel block
        f64 G_fd_norm = G_fd.norm();
        f64 A_fd_norm = G_fd.block<3, 3>(3, 0).norm();
        if (G_fd_norm < 1.0) G_fd_norm = 1.0;
        if (A_fd_norm < 1.0) A_fd_norm = 1.0;

        // consider fd as "truth"
        std::println("{}", name);
        std::println("G error = {}", G_err.norm());
        std::println("G rel error = {}", G_err.norm() / G_fd_norm);
        std::println("A error = {}", A_err.norm());
        std::println("A rel error = {}", A_err.norm() / A_fd_norm);
        std::println();
    };

    print_diag_title("OD Zonal Jacobian Diagnostic");

    ODDynamicsConfig cfg;
    cfg.mu = earth->mu;
    cfg.R_cb_ref = earth->semimajor_axis;
    cfg.J = earth->J;
    cfg.q_cb0 = earth->x_att.q;
    cfg.w_cb = earth->x_att.w;
    cfg.t0 = 0.0;

    cfg.tr_model = ODTrDynamicsModel::two_body;
    print_case("Two Body", 0.0, cfg);

    for (i32 degree = 2; degree <= 6; ++degree) {
        cfg.tr_model = ODTrDynamicsModel::zonal;
        cfg.zonal_degree = degree;

        cfg.update_body_attitude = false;
        cfg.att_model = ODAnchorAttModel::fixed;
        cfg.q_cb0 = q_identity;
        print_case("J" + std::to_string(degree) + " Fixed Identity", 0.0, cfg);

        cfg.update_body_attitude = false;
        cfg.att_model = ODAnchorAttModel::fixed;
        cfg.q_cb0 = earth->x_att.q;
        print_case("J" + std::to_string(degree) + " Fixed Tilted", 0.0, cfg);

        cfg.update_body_attitude = true;
        cfg.att_model = ODAnchorAttModel::simple_spin;
        cfg.q_cb0 = earth->x_att.q;
        print_case("J" + std::to_string(degree) + " Simple Spin", 1000.0, cfg);
    }
}

void run_world_ekf_step_diag() {
    World world;

    // Diagnostic parameters
    GravityModel truth_model = GravityModel::zonal;
    i32 truth_deg_ord = 4;
    GravityModel est_model = GravityModel::pointmass;
    i32 est_degree = 4;

    // Earth
    EntityId earth_id = wgs84(world);
    Celestial* earth = world.celestial(earth_id);
    earth->name = "Earth";
    earth->gravity_model = truth_model;
    earth->degree = truth_deg_ord;
    earth->order = truth_deg_ord;
    bool gfc_ok = read_gfc(
        pwd + "/assets/earth/EGM2008.gfc.txt",
        earth->C,
        earth->S,
        earth->degree,
        earth->order
    );
    if (!gfc_ok) {
        std::println("GFC Load Failed");
        return;
    }
    earth->propagate_tr = true;
    earth->propagate_att = true;
    earth->attitude_model = CelestialAttitudeModel::simple_spin;
    earth->x_att.q = dcm_to_ep(rotX(23.44, UAngle::degree));
    earth->set_spin_rate(earth->spin_rate());

    // Stations
    vec3d llh1 = vec3d{0.0, 0.0, 0.0}; // [lat, lon, h] = [deg, deg, sim units]
    EntityId stat_id = world.spawn_station();
    world.set_stat_anchor_detic(stat_id, earth_id, llh1);

    // Satellite
    EntityId sat_id = world.spawn_satellite();
    Satellite* sat = world.satellite(sat_id);
    StateTr x0_truth;
    x0_truth.r = vec3d{7000.0, 1000.0, 1300.0};
    x0_truth.v = vec3d{-0.5, 7.2, 1.0};
    sat->x_tr = x0_truth;

    // Orbit Determination Config
    ODDynamicsConfig dyn_config = make_od_cfg_from_celestial(*earth);
    dyn_config.tr_model = worldtrmodel_to_odtrmodel(
        est_model
    ); // have od model and "truth" model discrepancy
    dyn_config.zonal_degree = est_degree;
    dyn_config.update_body_attitude = true;

    f64 t_meas = 0.0;
    f64 sigma_range = 1e-3;

    ODEKFState filter;
    filter.x = x0_truth + StateTr{.r = {1.0, 0.25, 0.1}, .v = {0.001, 0.0005, 0.0}};
    filter.P = mat6d1;
    filter.t = t_meas;

    ObservationType obsv_type = ObservationType::range;
    vecXd z = world_predict_measurement(world, obsv_type, stat_id, sat_id);

    Measurement meas;
    meas.t = t_meas;
    meas.type = obsv_type;
    meas.z = z;
    meas.observer_id = stat_id;
    meas.target_id = sat_id;
    meas.R.resize(1, 1);
    meas.R << sigma_range * sigma_range;

    ODWorldMeasurementEvent event;
    event.measurement = meas;
    event.observer_id = stat_id;
    event.target_id = sat_id;

    i32 prop_steps = 100;
    mat6d Q = mat6d1 * sigma_range;

    std::mt19937_64 rng(12345);
    MeasurementNoiseOptions noise_opts{.rng = rng, .enabled = true};
    ODEKFStepResult world_step_result
        = od_ekf_step_world(world, filter, event, dyn_config, prop_steps, Q);

    StateTr x_tr_observer = world.stat_x_tr_inertial(stat_id);
    ODEKFStepInput input;
    input.filter = filter;
    input.dyn_config = dyn_config;
    input.measurement = meas;
    input.prop_steps = prop_steps;
    input.Q = Q;
    input.tol_time = tol12;
    input.x_tr_observer = x_tr_observer;
    ODEKFStepResult step_result = od_ekf_step(input);

    print_diag_title("World EKF Step Diagnostic");
    std::println("Observer ID = {}", event.observer_id);
    std::println("Target ID = {}", event.target_id);
    std::println("Measurement Type = {}", observation_type_str(obsv_type));
    std::println();

    f64 initial_err = statetr_to_vec6d(filter.x - x0_truth).norm();
    std::println("Initial Error = {}", initial_err);

    std::println("World EKF Success = {}", od_status_success(world_step_result.status));
    std::println("World EKF Status: {}", status_string(world_step_result.status));
    if (od_status_success(step_result.status)) {
        StateTr world_err = world_step_result.filter.x - x0_truth;
        f64 world_state_err = statetr_to_vec6d(world_err).norm();
        std::println("World EKF Final Error = {}", world_state_err);
        std::println("World EKF Position Error = {}", world_err.r.norm());
        std::println("World EKF Velocity Error = {}", world_err.v.norm());
        std::println("World EKF Residual Norm = {}", world_step_result.residual_norm);
        std::println(
            "World EKF Raw Residual Norm = {}",
            world_step_result.raw_residual_norm
        );
        std::println();
    }

    std::println("Direct EKF Success = {}", od_status_success(step_result.status));
    std::println("Direct EKF Status: {}", status_string(step_result.status));
    if (od_status_success(step_result.status)) {
        StateTr direct_err = step_result.filter.x - x0_truth;
        f64 direct_state_err = statetr_to_vec6d(direct_err).norm();
        std::println("Direct EKF Final Error = {}", direct_state_err);
        std::println("Direct EKF Position Error = {}", direct_err.r.norm());
        std::println("Direct EKF Velocity Error = {}", direct_err.v.norm());
        std::println("Direct EKF Residual Norm = {}", step_result.residual_norm);
        std::println("Direct EKF Raw Residual Norm = {}", step_result.raw_residual_norm);
        std::println();
    }

    if (od_status_success(step_result.status)
        && od_status_success(world_step_result.status)) {
        StateTr world_direct_err = world_step_result.filter.x - step_result.filter.x;
        f64 world_direct_state_err = statetr_to_vec6d(world_direct_err).norm();
        f64 covariance_err = (world_step_result.filter.P - step_result.filter.P).norm();
        f64 residual_err = (world_step_result.residual - step_result.residual).norm();
        std::println("World/Direct State Error = {}", world_direct_state_err);
        std::println("World/Direct Position Error = {}", world_direct_err.r.norm());
        std::println("World/Direct Velocity Error = {}", world_direct_err.v.norm());
        std::println("World/Direct Covariance Error = {}", covariance_err);
        std::println("World/Direct Residual Error = {}", residual_err);
    }
}

void run_ekf_prediction_only_diag() {
    World world;

    // Earth
    EntityId earth_id = wgs84(world);
    Celestial* earth = world.celestial(earth_id);
    earth->name = "Earth";
    earth->gravity_model = GravityModel::pointmass;
    earth->degree = 0;
    earth->order = 0;

    // Satellite truth/estimate
    StateTr x0_truth;
    x0_truth.r = vec3d{7000.0, 1000.0, 1300.0};
    x0_truth.v = vec3d{-0.5, 7.2, 1.0};

    ODEKFState filter;
    filter.x = x0_truth + StateTr{.r = {1.0, 0.25, 0.1}, .v = {0.001, 0.0005, 0.0}};
    filter.P = mat6d1;
    filter.t = 0.0;

    ODDynamicsConfig dyn_config = make_od_cfg_from_celestial(*earth);
    dyn_config.tr_model = ODTrDynamicsModel::two_body;

    f64 t_target = 60.0;
    i32 prop_steps = 100;
    mat6d Q = mat6d1 * 1e-6;

    ODEKFStepResult step_result
        = od_ekf_predict_step(filter, t_target, dyn_config, prop_steps, Q);
    ODEKFPredictResult predict_result
        = od_ekf_predict(filter, t_target, dyn_config, prop_steps, Q);

    StateTr step_direct_err = step_result.filter.x - predict_result.y.x;
    f64 state_err = statetr_to_vec6d(step_direct_err).norm();
    f64 covariance_err = (step_result.filter.P - predict_result.P).norm();
    bool step_state_finite = statetr_to_vec6d(step_result.filter.x).allFinite();
    bool step_cov_finite = step_result.filter.P.allFinite();

    print_diag_title("EKF Prediction Only Diagnostic");
    std::println("Predict Target Time = {}", t_target);
    std::println("Prediction Only Success = {}", od_status_success(step_result.status));
    std::println("Prediction Only Status: {}", status_string(step_result.status));
    std::println("Raw Predict Status: {}", status_string(predict_result.status));
    std::println("Predicted Time = {}", step_result.filter.t);
    std::println("Expected Time = {}", t_target);
    std::println("Predicted State Finite = {}", step_state_finite);
    std::println("Predicted Covariance Finite = {}", step_cov_finite);
    std::println("Step/Raw State Error = {}", state_err);
    std::println("Step/Raw Position Error = {}", step_direct_err.r.norm());
    std::println("Step/Raw Velocity Error = {}", step_direct_err.v.norm());
    std::println("Step/Raw Covariance Error = {}", covariance_err);
    std::println("Prediction Covariance Norm = {}", step_result.filter.P.norm());
}

static StatusCode make_realtime_ekf_diag_schedule(
    const World& world,
    f64 t,
    i32 i,
    i32 i_meas,
    svec<EntityId> stat_ids,
    EntityId sat_id,
    svec<ODRealtimeScheduleItem>& schedule
) {
    schedule.clear();

    if (i % i_meas != 0) {
        // prediction only
        ODRealtimeScheduleItem item;
        item.t = t;
        item.has_measurement = false;

        schedule.push_back(item);
    } else {
        for (const EntityId stat_id : stat_ids) {
            const Station* stat = world.station(stat_id);
            if (stat == nullptr) {
                return StatusCode::observer_not_found;
            }

            // estimation
            ODRealtimeScheduleItem item;
            item.t = t;
            item.has_measurement = true;
            item.observer_id = stat_id;
            item.target_id = sat_id;

            svec<InstrumentId> ids = stat->enabled_instrument_ids;

            for (InstrumentId id : ids) {
                const StationInstrument& instrument = stat->instruments.at(id);
                if (instrument.enabled) {
                    item.instrument_id = id;
                    item.type = instrument.type;
                    schedule.push_back(item);
                }
            }
        }
    }

    return StatusCode::ok;
}

static StatusCode make_realtime_ekf_diag_events_from_schedule(
    const World& world,
    const svec<ODRealtimeScheduleItem>& schedule,
    svec<ODRealtimeEvent>& events,
    const MeasurementNoiseOptions& noise_opts
) {
    // materialize the schedule
    events.clear();

    for (const auto& item : schedule) {
        ODWorldMeasurementEvent event;

        if (item.has_measurement) {
            if (item.instrument_id == kInvalidInstrumentId) {
                return StatusCode::instrument_not_found;
            }

            StatusCode status;
            if (noise_opts.enabled) {
                status = make_noisy_world_measurement_event_instrument(
                    world,
                    item.instrument_id,
                    item.observer_id,
                    item.target_id,
                    item.t,
                    event,
                    noise_opts
                );
            } else {
                status = make_world_measurement_event_instrument(
                    world,
                    item.instrument_id,
                    item.observer_id,
                    item.target_id,
                    item.t,
                    event
                );
            }
            if (!od_status_success(status)) {
                return status;
            }

            ODRealtimeEvent realtime_event;
            realtime_event.event = event;
            realtime_event.has_measurement = item.has_measurement;
            realtime_event.t = item.t;
            events.push_back(realtime_event);
        } else {
            ODRealtimeEvent realtime_event;
            realtime_event.has_measurement = item.has_measurement;
            realtime_event.t = item.t;
            events.push_back(realtime_event);
        }
    }

    return StatusCode::ok;
}

void run_realtime_ekf_world_update_diag() {
    print_diag_title("Realtime World EKF Update Diagnostic");

    // Diagnostic parameters
    GravityModel truth_model = GravityModel::zonal;
    i32 truth_deg_ord = 6;
    GravityModel est_model = GravityModel::zonal;
    i32 est_degree = 4;

    auto scenario = make_earth_sats_stats_scenario(truth_model, truth_deg_ord);
    if (!scenario.success) {
        std::println("Scenario Build Failed");
        return;
    }
    World& world = scenario.world;
    EntityId earth_id = scenario.earth_id;
    EntityId sat_id = scenario.sat1_id;
    EntityId stat1_id = scenario.stat1_id;
    EntityId stat2_id = scenario.stat2_id;
    svec<EntityId> stat_ids = {stat1_id, stat2_id};
    Celestial* earth = world.celestial(earth_id);
    Satellite* sat = world.satellite(sat_id);
    Station* stat1 = world.station(stat1_id);
    Station* stat2 = world.station(stat2_id);
    if (earth == nullptr || sat == nullptr || stat1 == nullptr || stat2 == nullptr) {
        std::println("Failed to load scenario");
        return;
    }

    // Orbit Determination Config
    ODDynamicsConfig dyn_config = make_od_cfg_from_celestial(*earth);
    dyn_config.tr_model = worldtrmodel_to_odtrmodel(
        est_model
    ); // have od model and "truth" model discrepancy
    dyn_config.zonal_degree = est_degree;
    dyn_config.update_body_attitude = true;

    f64 sigma_rad = 1e-5;
    f64 sigma_range = 1e-3;
    f64 sigma_r = 1e-3;
    f64 sigma_v = 1e-6;

    mat6d R_pv = mat6d1;
    R_pv.block<3, 3>(0, 0) *= sigma_r * sigma_r;
    R_pv.block<3, 3>(3, 3) *= sigma_v * sigma_v;

    StatusCode radec_status
        = add_radec_instrument(*stat1, mat2d1 * sigma_rad * sigma_rad);
    if (!od_status_success(radec_status)) {
        std::println("Failed to add instrument");
        return;
    }

    StatusCode range_status
        = add_range_instrument(*stat1, matXd1<1> * sigma_range * sigma_range);
    if (!od_status_success(range_status)) {
        std::println("Failed to add instrument");
        return;
    }

    StatusCode azel_status = add_azel_instrument(*stat2, mat2d1 * sigma_rad * sigma_rad);
    if (!od_status_success(azel_status)) {
        std::println("Failed to add instrument");
        return;
    }

    InstrumentId rel_pv_id;
    StatusCode rel_posvel_status = add_rel_posvel_instrument(*stat2, R_pv, rel_pv_id);
    if (!od_status_success(rel_posvel_status)) {
        std::println("Failed to add instrument");
        return;
    }
    // disable_station_instrument(*stat2, rel_pv_id);

    WorldStepperConfig cfg;
    cfg.step_tr = true;
    cfg.step_att = true;
    cfg.substeps = 1;
    cfg.ticks = 10;
    cfg.dt_scale = 1.0 / cfg.ticks;
    cfg.integrator_tr = IntegratorTypeFixed::rk4;
    cfg.integrator_att = IntegratorTypeFixed::rk4;

    f64 t_span = 1000.0;
    i32 n_steps = 10000;
    f64 dt = t_span / n_steps;
    i32 prop_steps = 100;
    mat6d Q = mat6d1 * 1e-5; // process noise

    ODEKFState filter;
    filter.x = sat->x_tr + StateTr{.r = {1.0, 0.25, 0.1}, .v = {0.001, 0.0005, 0.0}};
    filter.P = mat6d1;
    filter.t = world.t_sim();

    ODEKFStepResult last_result;
    i32 measurement_updates = 0;
    i32 prediction_updates = 0;
    i32 failed_updates = 0;
    i32 processed_events = 0;
    i32 schedule_items_generated = 0;
    i32 materialized_events = 0;
    i32 measurement_events_generated = 0;
    i32 prediction_events_generated = 0;

    std::mt19937_64 rng(12345);
    MeasurementNoiseOptions noise_opts{.rng = rng, .enabled = true, .diagonal = false};

    for (i32 i = 0; i < n_steps; ++i) {
        f64 t_meas = world.t_sim();

        svec<ODRealtimeScheduleItem> schedule;
        StatusCode schedule_status = make_realtime_ekf_diag_schedule(
            world,
            t_meas,
            i,
            2,
            stat_ids,
            sat_id,
            schedule
        );

        if (schedule_status != StatusCode::ok) {
            last_result.status = schedule_status;
            ++failed_updates;
            break;
        }
        schedule_items_generated += schedule.size();

        svec<ODRealtimeEvent> events;
        StatusCode event_status = make_realtime_ekf_diag_events_from_schedule(
            world,
            schedule,
            events,
            noise_opts
        );
        if (event_status != StatusCode::ok) {
            last_result.status = event_status;
            ++failed_updates;
            break;
        }
        materialized_events += static_cast<i32>(events.size());
        for (const ODRealtimeEvent& event : events) {
            if (event.has_measurement) {
                ++measurement_events_generated;
            } else {
                ++prediction_events_generated;
            }
        }

        event_status = validate_realtime_ekf_events(events, t_meas);
        if (event_status != StatusCode::ok) {
            last_result.status = event_status;
            ++failed_updates;
            break;
        }

        ODRealtimeEKFResult realtime_result = od_ekf_update_world_events(
            world,
            filter,
            events,
            dyn_config,
            prop_steps,
            Q
        );
        if (!od_status_success(realtime_result.status)) {
            last_result.status = realtime_result.status;
            ++failed_updates;
            break;
        }

        if (failed_updates > 0) break;

        prediction_updates += realtime_result.prediction_updates;
        measurement_updates += realtime_result.measurement_updates;
        processed_events += realtime_result.processed_events;
        filter = realtime_result.filter;
        last_result
            = copy_od_ekf_result<ODRealtimeEKFResult, ODEKFStepResult>(realtime_result);

        WorldStepResult step_result = step_world(world, dt, cfg);
        if (step_result.status != StatusCode::ok) {
            last_result.status = step_result.status;
            ++failed_updates;
            break;
        }
    }

    // last ekf update after final world step
    if (failed_updates == 0 && std::abs(world.t_sim() - filter.t) > tol12) {
        ODRealtimeEKFInput input_predict_final{
            .world = &world,
            .filter = filter,
            .t_target = world.t_sim(),
            .dyn_config = dyn_config,
            .prop_steps = prop_steps,
            .Q = Q
        };

        last_result = od_ekf_update_world(input_predict_final);
        if (od_status_success(last_result.status)) {
            filter = last_result.filter;
            if (last_result.status == StatusCode::prediction_only) {
                ++prediction_updates;
                ++processed_events;
            }
        } else {
            ++failed_updates;
        }
    }

    StateTr final_err = filter.x - sat->x_tr;
    f64 final_state_err = statetr_to_vec6d(final_err).norm();

    std::println("Station 1 Instruments:");
    print_station_instruments(*stat1);
    std::println();

    std::println("Station 2 Instruments:");
    print_station_instruments(*stat2);
    std::println();

    std::println("Final Status: {}", status_string(last_result.status));
    std::println("Final Success = {}", od_status_success(last_result.status));
    std::println();

    std::println("Schedule Items Generated = {}", schedule_items_generated);
    std::println("Materialized Events = {}", materialized_events);
    std::println("Measurement Events Generated = {}", measurement_events_generated);
    std::println("Prediction Events Generated = {}", prediction_events_generated);
    std::println("Processed Events = {}", processed_events);
    std::println("Measurement Updates = {}", measurement_updates);
    std::println("Prediction Only Updates = {}", prediction_updates);
    std::println("Failed Updates = {}", failed_updates);
    std::println();

    std::println("Filter Time = {}", filter.t);
    std::println("World Time = {}", world.t_sim());
    std::println("Final State Error = {}", final_state_err);
    std::println("Final Position Error = {}", final_err.r.norm());
    std::println("Final Velocity Error = {}", final_err.v.norm());
    std::println("Final Covariance Norm = {}", filter.P.norm());
    print_diag_title();
}

void run_station_instrument_diag() {

    auto scenario = make_earth_sats_stats_scenario();

    World& world = scenario.world;
    Station* stat1 = world.station(scenario.stat1_id);
    if (stat1 == nullptr) {
        std::println("Station Instrument Diagnostic: station not found");
        return;
    }

    f64 sigma_rad = 1e-5;
    f64 sigma_range = 1e-3;

    StationInstrument radec_instr1;
    radec_instr1.type = ObservationType::radec;
    radec_instr1.enabled = true;
    radec_instr1.name = "Ra/Dec Observer";
    radec_instr1.R = mat2d1 * sigma_rad * sigma_rad;
    StatusCode radec_status1
        = add_station_instrument(*stat1, radec_instr1, radec_instr1.id);

    StationInstrument radec_instr2;
    radec_instr2.type = ObservationType::radec;
    radec_instr2.enabled = true;
    radec_instr2.name = "Ra/Dec Observer";
    radec_instr2.R = mat2d1 * sigma_rad * sigma_rad;
    StatusCode radec_status2
        = add_station_instrument(*stat1, radec_instr2, radec_instr2.id);

    StationInstrument range_instr;
    range_instr.type = ObservationType::range;
    range_instr.enabled = true;
    range_instr.name = "Range Observer";
    range_instr.R = matXd1<1> * sigma_range * sigma_range;
    StatusCode range_status = add_station_instrument(*stat1, range_instr, range_instr.id);

    matXd R_radec;
    StatusCode radec_id_status
        = station_measurement_covariance(*stat1, radec_instr1.id, R_radec);

    matXd R_range;
    StatusCode range_type_status
        = station_measurement_covariance(*stat1, ObservationType::range, R_range);

    matXd R_radec_type;
    StatusCode radec_type_status
        = station_measurement_covariance(*stat1, ObservationType::radec, R_radec_type);

    f64 t = world.t_sim();
    EntityId stat_id = scenario.stat1_id;
    EntityId sat_id = scenario.sat1_id;

    ODWorldMeasurementEvent range_type_event;
    StatusCode range_type_event_status = make_world_measurement_event(
        world,
        ObservationType::range,
        stat_id,
        sat_id,
        t,
        range_type_event
    );

    ODWorldMeasurementEvent radec_type_event;
    StatusCode radec_type_event_status = make_world_measurement_event(
        world,
        ObservationType::radec,
        stat_id,
        sat_id,
        t,
        radec_type_event
    );

    ODWorldMeasurementEvent radec_id_event;
    StatusCode radec_id_event_status = make_world_measurement_event_instrument(
        world,
        radec_instr1.id,
        stat_id,
        sat_id,
        t,
        radec_id_event
    );

    ODWorldMeasurementEvent invalid_id_event;
    StatusCode invalid_id_event_status = make_world_measurement_event_instrument(
        world,
        kInvalidInstrumentId,
        stat_id,
        sat_id,
        t,
        invalid_id_event
    );

    print_diag_title("Station Instrument Diagnostic");
    std::println("Add RA/Dec 1 Status: {}", status_string(radec_status1));
    std::println("Add RA/Dec 1 ID: {}", radec_instr1.id);
    std::println("Add RA/Dec 2 Status: {}", status_string(radec_status2));
    std::println("Add RA/Dec 2 ID: {}", radec_instr2.id);
    std::println("Add Range Status: {}", status_string(range_status));
    std::println("Add Range ID: {}", range_instr.id);
    std::println("Station Instrument Count: {}", stat1->instruments.size());
    std::println("Next Instrument ID: {}", stat1->next_instrument_id);
    std::println("RA/Dec Query By ID Status: {}", status_string(radec_id_status));
    std::println("RA/Dec Query By ID R Norm: {}", R_radec.norm());
    std::println("Range Query By Type Status: {}", status_string(range_type_status));
    std::println("Range Query By Type R Norm: {}", R_range.norm());
    std::println("RA/Dec Query By Type Status: {}", status_string(radec_type_status));
    std::println("RA/Dec Query By Type Expected Ambiguous = true");
    std::println(
        "Range Event By Type Status: {}",
        status_string(range_type_event_status)
    );
    std::println("Range Event By Type z Size: {}", range_type_event.measurement.z.size());
    std::println("Range Event By Type R Norm: {}", range_type_event.measurement.R.norm());
    std::println(
        "RA/Dec Event By Type Status: {}",
        status_string(radec_type_event_status)
    );
    std::println("RA/Dec Event By Type Expected Ambiguous = true");
    std::println("RA/Dec Event By ID Status: {}", status_string(radec_id_event_status));
    std::println(
        "RA/Dec Event By ID Type: {}",
        observation_type_str(radec_id_event.measurement.type)
    );
    std::println("RA/Dec Event By ID z Size: {}", radec_id_event.measurement.z.size());
    std::println("RA/Dec Event By ID R Norm: {}", radec_id_event.measurement.R.norm());
    std::println(
        "Invalid Instrument Event Status: {}",
        status_string(invalid_id_event_status)
    );
    print_diag_title();
}

void run_world_history_diag() {
    print_diag_title("World History Diagnostic");

    // Diagnostic parameters
    GravityModel truth_model = GravityModel::zonal;
    i32 truth_deg_ord = 6;
    GravityModel est_model = GravityModel::zonal;
    i32 est_degree = 4;

    auto scenario = make_earth_sats_stats_scenario(truth_model, truth_deg_ord);
    if (!scenario.success) {
        std::println("Scenario Build Failed");
        return;
    }
    World& world = scenario.world;
    EntityId earth_id = scenario.earth_id;
    EntityId sat_id = scenario.sat1_id;
    EntityId stat1_id = scenario.stat1_id;
    EntityId stat2_id = scenario.stat2_id;
    svec<EntityId> stat_ids = {stat1_id, stat2_id};
    Celestial* earth = world.celestial(earth_id);
    Satellite* sat = world.satellite(sat_id);
    Station* stat1 = world.station(stat1_id);
    Station* stat2 = world.station(stat2_id);
    if (earth == nullptr || sat == nullptr || stat1 == nullptr || stat2 == nullptr) {
        std::println("Failed to load scenario");
        return;
    }
    sat->propagate_att = true;

    WorldStepperConfig cfg;
    cfg.step_tr = true;
    cfg.step_att = true;
    cfg.substeps = 1;
    cfg.ticks = 10;
    cfg.dt_scale = 1.0 / cfg.ticks;
    cfg.integrator_tr = IntegratorTypeFixed::rk4;
    cfg.integrator_att = IntegratorTypeFixed::rk4;

    f64 t0 = 0.0;
    f64 t1 = 30.0;
    f64 t2 = 2.0 * t1;

    i32 n_steps1 = 250;
    i32 n_steps2 = 250;
    f64 dt1 = (t1 - t0) / n_steps1;
    f64 dt2 = (t2 - t1) / n_steps2;

    f64 dt;
    WorldHistorySample sample0 = capture_world_history_sample(world);
    WorldHistorySample sample1;
    WorldHistorySample sample2;

    StateTr x_tr_stat_t1;
    StateAtt x_att_stat_t1;
    vecXd z_radec_t1;
    vecXd z_azel_t1;

    for (i32 i = 0; i < n_steps1 + n_steps2; ++i) {
        dt = i < n_steps1 ? dt1 : dt2;

        WorldStepResult step_result = step_world(world, dt, cfg);
        if (step_result.status != StatusCode::ok) {
            std::println("World Step Status = {}", status_string(step_result.status));
            return;
        }

        if (i == (n_steps1 - 1)) {
            sample1 = capture_world_history_sample(world);
            x_tr_stat_t1 = world.stat_x_tr_inertial(stat1_id);
            x_att_stat_t1 = world.stat_x_att_inertial(stat1_id);
            z_radec_t1 = world_predict_measurement(
                world,
                ObservationType::radec,
                stat1_id,
                sat_id
            );
            z_azel_t1 = world_predict_measurement(
                world,
                ObservationType::azel,
                stat1_id,
                sat_id
            );
        } else if (i == n_steps1 + n_steps2 - 1) {
            sample2 = capture_world_history_sample(world);
        }
    }

    WorldHistory history;
    history.max_samples = 3;
    push_world_history_sample(history, sample0);
    // push_world_history_sample(history, sample1);
    push_world_history_sample(history, sample2);

    vecXd z_radec_t1_interp;
    StatusCode z_radec_status = world_predict_measurement_history(
        world,
        history,
        ObservationType::radec,
        stat1_id,
        sat_id,
        t1,
        z_radec_t1_interp
    );
    vecXd z_azel_t1_interp;
    StatusCode z_azel_status = world_predict_measurement_history(
        world,
        history,
        ObservationType::azel,
        stat1_id,
        sat_id,
        t1,
        z_azel_t1_interp
    );

    StateTr x_tr_sat_t1 = sample1.x_tr.at(sat_id);
    StateAtt x_att_sat_t1 = sample1.x_att.at(sat_id);

    StateTr x_tr_sat_t1_interp;
    StatusCode sat_tr_status
        = sample_tr_interp_linear(history, sat_id, t1, x_tr_sat_t1_interp);
    StateAtt x_att_sat_t1_interp;
    StatusCode sat_att_status
        = sample_att_interp_linear(history, sat_id, t1, x_att_sat_t1_interp);

    StateTr x_tr_stat_t1_interp;
    StatusCode stat_tr_status = sample_station_tr_interp_linear(
        world,
        history,
        stat1_id,
        t1,
        x_tr_stat_t1_interp
    );
    StateAtt x_att_stat_t1_interp;
    StatusCode stat_att_status = sample_station_att_interp_linear(
        world,
        history,
        stat1_id,
        t1,
        x_att_stat_t1_interp
    );

    StateTr x_tr_sat_t1_err = x_tr_sat_t1 - x_tr_sat_t1_interp;
    StateAtt x_att_sat_t1_err = x_att_sat_t1 - x_att_sat_t1_interp;

    StateTr x_tr_stat_t1_err = x_tr_stat_t1 - x_tr_stat_t1_interp;
    StateAtt x_att_stat_t1_err = x_att_stat_t1 - x_att_stat_t1_interp;

    std::println("Sampled at t1 = {}s", t1);
    std::println("Satellite Queries ------");
    if (!od_status_success(sat_tr_status)) {
        std::println("Translation Status: {}", status_string(sat_tr_status));
    }
    if (!od_status_success(sat_att_status)) {
        std::println("Attitude Status: {}", status_string(sat_att_status));
    }
    std::println("Sampled Position = {}", x_tr_sat_t1.r);
    std::println("Interpolated Position = {}", x_tr_sat_t1_interp.r);
    std::println(
        "Translational State Error = {}",
        statetr_to_vec6d(x_tr_sat_t1_err).norm()
    );
    std::println("Position Error = {}", x_tr_sat_t1_err.r.norm());
    std::println("Velocity Error = {}", x_tr_sat_t1_err.v.norm());
    std::println("Sampled Attitude = {}", x_att_sat_t1.q);
    std::println("Interpolated Attitude = {}", x_att_sat_t1_interp.q);
    std::println("Attitude Error = {}", x_att_sat_t1_err.q.norm());
    std::println("Angular Velocity Error = {}", x_att_sat_t1_err.w.norm());
    std::println();
    std::println("Station Queries ------");
    if (!od_status_success(stat_tr_status)) {
        std::println("Translation Status: {}", status_string(stat_tr_status));
    }
    if (!od_status_success(stat_att_status)) {
        std::println("Attitude Status: {}", status_string(stat_att_status));
    }
    std::println("Sampled Position = {}", x_tr_stat_t1.r);
    std::println("Interpolated Position = {}", x_tr_stat_t1_interp.r);
    std::println(
        "Translational State Error = {}",
        statetr_to_vec6d(x_tr_stat_t1_err).norm()
    );
    std::println("Position Error = {}", x_tr_stat_t1_err.r.norm());
    std::println("Velocity Error = {}", x_tr_stat_t1_err.v.norm());
    std::println("Sampled Attitude = {}", x_att_stat_t1.q);
    std::println("Interpolated Attitude = {}", x_att_stat_t1_interp.q);
    std::println("Attitude Error = {}", x_att_stat_t1_err.q.norm());
    std::println("Angular Velocity Error = {}", x_att_stat_t1_err.w.norm());
    std::println();

    std::println("Radec Measurements ------");
    if (!od_status_success(z_radec_status)) {
        std::println("Measurement Status: {}", status_string(z_radec_status));
    }
    std::println("Sampled Measurement = {}", z_radec_t1);
    std::println("Interpolated Measurment = {}", z_radec_t1_interp);
    if (z_radec_t1.size() == 2 && z_radec_t1_interp.size() == 2) {
        std::println("Measurement Error = {}", (z_radec_t1 - z_radec_t1_interp).norm());
    } else {
        std::println("Error: measurement sizes don't match");
    }
    std::println();

    std::println("Az/El Measurements ------");
    if (!od_status_success(z_azel_status)) {
        std::println("Measurement Status: {}", status_string(z_azel_status));
    }
    std::println("Sampled Measurement = {}", z_azel_t1);
    std::println("Interpolated Measurment = {}", z_azel_t1_interp);
    if (z_azel_t1.size() == 2 && z_azel_t1_interp.size() == 2) {
        std::println("Measurement Error = {}", (z_azel_t1 - z_azel_t1_interp).norm());
    } else {
        std::println("Error: measurement sizes don't match");
    }
    std::println();

    // bad queries
    StateTr x_tr_bad;
    StatusCode status;
    status = sample_tr_interp_linear(history, sat_id, t0 - 1.0, x_tr_bad);
    std::println("Query t < t0, status: {}", status_string(status));
    status = sample_tr_interp_linear(history, sat_id, t2 + 1.0, x_tr_bad);
    std::println("Query t > t2, status: {}", status_string(status));
    status = sample_tr_interp_linear(history, sat_id + 100, t1, x_tr_bad);
    std::println("Query missing body, status: {}", status_string(status));
    WorldHistory history_empty;
    status = sample_tr_interp_linear(history_empty, sat_id, t1, x_tr_bad);
    std::println("Query empty history, status: {}", status_string(status));
}

static StatusCode make_realtime_ekf_diag_history_events_from_schedule(
    const World& world,
    const WorldHistory& history,
    const svec<ODRealtimeScheduleItem>& schedule,
    svec<ODRealtimeEvent>& events,
    const MeasurementNoiseOptions& noise_opts,
    const HistorySampleOptions& sample_opts
) {
    // materialize the schedule
    events.clear();

    for (const auto& item : schedule) {
        ODWorldMeasurementEvent event;

        if (item.has_measurement) {
            if (item.instrument_id == kInvalidInstrumentId) {
                return StatusCode::instrument_not_found;
            }

            StatusCode status;
            if (noise_opts.enabled) {
                status = make_noisy_world_measurement_event_history_instrument(
                    world,
                    history,
                    item.instrument_id,
                    item.observer_id,
                    item.target_id,
                    item.t,
                    event,
                    noise_opts,
                    sample_opts
                );
            } else {
                status = make_world_measurement_event_history_instrument(
                    world,
                    history,
                    item.instrument_id,
                    item.observer_id,
                    item.target_id,
                    item.t,
                    event,
                    sample_opts
                );
            }
            if (!od_status_success(status)) {
                return status;
            }

            ODRealtimeEvent realtime_event;
            realtime_event.event = event;
            realtime_event.has_measurement = item.has_measurement;
            realtime_event.t = item.t;
            events.push_back(realtime_event);
        } else {
            ODRealtimeEvent realtime_event;
            realtime_event.has_measurement = item.has_measurement;
            realtime_event.t = item.t;
            events.push_back(realtime_event);
        }
    }

    return StatusCode::ok;
}

void run_world_history_ekf_diag() {
    print_diag_title("Realtime World EKF Update Diagnostic (History-Backed)");

    // Diagnostic parameters
    GravityModel truth_model = GravityModel::zonal;
    i32 truth_deg_ord = 6;
    GravityModel est_model = GravityModel::zonal;
    i32 est_degree = 4;

    auto scenario = make_earth_sats_stats_scenario(truth_model, truth_deg_ord);
    if (!scenario.success) {
        std::println("Scenario Build Failed");
        return;
    }
    World& world = scenario.world;
    EntityId earth_id = scenario.earth_id;
    EntityId sat_id = scenario.sat1_id;
    EntityId stat1_id = scenario.stat1_id;
    EntityId stat2_id = scenario.stat2_id;
    svec<EntityId> stat_ids = {stat1_id, stat2_id};
    Celestial* earth = world.celestial(earth_id);
    Satellite* sat = world.satellite(sat_id);
    Station* stat1 = world.station(stat1_id);
    Station* stat2 = world.station(stat2_id);
    if (earth == nullptr || sat == nullptr || stat1 == nullptr || stat2 == nullptr) {
        std::println("Failed to load scenario");
        return;
    }

    // Orbit Determination Config
    ODDynamicsConfig dyn_config = make_od_cfg_from_celestial(*earth);
    dyn_config.tr_model = worldtrmodel_to_odtrmodel(
        est_model
    ); // have od model and "truth" model discrepancy
    dyn_config.zonal_degree = est_degree;
    dyn_config.update_body_attitude = true;

    f64 sigma_rad = 1e-5;
    f64 sigma_range = 1e-3;
    f64 sigma_r = 1e-3;
    f64 sigma_v = 1e-6;

    mat6d R_pv = mat6d1;
    R_pv.block<3, 3>(0, 0) *= sigma_r * sigma_r;
    R_pv.block<3, 3>(3, 3) *= sigma_v * sigma_v;

    StatusCode radec_status
        = add_radec_instrument(*stat1, mat2d1 * sigma_rad * sigma_rad);
    if (!od_status_success(radec_status)) {
        std::println("Failed to add instrument");
        return;
    }

    StatusCode range_status
        = add_range_instrument(*stat1, matXd1<1> * sigma_range * sigma_range);
    if (!od_status_success(range_status)) {
        std::println("Failed to add instrument");
        return;
    }

    StatusCode azel_status = add_azel_instrument(*stat2, mat2d1 * sigma_rad * sigma_rad);
    if (!od_status_success(azel_status)) {
        std::println("Failed to add instrument");
        return;
    }

    InstrumentId rel_pv_id;
    StatusCode rel_posvel_status = add_rel_posvel_instrument(*stat2, R_pv, rel_pv_id);
    if (!od_status_success(rel_posvel_status)) {
        std::println("Failed to add instrument");
        return;
    }
    // disable_station_instrument(*stat2, rel_pv_id);

    WorldStepperConfig cfg;
    cfg.step_tr = true;
    cfg.step_att = true;
    cfg.substeps = 1;
    cfg.ticks = 10;
    cfg.dt_scale = 1.0 / cfg.ticks;
    cfg.integrator_tr = IntegratorTypeFixed::rk4;
    cfg.integrator_att = IntegratorTypeFixed::rk4;

    f64 t_span = 1000.0;
    i32 n_steps = 1000;
    f64 dt = t_span / n_steps;
    i32 prop_steps = 100;
    mat6d Q = mat6d1 * 1e-5; // process noise

    ODEKFState filter;
    filter.x = sat->x_tr + StateTr{.r = {1.0, 0.25, 0.1}, .v = {0.001, 0.0005, 0.0}};
    filter.P = mat6d1;
    filter.t = world.t_sim();

    ODEKFStepResult last_result;
    i32 measurement_updates = 0;
    i32 prediction_updates = 0;
    i32 failed_updates = 0;
    i32 processed_events = 0;
    i32 schedule_items_generated = 0;
    i32 materialized_events = 0;
    i32 measurement_events_generated = 0;
    i32 prediction_events_generated = 0;

    std::mt19937_64 rng(12345);
    MeasurementNoiseOptions noise_opts{.rng = rng, .enabled = true, .diagonal = false};

    HistorySampleOptions sample_opts;
    sample_opts.tr_interp = HistoryInterpolation::nearest;
    sample_opts.att_interp = HistoryInterpolation::nearest;

    WorldHistory history;
    history.max_samples = 10;

    WorldHistorySample sample = capture_world_history_sample(world);
    push_world_history_sample(history, sample);

    for (i32 i = 0; i < n_steps; ++i) {
        f64 t_meas = world.t_sim();

        svec<ODRealtimeScheduleItem> schedule;
        StatusCode schedule_status = make_realtime_ekf_diag_schedule(
            world,
            t_meas,
            i,
            2,
            stat_ids,
            sat_id,
            schedule
        );

        if (schedule_status != StatusCode::ok) {
            last_result.status = schedule_status;
            ++failed_updates;
            break;
        }
        schedule_items_generated += schedule.size();
        // for (auto& item : schedule) {
        //     item.t = item.t - dt / 4.0;
        // }
        // TODO: the above runs if you comment out
        // if (event.t < t_prev - tol_time) {
        //     return StatusCode::time_mismatch;
        // }
        // in estimation_world.cpp and
        // if (dt < -tol) {
        //     result.status = StatusCode::propagation_failed;
        //     return result;
        // } else if (std::abs(dt) <= tol) {
        //     // same epoch
        //     dt = 0.0;
        // }
        // in estimation_recursive.
        // will likely need to allow backwards propagation with the above comments, do the
        // prediction/measurement/ekf, then propagate forward again back to the world step
        // time?

        svec<ODRealtimeEvent> events;
        StatusCode event_status = make_realtime_ekf_diag_history_events_from_schedule(
            world,
            history,
            schedule,
            events,
            noise_opts,
            sample_opts
        );
        if (event_status != StatusCode::ok) {
            last_result.status = event_status;
            ++failed_updates;
            break;
        }
        materialized_events += static_cast<i32>(events.size());
        for (const ODRealtimeEvent& event : events) {
            if (event.has_measurement) {
                ++measurement_events_generated;
            } else {
                ++prediction_events_generated;
            }
        }

        event_status = validate_realtime_ekf_events(events, t_meas);
        if (event_status != StatusCode::ok) {
            last_result.status = event_status;
            ++failed_updates;
            break;
        }

        ODRealtimeEKFResult realtime_result = od_ekf_update_world_events(
            world,
            filter,
            events,
            dyn_config,
            prop_steps,
            Q
        );
        if (!od_status_success(realtime_result.status)) {
            last_result.status = realtime_result.status;
            ++failed_updates;
            break;
        }

        if (failed_updates > 0) break;

        prediction_updates += realtime_result.prediction_updates;
        measurement_updates += realtime_result.measurement_updates;
        processed_events += realtime_result.processed_events;
        filter = realtime_result.filter;
        last_result
            = copy_od_ekf_result<ODRealtimeEKFResult, ODEKFStepResult>(realtime_result);

        WorldStepResult step_result = step_world(world, dt, cfg);
        if (step_result.status != StatusCode::ok) {
            last_result.status = step_result.status;
            ++failed_updates;
            break;
        }

        sample = capture_world_history_sample(world);
        push_world_history_sample(history, sample);
    }

    // last ekf update after final world step
    if (failed_updates == 0 && std::abs(world.t_sim() - filter.t) > tol12) {
        ODRealtimeEKFInput input_predict_final{
            .world = &world,
            .filter = filter,
            .t_target = world.t_sim(),
            .dyn_config = dyn_config,
            .prop_steps = prop_steps,
            .Q = Q
        };

        last_result = od_ekf_update_world(input_predict_final);
        if (od_status_success(last_result.status)) {
            filter = last_result.filter;
            if (last_result.status == StatusCode::prediction_only) {
                ++prediction_updates;
                ++processed_events;
            }
        } else {
            ++failed_updates;
        }
    }

    StateTr final_err = filter.x - sat->x_tr;
    f64 final_state_err = statetr_to_vec6d(final_err).norm();

    std::println("Station 1 Instruments:");
    print_station_instruments(*stat1);
    std::println();

    std::println("Station 2 Instruments:");
    print_station_instruments(*stat2);
    std::println();

    std::println("Final Status: {}", status_string(last_result.status));
    std::println("Final Success = {}", od_status_success(last_result.status));
    std::println();

    std::println("Schedule Items Generated = {}", schedule_items_generated);
    std::println("Materialized Events = {}", materialized_events);
    std::println("Measurement Events Generated = {}", measurement_events_generated);
    std::println("Prediction Events Generated = {}", prediction_events_generated);
    std::println("Processed Events = {}", processed_events);
    std::println("Measurement Updates = {}", measurement_updates);
    std::println("Prediction Only Updates = {}", prediction_updates);
    std::println("Failed Updates = {}", failed_updates);
    std::println();

    std::println("Filter Time = {}", filter.t);
    std::println("World Time = {}", world.t_sim());
    std::println("Final State Error = {}", final_state_err);
    std::println("Final Position Error = {}", final_err.r.norm());
    std::println("Final Velocity Error = {}", final_err.v.norm());
    std::println("Final Covariance Norm = {}", filter.P.norm());
    print_diag_title();

    HistoryCSVExportOptions csv_opts;
    csv_opts.include_attitude = true;
    StatusCode status = write_world_history_csv(
        world,
        history,
        pwd + "/assets/output/sat_csv_out_test.csv",
        csv_opts
    );
    std::println("Write Status: {}", status_string(status));
}

static void print_state_tr_config(
    const ScenarioStateTrConfig& x_tr,
    const string& indent
) {
    std::print("{}State (tr, {}): ", indent, state_tr_type_str(x_tr.input_type));
    switch (x_tr.input_type) {
    case StateTrInputType::pos_vel: {
        std::println("r = {}, v = {}", vec_string(x_tr.r), vec_string(x_tr.v));
    } break;
    case StateTrInputType::classical: {
        std::println(
            "sma = {}, ecc: {}, inc = {}, raan = {}, aop = {}, ta = {}",
            x_tr.coes.sma,
            x_tr.coes.ecc,
            x_tr.coes.inc,
            x_tr.coes.raan,
            x_tr.coes.aop,
            x_tr.coes.ta
        );
    } break;
    }
}
static void print_state_att_config(
    const ScenarioStateAttConfig& x_att,
    const string& indent
) {
    std::print("{}State (att, {}): ", indent, state_att_type_str(x_att.input_type));
    switch (x_att.input_type) {
    case AttitudeType::quaternion: {
        std::print("q = {}, ", vec_string(x_att.q));
    } break;
    case AttitudeType::dcm: {
        std::print("dcm: R = {}, ", vec_string(x_att.dcm));
    } break;
    case AttitudeType::axis_angle: {
        std::print(
            "axis-angle: lambda = {}, angle: theta = {}, ",
            x_att.axis,
            x_att.angle
        );
    } break;
    case AttitudeType::euler_angles: {
        std::print(
            "euler angles: angles = {}, sequence: {}, ",
            vec_string(x_att.angles),
            vec_string(x_att.sequence)
        );
    } break;
    case AttitudeType::crp: {
        std::print("classical rodriguez parameters: rho = {}, ", vec_string(x_att.axis));
    } break;
    case AttitudeType::mrp: {
        std::print("modified rodriguez parameters: rho = {}, ", vec_string(x_att.axis));
    } break;
    }
    std::println("w = {}", vec_string(x_att.w));
}
static void print_mass_properties_config(
    const ScenarioMassPropertiesConfig& mp,
    const string& indent
) {
    std::println("{}Mass: {}", indent, mp.mass);
    std::println("{}Inertia Tensor: I = {}", indent, vec_string(mp.inertia));
}

void run_scenario_loader_diag() {
    print_diag_title("JSON Scenario Loader");

    // TODO: create function to print body properties (cel, sat, stat, etc.)

    ScenarioConfig config;
    StatusCode status
        = load_scenario_json(pwd + "/scenarios/parser_stress_demo.json", config);
    // = load_scenario_json(pwd + "/scenarios/earth_moon_sat_demo.json", config);
    if (status == StatusCode::ok)
        std::println("JSON scenario loaded successfully");
    else
        std::println("{}", status_string(status));

    status = validate_scenario_config(config);
    if (status != StatusCode::ok) {
        std::println("{}", status_string(status));
    } else {
        std::println("Scenario configuration validated");

        std::println();
        std::println(
            "Schema: name = {}, version = {}",
            config.schema.name,
            config.schema.version
        );

        std::println(
            "Metadata: name = {}, rng seed = {}",
            config.metadata.name,
            config.metadata.rng_seed
        );
        std::println();

        std::print(
            "Time: type = {}, time scale = {}, ",
            date_type_str(config.time.date_type),
            time_scale_str(config.time.time_scale)
        );
        switch (config.time.date_type) {
        case DateType::cal: {
            print_cal(config.time.cal, CalendarPrintStyle::string);
        } break;
        case DateType::jd: {
            std::println("jd: {}", jd_to_scalar(config.time.jd));
        } break;
        case DateType::mjd: {
            std::println("mjd: {}", mjd_to_scalar(config.time.mjd));
        } break;
        }
        std::println();

        string indent = "    ";
        string indent2 = indent + indent;

        std::println("Providers count: {}", config.gravity_providers.size());
        for (const auto& gravity_provider : config.gravity_providers) {
            std::println("--- Gravity Provider ID: {}", gravity_provider.id);
            std::println("{}Type: {}", indent, gravity_provider.format);
            std::println("{}Filepath: {}", indent, gravity_provider.filepath);
            std::println("{}Line Skips: {}", indent, gravity_provider.lineskips);
            std::println("{}Normalized: {}", indent, gravity_provider.normalized);
        }
        std::println();

        std::println("Celestials count: {}", config.celestials.size());
        for (const auto& cel : config.celestials) {
            std::println("--- Celestial ID: {}", cel.id);
            std::println("{}Name: {}", indent, cel.name);
            std::println("{}Model: {}", indent, cel.model.id);
            std::println(
                "{}Gravity Model: {}",
                indent,
                gravity_model_str(cel.model.gravity_model.model)
            );
            std::println(
                "{}Gravity Coefficient Source: {}",
                indent,
                gravity_coefficient_source_str(cel.model.gravity_model.coefficient_source)
            );
            if (cel.model.gravity_model.coefficient_source
                == GravityCoefficientSource::direct_zonal) {
                std::println("{}J: {}", indent, cel.model.gravity_model.J);
            }
            print_state_tr_config(cel.x_tr, indent);
            print_state_att_config(cel.x_att, indent);
        }
        std::println();

        std::println("Satellites count: {}", config.satellites.size());
        for (const auto& sat : config.satellites) {
            std::println("--- Satellite ID: {}", sat.id);
            std::println("{}Name: {}", indent, sat.name);
            print_state_tr_config(sat.x_tr, indent);
            print_state_att_config(sat.x_att, indent);
            print_mass_properties_config(sat.mass_properties, indent);
        }

        std::println();
        std::println("Stations count: {}", config.stations.size());
        for (const auto& stat : config.stations) {
            std::println("--- Station ID: {}", stat.id);
            std::println("Name: {}", stat.name);
            if (stat.anchored) {
                std::println("{}Anchor name: {}", indent, stat.anchor);
                std::println("{}LLH (BCBF): {}", indent, stat.llh);
                std::println("{}Local Frame: {}", indent, stat.local_frame);
            } else {
                print_mass_properties_config(stat.mass_properties, indent);
                print_state_tr_config(stat.x_tr, indent);
                print_state_att_config(stat.x_att, indent);
            }

            std::println("{}Instruments count: {}", indent, stat.instruments.size());
            for (const auto& instrument : stat.instruments) {
                std::println("{}--- Instrument ID: {}", indent2, instrument.id);
                std::println(
                    "{}Instrument Type: {}",
                    indent2,
                    observation_type_str(instrument.type)
                );
                std::println(
                    "{}Measurement Covariance Dim: ({}, {})",
                    indent2,
                    instrument.covariance_cfg.covariance.rows(),
                    instrument.covariance_cfg.covariance.cols()
                );
                std::println(
                    "{}Measurement Covariance: R = {}",
                    indent2,
                    vec_string(instrument.covariance_cfg.covariance)
                );
            }
        }
    }

    print_diag_title();
}

void run_build_world_from_scenario() {
    print_diag_title("JSON Scenario World Builder");

    ScenarioConfig cfg;
    StatusCode status
        = load_scenario_json(pwd + "/scenarios/parser_stress_demo.json", cfg);
    // = load_scenario_json(pwd + "/scenarios/earth_moon_sat_demo.json", config);
    if (status == StatusCode::ok) {
        std::println("Parsing: Successful");
    } else {
        std::println("Parsing Error {}", status_string(status));
        return;
    }

    status = validate_scenario_config(cfg);
    if (status == StatusCode::ok) {
        std::println("Validation: Successful");
    } else {
        std::println("Validation Error: {}", status_string(status));
        return;
    }

    std::println("Scenario configuration validated");

    std::println();
    std::println("Schema: name = {}, version = {}", cfg.schema.name, cfg.schema.version);

    std::println(
        "Metadata: name = {}, rng seed = {}",
        cfg.metadata.name,
        cfg.metadata.rng_seed
    );
    std::println();

    std::print(
        "Time: type = {}, time scale = {}, ",
        date_type_str(cfg.time.date_type),
        time_scale_str(cfg.time.time_scale)
    );
    switch (cfg.time.date_type) {
    case DateType::cal: {
        print_cal(cfg.time.cal, CalendarPrintStyle::string);
    } break;
    case DateType::jd: {
        std::println("jd: {}", jd_to_scalar(cfg.time.jd));
    } break;
    case DateType::mjd: {
        std::println("mjd: {}", mjd_to_scalar(cfg.time.mjd));
    } break;
    }
    std::println();

    World world;
    WorldStepperConfig stepper_cfg;
    ScenarioBuildResult result;
    status = build_world_from_scenario_config(cfg, world, result, stepper_cfg);
    if (status == StatusCode::ok) {
        std::println("World Build: Successful");
    } else {
        std::println("World Build Error: {}", status_string(status));
        return;
    }

    std::println();
    std::println("World Counts");
    std::println("    Celestials: {}", world.num_active_celestials());
    std::println("    Satellites: {}", world.num_active_satellites());
    std::println("    Stations: {}", world.num_active_stations());
    std::println();

    std::println("Build Result Maps");
    std::println("    Body IDs: {}", result.body_ids.size());
    std::println("    Celestial IDs: {}", result.celestial_ids.size());
    std::println("    Satellite IDs: {}", result.satellite_ids.size());
    std::println("    Station IDs: {}", result.station_ids.size());
    std::println();

    std::println("Built Celestials");
    for (EntityId id : world.active_celestial_ids()) {
        const Celestial* cel = world.celestial(id);
        if (cel == nullptr) continue;
        print_celestial(*cel);
    }
    std::println();

    std::println("Built Satellites");
    for (EntityId id : world.active_satellite_ids()) {
        const Satellite* sat = world.satellite(id);
        if (sat == nullptr) continue;
        print_satellite(*sat);
    }
    std::println();

    std::println("Built Stations");
    for (EntityId id : world.active_station_ids()) {
        const Station* stat = world.station(id);
        if (stat == nullptr) continue;
        print_station(*stat);
    }

    print_diag_title();
}

void run_scenario_render_diag() {
    print_diag_title("Scenario Render Diag");

    const string scenario_filepath = pwd + "/scenarios/main_demo.json";
    ScenarioConfig scenario_cfg;
    StatusCode status = load_scenario_json(scenario_filepath, scenario_cfg);
    if (status != StatusCode::ok) {
        std::println("Parsing Error: {}", status_string(status));
        return;
    }

    status = validate_scenario_config(scenario_cfg);
    if (status != StatusCode::ok) {
        std::println("Validation Error: {}", status_string(status));
        return;
    }

    World world;
    WorldStepperConfig stepper_cfg;
    ScenarioBuildResult build_result;
    status = build_world_from_scenario_config(
        scenario_cfg,
        world,
        build_result,
        stepper_cfg
    );
    if (status != StatusCode::ok) {
        std::println("World Build Error: {}", status_string(status));
        return;
    }

    std::println(
        "World Build: celestials = {}, satellites = {}, stations = {}",
        world.num_active_celestials(),
        world.num_active_satellites(),
        world.num_active_stations()
    );

    // integrator
    f64 t0 = 0.0;
    f64 dt = 1.0;
    world.reset_time(t0);

    // windowing and graphics
    // TODO: make build_render_config
    RenderLoopConfig render_cfg{};
    render_cfg.window_title = "Scenario Render Diag";
    render_cfg.stepper_cfg = stepper_cfg;
    render_cfg.camera.position = vec3f{1.0f, 1.0f, 1.0f} * 75000.0f;
    render_cfg.camera.fovy = 45.0f;
    render_cfg.camera.projection = CAMERA_PERSPECTIVE;
    render_cfg.camera.up = axis_zf;
    render_cfg.camera.target = originf;
    rlSetClipPlanes(1.0e3, 1.0e6);

    SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_MSAA_4X_HINT | FLAG_WINDOW_HIGHDPI);
    ScenarioSession scenario;
    scenario.config = std::move(scenario_cfg);
    scenario.build_result = std::move(build_result);
    scenario.filepath = scenario_filepath;
    scenario.has_filepath = true;
    run_world_render_loop(world, render_cfg, dt, std::move(scenario));

    print_diag_title();
}

void run_world_dopri54_rk4_diag() {
    print_diag_title("World DOPRI5(4) vs RK4");

    ScenarioConfig rk4_scenario_cfg;
    StatusCode status = load_scenario_json(
        pwd + "/scenarios/earth_moon_sat_demo.json",
        rk4_scenario_cfg
    );
    if (status != StatusCode::ok) {
        std::println("RK4 Scenario Load Status = {}", status_string(status));
        return;
    }

    status = validate_scenario_config(rk4_scenario_cfg);
    if (status != StatusCode::ok) {
        std::println("RK4 Scenario Validation Status = {}", status_string(status));
        return;
    }

    ScenarioConfig adaptive_scenario_cfg;
    status = load_scenario_json(
        pwd + "/scenarios/earth_moon_sat_adaptive_demo.json",
        adaptive_scenario_cfg
    );
    if (status != StatusCode::ok) {
        std::println("Adaptive Scenario Load Status = {}", status_string(status));
        return;
    }

    status = validate_scenario_config(adaptive_scenario_cfg);
    if (status != StatusCode::ok) {
        std::println("Adaptive Scenario Validation Status = {}", status_string(status));
        return;
    }

    const std::filesystem::path roundtrip_filepath
        = std::filesystem::temp_directory_path()
        / "astrolib_earth_moon_sat_adaptive_roundtrip.json";

    status = save_scenario_json(roundtrip_filepath.string(), adaptive_scenario_cfg);
    if (status != StatusCode::ok) {
        std::println("Adaptive Scenario Save Status = {}", status_string(status));
        return;
    }

    ScenarioConfig adaptive_roundtrip_cfg;
    status = load_scenario_json(roundtrip_filepath.string(), adaptive_roundtrip_cfg);
    std::error_code remove_error;
    std::filesystem::remove(roundtrip_filepath, remove_error);
    if (status != StatusCode::ok) {
        std::println("Adaptive Scenario Reload Status = {}", status_string(status));
        return;
    }

    auto same_adaptive_opts = [](
                                  const AdaptiveIntegratorConfig& lhs,
                                  const AdaptiveIntegratorConfig& rhs
                              ) {
        return lhs.rel_tol == rhs.rel_tol && lhs.abs_tol_r == rhs.abs_tol_r
            && lhs.abs_tol_v == rhs.abs_tol_v
            && lhs.abs_tol_angle == rhs.abs_tol_angle && lhs.abs_tol_w == rhs.abs_tol_w
            && lhs.dt_initial == rhs.dt_initial && lhs.dt_min == rhs.dt_min
            && lhs.dt_max == rhs.dt_max && lhs.safety == rhs.safety
            && lhs.scale_min == rhs.scale_min && lhs.scale_max == rhs.scale_max
            && lhs.max_attempts == rhs.max_attempts
            && lhs.max_rejections == rhs.max_rejections;
    };

    auto same_stepper_config = [&](
                                   const ScenarioWorldStepperConfig& lhs,
                                   const ScenarioWorldStepperConfig& rhs
                               ) {
        return lhs.integrator_tr == rhs.integrator_tr
            && lhs.integrator_att == rhs.integrator_att
            && lhs.adaptive.use_substeps == rhs.adaptive.use_substeps
            && same_adaptive_opts(lhs.adaptive.opts, rhs.adaptive.opts)
            && lhs.paused == rhs.paused && lhs.step_tr == rhs.step_tr
            && lhs.step_att == rhs.step_att && lhs.substeps == rhs.substeps
            && lhs.ticks == rhs.ticks && lhs.dt_scale == rhs.dt_scale;
    };

    bool roundtrip_match = same_stepper_config(
        adaptive_scenario_cfg.world_stepper,
        adaptive_roundtrip_cfg.world_stepper
    );
    if (!roundtrip_match) {
        std::println("Adaptive Scenario Stepper Round Trip Match = false");
        return;
    }

    World rk4_world;
    WorldStepperConfig rk4_cfg;
    ScenarioBuildResult rk4_build;
    status
        = build_world_from_scenario_config(rk4_scenario_cfg, rk4_world, rk4_build, rk4_cfg);
    if (status != StatusCode::ok) {
        std::println("RK4 World Build Status = {}", status_string(status));
        return;
    }

    World dopri_world;
    WorldStepperConfig dopri_stepper_cfg;
    ScenarioBuildResult dopri_build;
    status = build_world_from_scenario_config(
        adaptive_roundtrip_cfg,
        dopri_world,
        dopri_build,
        dopri_stepper_cfg
    );
    if (status != StatusCode::ok) {
        std::println("DOPRI5(4) World Build Status = {}", status_string(status));
        return;
    }

    rk4_cfg.integrator_tr = IntegratorTypeFixed::rk4;
    rk4_cfg.integrator_att = IntegratorTypeFixed::rk4;
    f64 t_span = 1000.0;
    i32 n_steps = 10000;
    f64 dt = t_span / n_steps;

    bool runtime_config_match
        = dopri_stepper_cfg.integrator_tr
              == adaptive_roundtrip_cfg.world_stepper.integrator_tr
        && dopri_stepper_cfg.integrator_att
              == adaptive_roundtrip_cfg.world_stepper.integrator_att
        && dopri_stepper_cfg.adaptive.use_substeps
              == adaptive_roundtrip_cfg.world_stepper.adaptive.use_substeps
        && same_adaptive_opts(
            dopri_stepper_cfg.adaptive.opts,
            adaptive_roundtrip_cfg.world_stepper.adaptive.opts
        )
        && dopri_stepper_cfg.step_tr == adaptive_roundtrip_cfg.world_stepper.step_tr
        && dopri_stepper_cfg.step_att == adaptive_roundtrip_cfg.world_stepper.step_att
        && dopri_stepper_cfg.substeps == adaptive_roundtrip_cfg.world_stepper.substeps
        && dopri_stepper_cfg.ticks == adaptive_roundtrip_cfg.world_stepper.ticks
        && dopri_stepper_cfg.dt_scale == adaptive_roundtrip_cfg.world_stepper.dt_scale
        && dopri_stepper_cfg.paused == adaptive_roundtrip_cfg.world_stepper.paused;
    if (!runtime_config_match) {
        std::println("Adaptive Runtime Config Match = false");
        return;
    }

    WorldStepperWorkspace rk4_wksp;
    WorldStepperStats rk4_stats;
    StatusCode rk4_status = StatusCode::ok;
    auto rk4_start = std::chrono::high_resolution_clock::now();
    for (i32 i = 0; i < n_steps; ++i) {
        WorldStepResult step_result = step_world(rk4_world, dt, rk4_cfg, rk4_wksp);
        rk4_stats += step_result.stats;
        rk4_status = step_result.status;
        if (rk4_status != StatusCode::ok) break;
    }
    auto rk4_stop = std::chrono::high_resolution_clock::now();
    f64 rk4_runtime_ms
        = std::chrono::duration<f64, std::milli>(rk4_stop - rk4_start).count();

    WorldStepperWorkspace dopri_wksp;
    auto dopri_start = std::chrono::high_resolution_clock::now();
    WorldStepResult dopri_result = step_world(
        dopri_world,
        t_span,
        dopri_stepper_cfg,
        dopri_wksp
    );
    auto dopri_stop = std::chrono::high_resolution_clock::now();
    f64 dopri_runtime_ms
        = std::chrono::duration<f64, std::milli>(dopri_stop - dopri_start).count();

    if (rk4_status != StatusCode::ok || dopri_result.status != StatusCode::ok) {
        std::println("RK4 Status = {}", status_string(rk4_status));
        std::println("DOPRI5(4) Status = {}", status_string(dopri_result.status));
        std::println("DOPRI5(4) Final Time = {}", dopri_result.t);
        std::println(
            "DOPRI5(4) Attempts = {}",
            dopri_result.stats.adaptive.attempted_steps
        );
        std::println(
            "DOPRI5(4) Derivative Evaluations = {}",
            dopri_result.stats.adaptive.deriv_evals
        );
        return;
    }

    auto quaternion_error = [](const vec4d& q1, const vec4d& q2) {
        return std::min((q1 - q2).norm(), (q1 + q2).norm());
    };

    auto print_body_errors = [&](const string& config_id) {
        auto rk4_id_it = rk4_build.body_ids.find(config_id);
        auto dopri_id_it = dopri_build.body_ids.find(config_id);
        if (rk4_id_it == rk4_build.body_ids.end()
            || dopri_id_it == dopri_build.body_ids.end()) {
            std::println("{} Body Mapping Not Found", config_id);
            return;
        }

        const Body* rk4_body = rk4_world.body(rk4_id_it->second);
        const Body* dopri_body = dopri_world.body(dopri_id_it->second);
        if (rk4_body == nullptr || dopri_body == nullptr) {
            std::println("{} Body Lookup Failed", config_id);
            return;
        }

        std::println("{} State Errors", config_id);
        std::println("    Position = {}", (rk4_body->x_tr.r - dopri_body->x_tr.r).norm());
        std::println("    Velocity = {}", (rk4_body->x_tr.v - dopri_body->x_tr.v).norm());
        std::println(
            "    Quaternion = {}",
            quaternion_error(rk4_body->x_att.q, dopri_body->x_att.q)
        );
        std::println(
            "    Angular Velocity = {}",
            (rk4_body->x_att.w - dopri_body->x_att.w).norm()
        );
    };

    print_body_errors("earth");
    print_body_errors("moon");
    print_body_errors("sat1");

    std::println("Adaptive Scenario Round Trip");
    std::println("    Stepper Config Match = {}", roundtrip_match);
    std::println("    Runtime Config Match = {}", runtime_config_match);

    std::println("RK4 Stats");
    std::println("    Status = {}", status_string(rk4_status));
    std::println("    Final Time = {}", rk4_world.t_sim());
    std::println("    Simulated Time Advanced = {}", rk4_stats.dt_sim_advanced);
    std::println("    Ticks Completed = {}", rk4_stats.ticks_completed);
    std::println("    Substeps Completed = {}", rk4_stats.substeps_completed);
    i64 rk4_deriv_evals = i64(rk4_stats.substeps_completed) * i64(rk4_tableau.c.size());
    std::println("    Derivative Evaluations = {}", rk4_deriv_evals);
    std::println("    Runtime (ms) = {}", rk4_runtime_ms);

    std::println("DOPRI5(4) Stats");
    std::println("    Status = {}", status_string(dopri_result.status));
    std::println("    Final Time = {}", dopri_result.t);
    std::println("    Simulated Time Advanced = {}", dopri_result.stats.dt_sim_advanced);
    std::println("    Ticks Completed = {}", dopri_result.stats.ticks_completed);
    std::println("    Substeps Completed = {}", dopri_result.stats.substeps_completed);
    std::println("    Attempts = {}", dopri_result.stats.adaptive.attempted_steps);
    std::println("    Accepted Steps = {}", dopri_result.stats.adaptive.accepted_steps);
    std::println("    Rejected Steps = {}", dopri_result.stats.adaptive.rejected_steps);
    std::println(
        "    Derivative Evaluations = {}",
        dopri_result.stats.adaptive.deriv_evals
    );
    std::println("    Final Error Norm = {}", dopri_result.final_error_norm);
    std::println(
        "    Minimum Accepted dt = {}",
        dopri_result.stats.adaptive.min_accepted_dt
    );
    std::println(
        "    Maximum Accepted dt = {}",
        dopri_result.stats.adaptive.max_accepted_dt
    );
    std::println(
        "    Final Accepted dt = {}",
        dopri_result.stats.adaptive.final_accepted_dt
    );
    std::println("    Runtime (ms) = {}", dopri_runtime_ms);
}

void run_world_tableau_rk4_diag() {
    print_diag_title("World Tableau RK4 vs Legacy RK4");

    ScenarioConfig scenario_cfg;
    StatusCode status
        = load_scenario_json(pwd + "/scenarios/earth_moon_sat_demo.json", scenario_cfg);
    if (status != StatusCode::ok) {
        std::println("Scenario Load Status = {}", status_string(status));
        return;
    }

    status = validate_scenario_config(scenario_cfg);
    if (status != StatusCode::ok) {
        std::println("Scenario Validation Status = {}", status_string(status));
        return;
    }

    World legacy_world;
    WorldStepperConfig legacy_cfg;
    ScenarioBuildResult legacy_build;
    status = build_world_from_scenario_config(
        scenario_cfg,
        legacy_world,
        legacy_build,
        legacy_cfg
    );
    if (status != StatusCode::ok) {
        std::println("Legacy RK4 World Build Status = {}", status_string(status));
        return;
    }

    World tableau_world;
    WorldStepperConfig tableau_cfg;
    ScenarioBuildResult tableau_build;
    status = build_world_from_scenario_config(
        scenario_cfg,
        tableau_world,
        tableau_build,
        tableau_cfg
    );
    if (status != StatusCode::ok) {
        std::println("Tableau RK4 World Build Status = {}", status_string(status));
        return;
    }

    legacy_cfg.integrator_tr = IntegratorTypeFixed::rk4;
    legacy_cfg.integrator_att = IntegratorTypeFixed::rk4;
    tableau_cfg.integrator_tr = IntegratorTypeFixed::rk4;
    tableau_cfg.integrator_att = IntegratorTypeFixed::rk4;

    f64 t_span = 1000.0;
    i32 n_steps = 10000;
    f64 dt = t_span / n_steps;

    WorldStepperWorkspace legacy_wksp;
    WorldStepperStats legacy_stats;
    StatusCode legacy_status = StatusCode::ok;
    auto legacy_start = std::chrono::high_resolution_clock::now();
    for (i32 i = 0; i < n_steps; ++i) {
        WorldStepResult step_result
            = step_world_legacy(legacy_world, dt, legacy_cfg, legacy_wksp);
        legacy_stats += step_result.stats;
        legacy_status = step_result.status;
        if (legacy_status != StatusCode::ok) break;
    }
    auto legacy_stop = std::chrono::high_resolution_clock::now();
    f64 legacy_runtime_ms
        = std::chrono::duration<f64, std::milli>(legacy_stop - legacy_start).count();

    WorldStepperWorkspace tableau_wksp;
    WorldStepperStats tableau_stats;
    StatusCode tableau_status = StatusCode::ok;
    auto tableau_start = std::chrono::high_resolution_clock::now();
    for (i32 i = 0; i < n_steps; ++i) {
        WorldStepResult step_result
            = step_world(tableau_world, dt, tableau_cfg, tableau_wksp);
        tableau_stats += step_result.stats;
        tableau_status = step_result.status;
        if (tableau_status != StatusCode::ok) break;
    }
    auto tableau_stop = std::chrono::high_resolution_clock::now();
    f64 tableau_runtime_ms
        = std::chrono::duration<f64, std::milli>(tableau_stop - tableau_start).count();

    if (legacy_status != StatusCode::ok || tableau_status != StatusCode::ok) {
        std::println("Legacy RK4 Status = {}", status_string(legacy_status));
        std::println("Tableau RK4 Status = {}", status_string(tableau_status));
        return;
    }

    auto quaternion_error = [](const vec4d& q1, const vec4d& q2) {
        return std::min((q1 - q2).norm(), (q1 + q2).norm());
    };

    auto print_body_errors = [&](const string& config_id) {
        auto legacy_id_it = legacy_build.body_ids.find(config_id);
        auto tableau_id_it = tableau_build.body_ids.find(config_id);
        if (legacy_id_it == legacy_build.body_ids.end()
            || tableau_id_it == tableau_build.body_ids.end()) {
            std::println("{} Body Mapping Not Found", config_id);
            return;
        }

        const Body* legacy_body = legacy_world.body(legacy_id_it->second);
        const Body* tableau_body = tableau_world.body(tableau_id_it->second);
        if (legacy_body == nullptr || tableau_body == nullptr) {
            std::println("{} Body Lookup Failed", config_id);
            return;
        }

        std::println("{} State Errors", config_id);
        std::println(
            "    Position = {}",
            (legacy_body->x_tr.r - tableau_body->x_tr.r).norm()
        );
        std::println(
            "    Velocity = {}",
            (legacy_body->x_tr.v - tableau_body->x_tr.v).norm()
        );
        std::println(
            "    Quaternion = {}",
            quaternion_error(legacy_body->x_att.q, tableau_body->x_att.q)
        );
        std::println(
            "    Angular Velocity = {}",
            (legacy_body->x_att.w - tableau_body->x_att.w).norm()
        );
    };

    print_body_errors("earth");
    print_body_errors("moon");
    print_body_errors("sat1");

    i64 legacy_deriv_evals
        = i64(legacy_stats.substeps_completed) * i64(rk4_tableau.c.size());
    i64 tableau_deriv_evals
        = i64(tableau_stats.substeps_completed) * i64(rk4_tableau.c.size());

    std::println("Legacy RK4 Stats");
    std::println("    Status = {}", status_string(legacy_status));
    std::println("    Final Time = {}", legacy_world.t_sim());
    std::println("    Simulated Time Advanced = {}", legacy_stats.dt_sim_advanced);
    std::println("    Ticks Completed = {}", legacy_stats.ticks_completed);
    std::println("    Substeps Completed = {}", legacy_stats.substeps_completed);
    std::println("    Derivative Evaluations = {}", legacy_deriv_evals);
    std::println("    Runtime (ms) = {}", legacy_runtime_ms);

    std::println("Tableau RK4 Stats");
    std::println("    Status = {}", status_string(tableau_status));
    std::println("    Final Time = {}", tableau_world.t_sim());
    std::println("    Simulated Time Advanced = {}", tableau_stats.dt_sim_advanced);
    std::println("    Ticks Completed = {}", tableau_stats.ticks_completed);
    std::println("    Substeps Completed = {}", tableau_stats.substeps_completed);
    std::println("    Derivative Evaluations = {}", tableau_deriv_evals);
    std::println("    Runtime (ms) = {}", tableau_runtime_ms);
}

void run_world_embedded_rk_diag() {
    print_diag_title("World Tableau RK4 vs Legacy RK4");

    ScenarioConfig scenario_cfg;
    StatusCode status
        = load_scenario_json(pwd + "/scenarios/earth_moon_sat_demo.json", scenario_cfg);
    if (status != StatusCode::ok) {
        std::println("Scenario Load Status = {}", status_string(status));
        return;
    }

    status = validate_scenario_config(scenario_cfg);
    if (status != StatusCode::ok) {
        std::println("Scenario Validation Status = {}", status_string(status));
        return;
    }

    World reference_world;
    WorldStepperConfig reference_cfg;
    ScenarioBuildResult reference_build;
    status = build_world_from_scenario_config(
        scenario_cfg,
        reference_world,
        reference_build,
        reference_cfg
    );
    if (status != StatusCode::ok) {
        std::println("Tableau RK4 World Build Status = {}", status_string(status));
        return;
    }

    reference_cfg.integrator_tr = IntegratorTypeFixed::rk4;
    reference_cfg.integrator_att = IntegratorTypeFixed::rk4;

    f64 t_span = 1000.0;
    i32 n_steps = 10000;
    f64 dt = t_span / n_steps;

    // set up and propagate with reference: rk4
    WorldStepperWorkspace reference_wksp;
    WorldStepperStats reference_stats;
    StatusCode reference_status = StatusCode::ok;
    auto reference_start = std::chrono::high_resolution_clock::now();
    for (i32 i = 0; i < n_steps; ++i) {
        WorldStepResult step_result
            = step_world(reference_world, dt, reference_cfg, reference_wksp);
        reference_stats += step_result.stats;
        reference_status = step_result.status;
        if (reference_status != StatusCode::ok) break;
    }
    auto reference_stop = std::chrono::high_resolution_clock::now();
    f64 reference_runtime_ms
        = std::chrono::duration<f64, std::milli>(reference_stop - reference_start)
              .count();

    const array<IntegratorTypeAdaptive, 6> methods_embedded{
        IntegratorTypeAdaptive::rkf12,
        IntegratorTypeAdaptive::heuneuler21,
        IntegratorTypeAdaptive::bosha32,
        IntegratorTypeAdaptive::rkf54,
        IntegratorTypeAdaptive::cashkarp54,
        IntegratorTypeAdaptive::dopri54
    };

    const array<IntegratorTypeFixed, 1> methods_fixed{IntegratorTypeFixed::rk5_nystrom};

    auto make_case_world = [&](IntegratorTypeAdaptive type) {
        World case_world;
        WorldStepperConfig case_cfg;
        ScenarioBuildResult case_build;
        StatusCode status = build_world_from_scenario_config(
            scenario_cfg,
            case_world,
            case_build,
            case_cfg
        );
        if (status != StatusCode::ok) {
            std::println(
                "Tableau {} World Build Status = {}",
                integrator_name(type),
                status_string(status)
            );
        }
    };
}
