#include "core/diagnostics.hpp"
#include "core/body.hpp"
#include "core/dynamics_rotational.hpp"
#include "core/entity.hpp"
#include "core/estimation_batch.hpp"
#include "core/estimation_recursive.hpp"
#include "core/ingest.hpp"
#include "core/integrator.hpp"
#include "core/measurement.hpp"
#include "core/measurement_world.hpp"
#include "core/observations.hpp"
#include "core/od_dynamics.hpp"
#include "core/orbit_determination.hpp"
#include "core/orbital_elements.hpp"
#include "core/planets.hpp"
#include "core/state.hpp"
#include "core/station_geometry.hpp"
#include "core/time.hpp"
#include "core/transform.hpp"
#include "core/world.hpp"
#include "core/world_stepper.hpp"

#include "graphics/renderer.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include "graphics/raygen.hpp"
#include "graphics/rdraw.hpp"
#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"

#include <memory>
#include <random>

const std::string pwd = std::string(PROJECT_ROOT);

void print_diag_title(const std::string& title) {
    std::string line = "-----------------------------------------------------------";
    size_t width = line.size();

    std::string withSpace = title + " ";
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
    std::println("------------------------------------------------------------");
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
    std::println("------------------------------------------------------------");
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
    std::println("------------------------------------------------------------");
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
    std::println("------------------------------------------------------------");
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
    std::println("------------------------------------------------------------");
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

    ODDynamicsConfig cfg{.model = ODDynamicsModel::two_body, .mu = body.mu};
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

    vec6d x0_vec = statetr_to_vec6(x0);
    vec6d xf_vec = statetr_to_vec6(xf);

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
        StateTr x0_pert = vec6_to_statetr(x0_pert_vec);
        StateTr xf_pert = propagate_tr_od(0.0, x0_pert, dt, n_steps, cfg);
        vec6d xf_pert_vec = statetr_to_vec6(xf_pert);
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

void run_measurement_jacobian_diag() {
    MeasurementContext ctx;
    // zero observer state
    ctx.x_observer.r = vec3d{0.0, 0.0, 0.0};
    ctx.x_observer.v = vec3d{0.0, 0.0, 0.0};

    // arbitrary target state
    ctx.x_target.r = vec3d{7000.0, 1000.0, 1300.0};
    ctx.x_target.v = vec3d{-0.5, 7.2, 1.0};

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

struct EarthStationSatScenario {
    World world;
    EntityId earth_id = kInvalidEntityId;
    EntityId stat_id = kInvalidEntityId;
    EntityId sat_id = kInvalidEntityId;
};
static EarthStationSatScenario make_earth_station_sat_scenario(const vec3d& station_llh) {
    EarthStationSatScenario scenario;
    scenario.earth_id = wgs84(scenario.world);
    scenario.stat_id = scenario.world.spawn_station();
    scenario.sat_id = scenario.world.spawn_satellite();

    Celestial* earth = scenario.world.celestial(scenario.earth_id);

    // earth
    earth->x_att.q = vec4d{0.0, 0.0, 0.0, 1.0};
    earth->x_att.w = vec3d{0.0, 0.0, 7.292115000000000e-05};

    // station
    scenario.world
        .set_stat_anchor_detic(scenario.stat_id, scenario.earth_id, station_llh);

    return scenario;
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
        input.dyn_config.model = ODDynamicsModel::two_body;
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
            StateTr x_obs = world.stat_x_tr_inertial(stat_id);

            Measurement meas;

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
                meas.R = matXd::Identity(2, 2) * sigma_rad * sigma_rad;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_obs);
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
                meas.R = matXd::Identity(1, 1) * sigma_range * sigma_range;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_obs);
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
                meas.R = matXd::Identity(1, 1) * sigma_range_rate * sigma_range_rate;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_obs);
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
                input.observer_states.push_back(x_obs);
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
                input.observer_states.push_back(x_obs);
            }
        }

        // solve
        ODBatchResult result = od_batch_lumve(input);
        f64 initial_err = (statetr_to_vec6(input.x0_guess - x0_truth)).norm();
        f64 final_err = (statetr_to_vec6(result.x0_est - x0_truth)).norm();
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
        std::println("LUMVE Success = {}", result.success);
        std::println("LUMVE Status: {}", od_status_string(result.status));
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
        ODEKFInput input;
        input.initial_filter.x = x0_truth;
        input.initial_filter.t = 0.0;
        input.initial_filter.P = mat6d1;
        input.dyn_config.model = ODDynamicsModel::two_body;
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
            StateTr x_obs = world.stat_x_tr_inertial(stat_id);

            Measurement meas;

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
                meas.R = matXd::Identity(2, 2) * sigma_rad * sigma_rad;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_obs);
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
                meas.R = matXd::Identity(1, 1) * sigma_range * sigma_range;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_obs);
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
                meas.R = matXd::Identity(1, 1) * sigma_range_rate * sigma_range_rate;
                input.measurements.push_back(meas);
                input.observer_states.push_back(x_obs);
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
                input.observer_states.push_back(x_obs);
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
                input.observer_states.push_back(x_obs);
            }
        }

        // run ekf
        ODEKFResult result = od_ekf_offline(input);

        f64 initial_err = (statetr_to_vec6(input.initial_filter.x - x0_truth)).norm();
        f64 final_err = (statetr_to_vec6(result.filter.x - sat->x_tr)).norm();
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
        std::println("MIXED EKF Success = {}", result.success);
        std::println("MIXED EKF Status: {}", od_status_string(result.status));
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

    auto integrator_name = [](IntegratorType integrator) -> std::string {
        switch (integrator) {
        case IntegratorType::rk1: return "RK1";
        case IntegratorType::rk2: return "RK2";
        case IntegratorType::rk2_heun: return "RK2 Heun";
        case IntegratorType::rk2_ralston: return "RK2 Ralston";
        case IntegratorType::rk3: return "RK3";
        case IntegratorType::rk4: return "RK4";
        }
        return "Unknown";
    };

    auto run_case = [&](IntegratorType integrator_tr) {
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
        cfg.step_translation = true;
        cfg.step_attitude = true;
        cfg.substeps = 1;
        cfg.ticks = 10;
        cfg.time_scale = 1.0 / cfg.ticks;
        cfg.integrator_tr = integrator_tr;
        cfg.integrator_att = IntegratorType::rk4;

        f64 t_span = 100.0;
        f64 t0 = 0.0;
        i32 n_steps = 500;
        f64 dt = t_span / n_steps;
        world.reset_time(t0);

        ODDynamicsConfig ODcfg;
        ODcfg.body_radius = earth->mean_radius;
        ODcfg.integrator = integrator_tr;
        ODcfg.mu = earth->mu;
        ODcfg.model = ODDynamicsModel::two_body;

        // "truth"
        i32 n_ref_steps = n_steps * cfg.ticks * cfg.substeps;
        StateTr x_OD = propagate_tr_od(t0, x0, t_span, n_ref_steps, ODcfg);

        // world stepper
        WorldStepperStats stats;
        stats.success = true;
        for (i32 i = 0; i < n_steps; ++i) {
            stats += step_world(world, dt, cfg);
            if (!stats.success) break;
        }

        StateTr x_err = x_OD - sat->x_tr;
        f64 pos_err = x_err.r.norm();
        f64 vel_err = x_err.v.norm();
        f64 pos_rel_err = pos_err / x_OD.r.norm();
        f64 vel_rel_err = vel_err / x_OD.v.norm();

        std::println("{} Success = {}", integrator_name(integrator_tr), stats.success);
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

    run_case(IntegratorType::rk1);
    run_case(IntegratorType::rk2);
    run_case(IntegratorType::rk2_heun);
    run_case(IntegratorType::rk2_ralston);
    run_case(IntegratorType::rk3);
    run_case(IntegratorType::rk4);
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
        cfg.step_translation = false;
        cfg.step_attitude = true;
        cfg.integrator_tr = IntegratorType::rk1;
        cfg.integrator_att = IntegratorType::rk4;
        cfg.substeps = 1;
        cfg.ticks = 1;
        cfg.time_scale = 1.0;

        WorldStepperStats stats;
        stats.success = true;
        for (i32 i = 0; i < n_steps; ++i) {
            stats += step_world(world, dt, cfg);
            if (!stats.success) break;
        }

        f64 theta = spin_rate * world.t_sim();
        vec4d q_expected{0.0, 0.0, std::sin(theta / 2.0), std::cos(theta / 2.0)};

        std::println("Simple Spin Success = {}", stats.success);
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
        cfg.step_translation = true;
        cfg.step_attitude = true;
        cfg.substeps = 1;
        cfg.ticks = 10;
        cfg.time_scale = 1.0 / cfg.ticks;
        cfg.integrator_tr = IntegratorType::rk1;
        cfg.integrator_att = IntegratorType::rk4;

        f64 t_span = 100.0;
        f64 t0 = 0.0;
        i32 n_steps = 1000;
        f64 dt = t_span / n_steps;
        world.reset_time(t0);

        WorldStateSnapshot world_snapshot = world.capture_checkpoint();

        WorldStepperStats stats;
        auto prop = [&]() {
            stats.success = true;
            for (i32 i = 0; i < n_steps; ++i) {
                stats += step_world(world, dt, cfg);
                if (!stats.success) break;
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
    cfg.step_translation = true;
    cfg.step_attitude = true;
    cfg.substeps = 10;
    cfg.ticks = 1;
    cfg.time_scale = 1.0;
    cfg.integrator_tr = IntegratorType::rk4;
    cfg.integrator_att = IntegratorType::rk4;

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
    auto prop = [&](bool write_csv, bool print_progress) {
        stats.success = true;
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
            stats += step_world(world, dt, cfg);
            if (!stats.success) break;
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

    std::println("Moving Source Success = {}", stats_moving.success);
    std::println("Fixed Source Success = {}", stats_fixed.success);
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
        pwd + "/assets/EGM2008.gfc.txt",
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
    cfg.step_translation = true;
    cfg.step_attitude = true;
    cfg.substeps = 1;
    cfg.ticks = 1;
    cfg.time_scale = 1.0;
    cfg.integrator_tr = IntegratorType::rk4;
    cfg.integrator_att = IntegratorType::rk4;

    f64 t_span = 10000;
    f64 t0 = 0.0;
    i32 n_steps = 10000;
    f64 dt = t_span / n_steps;
    world.reset_time(t0);

    WorldStateSnapshot world_snapshot = world.capture_checkpoint();

    WorldStepperStats stats;
    auto prop = [&]() {
        stats.success = true;
        for (i32 i = 0; i < n_steps; ++i) {
            stats += step_world(world, dt, cfg);
            if (!stats.success) break;
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

        world.restore_checkpoint_state(world_snapshot);
        earth->gravity_model = gravity_model;
        earth->attitude_model = CelestialAttitudeModel::simple_spin;
        earth->propagate_att = true;
        prop();
        StateTr xf_rot = sat->x_tr;
        WorldStepperStats stats_rot = stats;

        StateTr x_err = xf_fixed - xf_rot;

        std::println("{} Fixed Success = {}", name, stats_fixed.success);
        std::println("{} Rotating Success = {}", name, stats_rot.success);
        std::println("{} Position Difference = {}", name, x_err.r.norm());
        std::println("{} Velocity Difference = {}", name, x_err.v.norm());
        std::println();
    };

    run_case(GravityModel::pointmass, "Pointmass");
    run_case(GravityModel::spherical_harmonics, "Spherical Harmonics");
    run_case(GravityModel::zonal, "Zonal");
}

void run_render_pipeline_diag() {
    World world;

    // Earth
    EntityId earth_id = wgs84(world);
    Celestial* earth = world.celestial(earth_id);
    earth->name = "Earth";
    earth->gravity_model = GravityModel::spherical_harmonics;
    earth->degree = 4;
    earth->order = 4;
    bool gfc_ok = read_gfc(
        pwd + "/assets/EGM2008.gfc.txt",
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
        svec<std::unique_ptr<Satellite>> tle_sats;
        svec<i32> sat_nums = {25544, 11, 5};
        TLEStatus tle_status = read_tle_data_satnums(
            pwd + "/assets/tle_all.txt",
            tles,
            sat_nums,
            tle_opts
        );
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
    cfg.step_translation = true;
    cfg.step_attitude = true;
    cfg.substeps = 1;
    cfg.ticks = 100;
    cfg.time_scale = 2.0;
    cfg.integrator_tr = IntegratorType::rk4;
    cfg.integrator_att = IntegratorType::rk4;

    f64 t0 = 0.0;
    f64 dt = 1.0 / cfg.ticks * cfg.time_scale;
    world.reset_time(t0);
    WorldStepperStats stats;
    stats.success = true;
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
        BeginDrawing();
        BeginMode3D(camera);
        ClearBackground(Color({30, 30, 30, 255}));

        RenderSceneSnapshot scene = build_render_scene_snapshot(world, assets);
        render_scene_snapshot(scene, sphere_model, cube_model, cylinder_model);

        if (IsKeyPressed(KEY_SPACE)) {
            paused = !paused;
        }
        if (!paused) {
            stats += step_world(world, dt, cfg, wksp);
            if (!stats.success) {
                std::println("Simulation Failure");
                break;
            }
        }
        svec<EntityId> stat_ids = world.station_ids();
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
        pwd + "/assets/EGM2008.gfc.txt",
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
    cfg.step_translation = true;
    cfg.step_attitude = true;
    cfg.substeps = 1;
    cfg.ticks = 1;
    cfg.time_scale = 1.0;
    cfg.integrator_tr = IntegratorType::rk4;
    cfg.integrator_att = IntegratorType::rk4;

    f64 t_span = 10000;
    f64 t0 = 0.0;
    i32 n_steps = 10000;
    f64 dt = t_span / n_steps;
    world.reset_time(t0);

    WorldStepperStats stats;
    stats.success = true;

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

    for (i32 i = 0; i < n_steps; ++i) {
        stats += step_world(world, dt, cfg, wksp);
        if (!stats.success) {
            std::println("Integration Failed");
            break;
        }
    }

    print_workspace("WORLD WKSP Final");
    std::println("WORLD WKSP Success = {}", stats.success);
    std::println("WORLD WKSP Ticks Completed = {}", stats.ticks_completed);
    std::println("WORLD WKSP Substeps Completed = {}", stats.substeps_completed);
    std::println("WORLD WKSP Time Advanced = {}", stats.dt_sim_advanced);
    std::println("WORLD WKSP Final Time = {}", world.t_sim());
    std::println("WORLD WKSP Expected Time = {}", t0 + t_span);
}
