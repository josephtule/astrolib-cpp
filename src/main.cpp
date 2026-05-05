#include "core/body.hpp"
#include "core/dynamics_rotational.hpp"
#include "core/earth_orientation.hpp"
#include "core/entity.hpp"
#include "core/estimation_batch.hpp"
#include "core/measurement.hpp"
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
#include "util/constants.hpp"
#include "util/math.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include <chrono>
#include <cmath>
#include <limits>
#include <print>

void run_batch_od_radec_diag(const Celestial& body);
void run_batch_od_pos_diag(const Celestial& body);
void run_batch_od_posvel_diag(const Celestial& body);
void run_gravity_diag(
    World& world,
    EntityId earth_id,
    EntityId urath_id,
    EntityId sat_id,
    EntityId stat_id
);
void run_station_geo_diag(
    World& world,
    EntityId earth_id,
    EntityId urath_id,
    EntityId sat_id,
    EntityId stat_id
);
void run_zonal_orientation_diag(
    World& world,
    EntityId earth_id,
    EntityId urath_id,
    EntityId sat_id,
    EntityId stat_id
);
void run_spherical_harmonics_diag(
    World& world,
    EntityId earth_id,
    EntityId urath_id,
    EntityId sat_id,
    EntityId stat_id
);
void run_sphh_longitude_diag(
    World& world,
    EntityId earth_id,
    EntityId urath_id,
    EntityId sat_id,
    EntityId stat_id
);
void run_epkde(
    World& world,
    EntityId earth_id,
    EntityId urath_id,
    EntityId sat_id,
    EntityId stat_id
);
void run_time_diag();
void run_earth_rot_diag(const JulianDate& jd, const EarthOrientationParams& eop);
void run_station_obs_geom_diag(
    Celestial& earth,
    const JulianDate& jd,
    const EarthOrientationParams& eop
);
void run_rv_coe_diag(const Celestial& body);
void run_radec_diag();
void run_iod_diag(const Celestial& body);
void run_od_prop_diag(const Celestial& body);
void run_measurement_jacobian_diag();

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int main() {
    std::string pwd = std::string(PROJECT_ROOT);
    World world;

    // Earth (at origin, tilted, spinning, zonal)
    EntityId earth_id = wgs84(world);
    Celestial* earth = world.celestial(earth_id);
    earth->gravity_model = GravityModel::zonal;
    earth->degree = 2;
    earth->use_simple_spin = true;
    earth->set_spin_rate(earth->spin_rate); // spin rate in rad/s
    mat3d earth_tilt_dcm = rotX(23.44 * deg_to_rad);
    earth->x_att.q = dcm_to_ep(earth_tilt_dcm);

    // Urath (offset, inertial aligned, fixed, pointmass)
    EntityId urath_id = wgs84(world);
    Celestial* urath = world.celestial(urath_id);
    urath->x_tr.r = {10000.0, 10000.0, 10000.0};
    urath->gravity_model = GravityModel::pointmass;
    urath->degree = 0;

    // Satellite
    EntityId sat_id = world.spawn_satellite();
    Satellite* sat = world.satellite(sat_id);
    sat->x_tr.r = {7000.0, 0.0, 0.0};

    // Station
    EntityId stat_id = world.spawn_station();
    Station* stat = world.station(stat_id);
    stat->anchored = true;
    stat->anchor_id = earth_id;
    stat->r_body = {earth->mean_radius, 0, 0};
    // earth rotated only about x-axis, station still on x-axis

    // run_gravity_diag(world, earth_id, urath_id, sat_id, stat_id);
    // run_station_geo_diag(world, earth_id, urath_id, sat_id, stat_id);
    // run_zonal_orientation_diag(world, earth_id, urath_id, sat_id, stat_id);
    // run_spherical_harmonics_diag(world, earth_id, urath_id, sat_id, stat_id);
    // run_sphh_longitude_diag(world, earth_id, urath_id, sat_id, stat_id);
    // run_epkde(world, earth_id, urath_id, sat_id, stat_id);

    LeapSecondParams lsp{
        .filename = pwd + "/assets/leap-seconds.list.txt",
        .lineskips = 85
    };
    EarthPolarMotionParams pmp{
        .filename = pwd + "/assets/EOP_20u24_C04_one_file_1962-now.txt",
        .lineskips = 6,
        .model = EarthPolarMotionModel::IAU2000A,
        .approx = false
    };
    EarthNutationParams enp{
        .filename = pwd + "/assets/nut_IAU1980.dat.txt",
        .lineskips = 3,
        .precision = 106,
        .approx = false,
        .model = EarthNutationModel::IAU1980
    };
    EarthOrientationParams eop{.leap_seconds = lsp, .nutation = enp, .polar_motion = pmp};
    bool eop_ok = load_all_eop(eop);
    std::println("EOP Loaded: {}", eop_ok);

    JulianDate jd; // j2000 utc
    get_time_offsets(jd, eop);

    // run_earth_rot_diag(jd, eop);
    // run_station_obs_geom_diag(*earth, jd, eop);
    // run_rv_coe_diag(*earth);
    // run_radec_diag();
    // run_iod_diag(*earth);
    // run_od_prop_diag(*earth);
    // run_measurement_jacobian_diag();

    std::println("-----------------------------------------------------------");
    auto start = std::chrono::high_resolution_clock::now();
    run_batch_od_radec_diag(*earth);
    run_batch_od_pos_diag(*earth);
    run_batch_od_posvel_diag(*earth);
    auto stop = std::chrono::high_resolution_clock::now();
    std::println("-----------------------------------------------------------");

    auto duration = (stop - start);
    print_chrono(duration, UTime::millisecond);

    // NOTE: quick todo list
    // - get measurement model started
    // - get IOD running
    // - get 2 body propagator (for estimation and perturbations, sim will use nbody
    // integrator (coupled simulation))
    // - create jacobians
    // - get batch estimator running (LUMVE)
    // - get ekf running

    return 0;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

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
    earth->x_att.q = q_default;
    earth->gravity_model = GravityModel::pointmass;
    vec3d a_point = world.gravity_accel_from(sat_id, earth_id);

    // pointmass (tilted)
    earth->x_att.q = earth_q0;
    vec3d a_point_tilt = world.gravity_accel_from(sat_id, earth_id);

    // zonal (no orientation)
    earth->gravity_model = GravityModel::zonal;
    earth->x_att.q = q_default;
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
    read_egm2008(
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
    read_egm2008(
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

    vec4d q = q_default;
    vec3d w = earth->spin_rate * axis_z;

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
    f64 tol = tol_strict
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
    StateTr x2_prop = propagate_tr_od_rk4(0.0, x1, dt, n_steps, cfg);

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

    VarStateTr yf = propagate_var_tr_od_rk4(0.0, y0, dt, n_steps, cfg);
    StateTr xf = propagate_tr_od_rk4(0.0, x0, dt, n_steps, cfg);

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
        StateTr x0_pert = vec6_to_statetr(x0_pert_vec);
        StateTr xf_pert = propagate_tr_od_rk4(0.0, x0_pert, dt, n_steps, cfg);
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

void run_batch_od_radec_diag(const Celestial& body) {
    StateTr x0_truth;
    x0_truth.r = vec3d{7000.0, 1000.0, 1300.0};
    x0_truth.v = vec3d{-0.5, 7.2, 1.0};

    ODBatchInput input;
    input.x0_guess = x0_truth;
    input.dyn_config.model = ODDynamicsModel::two_body;
    input.dyn_config.mu = body.mu;
    input.t0 = 0.0;
    input.max_iters = 10;
    input.prop_steps = 200;
    input.tol_dx = 1e-6;
    input.tol_residual = 1e-8;

    // add perturbations
    input.x0_guess.r += vec3d{1.0, -1.0, 0.5};
    input.x0_guess.v += vec3d{1e-3, -1e-3, 5e-4};

    f64 t_meas0 = 0.0;
    f64 t_measf = 600.0;
    f64 N_meas = 30;
    vecXd t_meas = vecXd::LinSpaced(N_meas, t_meas0, t_measf);
    f64 sigma_rad = 1e-5;
    for (i32 i = 0; i < N_meas; ++i) {
        // propagate to get measurements
        StateTr x_truth_i = propagate_tr_od_rk4(
            0.0,
            x0_truth,
            t_meas(i),
            input.prop_steps,
            input.dyn_config
        );
        StateTr x_obs;
        // use stationary observer for now
        x_obs.r = vec3d{body.mean_radius, 0.0, 0.0};
        x_obs.v = vec3d0;
        input.observer_states.push_back(x_obs);

        // get measurement
        MeasurementContext ctx;
        ctx.x_target = x_truth_i;
        ctx.x_observer = x_obs;
        Measurement meas;
        meas.t = t_meas(i);
        meas.type = ObservationType::radec;
        meas.z = predict_measurement(ObservationType::radec, ctx);
        meas.R = matXd::Identity(2, 2) * sigma_rad * sigma_rad;
        input.measurements.push_back(meas);
    }

    // solve
    ODBatchResult result = od_batch_lumve(input);
    f64 initial_err = (statetr_to_vec6d(input.x0_guess - x0_truth)).norm();
    f64 final_err = (statetr_to_vec6d(result.x0_est - x0_truth)).norm();
    std::println("RADEC Initial Error = {}", initial_err);
    std::println("RADEC Final Error = {}", final_err);
    std::println("RADEC LUMVE Success = {}", result.success);
    std::println("RADEC LUMVE Status: {}", od_batch_status_string(result.status));
    std::println("RADEC iterations: {}", result.iterations);
    std::println("RADEC Residual Norm = {}", result.residual_norm);
    std::println("RADEC Raw Residual Norm = {}", result.raw_residual_norm);
    std::println("RADEC Delta x Norm = {}", result.dx_norm);
}

void run_batch_od_pos_diag(const Celestial& body) {
    StateTr x0_truth;
    x0_truth.r = vec3d{7000.0, 1000.0, 1300.0};
    x0_truth.v = vec3d{-0.5, 7.2, 1.0};

    ODBatchInput input;
    input.x0_guess = x0_truth;
    input.dyn_config.model = ODDynamicsModel::two_body;
    input.dyn_config.mu = body.mu;
    input.t0 = 0.0;
    input.max_iters = 10;
    input.prop_steps = 200;
    input.tol_dx = 1e-6;
    input.tol_residual = 1e-8;

    // add perturbations
    input.x0_guess.r += vec3d{1.0, -1.0, 0.5};
    input.x0_guess.v += vec3d{1e-3, -1e-3, 5e-4};

    f64 t_meas0 = 0.0;
    f64 t_measf = 600.0;
    f64 N_meas = 30;
    vecXd t_meas = vecXd::LinSpaced(N_meas, t_meas0, t_measf);
    f64 sigma_r = 1e-3;
    for (i32 i = 0; i < N_meas; ++i) {
        // propagate to get measurements
        StateTr x_truth_i = propagate_tr_od_rk4(
            0.0,
            x0_truth,
            t_meas(i),
            input.prop_steps,
            input.dyn_config
        );
        StateTr x_obs;
        // use stationary observer for now
        x_obs.r = vec3d{body.mean_radius, 0.0, 0.0};
        x_obs.v = vec3d0;
        input.observer_states.push_back(x_obs);

        // get measurement
        MeasurementContext ctx;
        ctx.x_target = x_truth_i;
        ctx.x_observer = x_obs;
        Measurement meas;
        meas.t = t_meas(i);
        meas.type = ObservationType::pos;
        meas.z = predict_measurement(ObservationType::pos, ctx);
        meas.R = matXd::Identity(3, 3) * sigma_r * sigma_r;
        input.measurements.push_back(meas);
    }

    // solve
    ODBatchResult result = od_batch_lumve(input);
    f64 initial_err = (statetr_to_vec6d(input.x0_guess - x0_truth)).norm();
    f64 final_err = (statetr_to_vec6d(result.x0_est - x0_truth)).norm();
    std::println("POS Initial Error = {}", initial_err);
    std::println("POS Final Error = {}", final_err);
    std::println("POS LUMVE Success = {}", result.success);
    std::println("POS LUMVE Status: {}", od_batch_status_string(result.status));
    std::println("POS iterations: {}", result.iterations);
    std::println("POS Residual Norm = {}", result.residual_norm);
    std::println("POS Raw Residual Norm = {}", result.raw_residual_norm);
    std::println("POS Delta x Norm = {}", result.dx_norm);
}

void run_batch_od_posvel_diag(const Celestial& body) {
    StateTr x0_truth;
    x0_truth.r = vec3d{7000.0, 1000.0, 1300.0};
    x0_truth.v = vec3d{-0.5, 7.2, 1.0};

    ODBatchInput input;
    input.x0_guess = x0_truth;
    input.dyn_config.model = ODDynamicsModel::two_body;
    input.dyn_config.mu = body.mu;
    input.t0 = 0.0;
    input.max_iters = 10;
    input.prop_steps = 200;
    input.tol_dx = 1e-6;
    input.tol_residual = 1e-8;

    // add perturbations
    input.x0_guess.r += vec3d{1.0, -1.0, 0.5};
    input.x0_guess.v += vec3d{1e-3, -1e-3, 5e-4};

    f64 t_meas0 = 0.0;
    f64 t_measf = 600.0;
    f64 N_meas = 30;
    vecXd t_meas = vecXd::LinSpaced(N_meas, t_meas0, t_measf);
    f64 sigma_r = 1e-3;
    f64 sigma_v = 1e-6;
    for (i32 i = 0; i < N_meas; ++i) {
        // propagate to get measurements
        StateTr x_truth_i = propagate_tr_od_rk4(
            0.0,
            x0_truth,
            t_meas(i),
            input.prop_steps,
            input.dyn_config
        );
        StateTr x_obs;
        // use stationary observer for now
        x_obs.r = vec3d{body.mean_radius, 0.0, 0.0};
        x_obs.v = vec3d0;
        input.observer_states.push_back(x_obs);

        // get measurement
        MeasurementContext ctx;
        ctx.x_target = x_truth_i;
        ctx.x_observer = x_obs;
        Measurement meas;
        meas.t = t_meas(i);
        meas.type = ObservationType::pos_vel;
        meas.z = predict_measurement(ObservationType::pos_vel, ctx);
        meas.R = matXd::Zero(6, 6);
        meas.R.block(0, 0, 3, 3) = mat3d::Identity() * sigma_r * sigma_r;
        meas.R.block(3, 3, 3, 3) = mat3d::Identity() * sigma_v * sigma_v;
        input.measurements.push_back(meas);
    }

    // solve
    ODBatchResult result = od_batch_lumve(input);
    f64 initial_err = (statetr_to_vec6d(input.x0_guess - x0_truth)).norm();
    f64 final_err = (statetr_to_vec6d(result.x0_est - x0_truth)).norm();
    std::println("POSVEL Initial Error = {}", initial_err);
    std::println("POSVEL Final Error = {}", final_err);
    std::println("POSVEL LUMVE Success = {}", result.success);
    std::println("POSVEL LUMVE Status: {}", od_batch_status_string(result.status));
    std::println("POSVEL iterations: {}", result.iterations);
    std::println("POSVEL Residual Norm = {}", result.residual_norm);
    std::println("POSVEL Raw Residual Norm = {}", result.raw_residual_norm);
    std::println("POSVEL Delta x Norm = {}", result.dx_norm);
}

// TODO: consolidate batch OD diags
