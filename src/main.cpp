#include "core/body.hpp"
#include "core/dynamics_rotational.hpp"
#include "core/earth_orientation.hpp"
#include "core/entity.hpp"
#include "core/planets.hpp"
#include "core/state.hpp"
#include "core/time.hpp"
#include "core/transform.hpp"
#include "core/world.hpp"
#include "util/constants.hpp"
#include "util/math.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include <chrono>
#include <cmath>
#include <print>

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
        .filename = pwd + "/scratch/assets/leap-seconds.list.txt",
        .lineskips = 85
    };
    EarthPolarMotionParams pmp{
        .filename = pwd + "/scratch/assets/EOP_20u24_C04_one_file_1962-now.txt",
        .lineskips = 6,
        .model = EarthPolarMotionModel::IAU2000A,
        .approx = false
    };
    EarthNutationParams enp{
        .filename = pwd + "/scratch/assets/nut_IAU1980.dat.txt",
        .lineskips = 3,
        .precision = 106,
        .approx = false,
        .model = EarthNutationModel::IAU1980
    };
    EarthOrientationParams eop{.leap_seconds = lsp, .nutation = enp, .polar_motion = pmp};
    bool eop_ok = load_all_eop(eop);

    JulianDate jd; // j2000 utc
    get_time_offsets(jd, eop);

    auto start = std::chrono::high_resolution_clock::now();
    mat3d R = mat3d1;
    R = rot_eci_to_ecef(
        jd,
        EarthFrame(1),
        EarthFrame(7),
        eop,
        eop.offsets,
        TimeScale::utc
    );
    auto stop = std::chrono::high_resolution_clock::now();

    std::println("EOP Loaded: {}", eop_ok);
    std::cout << "ECI to ECEF (prolly wrong):\n" << R << std::endl;
    auto duration = (stop - start);
    print_chrono(duration, UTime::millisecond);
    return 0;
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
        std::string(PROJECT_ROOT) + "/scratch/assets/egm2008_120.txt",
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
        std::string(PROJECT_ROOT) + "/scratch/assets/egm2008_120.txt",
        earth->C,
        earth->S,
        earth->degree,
        earth->order
    );

    // Spherical Harmonics, zonal only
    vec3d llh = vec3d{-45.0, 0, 7000.0};       // in body fixed frame
    vec3d r_body = detic_to_body(llh, *earth); // in body fixed frame
    sat->x_tr.r
        = ep_rotate_fast_passive(ep_conj(earth->x_att.q), r_body); // back to inertial
    earth->degree = degree;
    earth->order = 0;
    vec3d a_sphh_zonal_long1 = world.gravity_accel_from(sat_id, earth_id);

    // Spherical Harmonics, zonal only, different longitude
    llh = vec3d{-45.0, 57.0, 7000.0};
    r_body = detic_to_body(llh, *earth);
    sat->x_tr.r = ep_rotate_fast_passive(ep_conj(earth->x_att.q), r_body);
    vec3d a_sphh_zonal_long2 = world.gravity_accel_from(sat_id, earth_id);

    // Spherical Harmonics
    llh = vec3d{-45.0, 0.0, 7000.0};
    r_body = detic_to_body(llh, *earth);
    sat->x_tr.r = ep_rotate_fast_passive(ep_conj(earth->x_att.q), r_body);
    earth->order = order;
    vec3d a_sphh_long1 = world.gravity_accel_from(sat_id, earth_id);

    // Spherical Harmonics, different longitude
    llh = vec3d{-45.0, 57.0, 7000.0};
    r_body = detic_to_body(llh, *earth);
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
