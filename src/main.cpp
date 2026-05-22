#include "core/body.hpp"
#include "core/diagnostics.hpp"
#include "core/earth_orientation.hpp"
#include "core/entity.hpp"
#include "core/planets.hpp"
#include "core/state.hpp"
#include "core/time.hpp"
#include "core/transform.hpp"
#include "core/world.hpp"
#include "graphics/raygen.hpp"
#include "util/constants.hpp"
#include "util/math.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include <chrono>
#include <iostream>
#include <print>

const std::string pwd = std::string(PROJECT_ROOT);

int main() {
    World world;

    // Earth (at origin, tilted, spinning, zonal)
    EntityId earth_id = wgs84(world);
    Celestial* earth = world.celestial(earth_id);
    earth->gravity_model = GravityModel::zonal;
    earth->degree = 2;
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
    stat->r_body_BCBF = {earth->mean_radius, 0, 0};
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
    // std::println("EOP Loaded: {}", eop_ok);

    JulianDate jd; // j2000 utc
    get_time_offsets(jd, eop);

    // run_earth_rot_diag(jd, eop);
    // run_station_obs_geom_diag(*earth, jd, eop);
    // run_rv_coe_diag(*earth);
    // run_radec_diag();
    // run_iod_diag(*earth);
    // run_od_prop_diag(*earth);
    // run_measurement_jacobian_diag();
    // run_batch_od_diag();
    // run_checkpoint_diag();
    // run_station_anchor_diag();
    // run_world_measurement_diag();
    // run_world_stepper_diag();
    // run_body_fixed_gravity_timing_diag();
    // run_moving_source_world_diag();
    // run_staged_attitude_gravity_diag();
    // run_world_workspace_diag();
    // run_tle_status_reader_diag();
    // run_world_measurement_context_diag();
    // run_batch_od_diag();
    // run_ekf_world_diag();
    // run_iod_lumve_ekf_init_diag();
    // run_od_zonal_jacobian_diag();

    // Current diagnostic(s)
    std::println("-----------------------------------------------------------");
    auto start = std::chrono::high_resolution_clock::now();
    // run_render_pipeline_diag();
    // run_world_ekf_step_diag();
    // print_diag_title("");
    // run_ekf_prediction_only_diag();
    // run_realtime_ekf_world_update_diag();
    // run_station_instrument_diag();
    // run_world_history_diag();
    run_world_history_ekf_diag();
    auto stop = std::chrono::high_resolution_clock::now();
    std::println("-----------------------------------------------------------");

    auto duration = (stop - start);
    print_chrono(duration, UTime::millisecond);

    return 0;
}
