// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/entity.hpp"

#include <string>

class World;
struct Celestial;
struct EarthOrientationParams;
struct JulianDate;

void print_diag_title(const std::string& title);
void run_batch_od_diag();
void run_gravity_diag(
    World& world,
    EntityId earth_id,
    EntityId urath_id,
    EntityId sat_id,
    EntityId id
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
void run_od_prop_adaptive_diag(const Celestial& body);
void run_measurement_jacobian_diag();
void run_checkpoint_diag();
void run_station_anchor_diag();
void run_world_measurement_diag();
void run_ekf_world_diag();
void run_world_stepper_diag();
void run_world_dopri54_rk4_diag();
void run_world_tableau_rk4_diag();
void run_body_fixed_gravity_timing_diag();
void run_moving_source_world_diag();
void run_tle_status_reader_diag();
void run_staged_attitude_gravity_diag();
void run_render_pipeline_diag();
void run_world_workspace_diag();
void run_world_measurement_context_diag();
void run_iod_lumve_ekf_init_diag();
void run_od_zonal_jacobian_diag();
void run_world_ekf_step_diag();
void run_ekf_prediction_only_diag();
void run_realtime_ekf_world_update_diag();
void run_station_instrument_diag();
void run_world_history_diag();
void run_world_history_ekf_diag();
void run_scenario_loader_diag();
void run_build_world_from_scenario();
void run_scenario_render_diag();
void run_world_embedded_rk_diag();
void run_ephemeris_io_diag();
void run_ephemeris_sampling_diag();
void run_ephemeris_provider_diag();
void run_world_provider_diag();
void run_cspice_provider_diag();
