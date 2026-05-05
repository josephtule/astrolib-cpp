#pragma once

#include "core/earth_orientation.hpp"
#include "core/entity.hpp"
#include "core/time.hpp"
#include "core/world.hpp"

void run_batch_od_diag(const Celestial& body);
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
void run_ekf_mixed_measurement_diag(const Celestial& body);
void run_checkpoint_diag();