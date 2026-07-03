#pragma once

#include "core/orbital_elements.hpp"
#include "core/planets.hpp"
#include "core/state.hpp"
#include "core/transform.hpp"
#include "core/world.hpp"

enum struct WorldScenario {};

struct EarthStationSatScenario {
    World world;
    EntityId earth_id = kInvalidEntityId;
    EntityId stat_id = kInvalidEntityId;
    EntityId sat_id = kInvalidEntityId;
    bool success = false;
};
inline EarthStationSatScenario make_earth_station_sat_scenario(const vec3d& station_llh) {
    EarthStationSatScenario scenario;
    scenario.earth_id = wgs84(scenario.world);
    scenario.stat_id = scenario.world.spawn_station();
    scenario.sat_id = scenario.world.spawn_satellite();

    Celestial* earth = scenario.world.celestial(scenario.earth_id);
    Satellite* sat = scenario.world.satellite(scenario.sat_id);
    Station* stat = scenario.world.station(scenario.stat_id);
    if (earth == nullptr || sat == nullptr || stat == nullptr) {
        return scenario;
    }

    // earth
    earth->x_att.q = vec4d{0.0, 0.0, 0.0, 1.0};
    earth->x_att.w = vec3d{0.0, 0.0, 7.292115000000000e-05};

    // station
    bool stat_set = scenario.world.set_stat_anchor_detic(
        scenario.stat_id,
        scenario.earth_id,
        station_llh
    );
    if (!stat_set) {
        return scenario;
    }

    scenario.success = true;
    return scenario;
}

struct EarthSatsStatsScenario {
    World world;
    EntityId earth_id = kInvalidEntityId;
    EntityId sat1_id = kInvalidEntityId;
    EntityId sat2_id = kInvalidEntityId;
    EntityId stat1_id = kInvalidEntityId;
    EntityId stat2_id = kInvalidEntityId;
    bool success = false;
};

inline EarthSatsStatsScenario make_earth_sats_stats_scenario(
    GravityModel model = GravityModel::zonal,
    i32 degree_order = 6
) {
    EarthSatsStatsScenario scenario;
    World& world = scenario.world;

    // Earth
    scenario.earth_id = wgs84(world);
    Celestial* earth = world.celestial(scenario.earth_id);
    if (earth == nullptr) {
        return scenario;
    }
    earth->name = "Earth";
    earth->gravity_model = model;
    earth->degree = degree_order;
    earth->order = degree_order;
    if (model == GravityModel::spherical_harmonics) {
        bool gfc_ok = read_gfc(
            pwd + "/assets/EGM2008.gfc.txt",
            earth->C,
            earth->S,
            earth->degree,
            earth->order
        );
        if (!gfc_ok) {
            std::println("GFC Load Failed");
            return scenario;
        }
    }
    earth->propagate_tr = true;
    earth->propagate_att = true;
    earth->attitude_model = CelestialAttitudeModel::simple_spin;
    earth->x_att.q = dcm_to_ep(rotX(23.44, UAngle::degree));
    earth->set_spin_rate(earth->spin_rate());

    // Stations
    vec3d llh1 = vec3d{0.0, 0.0, 0.0}; // [lat, lon, h] = [deg, deg, sim units]
    scenario.stat1_id = world.spawn_station();
    Station* stat1 = world.station(scenario.stat1_id);
    bool stat1_set
        = world.set_stat_anchor_detic(scenario.stat1_id, scenario.earth_id, llh1);
    if (stat1 == nullptr || !stat1_set) {
        return scenario;
    }

    vec3d llh2 = vec3d{45.0, 30.0, 0.0};
    scenario.stat2_id = world.spawn_station();
    Station* stat2 = world.station(scenario.stat2_id);
    bool stat2_set
        = world.set_stat_anchor_detic(scenario.stat2_id, scenario.earth_id, llh2);
    if (stat2 == nullptr || !stat2_set) {
        return scenario;
    }

    // Satellites
    mat3d I_sat = vec3d{100.0, 200.0, 300.0}.asDiagonal();

    scenario.sat1_id = world.spawn_satellite();
    Satellite* sat1 = world.satellite(scenario.sat1_id);
    if (sat1 == nullptr) {
        return scenario;
    }
    sat1->x_tr.r = vec3d{7000.0, 1000.0, 1300.0};
    sat1->x_tr.v = vec3d{-0.5, 7.2, 1.0};
    sat1->x_att.q = q_identity;
    sat1->x_att.w = vec3d{0.000001, 0.0025, 0.000001};
    sat1->set_I(I_sat);

    scenario.sat2_id = world.spawn_satellite();
    Satellite* sat2 = world.satellite(scenario.sat2_id);
    if (sat2 == nullptr) {
        return scenario;
    }
    sat2->x_tr = classical_to_rv(
        26560.0,
        0.01,
        55.0,
        30.0,
        0.0,
        45.0,
        earth->mu,
        UAngle::degree
    );
    sat2->x_att.q = q_identity;
    sat2->x_att.w = vec3d{-0.0005, -0.025, -0.00025};
    sat2->set_I(I_sat);

    scenario.success = true;
    return scenario;
}
