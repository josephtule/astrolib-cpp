#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/planets.hpp"
#include "core/state.hpp"
#include "core/world.hpp"

#include <print>

int main() {

    World world;

    EntityId sat_id = world.spawn_satellite();
    world.satellite(sat_id)->x_tr.r = {7000.0, 0.0, 0.0};

    EntityId earth_id = wgs84(world);
    EntityId earth2_id = wgs84(world);
    world.celestial(earth2_id)->x_tr.r = {10000.0, 10000.0, 10000.0};

    auto a_from1 = world.gravity_accel_from(sat_id, earth_id);
    auto a_from2 = world.gravity_accel_from(sat_id, earth2_id);
    auto a_on = world.gravity_accel_on(sat_id);
    auto earth = world.celestial(earth_id);
    auto earth2 = world.celestial(earth2_id);

    std::println(
        "from earth1: {} (mag = {})\nfrom earth2: {} (mag = {})\non: {} (mag = {})",
        a_from1,
        a_from1.norm(),
        a_from2,
        a_from2.norm(),
        a_on,
        a_on.norm()
    );
    std::println("error: {}", a_from1 + a_from2 - a_on);

    std::println();
    std::println("{}", world.gravity_accel_from(earth_id, earth_id));

    earth->use_simple_spin = true;
    earth->spin_rate = 100;
    earth->x_att.w[2] = earth->spin_rate;

    std::println("{}", world.celestial(earth_id)->x_att.w);

    EntityId stat_id = world.spawn_station();
    auto station = world.station(stat_id);
    station->anchor_id = earth_id;
    station->r_body = {earth->mean_radius, 0, 0};

    std::println("stat_r_inertial: {}", world.stat_r_inertial(stat_id));
    std::println("stat_v_inertial: {}", world.stat_v_inertial(stat_id));
    std::println("stat_r_inertial: {}", world.stat_x_tr_inertial(stat_id).r);
    std::println("stat_v_inertial: {}", world.stat_x_tr_inertial(stat_id).v);

    std::println();

    std::println("z inertial earth: {}", world.body_z_inertial(earth_id));
    std::println("w inertial earth: {}", world.body_w_inertial(earth_id));
    vec4d q_earth0 = earth->x_att.q;

    earth->x_att.q = vec4d{1, 2, 3, 4}.normalized();
    auto a_from1_rot = world.gravity_accel_from(sat_id, earth_id);

    earth->x_att.q = q_default;
    earth->gravity_model = GravityModel::zonal;
    earth->degree = 2;
    auto a_from1_j2 = world.gravity_accel_from(sat_id, earth_id);

    earth->x_att.q = vec4d{1, 2, 3, 4}.normalized();
    auto a_from1_j2rot = world.gravity_accel_from(sat_id, earth_id);

    std::println(
        "from earth (point): {} (mag = {})\nfrom earth (point rot): {} (mag = {})\nfrom "
        "earth1 (j2): {} (mag = {})\nfrom earth1 "
        "(j2rot): {} (mag = {})",
        a_from1,
        a_from1.norm(),
        a_from1_rot,
        a_from1_rot.norm(),
        a_from1_j2,
        a_from1_j2.norm(),
        a_from1_j2rot,
        a_from1_j2rot.norm()
    );
    std::println("J2 pertrubation mag: {}", (a_from1 - a_from1_j2rot).norm());
    std::println("Earth J2: {}", earth->J[2]);

    earth->gravity_model = GravityModel::pointmass;
    earth->x_att.q = q_earth0;
    vec3d a_point = world.gravity_accel_from(sat_id, earth_id);

    earth->x_att.q = vec4d{1, 2, 3, 4}.normalized();
    vec3d a_point_rot = world.gravity_accel_from(sat_id, earth_id);

    earth->gravity_model = GravityModel::zonal;
    earth->degree = 2;
    earth->x_att.q = q_earth0;
    vec3d a_zonal = world.gravity_accel_from(sat_id, earth_id);

    earth->x_att.q = vec4d{1, 2, 3, 4}.normalized();
    vec3d a_zonal_rot = world.gravity_accel_from(sat_id, earth_id);

    std::println("zonal - point: {}", a_zonal - a_point);
    std::println("zonal_rot - point: {}", a_zonal_rot - a_point);
    std::println("zonal_rot - zonal: {}", a_zonal_rot - a_zonal);

    earth->x_att.q = q_earth0;

    return 0;
}