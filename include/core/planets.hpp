#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/world.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"
#include <memory>

inline Celestial wgs84(ULength u_len = ULength::kilometer) {

    Celestial earth;
    earth.gravity_model = GravityModel::pointmass;

    f64 omega = 7.292115000000000e-05;
    switch (u_len) {
    case ULength::meter:
        earth.mu = 3.986004418000000e+14;
        earth.semimajor_axis = 6378137.;
        earth.semiminor_axis = 6.356752314245179e+06;
        earth.eccentricity = 0.081819190842621;
        earth.flattening = 0.003352810664747;
        // inverse_flattening = 2.982572235630000e+02;
        // third_flattening   = 0.001679220386384;
        earth.mean_radius = 6.371008771415059e+06;
        // surface_area   = 5.100656217240886e+14;
        // volume         = 1.083207319801408e+21;
        // earth.mass = earth.mu / G_m;
        break;
    case ULength::kilometer:
        earth.mu = 3.986004418000000e+05;
        earth.semimajor_axis = 6378.137;
        earth.semiminor_axis = 6.356752314245179e+03;
        earth.eccentricity = 0.081819190842621;
        earth.flattening = 0.003352810664747;
        // inverse_flattening = 2.982572235630000e+02;
        // third_flattening   = 0.001679220386384;
        earth.mean_radius = 6.371008771415059e+03;
        // surface_area   = 5.100656217240886e+08;
        // volume         = 1.083207319801408e+12;
        // earth.mass = earth.mu / G_km;
        break;
    default:
    }

    earth.x_att.q = {0.0, 0.0, 0.0, 1.0};
    earth.x_att.w = {0.0, 0.0, omega};
    earth.use_simple_spin = true;
    earth.spin_rate = omega;

    earth.J = {
        0,
        0,
        0.001082626173852,
        -0.000002532410519,
        -0.000001619897600,
        -0.000000227753591,
        0.000000540666576,
    };

    earth.name = "Earth (WGS84)";
    return earth;
}

inline EntityId wgs84(World& world, ULength u_len = ULength::kilometer) {
    auto earth = std::make_unique<Celestial>(wgs84(u_len));
    EntityId id = world.insert_celestial(std::move(earth));
    return id;
}