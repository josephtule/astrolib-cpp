#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/world.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

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
    earth.attitude_model = CelestialAttitudeModel::simple_spin;
    earth.set_spin_rate(omega);

    earth.J = -vec7d{
        0.000000000000000000,
        0.000000000000000000,
        -0.001082626173852223,
        0.000002532410518568,
        0.000001619897599917,
        0.000000227753590731,
        -0.000000540666576284,
    };

    earth.name = "Earth (WGS84)";
    return earth;
}

inline EntityId wgs84(World& world, ULength u_len = ULength::kilometer) {
    auto earth = std::make_unique<Celestial>(wgs84(u_len));
    EntityId id = world.insert_celestial(std::move(earth));
    return id;
}

inline bool read_egm2008(
    const std::string& filename,
    matXd& C,
    matXd& S,
    i32 degree,
    i32 order,
    i32 lineskips = 0
) {
    order = std::min(order, degree);
    C = matXd::Zero(degree + 1, order + 1);
    S = matXd::Zero(degree + 1, order + 1);

    std::ifstream file(filename);
    if (!file) return false;

    std::string line;

    for (i32 i = 0; i < lineskips && std::getline(file, line); ++i) {}

    while (getline(file, line)) {

        std::replace(line.begin(), line.end(), 'D', 'E');
        std::replace(line.begin(), line.end(), 'd', 'e');

        std::istringstream iss(line);

        i32 n, m;
        f64 c, s;
        if (!(iss >> n >> m >> c >> s)) continue;

        if (n > degree) break;
        if (m > order) continue;

        C(n, m) = c;
        S(n, m) = s;
    }
    return true;
}

inline f64 kronecker_delta(i32 n, i32 m) { return n == m ? 1 : 0; }
inline i64 factorial(int n) {
    if (n < 0) return 0;

    i64 result = 1;
    for (int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}
inline f64 norm_factor(i32 n, i32 m) {
    // multiply to denormalize, divide to normalize
    f64 factor = std::sqrt(
        (2.0 - kronecker_delta(0, m)) * (2.0 * n + 1) * static_cast<f64>(factorial(n - m))
        / static_cast<f64>(factorial(n + m))
    );
    return factor;
}
