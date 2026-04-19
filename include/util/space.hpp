#pragma once

#include "core/body.hpp"
#include "core/state.hpp"
#include "util/constants.hpp"
#include "util/transform.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"
#include <cmath>
#include <complex>

struct OEClassical {
    f64 sma = 0.0;  // semimajor axis
    f64 ecc = 0.0;  // eccentricity
    f64 inc = 0.0;  // inclination
    f64 raan = 0.0; // right ascension of the ascending node
    f64 aop = 0.0;  // argument of perigee
    f64 ta = 0.0;   // true anomaly
};

struct OEEquinoctial {};
struct OEFlight {}; // AKA ADBARV
struct OEPoincare {};
struct OEDelaunay {};

inline OEClassical rv_to_coe(
    vec3d r,
    vec3d v,
    f64 mu,
    UAngle uangle_out = UAngle::degree,
    f64 tol = tol_strict
) {
    // r and v in inertial frame
    OEClassical coe;

    f64 r_mag = r.norm();
    f64 v_mag2 = v.squaredNorm();
    f64 v_mag = std::sqrt(v_mag2);

    // Specific angular momentum
    vec3d h = r.cross(v);
    f64 h_mag = h.norm();

    coe.inc = std::acos(h(2) / h_mag);

    vec3d N = axis_z.cross(h);
    f64 N_mag = N.norm();

    coe.raan = std::acos(N(0) / N_mag);
    if (N(1) < 0.0) coe.raan = 2.0 * pi - coe.raan;

    f64 v_r = r.dot(v) / r_mag;
    vec3d e = 1.0 / mu * ((v_mag2 - mu / r_mag) * r - r.dot(v) * v);
    coe.ecc = e.norm();

    coe.aop = std::acos(N.dot(e) / N_mag / coe.ecc);
    if (e(2) < 0.0) coe.aop = 2.0 * pi - coe.aop;

    coe.ta = std::real(std::acos(e.dot(r) / coe.ecc / r_mag));
    if (v_r < 0.0) coe.ta = 2.0 * pi - coe.ta;

    f64 E = v_mag2 / 2 - mu / r_mag;
    coe.sma = -mu / E / 2.0;

    if (std::abs(coe.raan - 2.0 * pi) < tol) coe.raan = 0.0;
    if (std::abs(coe.aop - 2.0 * pi) < tol) coe.aop = 0.0;
    if (uangle_out != UAngle::radian) {
        coe.inc = convert_angle(coe.inc, UAngle::radian, uangle_out);
        coe.raan = convert_angle(coe.raan, UAngle::radian, uangle_out);
        coe.aop = convert_angle(coe.aop, UAngle::radian, uangle_out);
        coe.ta = convert_angle(coe.ta, UAngle::radian, uangle_out);
    }

    return coe;
}

inline OEClassical rv_to_coe(
    StateTr x_tr,
    f64 mu,
    UAngle uangle_out = UAngle::degree,
    f64 tol = tol_strict
) {
    OEClassical coe;
    rv_to_coe(x_tr.r, x_tr.v, mu, uangle_out, tol);
    return coe;
}

inline StateTr coe_to_rv(
    f64 sma,
    f64 ecc,
    f64 inc,
    f64 raan,
    f64 aop,
    f64 ta,
    f64 mu,
    UAngle uangle_in = UAngle::degree,
    f64 tol = tol_strict
) {
    StateTr rv;

    if (uangle_in != UAngle::radian) {
        inc = convert_angle(inc, uangle_in, UAngle::radian);
        raan = convert_angle(raan, uangle_in, UAngle::radian);
        aop = convert_angle(aop, uangle_in, UAngle::radian);
        ta = convert_angle(ta, uangle_in, UAngle::radian);
    }

    // Orbit type and anomaly
    f64 anom = 0.0;
    if (ecc < tol && inc < tol) {
        // Circular equatorial
        anom = ta + raan + aop;
        aop = 0.0;
        raan = 0.0;
    } else if (ecc < tol && inc > tol) {
        // Circular inclined
        anom = ta + aop;
        aop = 0.0;
    } else if (ecc > tol && inc < tol) {
        // Elliptic equatorial
        anom = raan + aop;
        raan = 0;
    } else {
        anom = ta;
    }

    // Orbit parameter
    f64 p = sma * (1.0 - ecc * ecc);

    // State in orbital plane wrt periapsis
    f64 sanom = std::sin(anom);
    f64 canom = std::cos(anom);
    vec3d r_orbit
        = vec3d{p * canom / (1.0 + ecc * canom), p * sanom / (1.0 + ecc * canom), 0.0};
    vec3d v_orbit
        = vec3d{-std::sqrt(mu / p) * sanom, std::sqrt(mu / p) * (ecc + canom), 0.0};

    // Rotate into inertial frame
    std::array<RotAxis, 3> seq = {RotAxis::z, RotAxis::x, RotAxis::z};
    mat3d R = ea_to_dcm(vec3d{-aop, -inc, -raan}, seq);

    rv.r = R * r_orbit;
    rv.v = R * v_orbit;

    return rv;
}

inline StateTr coe_to_rv(
    OEClassical coe,
    f64 mu,
    UAngle uangle_in = UAngle::degree,
    f64 tol = tol_strict
) {
    StateTr rv;
    rv = coe_to_rv(
        coe.sma,
        coe.ecc,
        coe.inc,
        coe.raan,
        coe.aop,
        coe.ta,
        mu,
        uangle_in,
        tol
    );
    return rv;
}