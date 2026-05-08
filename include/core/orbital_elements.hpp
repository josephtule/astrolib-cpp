#pragma once

#include "core/state.hpp"
#include "core/transform.hpp"
#include "util/constants.hpp"
#include "util/math.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include <cmath>
#include <print>

struct OEClassical {
    // angles stored as radians internally
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

inline OEClassical rv_to_classical(
    vec3d r,
    vec3d v,
    f64 mu,
    UAngle uangle_out = UAngle::radian,
    f64 tol = tol12
) {
    // r and v in inertial frame
    OEClassical coe;

    f64 r_mag = r.norm();
    f64 v_mag2 = v.squaredNorm();
    f64 v_mag = std::sqrt(v_mag2);

    // Specific angular momentum
    vec3d h = r.cross(v);
    f64 h_mag = h.norm();
    if (r_mag < tol || h_mag < tol || mu <= 0.0) return coe;

    coe.inc = std::acos(clamp_unit(h(2) / h_mag));

    vec3d N = axis_z.cross(h);
    f64 N_mag = N.norm();

    f64 v_r = r.dot(v) / r_mag;
    vec3d e = 1.0 / mu * ((v_mag2 - mu / r_mag) * r - r.dot(v) * v);
    coe.ecc = e.norm();

    f64 E = v_mag2 / 2 - mu / r_mag;
    coe.sma = -mu / E / 2.0;

    bool equatorial = N_mag < tol;
    bool circular = coe.ecc < tol;

    if (!equatorial) {
        coe.raan = std::acos(clamp_unit(N(0) / N_mag));
        if (N(1) < 0.0) coe.raan = 2.0 * pi - coe.raan;
    }

    if (!circular && !equatorial) {
        // non-circular non-equatorial
        coe.aop = std::acos(clamp_unit(N.dot(e) / (N_mag * coe.ecc)));
        if (e(2) < 0.0) coe.aop = 2.0 * pi - coe.aop;

        coe.ta = std::acos(clamp_unit(e.dot(r) / (coe.ecc * r_mag)));
        if (v_r < 0.0) coe.ta = 2.0 * pi - coe.ta;
    } else if (!circular && equatorial) {
        // non-circular equatorial
        // RAAN undefined, store longitude of periapsis in aop
        coe.raan = 0.0;
        coe.aop = std::acos(clamp_unit(e(0) / coe.ecc));
        if (e(1) < 0.0) coe.aop = 2.0 * pi - coe.aop;

        coe.ta = std::acos(clamp_unit(e.dot(r) / (coe.ecc * r_mag)));
        if (v_r < 0.0) coe.ta = 2.0 * pi - coe.ta;
    } else if (circular && !equatorial) {
        // circular non-equatorial
        // AOP undefined, store argument of latitude in ta
        coe.aop = 0.0;
        coe.ta = std::acos(clamp_unit(N.dot(r) / (N_mag * r_mag)));
        if (r(2) < 0.0) coe.ta = 2.0 * pi - coe.ta;
    } else {
        // circular equatorial
        // RAAN/AOP undefined, store true longitude in ta
        coe.raan = 0.0;
        coe.aop = 0.0;
        coe.ta = std::acos(clamp_unit(r(0) / r_mag));
        if (r(1) < 0.0) coe.ta = 2.0 * pi - coe.ta;
    }

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

inline OEClassical rv_to_classical(
    StateTr x_tr,
    f64 mu,
    UAngle uangle_out = UAngle::radian,
    f64 tol = tol12
) {
    return rv_to_classical(x_tr.r, x_tr.v, mu, uangle_out, tol);
}

inline StateTr classical_to_rv(
    f64 sma,
    f64 ecc,
    f64 inc,
    f64 raan,
    f64 aop,
    f64 ta,
    f64 mu,
    UAngle uangle_in = UAngle::radian,
    f64 tol = tol12
) {
    StateTr rv;

    if (uangle_in != UAngle::radian) {
        inc = convert_angle(inc, uangle_in, UAngle::radian);
        raan = convert_angle(raan, uangle_in, UAngle::radian);
        aop = convert_angle(aop, uangle_in, UAngle::radian);
        ta = convert_angle(ta, uangle_in, UAngle::radian);
    }

    bool equatorial = inc < tol;
    bool circular = ecc < tol;

    // Orbit type and anomaly
    f64 anom = 0.0;
    if (circular && equatorial) {
        // RAAN/AOP undefined, ta stores true longitude
        anom = ta + raan + aop;
        aop = 0.0;
        raan = 0.0;
    } else if (circular && !equatorial) {
        // AOP undefined, ta stores argument of latitude
        anom = ta + aop;
        aop = 0.0;
    } else if (!circular && equatorial) {
        // RAAN undefined, aop stores longitude of periapsis
        aop = raan + aop;
        raan = 0.0;
        anom = ta;
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

inline StateTr classical_to_rv(
    OEClassical coe,
    f64 mu,
    UAngle uangle_in = UAngle::radian,
    f64 tol = tol12
) {
    return classical_to_rv(
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
}

inline void print_coe(const OEClassical& coe) {
    std::println("Semimajor Axis: {}", coe.sma);
    std::println("Eccentricity: {}", coe.ecc);
    std::println("Inclination: {}", coe.inc);
    std::println("RA of Asc Node: {}", coe.raan);
    std::println("Arg of Perigee: {}", coe.aop);
    std::println("True Anomaly: {}", coe.ta);
}