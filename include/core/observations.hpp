#pragma once

#include "core/state.hpp"
#include "util/constants.hpp"
#include "util/math.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"
#include <cmath>

inline vec3d radec_from_rel(
    const vec3d& r_rel,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol_strict
) {
    vec3d radec = vec3d0;

    // Topocentric (observer relative position of obj)
    f64 x = r_rel(0);
    f64 y = r_rel(1);
    f64 z = r_rel(2);
    f64 rho = r_rel.norm();
    if (rho <= tol) return radec;

    // Compute right ascension and declination
    f64 ra = wrap_angle(atan2(y, x), 0.0, twopi);
    f64 dec = atan2(z, std::sqrt(x * x + y * y));

    if (angle_out != UAngle::radian) {
        ra = convert_angle(ra, UAngle::radian, angle_out);
        dec = convert_angle(dec, UAngle::radian, angle_out);
    }

    radec = vec3d{ra, dec, rho};

    return radec;
}

inline vec3d radec_from_pos(
    const vec3d& r_target,
    const vec3d& r_observer,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol_strict
) {
    // r_target and r_observer must be in the same frame
    return radec_from_rel(r_target - r_observer, angle_out, tol);
}

inline vec3d radec_rates_from_rel(
    const vec3d& r_rel,
    const vec3d& v_rel,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol_strict
) {
    vec3d radec_dot = vec3d0;

    f64 x = r_rel(0);
    f64 y = r_rel(1);
    f64 z = r_rel(2);
    f64 x_dot = v_rel(0);
    f64 y_dot = v_rel(1);
    f64 z_dot = v_rel(2);

    f64 rho2 = r_rel.squaredNorm();
    f64 rho = std::sqrt(rho2);
    f64 rho_xy2 = (r_rel.segment<2>(0)).squaredNorm();
    f64 rho_xy = std::sqrt(rho_xy2);
    if (rho <= tol || rho_xy <= tol) return radec_dot;

    f64 rho_dot = r_rel.dot(v_rel) / rho;
    f64 ra_dot = (x * y_dot - y * x_dot) / rho_xy2;
    f64 dec_dot = 1.0 / std::sqrt(1.0 - std::pow(z / rho, 2))
                  * (z_dot * rho - z * rho_dot) / rho2;

    if (angle_out != UAngle::radian) {
        ra_dot = convert_angle(ra_dot, UAngle::radian, angle_out);
        dec_dot = convert_angle(dec_dot, UAngle::radian, angle_out);
    }

    radec_dot = vec3d{ra_dot, dec_dot, rho_dot};
    return radec_dot;
}

inline vec3d radec_rates_from_state(
    const StateTr& x_target,
    const StateTr& x_observer,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol_strict
) {
    return radec_rates_from_rel(
        x_target.r - x_observer.r,
        x_target.v - x_observer.v,
        angle_out,
        tol
    );
}

inline vec3d los_from_radec(const vec2d& radec, UAngle angle_in = UAngle::radian) {
    f64 ra = radec(0);
    f64 dec = radec(1);
    if (angle_in != UAngle::radian) {
        ra = convert_angle(ra, angle_in, UAngle::radian);
        dec = convert_angle(dec, angle_in, UAngle::radian);
    }
    return {std::cos(dec) * std::cos(ra), std::cos(dec) * std::sin(ra), std::sin(dec)};
}

inline vec3d los_from_radec(const vec3d& radec, UAngle angle_in = UAngle::radian) {
    return los_from_radec(vec2d{radec(0), radec(1)}, angle_in);
}

// inline vec3d radec_to_azel(vec3d radec, vec3d stat_geod, f64 theta) {
// TODO: do later
// }