#pragma once

#include "util/constants.hpp"
#include "util/vecdefs.hpp"
#include <cmath>

inline f64 inv_r3_safe(const ecref<vec3d> r, f64 epsilon = tol_strict) {
    f64 r_mag = std::sqrt(r.squaredNorm());
    if (r_mag < epsilon) return 0.0;
    f64 inv_r = 1.0 / r_mag;
    return inv_r * inv_r * inv_r;
}
inline vec3d accel_gravity_pointmass(
    ecref<vec3d> r_rel,
    f64 mu,
    f64 epsilon = tol_strict
) {
    return -mu * r_rel * inv_r3_safe(r_rel, epsilon);
}
inline vec3d accel_gravity_zonal(
    ecref<vec3d> r_rel,
    f64 mu,   // target body gravitational parameter
    f64 R_cb, // target body radius
    i32 degree,
    vec7d J, // zonal gravity parameters [J0, J1, J2, J3, J4, J5, J6]
    f64 epsilon = tol_strict
) {
    vec3d a = accel_gravity_pointmass(r_rel, mu);

    f64 r_mag2 = r_rel.squaredNorm();
    f64 r_mag = std::sqrt(r_mag2);
    f64 Rr = R_cb / r_mag;
    f64 mur2 = mu / r_mag2;
    f64 r0r = r_rel[0] / r_mag;
    f64 r1r = r_rel[1] / r_mag;
    f64 r2r = r_rel[2] / r_mag;

    f64 zr2 = r2r * r2r;
    f64 zr4 = zr2 * zr2;
    f64 Rr2 = Rr * Rr;

    switch (degree) {
    case 6: {
        f64 Rr6 = Rr2 * Rr2 * Rr2;
        f64 zr6 = zr4 * zr2;
        f64 coef = -7. / 16. * J[6] * mur2 * Rr6;
        vec3d dir = {
            (5. - 135. * zr2 + 495. * zr4 - 429. * zr6) * r0r,  //
            (5. - 135. * zr2 + 495. * zr4 - 429. * zr6) * r1r,  //
            (35. - 315. * zr2 + 693. * zr4 - 429. * zr6) * r2r, //
        };
        a += coef * dir;
    }
        [[fallthrough]];
    case 5: {
        f64 Rr5 = Rr2 * Rr2 * Rr;
        f64 zr5 = zr4 * r2r;
        f64 coef = 3. / 8. * J[5] * mur2 * Rr5;
        vec3d dir = {
            7. * (5. * r2r - 30. * zr4 + 33. * zr5) * r0r,
            7. * (5. * r2r - 30. * zr4 + 33. * zr5) * r1r,
            -(5. - 105. * zr2 + 315. * zr4 - 231. * zr5),
        };
        a += coef * dir;
    }
        [[fallthrough]];
    case 4: {
        f64 Rr4 = Rr2 * Rr2;
        f64 coef = 5. / 8. * J[4] * mur2 * Rr4;
        vec3d dir = {
            3. * (1. - 14. * zr2 + 21. * zr4) * r0r,
            3. * (1. - 14. * zr2 + 21. * zr4) * r1r,
            (15. - 70. * zr2 + 63. * zr4) * r2r,
        };
        a += coef * dir;
    }
        [[fallthrough]];
    case 3: {
        f64 Rr3 = Rr2 * Rr;
        f64 zr3 = zr2 * r2r;
        f64 coef = -1. / 2. * J[3] * mur2 * Rr3;
        vec3d dir = {
            5. * (3. * r2r - 7. * zr3) * r0r,
            5. * (3. * r2r - 7. * zr3) * r1r,
            -(3. - 30. * zr2 + 35. * zr3 * r2r),
        };
        a += coef * dir;
    }
        [[fallthrough]];
    case 2: {
        f64 coef = -3. / 2. * J[2] * mur2 * Rr2;
        vec3d dir = {
            (1. - 5. * zr2) * r0r, //
            (1. - 5. * zr2) * r1r, //
            (3. - 5. * zr2) * r2r, //
        };
        a += coef * dir;
        break;
    }
    default:
    }
    return a;
}

inline vec3d accel_gravity_zonal_single(
    ecref<vec3d> r_rel,
    f64 mu,   // target body gravitational parameter
    f64 R_cb, // target body radius
    i32 degree,
    vec7d J, // zonal gravity parameters [J0, J1, J2, J3, J4, J5, J6]
    f64 epsilon = tol_strict
) {
    vec3d a = accel_gravity_pointmass(r_rel, mu);

    f64 r_mag2 = r_rel.squaredNorm();
    f64 r_mag = std::sqrt(r_mag2);
    f64 Rr = R_cb / r_mag;
    f64 mur2 = mu / r_mag2;
    f64 r0r = r_rel[0] / r_mag;
    f64 r1r = r_rel[1] / r_mag;
    f64 r2r = r_rel[2] / r_mag;

    f64 zr2 = r2r * r2r;
    f64 zr4 = zr2 * zr2;
    f64 Rr2 = Rr * Rr;

    switch (degree) {
    case 6: {
        f64 Rr6 = Rr2 * Rr2 * Rr2;
        f64 zr6 = zr4 * zr2;
        f64 coef = -7. / 16. * J[6] * mur2 * Rr6;
        vec3d dir = {
            (5. - 135. * zr2 + 495. * zr4 - 429. * zr6) * r0r,  //
            (5. - 135. * zr2 + 495. * zr4 - 429. * zr6) * r1r,  //
            (35. - 315. * zr2 + 693. * zr4 - 429. * zr6) * r2r, //
        };
        a += coef * dir;
    } break;
    case 5: {
        f64 Rr5 = Rr2 * Rr2 * Rr;
        f64 zr5 = zr4 * r2r;
        f64 coef = 3. / 8. * J[5] * mur2 * Rr5;
        vec3d dir = {
            7. * (5. * r2r - 30. * zr4 + 33. * zr5) * r0r,
            7. * (5. * r2r - 30. * zr4 + 33. * zr5) * r1r,
            -(5. - 105. * zr2 + 315. * zr4 - 231. * zr5),
        };
        a += coef * dir;
    } break;
    case 4: {
        f64 Rr4 = Rr2 * Rr2;
        f64 coef = 5. / 8. * J[4] * mur2 * Rr4;
        vec3d dir = {
            3. * (1. - 14. * zr2 + 21. * zr4) * r0r,
            3. * (1. - 14. * zr2 + 21. * zr4) * r1r,
            (15. - 70. * zr2 + 63. * zr4) * r2r,
        };
        a += coef * dir;
    } break;
    case 3: {
        f64 Rr3 = Rr2 * Rr;
        f64 zr3 = zr2 * r2r;
        f64 coef = -1. / 2. * J[3] * mur2 * Rr3;
        vec3d dir = {
            5. * (3. * r2r - 7. * zr3) * r0r,
            5. * (3. * r2r - 7. * zr3) * r1r,
            -(3. - 30. * zr2 + 35. * zr3 * r2r),
        };
        a += coef * dir;
    } break;
    case 2: {
        f64 coef = -3. / 2. * J[2] * mur2 * Rr2;
        vec3d dir = {
            (1. - 5. * zr2) * r0r, //
            (1. - 5. * zr2) * r1r, //
            (3. - 5. * zr2) * r2r, //
        };
        a += coef * dir;
        break;
    }
    default:
    }
    return a;
}