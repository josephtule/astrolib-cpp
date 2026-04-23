#pragma once

#include "util/constants.hpp"
#include "util/vecdefs.hpp"
#include <cmath>
#include <cstddef>


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
    ecref<vec3d> r_rel, // relative pos of obj wrt to gravity source in body frame
    f64 mu,             // target body gravitational parameter
    f64 R_cb,           // target body radius
    i32 degree,
    vec7d J, // zonal gravity parameters [J0, J1, J2, J3, J4, J5, J6]
    f64 epsilon = tol_strict
) {
    vec3d a = accel_gravity_pointmass(r_rel, mu);
    degree = std::min(degree, 6);

    f64 r_mag2 = r_rel.squaredNorm();
    f64 r_mag = std::sqrt(r_mag2);
    if (r_mag < epsilon) return a;

    f64 Rr = R_cb / r_mag;
    f64 mur2 = mu / r_mag2;
    f64 r0r = r_rel(0) / r_mag;
    f64 r1r = r_rel(1) / r_mag;
    f64 r2r = r_rel(2) / r_mag;

    f64 zr2 = r2r * r2r;
    f64 zr4 = zr2 * zr2;
    f64 Rr2 = Rr * Rr;

    switch (degree) {
    case 6: {
        f64 Rr6 = Rr2 * Rr2 * Rr2;
        f64 zr6 = zr4 * zr2;
        f64 coef = -7. / 16. * J(6) * mur2 * Rr6;
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
        f64 zr3 = zr2 * r2r;
        f64 zr5 = zr4 * r2r;
        f64 coef = 3. / 8. * J(5) * mur2 * Rr5;
        vec3d dir = {
            7. * (5. * r2r - 30. * zr3 + 33. * zr5) * r0r,
            7. * (5. * r2r - 30. * zr3 + 33. * zr5) * r1r,
            -(5. - 105. * zr2 + 315. * zr4 - 231. * zr5 * r2r),
        };
        a += coef * dir;
    }
        [[fallthrough]];
    case 4: {
        f64 Rr4 = Rr2 * Rr2;
        f64 coef = 5. / 8. * J(4) * mur2 * Rr4;
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
        f64 coef = -1. / 2. * J(3) * mur2 * Rr3;
        vec3d dir = {
            5. * (3. * r2r - 7. * zr3) * r0r,
            5. * (3. * r2r - 7. * zr3) * r1r,
            -(3. - 30. * zr2 + 35. * zr3 * r2r),
        };
        a += coef * dir;
    }
        [[fallthrough]];
    case 2: {
        f64 coef = -3. / 2. * J(2) * mur2 * Rr2;
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

inline void norm_legendre(
    f64 phi,
    i32 degree,
    i32 order,
    eref<matXd> P,
    eref<matXd> scales,
    f64 tol = tol_strict
) {
    f64 cphi = std::cos(pio2 - phi);
    f64 sphi = std::sin(pio2 - phi);

    if (std::abs(cphi) <= tol) cphi = 0;
    if (std::abs(sphi) <= tol) sphi = 0;

    f64 sqrt3 = std::sqrt(3);
    P(0, 0) = 1.0;
    P(1, 0) = sqrt3 * cphi;
    P(1, 1) = sqrt3 * sphi;
    scales(0, 0) = 0.0;
    scales(1, 0) = 1.0;
    scales(1, 1) = 0.0;

    for (i32 n = 2; n < degree + 1; ++n) {
        f64 sqrt2np1 = std::sqrt(2.0 * n + 1);

        // Compute through order+1 because dU/dphi uses P(n, m+1).
        i32 m_max = std::min(n, order + 1) + 1;
        for (i32 m = 0; m < m_max; ++m) {
            if (n == m) { // Diagonals
                P(n, n) = sqrt2np1 / std::sqrt(2.0 * n) * sphi * P(n - 1, n - 1);
                scales(n, n) = 0;
            } else if (m == 0) { // Zonals
                P(n, m) = sqrt2np1 / n
                          * (std::sqrt(2.0 * n - 1) * cphi * P(n - 1, m)
                             - (n - 1.0) / std::sqrt(2.0 * n - 3.0) * P(n - 2, m));
                scales(n, m) = std::sqrt((n + 1.0) * n / 2.0);
            } else {
                P(n, m) = sqrt2np1 / (std::sqrt(n + m) * std::sqrt(n - m))
                          * (std::sqrt(2.0 * n - 1.0) * cphi * P(n - 1, m)
                             - std::sqrt(n + m - 1.0) * std::sqrt(n - m - 1.0)
                                   / std::sqrt(2.0 * n - 3.0) * P(n - 2, m));
                scales(n, m) = std::sqrt((n + m + 1.0) * (n - m));
            }
        }
    }
}

inline vec3d accel_gravity_spherical_harmonics(
    ecref<vec3d> r_rel, // relative pos of obj wrt to gravity source in body frame
    f64 mu,
    f64 R_cb,
    i32 degree,
    i32 order,
    ecref<matXd> C,
    ecref<matXd> S,
    f64 tol = tol_strict
) {
    vec3d a = vec3d::Zero();

    // Order and degree guards
    if (degree < 2 || order < 0) return a;
    order = std::min(order, degree);

    // Position and projected norms
    f64 x = r_rel(0);
    f64 y = r_rel(1);
    f64 z = r_rel(2);
    f64 r_mag2 = r_rel.squaredNorm();
    f64 r_mag = std::sqrt(r_mag2);
    f64 xxyy = x * x + y * y;
    f64 sqrtxxyy = std::sqrt(xxyy);
    bool at_pole = sqrtxxyy < tol;

    //  Geocentric lattitude
    if (r_mag < tol) return a;
    f64 phi = std::asin(z / r_mag);

    // Legendre polynomials
    matXd P = matXd::Zero(degree + 3, order + 3);
    matXd scales = matXd::Zero(degree + 3, order + 3);
    norm_legendre(phi, degree, order, P, scales);

    // Longitude and trig
    f64 lambda = std::atan2(y, x);
    f64 slam = std::sin(lambda);
    f64 clam = std::cos(lambda);

    vecXd slams = vecXd::Zero(order + 1);
    vecXd clams = vecXd::Zero(order + 1);
    slams(0) = 0.0;
    clams(0) = 1.0;
    if (order >= 1) {
        slams(1) = slam;
        clams(1) = clam;
    }
    for (i32 m = 2; m < order + 1; ++m) {
        slams(m) = 2.0 * clam * slams(m - 1) - slams(m - 2);
        clams(m) = 2.0 * clam * clams(m - 1) - clams(m - 2);
    }

    f64 r_ratio = R_cb / r_mag;
    f64 r_ratio_n = r_ratio;

    // Initialize gravity summation (radial; over n)
    f64 dUdr_sumN = 1.0;
    f64 dUdphi_sumN = 0.0;
    f64 dUdlambda_sumN = 0.0;

    // Summation loop
    for (i32 n = 2; n < degree + 1; ++n) {
        r_ratio_n = r_ratio_n * r_ratio;

        // Initialize gravity summation (radial; over m)
        f64 dUdr_sumM = 0.0;
        f64 dUdphi_sumM = 0.0;
        f64 dUdlambda_sumM = 0.0;

        i32 m_max = std::min(n, order) + 1;
        for (i32 m = 0; m < m_max; ++m) {
            f64 phi_factor = at_pole ? 0.0 : z / sqrtxxyy;
            dUdr_sumM += P(n, m) * (C(n, m) * clams(m) + S(n, m) * slams(m));
            dUdphi_sumM += (P(n, m + 1) * scales(n, m) - phi_factor * m * P(n, m))
                           * (C(n, m) * clams(m) + S(n, m) * slams(m));
            dUdlambda_sumM += m * P(n, m) * (S(n, m) * clams(m) - C(n, m) * slams(m));
        }
        dUdr_sumN += dUdr_sumM * r_ratio_n * (n + 1);
        dUdphi_sumN += dUdphi_sumM * r_ratio_n;
        dUdlambda_sumN += dUdlambda_sumM * r_ratio_n;
    }

    // Acceleration in spherical coordinates
    f64 dUdr = -mu / r_mag2 * dUdr_sumN;
    f64 dUdphi = mu / r_mag * dUdphi_sumN;
    f64 dUdlambda = mu / r_mag * dUdlambda_sumN;

    // Acceleration in body-fixed Cartesian coordinates
    if (at_pole) {
        // if (std::abs(std::atan2(z, sqrtxxyy)) == pio2) {
        // Special case for poles
        a = {0.0, 0.0, 1.0 / r_mag * dUdr * z};
    } else {
        a(0) = (1.0 / r_mag * dUdr - z / (r_mag2 * sqrtxxyy) * dUdphi) * x
               - (dUdlambda / xxyy) * y;
        a(1) = (1.0 / r_mag * dUdr - z / (r_mag2 * sqrtxxyy) * dUdphi) * y
               + (dUdlambda / xxyy) * x;
        a(2) = 1.0 / r_mag * dUdr * z + sqrtxxyy / r_mag2 * dUdphi;
    }

    return a;
}

