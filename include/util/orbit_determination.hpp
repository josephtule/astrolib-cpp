#pragma once

#include "core/state.hpp"
#include "util/constants.hpp"
#include "util/vecdefs.hpp"

inline StateTr iod_gauss(
    f64 t1,
    f64 t2,
    f64 t3,
    ecref<vec3d> L1,
    ecref<vec3d> L2,
    ecref<vec3d> L3,
    ecref<vec3d> R1,
    ecref<vec3d> R2,
    ecref<vec3d> R3,
    f64 mu,
    f64 tol = tol_med
) {
    // Time intervals
    f64 tau1 = t1 - t2;
    f64 tau3 = t3 - t2;
    f64 tau13 = t3 - t1;
    f64 tau1_2 = tau1 * tau1;
    f64 tau3_2 = tau3 * tau3;
    f64 tau13_2 = tau13 * tau13;

    // Compute D values
    f64 D0 = L1.dot(L2.cross(L3));
    mat3d tempmat1;
    tempmat1.row(0) = R1;
    tempmat1.row(1) = R2;
    tempmat1.row(2) = R3;
    mat3d tempmat2;
    tempmat2.col(0) = L2.cross(L3);
    tempmat2.col(1) = L1.cross(L3);
    tempmat2.col(2) = L1.cross(L2);
    mat3d D = tempmat1 * tempmat2;

    // Compute Alphabet coefs
    f64 A = 1.0 / D0 * (-tau3 / tau13 * D(0, 1) + D(1, 1) + tau1 / tau13 * D(2, 1));
    f64 B = 1.0 / (6.0 * D0)
            * (-(tau13_2 - tau3_2) * tau3 / tau13 * D(0, 1) * (tau13_2 - tau1_2) * tau1
               / tau13 * D(2, 1));

    f64 a = -A * A - 2.0 * A * L2.dot(R2) - R2.squaredNorm();
    f64 b = -2.0 * mu * B * (A + L2.dot(R2));
    f64 c = -mu * mu * B * B;

    // Root solve for ||r2||
    // Using Newton's method for root finding, simple version w/o no heavy setup with
    // optimization library

    // arbitrary, larger than ranges
    f64 x0 = std::max(R1.norm(), std::max(R2.norm(), R3.norm())) * 5.0;

    auto fun = [a, b, c](f64 x) -> f64 {
        return std::pow(x, 8) + a * std::pow(x, 6) + b * std::pow(x, 3) + c;
    };
    auto dfun = [a, b](f64 x) -> f64 {
        return 8.0 * std::pow(x, 7) + 6.0 * a * std::pow(x, 5) + 3.0 * b * std::pow(x, 2);
    };
    i32 max_iter = 10000;
    i32 iter = 0;
    f64 x_iter = x0;
    while (iter < max_iter) {
        f64 dfx = dfun(x_iter);
        if (std::abs(dfun(x_iter)) < tol) { // gradient convergence
            break;
        }

        f64 x_next = x_iter - fun(x_iter) / dfx;

        if (std::abs(x_next - x_iter) < tol) { // delta x convergence
            break;
        }

        x_iter = x_next;
        ++iter;
    }
    // TODO: handle max iterations
    f64 r2_mag = x_iter;
    f64 r2_mag3 = r2_mag * r2_mag * r2_mag;

    // Compute ranges
    f64 rho1 = 1.0 / D0
               * ((6 * (D(2, 0) * tau1 / tau3 + D(1, 0) * tau13 / tau3) * r2_mag3
                   + mu * D(2, 0) * (tau13_2 - tau1_2) * tau1 / tau3)
                      / (6.0 * r2_mag3 + mu * (tau13_2 - tau3_2))
                  - D(0, 0));
    f64 rho2 = A + mu * B / r2_mag3;
    f64 rho3 = 1 / D0
               * ((6 * (D(0, 2) * tau3 / tau1 - D(1, 2) * tau13 / tau1) * r2_mag3
                   + mu * D(0, 2) * (tau13_2 - tau3_2) * tau3 / tau1)
                      / (6 * r2_mag3 + mu * (tau13_2 - tau1_2))
                  - D(2, 2));

    // Compute position vectors
    vec3d r1 = R1 + rho1 * L1;
    vec3d r2 = R2 + rho2 * L2;
    vec3d r3 = R3 + rho3 * L3;

    r2_mag = r2.norm();
    r2_mag3 = r2_mag * r2_mag * r2_mag;
    f64 f1 = 1.0 - 0.5 * mu / r2_mag3 * tau1_2;
    f64 f3 = 1.0 - 0.5 * mu / r2_mag3 * tau3_2;
    f64 g1 = tau1 - 1.0 / 6.0 * mu / r2_mag3 * tau1_2 * tau1;
    f64 g3 = tau3 - 1.0 / 6.0 * mu / r2_mag3 * tau3_2 * tau3;

    vec3d v2 = 1.0 / (f1 * g3 - f3 - g1) * (f1 * r3 - f3 * r1);

    StateTr x2 = StateTr{.r = r2, .v = v2};

    return x2;
}

inline StateTr iod_gauss(
    std::array<f64, 3> t,
    const std::array<ecref<vec3d>, 3> L,
    const std::array<ecref<vec3d>, 3> R,
    f64 mu
) {
    return iod_gauss(t[0], t[1], t[2], L[0], L[1], L[2], R[0], R[1], R[2], mu);
}

inline StateTr iod_gauss(svec<f64> t, const svec<vec3d> L, const svec<vec3d> R, f64 mu) {
    return iod_gauss(t[0], t[1], t[2], L[0], L[1], L[2], R[0], R[1], R[2], mu);
}

inline StateTr iod_gauss(ecref<vec3d> t, ecref<mat3d> L, ecref<mat3d> R, f64 mu) {
    const vec3d L1 = L.col(0);
    const vec3d L2 = L.col(1);
    const vec3d L3 = L.col(2);
    const vec3d R1 = R.col(0);
    const vec3d R2 = R.col(1);
    const vec3d R3 = R.col(2);

    return iod_gauss(t(0), t(1), t(2), L1, L2, L3, R1, R2, R3, mu);
}

inline StateTr iod_gibbs(
    ecref<vec3d> r1,
    ecref<vec3d> r2,
    ecref<vec3d> r3,
    f64 mu,
    f64 tol = tol_med
) {

    // Precompute cross prodcuts
    vec3d r2x3 = r2.cross(r3);
    vec3d r3x1 = r3.cross(r1);
    vec3d r1x2 = r1.cross(r2);

    // Magnitudes
    f64 r1_mag = r1.norm();
    f64 r2_mag = r2.norm();
    f64 r3_mag = r3.norm();

    // Plane condition
    vec3d r1_hat = r1 / r1_mag;
    vec3d r2x3_hat = r2x3.normalized();
    if (r1_hat.dot(r2x3_hat) > tol * r2_mag) {
        // TODO: handle vectors not in the same plane
    }

    vec3d n = r1_mag * r2x3 + r2_mag * r3x1 + r3_mag * r1x2;
    vec3d d = r1x2 + r2x3 + r3x1;
    vec3d s = r1 * (r2_mag - r3_mag) + r2 * (r3_mag - r1_mag) + r3 * (r1_mag - r2_mag);

    vec3d v2 = std::sqrt(mu / (n.norm() * d.norm())) * (d.cross(r2) / r2_mag + s);

    return StateTr{.r = r2, .v = v2};
}

inline StateTr iod_herrickgibbs(
    f64 t1,
    f64 t2,
    f64 t3,
    ecref<vec3d> r1,
    ecref<vec3d> r2,
    ecref<vec3d> r3,
    f64 mu,
    f64 tol = tol_med
) {
    // Precompute cross prodcuts
    vec3d r2x3 = r2.cross(r3);

    // Magnitudes
    f64 r1_mag = r1.norm();
    f64 r2_mag = r2.norm();
    f64 r3_mag = r3.norm();
    f64 r1_mag3 = std::pow(r1_mag, 3);
    f64 r2_mag3 = std::pow(r2_mag, 3);
    f64 r3_mag3 = std::pow(r3_mag, 3);

    // Plane condition
    vec3d r1_hat = r1 / r1_mag;
    vec3d r2x3_hat = r2x3.normalized();
    if (r1_hat.dot(r2x3_hat) > tol * r2_mag) {
        // TODO: handle vectors not in the same plane
    }

    // Delta times
    f64 dt31 = t3 - t1;
    f64 dt21 = t2 - t1;
    f64 dt32 = t3 - t2;

    vec3d v2 = -dt32 * (1.0 / (dt21 * dt31) + mu / (12.0 * r1_mag3)) * r1
               + (dt32 - dt21) * (1.0 / (dt21 * dt32) + mu / (12.0 * r2_mag3)) * r2
               + dt21 * (1.0 / (dt32 * dt31) + mu / (12.0 * r3_mag3)) * r3;

    return StateTr{.r = r2, .v = v2};
}

inline StateTr iod_laplace(
    f64 t1,
    f64 t2,
    f64 t3,
    ecref<vec3d> L1,
    ecref<vec3d> L2,
    ecref<vec3d> L3,
    ecref<vec3d> R1,
    ecref<vec3d> R2,
    ecref<vec3d> R3,
    f64 mu,
    vec3d w,
    f64 tol = tol_med
) {

    // Middle observer position derivatives
    vec3d R2_dot = w.cross(R2);
    vec3d R2_ddot = w.cross(R2_dot);

    // Delta times
    f64 dt23 = t2 - t3;
    f64 dt12 = t1 - t2;
    f64 dt13 = t1 - t3;
    f64 dt21 = -dt12;
    f64 dt31 = -dt13;
    f64 dt32 = -dt23;

    // Lagrange interpolation coefs
    f64 s1 = dt23 / (dt12 * dt13);
    f64 s2 = (dt21 + dt23) / (dt12 * dt23);
    f64 s3 = dt21 / (dt31 * dt32);
    f64 s4 = 2.0 / (dt12 * dt13);
    f64 s5 = 2.0 / (dt21 * dt23);
    f64 s6 = 2.0 / (dt31 * dt32);

    // Middle observer line of sight derivatives
    vec3d L2_dot = s1 * L2 + s2 * L2 + s3 * L3;
    vec3d L2_ddot = s5 * L1 + s5 * L2 + s6 * L3;

    // Determinants
    mat3d mat0;
    mat0.col(0) = L2;
    mat0.col(1) = L2_dot;
    mat0.col(2) = L2_ddot;
    f64 det0 = 2.0 * mat0.determinant();

    mat3d mat1;
    mat1.col(0) = L2;
    mat1.col(1) = L2_dot;
    mat1.col(2) = R2_ddot;
    f64 det1 = mat1.determinant();

    mat3d mat2;
    mat2.col(0) = L2;
    mat2.col(1) = L2_dot;
    mat2.col(2) = R2;
    f64 det2 = mat2.determinant();

    mat3d mat3;
    mat3.col(0) = L2;
    mat3.col(1) = R2_ddot;
    mat3.col(2) = L2_ddot;
    f64 det3 = mat3.determinant();

    mat3d mat4;
    mat4.col(0) = L2;
    mat4.col(1) = R2;
    mat4.col(2) = L2_ddot;
    f64 det4 = mat4.determinant();

    // Polynomial coefs
    f64 det10 = det1 / det0;
    f64 det20 = det2 / det0;
    f64 det30 = det3 / det0;
    f64 det40 = det4 / det0;

    f64 a = -4.0 * std::pow(det10, 2) + 4.0 * det10 * L2.dot(R2) - R2.squaredNorm();
    f64 b = -8.0 * mu * det10 * det20 + 4.0 * mu * det20 * L2.dot(R2);
    f64 c = -4.0 * mu * mu * std::pow(det20, 2);

    // Root solve for ||r2||
    // Using Newton's method for root finding, simple version w/o no heavy setup with
    // optimization library

    // arbitrary, larger than ranges
    f64 x0 = std::max(R1.norm(), std::max(R2.norm(), R3.norm())) * 5.0;

    auto fun = [a, b, c](f64 x) -> f64 {
        return std::pow(x, 8) + a * std::pow(x, 6) + b * std::pow(x, 3) + c;
    };
    auto dfun = [a, b](f64 x) -> f64 {
        return 8.0 * std::pow(x, 7) + 6.0 * a * std::pow(x, 5) + 3.0 * b * std::pow(x, 2);
    };
    i32 max_iter = 10000;
    i32 iter = 0;
    f64 x_iter = x0;
    while (iter < max_iter) {
        f64 dfx = dfun(x_iter);
        if (std::abs(dfun(x_iter)) < tol) { // gradient convergence
            break;
        }

        f64 x_next = x_iter - fun(x_iter) / dfx;

        if (std::abs(x_next - x_iter) < tol) { // delta x convergence
            break;
        }

        x_iter = x_next;
        ++iter;
    }
    // TODO: handle max iterations
    f64 r2_mag = x_iter;
    f64 r2_mag3 = std::pow(r2_mag, 3);

    // Range and range rate
    f64 rho = -2.0 * det10 - 2.0 * mu / r2_mag3 * det20;
    f64 rho_dot = -det30 - mu / r2_mag3 * det40;

    vec3d r2 = R2 + rho * L2;
    vec3d v2 = R2_dot + rho_dot * L2 + rho * L2_dot;

    return StateTr{.r = r2, .v = v2};
}
// TODO: overloads for previous iod

enum struct StatObservationModel {
    // Coordinates that ground stations record from space objects
    unknown,
    topo_radec,
    topo_azel,
    topo_inertial_cart, // inertial relative to station
    topo_fixed_cart,    // body fixed relative to station
    // TODO: create observation models and jacobians for these (to use in EKF)
};
// TODO: not sure if i should do structs and enum for observation types or variants
struct StatObservation {
    StatObservationModel obs_model = StatObservationModel::unknown;
};
struct StatRADecObsv : public StatObservation {
    f64 ra = 0.0;
    f64 dec = 0.0;
};
struct StatAzElObsv : public StatObservation {
    f64 az = 0.0;
    f64 el = 0.0;
};

enum struct VehicleObservationModel {};

inline void LUMVEstat(
    StateTr x_0, // initial guess
    StatObservationModel obs_model = StatObservationModel::topo_inertial_cart,
    i32 max_iter = 10
) {
    // Nonlinear Least Squares / Linear, Unbiased, Minimum Variance Estimate (NLS/LUMVE)
    // This is a nonlinear least squares iterative batch estimator for ground station
    // observers to perform orbit determination corrections
}