#pragma once

#include "Eigen/Geometry"
#include "util/constants.hpp"
#include "util/math.hpp"
#include "util/vecdefs.hpp"

inline vecXd interp_linear(f64 t, vecXd v1, f64 t1, vecXd v2, f64 t2) {
    // lerp
    if (v1.size() != v2.size()) {
        return vecXd{};
    }

    f64 alpha = (t - t1) / (t2 - t1);

    return v1 + (v2 - v1) * alpha;
}

inline vec4d interp_quat_linear(
    f64 t,
    vec4d q1,
    f64 t1,
    vec4d q2,
    f64 t2,
    f64 tol = tol12
) {
    // spherical linear interpolation (slerp)
    if (q1.norm() <= tol || q2.norm() <= tol) {
        return vec4d{};
    }
    f64 dot = q1.dot(q2); // cos(theta)
    if (dot < 0.0) {
        q1 = -q1;
        dot = -dot;
    }

    if (dot > tol * 1000.0) {
        vec4d qt = interp_linear(t, q1, t1, q2, t2);
        if (qt.norm() < tol) {
            return vec4d{};
        }
        return qt.normalized();
    }

    dot = std::clamp(dot, -1.0, 1.0);

    double theta = std::acos(dot);
    double theta_t = theta * t;

    double sin_theta = std::sin(theta);
    double sin_theta_t = std::sin(theta_t);

    double s1 = std::cos(theta_t) - dot * sin_theta_t / sin_theta;
    double s2 = sin_theta_t / sin_theta;

    // Eigen::Quaterniond qa(q1.w(), q1.x(), q1.y(), q1.z());
    // Eigen::Quaterniond qb(q2.w(), q2.x(), q2.y(), q2.z());
    // f64 alpha = (t - t1) / (t2 - t1);
    // Eigen::Quaterniond qt_eig = qa.slerp(alpha, qb);
    // return vec4d{qt_eig.x(), qt_eig.y(), qt_eig.z(), qt_eig.w()};

    return (s1 * q1) + (s2 * q2);
}