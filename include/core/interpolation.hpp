// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "util/constants.hpp"
#include "util/vecdefs.hpp"

#include <algorithm>
#include <cmath>

enum struct HistoryInterpolation {
    nearest,
    linear,
    slerp,
    // hermite,
    // lagrange,
    // chebyshev,
};

enum struct HistoryExtrapolation {
    reject,
    hold,
    constant_velocity, // both linear and angular
    // TODO: add more
};

struct HistorySampleOptions {
    HistoryInterpolation tr_interp = HistoryInterpolation::linear;
    HistoryInterpolation att_interp = HistoryInterpolation::slerp;
    HistoryExtrapolation tr_extrap = HistoryExtrapolation::reject;
    HistoryExtrapolation att_extrap = HistoryExtrapolation::reject;
    f64 tol = tol12;
};

inline vecXd interp_linear(f64 t, vecXd v1, f64 t1, vecXd v2, f64 t2, f64 tol = tol12) {
    // lerp
    if (v1.size() != v2.size()) {
        return vecXd{};
    }
    if (std::abs(t1 - t2) <= tol) {
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
        return vec4d0;
    }
    if (std::abs(t1 - t2) <= tol) {
        return vec4d0;
    }
    f64 alpha = (t - t1) / (t2 - t1);

    f64 dot = q1.dot(q2); // cos(theta)
    if (dot < 0.0) {
        q1 = -q1;
        dot = -dot;
    }

    if (dot > 1.0 - tol * 1000.0) {
        vec4d qt = interp_linear(t, q1, t1, q2, t2);
        if (qt.norm() < tol) {
            return vec4d0;
        }
        return qt.normalized();
    }

    dot = std::clamp(dot, -1.0, 1.0);

    f64 theta = std::acos(dot);
    f64 theta_t = theta * alpha;

    f64 sin_theta = std::sin(theta);
    f64 sin_theta_t = std::sin(theta_t);

    f64 s1 = std::cos(theta_t) - dot * sin_theta_t / sin_theta;
    f64 s2 = sin_theta_t / sin_theta;
    // Eigen::Quaterniond qa(q1.w(), q1.x(), q1.y(), q1.z());
    // Eigen::Quaterniond qb(q2.w(), q2.x(), q2.y(), q2.z());
    // Eigen::Quaterniond qt_eig = qa.slerp(alpha, qb);
    // return vec4d{qt_eig.x(), qt_eig.y(), qt_eig.z(), qt_eig.w()};

    vec4d qt = (s1 * q1) + (s2 * q2);
    if (qt.norm() <= tol) {
        return vec4d0;
    }

    return qt.normalized();
}
