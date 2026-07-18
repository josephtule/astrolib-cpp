// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0
/**
 * @file math.hpp
 */

#pragma once

#include "util/constants.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"
#include <cmath>

/**
 * @brief Computes the power of an input
 *
 * @tparam T
 * @tparam I
 * @param x
 * @param n
 * @return T
 */
template <typename T, typename I>
inline T pow_Ti(T x, I n) {
    if (n == 0) return T(1.0);
    if (n < 0) {
        x = T(1.0) / x;
        n = -n;
    }

    T result = T(1.0);
    while (n) {
        if (n & 1) result *= x;
        x *= x;
        n >>= 1;
    }
    return result;
}

template <typename T>
inline T sign(T x, T eps = T(0)) {
    if (std::abs(x) <= eps) return T(0);
    return (x > T(0)) ? T(1) : T(-1);
}

template <typename T>
inline T wrap_pi(T a) {
    return std::atan2(std::sin(a), std::cos(a));
};

template <typename T>
inline T deg(T val) {
    return rad_to_deg * val;
}
template <typename T>
inline T rad(T val) {
    return deg_to_rad * val;
}

template <typename T>
inline svec<T> vieta(const eref<const vecX<T>>& poles) {
    // Vieta's formula for real or complex poles (works for any scalar T)
    svec<T> poly = {T(1)}; // start with 1

    for (int i = 0; i < poles.size(); ++i) {
        svec<T> new_poly(poly.size() + 1, T(0));
        for (i32 j = 0; j < poly.size(); ++j) {
            new_poly(j) += poly(j);                 // coefficient without this root
            new_poly(j + 1) += -poles(i) * poly(j); // include this root
        }
        poly = new_poly;
    }

    return poly;
}

template <typename T>
vecXd conv(ecref<vecX<T>> a, ecref<vecX<T>> b) {
    static_assert(std::is_arithmetic_v<T>, "conv<T>: T must be an arithmetic type");

    if (a.size() == 0 || b.size() == 0) return vecX<T>{};

    vecX<T> y = vecX<T>::Zero(a.size() + b.size() - 1);

    for (i32 i = 0; i < a.size(); ++i) {
        for (i32 j = 0; j < b.size(); ++j) {
            y(i + j) += a(i) * b(j);
        }
    }

    return y;
}

template <typename T>
constexpr T eps(T x = 1.) {
    static_assert(
        std::is_floating_point<T>::value,
        "eps(x) requires a floating-point type"
    );

    return std::nextafter(x, std::numeric_limits<T>::infinity()) - x;
}

template <typename T>
inline T vector_angle(vecX<T> a, vecX<T> b, UAngle u_out = UAngle::radian) {
    T a_mag = a.norm();
    T b_mag = b.norm();

    T angle = std::acos(a.dot(b) / (a_mag * b_mag));

    if (u_out != UAngle::radian) {
        angle = convert_angle(angle, UAngle::radian, u_out);
    }

    return angle;
}

template <typename T>
inline T wrap_angle(
    T angle,
    T min = 0,
    T max = 2.0 * pi,
    UAngle u_in = UAngle::radian,
    UAngle u_out = UAngle::radian
) {
    if (u_in != UAngle::radian) {
        angle = convert_angle(angle, u_in, UAngle::radian);
        min = convert_angle(min, u_in, UAngle::radian);
        max = convert_angle(max, u_in, UAngle::radian);
    }

    T width = max - min;
    T wrapped = std::fmod(angle - min, width);
    if (wrapped < 0) {
        wrapped += width;
    }
    angle = wrapped + min;

    if (u_out != UAngle::radian) {
        angle = convert_angle(angle, UAngle::radian, u_out);
    }

    return angle;
}

struct DMSAngle {
    f64 degree = 0.0;
    f64 minute = 0.0;
    f64 second = 0.0;
};

inline f64 DMS_to_deg(f64 degrees, f64 minutes, f64 seconds) {
    f64 degs = degrees + minutes / 60.0 + seconds / 3600.0;
    return degs;
}

inline f64 DMS_to_deg(DMSAngle dms) {
    return DMS_to_deg(dms.degree, dms.minute, dms.second);
}

template <typename T>
inline T sind(T val) {
    return std::sin(val * deg_to_rad);
}

template <typename T>
inline T cosd(T val) {
    return std::cos(val * deg_to_rad);
}

inline f64 clamp_unit(f64 x) { return std::clamp(x, -1.0, 1.0); }

template <typename T>
inline void normalize_quaternion_inplace(vec4d& q, T tol = tol12) {
    T q_mag = q.norm();
    if (q_mag <= tol) return;
    q /= q_mag;
}

template <class T>
inline mat3<T> matrix_cross(const vec3<T>& v) {
    const T T0 = T(0);

    return mat3<T>{T0, -v(2), v(1), v(2), T0, -v(0), -v(1), v(0), T0};
}

// scalar checks

template <class T>
inline bool finite_pos(T val) {
    return std::isfinite(val) && val > T(0);
}

template <class T>
inline bool finite_nonneg(T val) {
    return std::isfinite(val) && val >= T(0);
}

template <class T>
inline bool finite_nonzero(T val, T tol = tol12) {
    return std::isfinite(val) && std::isfinite(tol) && tol >= T(0) && std::abs(val) > tol;
}

template <class T>
inline bool finite_inrange(
    T val,
    T low,
    T up,
    bool include_low = false,
    bool include_up = false
) {
    if (!std::isfinite(val)) return false;
    if (!std::isfinite(low)) return false;
    if (!std::isfinite(up)) return false;
    if (include_low && include_up) {
        if (low > up) return false;
    } else {
        if (low >= up) return false;
    }

    bool above_low = include_low ? val >= low : val > low;
    if (!above_low) return false;

    bool below_up = include_up ? val <= up : val < up;
    if (!below_up) return false;

    return true;
}

// matrix/vector checks

template <class Derived>
inline bool finite_dense(const eig::DenseBase<Derived>& x) {
    return x.allFinite();
}

template <class Derived>
inline bool finite_nonempty_dense(const eig::DenseBase<Derived>& x) {
    return x.rows() > 0 && x.cols() > 0 && x.allFinite();
}

template <class Derived>
inline bool finite_vec(const eig::DenseBase<Derived>& v) {
    return finite_dense(v);
}

template <class Derived>
inline bool finite_nonempty_vec(const eig::DenseBase<Derived>& v) {
    return v.size() > 0 && v.allFinite();
}

template <class Derived>
inline bool finite_mat(const eig::DenseBase<Derived>& A) {
    return finite_dense(A);
}

template <class Derived>
inline bool finite_nonempty_mat(const eig::DenseBase<Derived>& A) {
    return finite_nonempty_dense(A);
}

template <class Derived, class T>
inline bool finite_norm_nonzero(const eig::MatrixBase<Derived>& v, T tol = tol12) {
    return finite_dense(v) && finite_nonzero(v.norm(), tol);
}

template <class Derived>
inline bool finite_norm_pos(const eig::MatrixBase<Derived>& v) {
    return finite_dense(v) && finite_pos(v.norm());
}

// integer operations

inline i32 floor_div(i32 a, i32 b) {
    i32 q = a / b;
    i32 r = a % b;

    if (r != 0 && ((r > 0) != (b > 0))) {
        --q;
    }

    return q;
}

inline i32 floor_mod(i32 a, i32 b) { return a - b * floor_div(a, b); }

inline i32 count_digits(i32 n) {
    if (n == 0) return 1;

    n = abs(n); // Handle negative numbers
    int count = 0;

    while (n > 0) {
        n /= 10;
        count++;
    }

    return count;
}

template <class T>
inline T ecc_from_semiaxes(T a, T b) {
    // a : semimajor axis, b : semiminor axis
    static_assert(std::is_floating_point_v<T>, "T must be a floating-point type");
    if (a <= T(0) || b < T(0) || b > a) return T(0);
    return std::sqrt(T(1.0) - b * b / (a * a));
}

template <class T>
inline T flat_from_semiaxes(T a, T b) {
    static_assert(std::is_floating_point_v<T>, "T must be a floating-point type");
    if (a <= T(0) || b < T(0) || b > a) return T(0);
    return (a - b) / a;
}

template <class T>
inline T mean_from_semiaxes(T a, T b) {
    static_assert(std::is_floating_point_v<T>, "T must be a floating-point type");
    if (a <= T(0) || b < T(0) || b > a) return T(0);
    return std::pow(a * a * b, T(1) / T(3));
}

template <class T>
inline mat3<T> inertia_PAT(const mat3<T>& I_o, T mass, const vec3<T>& offset) {
    static_assert(std::is_floating_point_v<T>, "T must be a floating-point type");
    const T d2 = offset.dot(offset);
    mat3<T> D = d2 * eig::Matrix<T, 3, 3>::Identity() - offset * offset.transpose();
    return I_o - mass * D;
}
