#pragma once

#include "util/constants.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

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
inline T wrap_pi(f64 a) {
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
