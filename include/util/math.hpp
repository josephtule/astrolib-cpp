#pragma once

#include "util/constants.hpp"
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
svec<T> vieta(const eref<const vecX<T>>& poles) {
    // Vieta's formula for real or complex poles (works for any scalar T)
    svec<T> poly = {T(1)}; // start with 1

    for (int i = 0; i < poles.size(); ++i) {
        svec<T> new_poly(poly.size() + 1, T(0));
        for (i32 j = 0; j < poly.size(); ++j) {
            new_poly[j] += poly[j];                 // coefficient without this root
            new_poly[j + 1] += -poles(i) * poly[j]; // include this root
        }
        poly = new_poly;
    }

    return poly;
}

template <typename T>
vecX<T> conv(const eref<const vecX<T>>& a, const eref<const vecX<T>>& b) {
    static_assert(std::is_arithmetic_v<T>, "conv<T>: T must be an arithmetic type");

    if (a.empty() || b.empty()) return vecX<T>{};

    std::vector<T> y(a.size() + b.size() - 1, T{});

    for (i32 i = 0; i < a.size(); ++i) {
        for (i32 j = 0; j < b.size(); ++j) {
            y[i + j] += a[i] * b[j];
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

// template<typename T>
// inline vec4<T> inv_q(ecref<vec4<T>>){
// }