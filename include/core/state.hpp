#pragma once

#include "util/vecdefs.hpp"

struct StateTr {
    vec3d r = vec3d::Zero(); // position
    vec3d v = vec3d::Zero(); // velocity
};

inline vec4d q_default{0.0, 0.0, 0.0, 1.0};
struct StateAtt {
    vec4d q = q_default;     // euler-parameter/quaternion
    vec3d w = vec3d::Zero(); // angular velocity
};

// Derivatives
struct DerivTr {
    vec3d dr = vec3d::Zero();
    vec3d dv = vec3d::Zero();
};

struct DerivAtt {
    vec4d dq = vec4d::Zero();
    vec3d dw = vec3d::Zero();
};

// Translation Operations
inline StateTr operator+(const StateTr& x1, const StateTr& x2) {
    return StateTr{.r = x1.r + x2.r, .v = x1.v + x2.v};
}
inline StateTr operator-(const StateTr& x1, const StateTr& x2) {
    return StateTr{.r = x1.r - x2.r, .v = x1.v - x2.v};
}
inline StateTr operator*(const StateTr& x, f64 scalar) {
    return StateTr{.r = x.r * scalar, .v = x.v * scalar};
}
inline StateTr operator*(f64 scalar, const StateTr& x) { return x * scalar; }

inline StateTr operator+(const StateTr& x, const DerivTr& mdx) {
    return StateTr{.r = x.r + mdx.dr, .v = x.v + mdx.dv};
}
inline StateTr operator+(const DerivTr& mdx, const StateTr& x) { return x + mdx; }
inline StateTr& operator+=(StateTr& x1, const StateTr& x2) {
    x1.r += x2.r;
    x1.v += x2.v;
    return x1;
}
inline StateTr& operator+=(StateTr& x, const DerivTr& dx) {
    x.r += dx.dr;
    x.v += dx.dv;
    return x;
}
inline StateTr operator-(const StateTr& x) { return StateTr{.r = -x.r, .v = -x.v}; }

inline DerivTr operator+(const DerivTr& dx1, const DerivTr& dx2) {
    return DerivTr{.dr = dx1.dr + dx2.dr, .dv = dx1.dv + dx2.dv};
}
inline DerivTr operator-(const DerivTr& dx1, const DerivTr& dx2) {
    return DerivTr{.dr = dx1.dr - dx2.dr, .dv = dx1.dv - dx2.dv};
}
inline DerivTr operator*(const DerivTr& dx, f64 scalar) {
    return DerivTr{.dr = dx.dr * scalar, .dv = dx.dv * scalar};
}
inline DerivTr operator*(f64 scalar, const DerivTr& dx) { return dx * scalar; }
inline DerivTr operator-(const DerivTr& dx) {
    return DerivTr{.dr = -dx.dr, .dv = -dx.dv};
}

inline vec6d statetr_to_vec6d(const StateTr& x) {
    vec6d out;
    out << x.r, x.v;
    return out;
}
inline vec6d derivtr_to_vec6d(const DerivTr& x) {
    vec6d out;
    out << x.dr, x.dv;
    return out;
}
inline StateTr vec6_to_statetr(const vec6d& x) {
    StateTr out;
    out.r = x.segment<3>(0);
    out.v = x.segment<3>(3);
    return out;
};
inline DerivTr vec6_to_derivtr(const vec6d& dx) {
    DerivTr out;
    out.dr = dx.segment<3>(0);
    out.dv = dx.segment<3>(3);
    return out;
};

inline vec7d stateatt_to_vec7d(const StateAtt& x) {
    vec7d out;
    out << x.q, x.w;
    return out;
}
inline vec7d derivatt_to_vec7d(const DerivAtt& x) {
    vec7d out;
    out << x.dq, x.dw;
    return out;
}

// Attitude Operations
