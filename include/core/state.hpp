#pragma once

#include "util/vecdefs.hpp"

struct StateTr {
    vec3d r = vec3d0; // position, inertial
    vec3d v = vec3d0; // velocity, inertial
};

struct StateAtt {
    vec4d q = q_identity; // euler-parameter/quaternion, N -> B
    vec3d w = vec3d0;     // angular velocity, in body frame
};

// Derivatives
struct DerivTr {
    vec3d dr = vec3d0;
    vec3d dv = vec3d0;
};

struct DerivAtt {
    vec4d dq = vec4d0;
    vec3d dw = vec3d0;
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
inline StateTr operator/(const StateTr& x, f64 scalar) { return x * (1.0 / scalar); }
inline StateTr operator+(const StateTr& x, const DerivTr& dx) {
    return StateTr{.r = x.r + dx.dr, .v = x.v + dx.dv};
}
inline StateTr operator+(const DerivTr& dx, const StateTr& x) { return x + dx; }
inline StateTr operator-(const StateTr& x, const DerivTr& dx) {
    return {.r = x.r - dx.dr, .v = x.v - dx.dv};
}
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
inline DerivTr operator/(const DerivTr& dx, f64 scalar) { return dx * (1.0 / scalar); }
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
inline StateTr vec6d_to_statetr(const vec6d& x) {
    StateTr out;
    out.r = x.segment<3>(0);
    out.v = x.segment<3>(3);
    return out;
};
inline DerivTr vec6d_to_derivtr(const vec6d& dx) {
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
inline StateAtt operator+(const StateAtt& x1, const StateAtt& x2) {
    return StateAtt{.q = x1.q + x2.q, .w = x1.w + x2.w};
}
inline StateAtt operator-(const StateAtt& x1, const StateAtt& x2) {
    return StateAtt{.q = x1.q - x2.q, .w = x1.w - x2.w};
}
inline StateAtt operator*(const StateAtt& x, f64 scalar) {
    return StateAtt{.q = x.q * scalar, .w = x.w * scalar};
}
inline StateAtt operator*(f64 scalar, const StateAtt& x) { return x * scalar; }
inline StateAtt operator/(const StateAtt& x, f64 scalar) { return x * (1.0 / scalar); }
inline StateAtt operator+(const StateAtt& x, const DerivAtt& dx) {
    return StateAtt{.q = x.q + dx.dq, .w = x.w + dx.dw};
}
inline StateAtt operator+(const DerivAtt& dx, const StateAtt& x) { return x + dx; }
inline StateAtt operator-(const StateAtt& x, const DerivAtt& dx) {
    return StateAtt{.q = x.q - dx.dq, .w = x.w - dx.dw};
}
inline StateAtt& operator+=(StateAtt& x1, const StateAtt& x2) {
    x1.q += x2.q;
    x1.w += x2.w;
    return x1;
}
inline StateAtt& operator+=(StateAtt& x, const DerivAtt& dx) {
    x.q += dx.dq;
    x.w += dx.dw;
    return x;
}
inline StateAtt operator-(const StateAtt& x) { return StateAtt{.q = -x.q, .w = -x.w}; }

inline DerivAtt operator+(const DerivAtt& dx1, const DerivAtt& dx2) {
    return DerivAtt{.dq = dx1.dq + dx2.dq, .dw = dx1.dw + dx2.dw};
}
inline DerivAtt operator-(const DerivAtt& dx1, const DerivAtt& dx2) {
    return DerivAtt{.dq = dx1.dq - dx2.dq, .dw = dx1.dw - dx2.dw};
}
inline DerivAtt operator*(const DerivAtt& dx, f64 scalar) {
    return DerivAtt{.dq = dx.dq * scalar, .dw = dx.dw * scalar};
}
inline DerivAtt operator*(f64 scalar, const DerivAtt& dx) { return dx * scalar; }
inline DerivAtt operator/(const DerivAtt& dx, f64 scalar) { return dx * (1.0 / scalar); }
inline DerivAtt operator-(const DerivAtt& dx) {
    return DerivAtt{.dq = -dx.dq, .dw = -dx.dw};
}

inline StateTr operator*(const mat3d& R, const StateTr& x) {
    // NOTE: this is a pure rotation of coordinates, it does not return true velocity in
    // the new frame

    return StateTr{.r = R * x.r, .v = R * x.v};
}

inline DerivTr operator*(const mat3d& R, const DerivTr& dx) {
    // NOTE: this is a pure rotation of coordinates, it does not return true velocity in
    // the new frame

    return DerivTr{.dr = R * dx.dr, .dv = R * dx.dv};
}
