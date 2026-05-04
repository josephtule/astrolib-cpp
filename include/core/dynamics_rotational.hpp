#pragma once

#include "core/state.hpp"
#include "util/vecdefs.hpp"

inline vec4d k_eulerparams(ecref<vec4d> q, ecref<vec3d> w) {
    vec4d dqdt = vec4d::Zero();
    dqdt
        = {w(2) * q(1) - w(1) * q(2) + w(0) * q(3),
           -w(2) * q(0) + w(0) * q(2) + w(1) * q(3),
           w(1) * q(0) - w(0) * q(1) + w(2) * q(3),
           -w(0) * q(0) - w(1) * q(1) - w(2) * q(2)};
    dqdt *= 0.5;
    return dqdt;
}
inline vec4d k_eulerparams(StateAtt x) { return k_eulerparams(x.q, x.w); }

inline vec3d d_angularvelocity_PA(ecref<vec3d> w, ecref<mat3d> I) {
    // this assumes the inertia matrix is diagonal (body frame = principle axes)
    vec3d dwdt = {
        ((I(1, 1) - I(2, 2)) * w(1) * w(2)) / I(0, 0),
        ((I(2, 2) - I(0, 0)) * w(0) * w(2)) / I(1, 1),
        ((I(0, 0) - I(1, 1)) * w(0) * w(1)) / I(2, 2),
    };

    return dwdt;
}
inline vec7d d_rigidbody_PA(f64 t, StateAtt& x, ecref<mat3d> I) {
    // this assumes the inertia matrix is diagonal (body frame = principle axes)
    vec4d q = x.q;
    vec3d w = x.w;

    vec4d dqdt = k_eulerparams(q, w);
    vec3d dwdt = d_angularvelocity_PA(w, I);

    vec7d dxdt;
    dxdt << dqdt, dwdt;

    return dxdt;
}

// TODO: add override for rotation of celestial bodies (With Earth Orientation
// Parameters/etc.)