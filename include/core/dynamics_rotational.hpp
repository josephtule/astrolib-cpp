#pragma once

#include "core/state.hpp"
#include "util/vecdefs.hpp"

inline vec4d k_eulerparams(ecref<vec4d> q, ecref<vec3d> w) {
    vec4d dqdt = vec4d0;
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

inline vec3d d_angularvelocity_nPA(const vec3d& w, const mat3d& I, const mat3d& I_inv) {
    // does not assume inertia matrix is diagonal
    // specific angular momentum
    vec3d h = {
        I(0, 0) * w(0) + I(0, 1) * w(1) + I(0, 2) * w(2),
        I(1, 0) * w(0) + I(1, 1) * w(1) + I(1, 2) * w(2),
        I(2, 0) * w(0) + I(2, 1) * w(1) + I(2, 2) * w(2),
    };

    vec3d dhdt_tf = vec3d{
        h(2) * w(1) - h(1) * w(2),
        h(0) * w(2) - h(2) * w(0),
        h(1) * w(0) - h(0) * w(1),
    }; // torque-free

    // vec3d dhdt = torque - dhdt_tf;
    vec3d dhdt = -dhdt_tf;
    
    // vec3T<T> dwdt = I_inv * (torque - dhdt_tf);
    vec3d dwdt = vec3d{
        I_inv(0, 0) * dhdt(0) + I_inv(0, 1) * dhdt(1) + I_inv(0, 2) * dhdt(2),
        I_inv(1, 0) * dhdt(0) + I_inv(1, 1) * dhdt(1) + I_inv(1, 2) * dhdt(2),
        I_inv(2, 0) * dhdt(0) + I_inv(2, 1) * dhdt(1) + I_inv(2, 2) * dhdt(2),
    };

    return dwdt;
}

inline DerivAtt d_rigidbody(f64 t, const StateAtt& x, const mat3d& I) {
    // this assumes the inertia matrix is diagonal (body frame = principle axes)
    DerivAtt dx;
    vec4d q = x.q;
    vec3d w = x.w;

    dx.dq = k_eulerparams(q, w);
    dx.dw = d_angularvelocity_PA(w, I);

    return dx;
}
inline DerivAtt d_rigidbody(
    f64 t,
    const StateAtt& x,
    const mat3d& I,
    const mat3d& I_inv
) {
    // this assumes the inertia matrix is diagonal (body frame = principle axes)
    DerivAtt dx;
    vec4d q = x.q;
    vec3d w = x.w;

    dx.dq = k_eulerparams(q, w);
    dx.dw = d_angularvelocity_nPA(w, I, I_inv);

    return dx;
}

// TODO: add override for rotation of celestial bodies (With Earth Orientation
// Parameters/etc.)