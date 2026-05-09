#pragma once

#include "core/transform.hpp"
#include "raylib.h"
#include "raymath.h"
#include "util/vecdefs.hpp"

template <typename T>
inline vec3<T> rl_to_eig(const Vector3& v) {
    return vec3<T>{v.x, v.y, v.z};
}
template <typename T>
inline vec2<T> rl_to_eig(const Vector2& v) {
    return vec2<T>{v.x, v.y};
}
template <typename T>
inline vec4<T> rl_to_eig(const Vector4& v) {
    return vec4<T>{v.x, v.y, v.z, v.w};
}

template <typename T>
inline Vector2 eig_to_rl(const vec2<T>& v) {
    return Vector2(f32(v[0]), f32(v[1]));
}
template <typename T>
inline Vector3 eig_to_rl(const vec3<T>& v) {
    return Vector3{f32(v[0]), f32(v[1]), f32(v[2])};
}
template <typename T>
inline Vector4 eig_to_rl(const vec4<T>& v) {
    return Vector4(f32(v[0]), f32(v[1]), f32(v[2]), f32(v[3]));
}

template <class T>
inline void set_rotation(Matrix& M, const mat3<T>& R) {
    M.m0 = R(0, 0);
    M.m4 = R(0, 1);
    M.m8 = R(0, 2);

    M.m1 = R(1, 0);
    M.m5 = R(1, 1);
    M.m9 = R(1, 2);

    M.m2 = R(2, 0);
    M.m6 = R(2, 1);
    M.m10 = R(2, 2);
}

template <class T>
inline void set_rotation(Matrix& M, const vec4<T>& q) {
    mat3<T> R = ep_to_dcm(q);
    set_rotation(M, R);
}

template <class T>
inline void set_translation(Matrix& M, const vec3<T>& r) {
    M.m12 = r.x();
    M.m13 = r.y();
    M.m14 = r.z();
}

inline Matrix make_transform(
    const vec3f& r,
    const mat3f& R, // passive transform
    const vec3f& size,
    f32 scale = 1.0
) {
    mat3f S = size.asDiagonal() * scale; // scaled size matrix
    mat3f RS = R.transpose() * S;        // active/direct transform

    Matrix M = MatrixIdentity();
    set_rotation(M, RS);
    set_translation(M, r);

    return M;
}

inline Matrix make_transform(
    const StateTr& x_tr,
    const StateAtt& x_att,
    const vec3f& size,
    f32 scale = 1
) {
    mat3f S = size.asDiagonal() * scale;

    // q is passive N -> B.
    // ep_to_dcm(q) gives R_NB.
    // ep_to_dcm(conj(q)) gives R_BN.
    mat3f R = mf32(ep_to_dcm(ep_conj(x_att.q)));
    mat3f RS = R * S;

    vec3f r = vf32(x_tr.r);

    Matrix M = MatrixIdentity();
    set_rotation(M, RS);
    set_translation(M, r);

    return M;
}