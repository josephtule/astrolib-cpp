#pragma once

#include "core/body.hpp"
#include "util/constants.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

// NOTE: All rotations are passive (DCM) unless specified
// NOTE: Use passive (DCM) for frame transforms, active for rendering

template <typename T>
inline mat3<T> rotX(T angle, UAngle uin = UAngle::radian) {
    angle = convert_angle(angle, uin, UAngle::radian);

    T c = std::cos(angle), s = std::sin(angle);

    mat3<T> R;
    R << 1, 0, 0, 0, c, s, 0, -s, c;
    return R;
}
template <typename T>
inline mat3<T> rotX_active(T angle, UAngle uin = UAngle::radian) {
    return rotX(angle, uin).transpose();
}
template <typename T>
inline mat3<T> rotY(T angle, UAngle uin = UAngle::radian) {
    angle = convert_angle(angle, uin, UAngle::radian);

    T c = std::cos(angle), s = std::sin(angle);

    mat3<T> R;
    R << c, 0, -s, 0, 1, 0, s, 0, c;
    return R;
}
template <typename T>
inline mat3<T> rotY_active(T angle, UAngle uin = UAngle::radian) {
    return rotY(angle, uin).transpose();
}
template <typename T>
inline mat3<T> rotZ(T angle, UAngle uin = UAngle::radian) {
    angle = convert_angle(angle, uin, UAngle::radian);

    T c = std::cos(angle), s = std::sin(angle);

    mat3<T> R;
    R << c, s, 0, -s, c, 0, 0, 0, 1;
    return R;
}
template <typename T>
inline mat3<T> rotZ_active(T angle, UAngle uin = UAngle::radian) {
    return rotZ(angle, uin).transpose();
}

enum struct RotAxis : i32 { x = 1, y = 2, z = 3 };

template <typename T>
inline mat3<T> rot(T angle, RotAxis axis, UAngle uin = UAngle::radian) {
    mat3<T> R;
    switch (axis) {
    case RotAxis::x: R = rotX(angle, uin); break;
    case RotAxis::y: R = rotY(angle, uin); break;
    case RotAxis::z: R = rotZ(angle, uin); break;
    }
    return R;
}
template <typename T>
inline mat3<T> rot_active(T angle, RotAxis axis, UAngle uin = UAngle::radian) {
    mat3<T> R;
    switch (axis) {
    case RotAxis::x: R = rotX(angle, uin).transpose(); break;
    case RotAxis::y: R = rotY(angle, uin).transpose(); break;
    case RotAxis::z: R = rotZ(angle, uin).transpose(); break;
    }
    return R;
}

template <typename T>
inline mat3<T> ea_to_dcm(
    vec3<T> angles,
    std::array<RotAxis, 3> seq,
    UAngle units_in = UAngle::radian
) {
    for (int i = 0; i < 3; i++) {
        angles(i) = convert_angle(angles(i), units_in, UAngle::radian);
    }

    mat3<T> R = rot(angles(2), seq[2]) * rot(angles(1), seq[1]) * rot(angles(0), seq[0]);

    return R;
}

// Euler-Parameters / Unit-Quaternions
template <typename T>
inline mat3<T> ep_to_dcm(const vec4<T>& q) {
    mat3<T> R;
    T q1 = q(0);
    T q2 = q(1);
    T q3 = q(2);
    T q4 = q(3);

    T q11 = q1 * q1;
    T q22 = q2 * q2;
    T q33 = q3 * q3;
    T q44 = q4 * q4;
    T q12 = q1 * q2;
    T q13 = q1 * q3;
    T q14 = q1 * q4;
    T q23 = q2 * q3;
    T q24 = q2 * q4;
    T q34 = q3 * q4;

    R(0, 0) = +q11 - q22 - q33 + q44; // = 1.0 - 2.0 * q22 - 2.0 * q33;
    R(0, 1) = 2.0 * (q12 + q34);
    R(0, 2) = 2.0 * (q13 - q24);
    R(1, 0) = 2.0 * (q12 - q34);
    R(1, 1) = -q11 + q22 - q33 + q44; // = 1.0 - 2.0 * q11 - 2.0 * q33;
    R(1, 2) = 2.0 * (q23 + q14);
    R(2, 0) = 2.0 * (q13 + q24);
    R(2, 1) = 2.0 * (q23 - q14);
    R(2, 2) = -q11 - q22 + q33 + q44; // = 1.0 - 2.0 * q11 - 2.0 * q22;

    return R;
}

template <typename T>
inline vec4<T> ep_mult(const vec4<T>& a, const vec4<T>& b) {
    vec3<T> av = a.template segment<3>(0);
    vec3<T> bv = b.template segment<3>(0);
    T as = a(3);
    T bs = b(3);

    T qs = as * bs - av.dot(bv);
    vec3<T> qv = as * bv + bs * av + av.cross(bv);
    vec4<T> q;
    q << qv, qs;
    return q;
}

template <typename T>
inline vec4<T> ep_conj(const vec4<T>& q) {
    vec4<T> qc = {-q(0), -q(1), -q(2), q(3)};
    return qc;
}

template <typename T>
inline vec3<T> ep_rotate_active(const vec4<T>& q, const vec3<T>& v) {
    // this is an active rotation B -> N
    // assume q is normalized
    vec4<T> vq = {v(0), v(1), v(2), static_cast<T>(0)};
    vec4<T> qc = ep_conj(q);
    vec4<T> vp = ep_mult(ep_mult(q, vq), qc);

    return vec3<T>{vp(0), vp(1), vp(2)};
}
template <typename T>
inline vec3<T> ep_rotate_fast_active(const vec4<T>& q, const vec3<T>& v) {
    // this is an active rotation [NB]: B -> N
    // assume q is normalized
    vec3<T> u = q.template segment<3>(0);
    T s = q(3);
    vec3<T> uxv = u.cross(v);
    vec3<T> vp = v + 2 * s * uxv + 2 * (u.dot(v) * u - u.dot(u) * v);

    return vp;
}

template <typename T>
inline vec3<T> ep_rotate_passive(const vec4<T>& q, const vec3<T>& v) {
    // this is an passive rotation [BN]: N -> B
    // assume q is normalized
    return ep_rotate_active(ep_conj(q), v);
}
template <typename T>
inline vec3<T> ep_rotate_fast_passive(const vec4<T>& q, const vec3<T>& v) {
    // this is an passive rotation [BN]: N -> B
    // assume q is normalized
    return ep_rotate_fast_active(ep_conj(q), v);
}

template <typename T>
inline vec4<T> dcm_to_ep(const mat3<T>& R) {
    T Rtr = R.trace();
    T q1, q2, q3, q4;
    // Shepperd's Selection Algorithm
    if ((Rtr > R(0, 0)) && (Rtr > R(1, 1)) && (Rtr > R(2, 2))) {
        q4 = std::sqrt((1. + Rtr) / 4.);
        q1 = (R(1, 2) - R(2, 1)) / (4. * q4);
        q2 = (R(2, 0) - R(0, 2)) / (4. * q4);
        q3 = (R(0, 1) - R(1, 0)) / (4. * q4);
    } else if (R(0, 0) > R(1, 1) && R(0, 0) > R(2, 2)) {
        q1 = std::sqrt((1. + R(0, 0) - R(1, 1) - R(2, 2)) / 4.);
        q4 = (R(1, 2) - R(2, 1)) / (4. * q1);
        q2 = (R(0, 1) + R(1, 0)) / (4. * q1);
        q3 = (R(2, 0) + R(0, 2)) / (4. * q1);
    } else if (R(1, 1) > R(2, 2)) {
        q2 = std::sqrt((1. - R(0, 0) + R(1, 1) - R(2, 2)) / 4.);
        q4 = (R(2, 0) - R(0, 2)) / (4. * q2);
        q1 = (R(0, 1) + R(1, 0)) / (4. * q2);
        q3 = (R(1, 2) + R(2, 1)) / (4. * q2);
    } else {
        q3 = std::sqrt((1. - R(0, 0) - R(1, 1) + R(2, 2)) / 4.);
        q4 = (R(0, 1) - R(1, 0)) / (4. * q3);
        q1 = (R(2, 0) + R(0, 2)) / (4. * q3);
        q2 = (R(1, 2) + R(2, 1)) / (4. * q3);
    }
    vec4<T> q;
    q << q1, q2, q3, q4;

    return q;
}

inline vec3d bcbf_to_centric(ecref<vec3d> r_bcbf, UAngle u_out = UAngle::degree) {
    vec3d llr = vec3d0;

    f64 x = r_bcbf(0);
    f64 y = r_bcbf(1);
    f64 z = r_bcbf(2);

    f64 r = r_bcbf.norm();
    f64 r_tilde = std::sqrt(x * x + y * y);

    // Planetocentric latitude/longitude (radians)
    f64 longitude = (r_tilde > 0.0) ? std::atan2(y, x) : 0.0;
    f64 latitude = std::atan2(z, r_tilde);

    // Convert
    if (u_out == UAngle::degree) {
        longitude *= rad_to_deg;
        latitude *= rad_to_deg;
    }

    llr << latitude, longitude, r;
    return llr;
}

inline vec3d centric_to_bcbf(ecref<vec3d> llr, UAngle u_in = UAngle::degree) {
    vec3d r_bcbf = vec3d0;

    f64 latitude = llr(0);
    f64 longitude = llr(1);
    f64 r = llr(2);

    if (u_in != UAngle::radian) {
        latitude = convert_angle(latitude, u_in, UAngle::radian);
        longitude = convert_angle(longitude, u_in, UAngle::radian);
    }

    f64 slat = std::sin(latitude);
    f64 clat = std::cos(latitude);
    f64 slon = std::sin(longitude);
    f64 clon = std::cos(longitude);

    r_bcbf = vec3d{r * clat * clon, r * clat * slon, r * slat};

    return r_bcbf;
}

inline vec3d bcbf_to_detic(
    ecref<vec3d> r_bcbf,
    Celestial& body,
    UAngle u_out = UAngle::degree,
    f64 tol = tol_strict
) {
    vec3d llh = vec3d0;

    f64 x = r_bcbf(0);
    f64 y = r_bcbf(1);
    f64 z = r_bcbf(2);

    f64 r_tilde2 = x * x + y * y;
    f64 r_tilde = std::sqrt(r_tilde2);

    // Ellipsoid
    f64 a = body.semimajor_axis;
    f64 f = body.flattening;
    f64 b = a * (1.0 - f);

    f64 e2 = (a * a - b * b) / (a * a);
    f64 ep2 = (a * a - b * b) / (b * b);

    f64 longitude = (r_tilde > tol) ? std::atan2(y, x) : 0.0;

    f64 latitude, h;

    // Special handling near poles
    if (r_tilde < tol) {
        latitude = (z >= 0.0 ? pio2 : -pio2);
        h = std::abs(z) - b;
    } else {
        // Closed-form solution
        f64 F = 54.0 * b * b * z * z;
        f64 G = r_tilde2 + (1.0 - e2) * z * z - e2 * (a * a - b * b);

        f64 c = (std::abs(G) > tol) ? (e2 * e2 * F * r_tilde2) / (G * G * G) : 0.0;

        f64 s = std::cbrt(1.0 + c + std::sqrt(c * c + 2.0 * c));
        f64 P = F / (3.0 * (s + 1.0 / s + 1.0) * (s + 1.0 / s + 1.0) * G * G);
        f64 Q = std::sqrt(1.0 + 2.0 * e2 * e2 * P);

        f64 r0 = -P * e2 * r_tilde / (1.0 + Q)
                 + std::sqrt(
                     0.5 * a * a * (1.0 + 1.0 / Q)
                     - P * (1.0 - e2) * z * z / (Q * (1.0 + Q)) - 0.5 * P * r_tilde2
                 );

        f64 U = std::sqrt((r_tilde - e2 * r0) * (r_tilde - e2 * r0) + z * z);
        f64 V = std::sqrt((r_tilde - e2 * r0) * (r_tilde - e2 * r0) + (1.0 - e2) * z * z);
        f64 z0 = b * b * z / (a * V);

        latitude = std::atan2(z + ep2 * z0, r_tilde);
        h = U * (1.0 - b * b / (a * V));
    }

    if (u_out == UAngle::degree) {
        longitude *= rad_to_deg;
        latitude *= rad_to_deg;
    }

    llh << latitude, longitude, h;
    return llh;
}

inline vec3d detic_to_bcbf(
    ecref<vec3d> llh,
    const Celestial& body,
    UAngle u_in = UAngle::degree
) {
    vec3d r_bcbf = vec3d0;

    f64 latitude = llh(0);
    f64 longitude = llh(1);
    f64 h = llh(2);

    if (u_in != UAngle::radian) {
        latitude = convert_angle(latitude, u_in, UAngle::radian);
        longitude = convert_angle(longitude, u_in, UAngle::radian);
    }

    f64 slat = std::sin(latitude);
    f64 clat = std::cos(latitude);
    f64 slon = std::sin(longitude);
    f64 clon = std::cos(longitude);

    // Ellipsoid
    f64 a = body.semimajor_axis;
    f64 f = body.flattening;
    f64 e2 = 2.0 * f - f * f;

    f64 N = a / std::sqrt(1.0 - e2 * slat * slat);

    r_bcbf = vec3d{
        (N + h) * clat * clon,
        (N + h) * clat * slon,
        (N * (1.0 - e2) + h) * slat
    };

    return r_bcbf;
}

