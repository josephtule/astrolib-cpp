#pragma once

#include "util/typedefs.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

namespace eig = Eigen;

// clang-format off
template <typename M> using eref = eig::Ref<M>;
template <typename M> using ecref = const eig::Ref<const M>;
constexpr auto eUp = eig::Upper;
constexpr auto eSUp = eig::StrictlyUpper;
constexpr auto eLo = eig::Lower;
constexpr auto eSLo = eig::StrictlyLower;

template <typename T, int N, int M> using mat = eig::Matrix<T, N, M>;
template <int N, int M> using matd = eig::Matrix<f64, N, M>;
template <int N, int M> using matf = eig::Matrix<f32, N, M>;
template <typename T> using matX = eig::MatrixX<T>;
using matXd = eig::MatrixXd;
using matXf = eig::MatrixXf;

template <typename T, int N> using vec = eig::Vector<T, N>;
template <int N> using vecd = eig::Vector<f64, N>;
template <int N> using vecf = eig::Vector<f32, N>;
template <typename T> using vecX = eig::VectorX<T>;
using vecXd = eig::VectorXd;
using vecXf = eig::VectorXf;

template <typename T> using vec2 = eig::Vector<T, 2>;
template <typename T> using vec3 = eig::Vector<T, 3>;
template <typename T> using vec4 = eig::Vector<T, 4>;
template <typename T> using vec5 = eig::Vector<T, 5>;
template <typename T> using vec6 = eig::Vector<T, 6>;
template <typename T> using vec7 = eig::Vector<T, 7>;
template <typename T> using vec8 = eig::Vector<T, 8>;
template <typename T> using vec9 = eig::Vector<T, 9>;
template <typename T> using vec10 = eig::Vector<T, 10>;
template <typename T> using vec11 = eig::Vector<T, 11>;
template <typename T> using vec12 = eig::Vector<T, 12>;

template <typename T> using mat2 = eig::Matrix<T, 2, 2>;
template <typename T> using mat3 = eig::Matrix<T, 3, 3>;
template <typename T> using mat4 = eig::Matrix<T, 4, 4>;
template <typename T> using mat5 = eig::Matrix<T, 5, 5>;
template <typename T> using mat6 = eig::Matrix<T, 6, 6>;
template <typename T> using mat7 = eig::Matrix<T, 7, 7>;
template <typename T> using mat8 = eig::Matrix<T, 8, 8>;
template <typename T> using mat9 = eig::Matrix<T, 9, 9>;
template <typename T> using mat10 = eig::Matrix<T, 10, 10>;
template <typename T> using mat11 = eig::Matrix<T, 11, 11>;
template <typename T> using mat12 = eig::Matrix<T, 12, 12>;
// clang-format on

using vec2d = eig::Vector<f64, 2>;
using vec3d = eig::Vector<f64, 3>;
using vec4d = eig::Vector<f64, 4>;
using vec5d = eig::Vector<f64, 5>;
using vec6d = eig::Vector<f64, 6>;
using vec7d = eig::Vector<f64, 7>;
using vec8d = eig::Vector<f64, 8>;
using vec9d = eig::Vector<f64, 9>;
using vec10d = eig::Vector<f64, 10>;
using vec11d = eig::Vector<f64, 11>;
using vec12d = eig::Vector<f64, 12>;
using vec13d = eig::Vector<f64, 13>;
using vec14d = eig::Vector<f64, 14>;
using vec15d = eig::Vector<f64, 15>;
using vec16d = eig::Vector<f64, 16>;
using vec17d = eig::Vector<f64, 17>;
using vec18d = eig::Vector<f64, 18>;
using vec19d = eig::Vector<f64, 19>;
using vec20d = eig::Vector<f64, 20>;
using vec21d = eig::Vector<f64, 21>;
using vec22d = eig::Vector<f64, 22>;
using vec23d = eig::Vector<f64, 23>;
using vec24d = eig::Vector<f64, 24>;

using vec2f = eig::Vector<f32, 2>;
using vec3f = eig::Vector<f32, 3>;
using vec4f = eig::Vector<f32, 4>;
using vec5f = eig::Vector<f32, 5>;
using vec6f = eig::Vector<f32, 6>;
using vec7f = eig::Vector<f32, 7>;
using vec8f = eig::Vector<f32, 8>;
using vec9f = eig::Vector<f32, 9>;
using vec10f = eig::Vector<f32, 10>;
using vec11f = eig::Vector<f32, 11>;
using vec12f = eig::Vector<f32, 12>;

using mat2d = eig::Matrix<f64, 2, 2>;
using mat3d = eig::Matrix<f64, 3, 3>;
using mat4d = eig::Matrix<f64, 4, 4>;
using mat5d = eig::Matrix<f64, 5, 5>;
using mat6d = eig::Matrix<f64, 6, 6>;
using mat7d = eig::Matrix<f64, 7, 7>;
using mat8d = eig::Matrix<f64, 8, 8>;
using mat9d = eig::Matrix<f64, 9, 9>;
using mat10d = eig::Matrix<f64, 10, 10>;
using mat11d = eig::Matrix<f64, 11, 11>;
using mat12d = eig::Matrix<f64, 12, 12>;

using mat2f = eig::Matrix<f32, 2, 2>;
using mat3f = eig::Matrix<f32, 3, 3>;
using mat4f = eig::Matrix<f32, 4, 4>;
using mat5f = eig::Matrix<f32, 5, 5>;
using mat6f = eig::Matrix<f32, 6, 6>;
using mat7f = eig::Matrix<f32, 7, 7>;
using mat8f = eig::Matrix<f32, 8, 8>;
using mat9f = eig::Matrix<f32, 9, 9>;
using mat10f = eig::Matrix<f32, 10, 10>;
using mat11f = eig::Matrix<f32, 11, 11>;
using mat12f = eig::Matrix<f32, 12, 12>;

inline const vec2d vec2d0 = vec2d::Zero();
inline const vec3d vec3d0 = vec3d::Zero();
inline const vec4d vec4d0 = vec4d::Zero();
inline const vec5d vec5d0 = vec5d::Zero();
inline const vec6d vec6d0 = vec6d::Zero();
inline const vec7d vec7d0 = vec7d::Zero();
inline const vec8d vec8d0 = vec8d::Zero();
inline const vec9d vec9d0 = vec9d::Zero();
inline const vec10d vec10d0 = vec10d::Zero();
inline const vec11d vec11d0 = vec11d::Zero();
inline const vec12d vec12d0 = vec12d::Zero();

inline const mat2d mat2d0 = mat2d::Zero();
inline const mat3d mat3d0 = mat3d::Zero();
inline const mat4d mat4d0 = mat4d::Zero();
inline const mat5d mat5d0 = mat5d::Zero();
inline const mat6d mat6d0 = mat6d::Zero();
inline const mat7d mat7d0 = mat7d::Zero();
inline const mat8d mat8d0 = mat8d::Zero();
inline const mat9d mat9d0 = mat9d::Zero();
inline const mat10d mat10d0 = mat10d::Zero();
inline const mat11d mat11d0 = mat11d::Zero();
inline const mat12d mat12d0 = mat12d::Zero();

inline const vec2d vec2d1 = vec2d::Identity();
inline const vec3d vec3d1 = vec3d::Identity();
inline const vec4d vec4d1 = vec4d::Identity();
inline const vec5d vec5d1 = vec5d::Identity();
inline const vec6d vec6d1 = vec6d::Identity();
inline const vec7d vec7d1 = vec7d::Identity();
inline const vec8d vec8d1 = vec8d::Identity();
inline const vec9d vec9d1 = vec9d::Identity();
inline const vec10d vec10d1 = vec10d::Identity();
inline const vec11d vec11d1 = vec11d::Identity();
inline const vec12d vec12d1 = vec12d::Identity();

inline const mat2d mat2d1 = mat2d::Identity();
inline const mat3d mat3d1 = mat3d::Identity();
inline const mat4d mat4d1 = mat4d::Identity();
inline const mat5d mat5d1 = mat5d::Identity();
inline const mat6d mat6d1 = mat6d::Identity();
inline const mat7d mat7d1 = mat7d::Identity();
inline const mat8d mat8d1 = mat8d::Identity();
inline const mat9d mat9d1 = mat9d::Identity();
inline const mat10d mat10d1 = mat10d::Identity();
inline const mat11d mat11d1 = mat11d::Identity();
inline const mat12d mat12d1 = mat12d::Identity();

inline const vec2f vec2f1 = vec2f::Identity();
inline const vec3f vec3f1 = vec3f::Identity();
inline const vec4f vec4f1 = vec4f::Identity();
inline const vec5f vec5f1 = vec5f::Identity();
inline const vec6f vec6f1 = vec6f::Identity();
inline const vec7f vec7f1 = vec7f::Identity();
inline const vec8f vec8f1 = vec8f::Identity();
inline const vec9f vec9f1 = vec9f::Identity();
inline const vec10f vec10f1 = vec10f::Identity();
inline const vec11f vec11f1 = vec11f::Identity();
inline const vec12f vec12f1 = vec12f::Identity();

inline const mat2f mat2f1 = mat2f::Identity();
inline const mat3f mat3f1 = mat3f::Identity();
inline const mat4f mat4f1 = mat4f::Identity();
inline const mat5f mat5f1 = mat5f::Identity();
inline const mat6f mat6f1 = mat6f::Identity();
inline const mat7f mat7f1 = mat7f::Identity();
inline const mat8f mat8f1 = mat8f::Identity();
inline const mat9f mat9f1 = mat9f::Identity();
inline const mat10f mat10f1 = mat10f::Identity();
inline const mat11f mat11f1 = mat11f::Identity();
inline const mat12f mat12f1 = mat12f::Identity();

inline const vec2f vec2f0 = vec2f::Zero();
inline const vec3f vec3f0 = vec3f::Zero();
inline const vec4f vec4f0 = vec4f::Zero();
inline const vec5f vec5f0 = vec5f::Zero();
inline const vec6f vec6f0 = vec6f::Zero();
inline const vec7f vec7f0 = vec7f::Zero();
inline const vec8f vec8f0 = vec8f::Zero();
inline const vec9f vec9f0 = vec9f::Zero();
inline const vec10f vec10f0 = vec10f::Zero();
inline const vec11f vec11f0 = vec11f::Zero();
inline const vec12f vec12f0 = vec12f::Zero();

inline const mat2f mat2f0 = mat2f::Zero();
inline const mat3f mat3f0 = mat3f::Zero();
inline const mat4f mat4f0 = mat4f::Zero();
inline const mat5f mat5f0 = mat5f::Zero();
inline const mat6f mat6f0 = mat6f::Zero();
inline const mat7f mat7f0 = mat7f::Zero();
inline const mat8f mat8f0 = mat8f::Zero();
inline const mat9f mat9f0 = mat9f::Zero();
inline const mat10f mat10f0 = mat10f::Zero();
inline const mat11f mat11f0 = mat11f::Zero();
inline const mat12f mat12f0 = mat12f::Zero();

inline const vec3d origin = vec3d0;
inline const vec3f originf = vec3f0;
inline const vec3d axis_x = vec3d{1.0, 0.0, 0.0};
inline const vec3d axis_y = vec3d{0.0, 1.0, 0.0};
inline const vec3d axis_z = vec3d{0.0, 0.0, 1.0};
inline const vec3f axis_xf = vec3f{1.0f, 0.0f, 0.0f};
inline const vec3f axis_yf = vec3f{0.0f, 1.0f, 0.0f};
inline const vec3f axis_zf = vec3f{0.0f, 0.0f, 1.0f};

template <typename T, int N>
inline vec<f32, N> vf32(const vec<T, N>& v) {
    return v.template cast<f32>();
}
template <class T, int N, int M>
inline mat<f32, N, M> mf32(const mat<T, N, M>& m) {
    return m.template cast<f32>();
}

template <typename T, int N>
inline vec<f64, N> vf64(const vec<T, N>& v) {
    return v.template cast<f64>();
}
template <class T, int N, int M>
inline mat<f64, N, M> mf64(const mat<T, N, M>& m) {
    return m.template cast<f64>();
}
