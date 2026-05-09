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

const vec2d vec2d0 = vec2d::Zero();
const vec3d vec3d0 = vec3d::Zero();
const vec4d vec4d0 = vec4d::Zero();
const vec5d vec5d0 = vec5d::Zero();
const vec6d vec6d0 = vec6d::Zero();
const vec7d vec7d0 = vec7d::Zero();
const vec8d vec8d0 = vec8d::Zero();
const vec9d vec9d0 = vec9d::Zero();
const vec10d vec10d0 = vec10d::Zero();
const vec11d vec11d0 = vec11d::Zero();
const vec12d vec12d0 = vec12d::Zero();

const mat2d mat2d0 = mat2d::Zero();
const mat3d mat3d0 = mat3d::Zero();
const mat4d mat4d0 = mat4d::Zero();
const mat5d mat5d0 = mat5d::Zero();
const mat6d mat6d0 = mat6d::Zero();
const mat7d mat7d0 = mat7d::Zero();
const mat8d mat8d0 = mat8d::Zero();
const mat9d mat9d0 = mat9d::Zero();
const mat10d mat10d0 = mat10d::Zero();
const mat11d mat11d0 = mat11d::Zero();
const mat12d mat12d0 = mat12d::Zero();

const vec2d vec2d1 = vec2d::Identity();
const vec3d vec3d1 = vec3d::Identity();
const vec4d vec4d1 = vec4d::Identity();
const vec5d vec5d1 = vec5d::Identity();
const vec6d vec6d1 = vec6d::Identity();
const vec7d vec7d1 = vec7d::Identity();
const vec8d vec8d1 = vec8d::Identity();
const vec9d vec9d1 = vec9d::Identity();
const vec10d vec10d1 = vec10d::Identity();
const vec11d vec11d1 = vec11d::Identity();
const vec12d vec12d1 = vec12d::Identity();

const mat2d mat2d1 = mat2d::Identity();
const mat3d mat3d1 = mat3d::Identity();
const mat4d mat4d1 = mat4d::Identity();
const mat5d mat5d1 = mat5d::Identity();
const mat6d mat6d1 = mat6d::Identity();
const mat7d mat7d1 = mat7d::Identity();
const mat8d mat8d1 = mat8d::Identity();
const mat9d mat9d1 = mat9d::Identity();
const mat10d mat10d1 = mat10d::Identity();
const mat11d mat11d1 = mat11d::Identity();
const mat12d mat12d1 = mat12d::Identity();

const vec2f vec2f1 = vec2f::Identity();
const vec3f vec3f1 = vec3f::Identity();
const vec4f vec4f1 = vec4f::Identity();
const vec5f vec5f1 = vec5f::Identity();
const vec6f vec6f1 = vec6f::Identity();
const vec7f vec7f1 = vec7f::Identity();
const vec8f vec8f1 = vec8f::Identity();
const vec9f vec9f1 = vec9f::Identity();
const vec10f vec10f1 = vec10f::Identity();
const vec11f vec11f1 = vec11f::Identity();
const vec12f vec12f1 = vec12f::Identity();

const mat2f mat2f1 = mat2f::Identity();
const mat3f mat3f1 = mat3f::Identity();
const mat4f mat4f1 = mat4f::Identity();
const mat5f mat5f1 = mat5f::Identity();
const mat6f mat6f1 = mat6f::Identity();
const mat7f mat7f1 = mat7f::Identity();
const mat8f mat8f1 = mat8f::Identity();
const mat9f mat9f1 = mat9f::Identity();
const mat10f mat10f1 = mat10f::Identity();
const mat11f mat11f1 = mat11f::Identity();
const mat12f mat12f1 = mat12f::Identity();

const vec2f vec2f0 = vec2f::Zero();
const vec3f vec3f0 = vec3f::Zero();
const vec4f vec4f0 = vec4f::Zero();
const vec5f vec5f0 = vec5f::Zero();
const vec6f vec6f0 = vec6f::Zero();
const vec7f vec7f0 = vec7f::Zero();
const vec8f vec8f0 = vec8f::Zero();
const vec9f vec9f0 = vec9f::Zero();
const vec10f vec10f0 = vec10f::Zero();
const vec11f vec11f0 = vec11f::Zero();
const vec12f vec12f0 = vec12f::Zero();

const mat2f mat2f0 = mat2f::Zero();
const mat3f mat3f0 = mat3f::Zero();
const mat4f mat4f0 = mat4f::Zero();
const mat5f mat5f0 = mat5f::Zero();
const mat6f mat6f0 = mat6f::Zero();
const mat7f mat7f0 = mat7f::Zero();
const mat8f mat8f0 = mat8f::Zero();
const mat9f mat9f0 = mat9f::Zero();
const mat10f mat10f0 = mat10f::Zero();
const mat11f mat11f0 = mat11f::Zero();
const mat12f mat12f0 = mat12f::Zero();

const vec3d origin = vec3d0;
const vec3d axis_x = vec3d{1.0, 0.0, 0.0};
const vec3d axis_y = vec3d{0.0, 1.0, 0.0};
const vec3d axis_z = vec3d{0.0, 0.0, 1.0};

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
    return m.template cast<64>();
}
