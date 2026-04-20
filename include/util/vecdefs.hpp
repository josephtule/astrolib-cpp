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

const vec3d axis_x = vec3d{1.0, 0.0, 0.0};
const vec3d axis_y = vec3d{0.0, 1.0, 0.0};
const vec3d axis_z = vec3d{0.0, 0.0, 1.0};

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
