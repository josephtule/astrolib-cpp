#pragma once

#include <complex>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Type aliases
using f32 = float;
using f64 = double;
using c64 = std::complex<f32>;
using c128 = std::complex<f64>;

using size_t = std::size_t;
using i8 = int8_t;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;
using uint = unsigned int;
using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

template <class T>
using svec = std::vector<T>;
template <class T1, class T2>
using umap = std::unordered_map<T1, T2>;
template <class T1>
using uset = std::unordered_set<T1>;
template <class C>
using uptr = std::unique_ptr<C>;

template <typename T>
inline constexpr T qNaN = std::numeric_limits<T>::quiet_NaN();
template <typename T>
inline constexpr T inf = std::numeric_limits<T>::infinity();
