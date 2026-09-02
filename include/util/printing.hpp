// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/state.hpp"
#include "util/vecdefs.hpp"
#include <print>
#include <sstream>
#include <string>

inline eig::IOFormat eigen_vec_format(
    int precision = eig::FullPrecision,
    int flags = eig::DontAlignCols,
    const string& coeff_sep = ", ",
    const string& row_sep = "; ",
    const string& row_prefix = "",
    const string& row_suffix = "",
    const string& mat_prefix = "[",
    const string& mat_suffix = "]"
) {
    return eig::IOFormat(
        precision,
        flags,
        coeff_sep,
        row_sep,
        row_prefix,
        row_suffix,
        mat_prefix,
        mat_suffix
    );
}

template <class Derived>
inline string vec_string(
    const eig::MatrixBase<Derived>& x,
    const eig::IOFormat& fmt = eigen_vec_format()
) {
    std::ostringstream ss;
    ss << x.format(fmt);
    return ss.str();
}

template <class T>
inline string vec_string(
    const svec<T>& x,
    const eig::IOFormat& fmt = eigen_vec_format()
) {
    if (!std::is_arithmetic<T>::value) return "Type is not numeric";
    if (x.empty()) return "[]";

    eig::Map<const eig::Matrix<T, eig::Dynamic, 1>> x_eig(
        x.data(),
        static_cast<eig::Index>(x.size())
    );

    return vec_string(x_eig, fmt);
}

template <class T, std::size_t N>
inline string vec_string(
    const array<T, N>& x,
    const eig::IOFormat& fmt = eigen_vec_format()
) {
    if (!std::is_arithmetic<T>::value) return "Type is not numeric";
    eig::Map<const eig::Matrix<T, eig::Dynamic, 1>> x_eig(
        x.data(),
        static_cast<eig::Index>(N)
    );

    return vec_string(x_eig, fmt);
}

inline string state_string(
    const StateTr& x,
    const eig::IOFormat& fmt = eigen_vec_format()
) {
    std::ostringstream ss;
    ss << "r = " << vec_string(x.r, fmt) << ", v = " << vec_string(x.v, fmt);
    return ss.str();
}

inline string state_string(
    const StateAtt& x,
    const eig::IOFormat& fmt = eigen_vec_format()
) {
    std::ostringstream ss;
    ss << "q = " << vec_string(x.q, fmt) << ", w = " << vec_string(x.w, fmt);
    return ss.str();
}

inline string deriv_string(
    const DerivTr& dx,
    const eig::IOFormat& fmt = eigen_vec_format()
) {
    std::ostringstream ss;
    ss << "dr = " << vec_string(dx.dr, fmt) << ", dv = " << vec_string(dx.dv, fmt);
    return ss.str();
}

inline string deriv_string(
    const DerivAtt& dx,
    const eig::IOFormat& fmt = eigen_vec_format()
) {
    std::ostringstream ss;
    ss << "dq = " << vec_string(dx.dq, fmt) << ", dw = " << vec_string(dx.dw, fmt);
    return ss.str();
}

inline void print_state_tr(const StateTr& x_tr, const string& indent = "") {
    std::println(
        "{}State (tr): r = {}, v = {}",
        indent,
        vec_string(x_tr.r),
        vec_string(x_tr.v)
    );
}

inline void print_state_att(const StateAtt& x_att, const string& indent = "") {
    std::println(
        "{}State (att): q = {}, w = {}",
        indent,
        vec_string(x_att.q),
        vec_string(x_att.w)
    );
}
