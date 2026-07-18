// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/body.hpp"
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

inline string celestial_attitude_model_str(const CelestialAttitudeModel model) {
    switch (model) {
    case CelestialAttitudeModel::fixed: return "fixed";
    case CelestialAttitudeModel::simple_spin: return "simple spin";
    case CelestialAttitudeModel::provider: return "provider";
    }

    return "unknown";
}

inline string radiation_model_str(const RadiationModel model) {
    switch (model) {
    case RadiationModel::none: return "none";
    case RadiationModel::isotropic: return "isotropic";
    }

    return "unknown";
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

inline void print_mass_properties(
    const MassProperties& mp,
    const string& indent = ""
) {
    std::println("{}Mass: {}", indent, mp.mass);
    std::println("{}Inertia Tensor: I = {}", indent, vec_string(mp.I));
    std::println("{}Inverse Inertia Tensor: I_inv = {}", indent, vec_string(mp.I_inv));
    std::println("{}Principle Axes: {}", indent, mp.principal_axes);
    std::println("{}Active: {}", indent, mp.active);
}

inline void print_body_common(const Body& body, const string& indent = "") {
    std::println("{}ID: {}", indent, body.id);
    std::println("{}Name: {}", indent, body.name);
    std::println("{}Body Type: {}", indent, body_type_str(body.body_type));
    std::println("{}Propagate Translation: {}", indent, body.propagate_tr);
    std::println("{}Propagate Attitude: {}", indent, body.propagate_att);
    std::println("{}Emits Gravity: {}", indent, body.emits_gravity);
    std::println("{}Emits Radiation: {}", indent, body.emits_radiation);
    std::println("{}Has Atmosphere: {}", indent, body.has_atmosphere);
    print_state_tr(body.x_tr, indent);
    print_state_att(body.x_att, indent);
}

inline void print_celestial(const Celestial& cel, const string& indent = "") {
    const string indent2 = indent + "    ";

    std::println("{}--- Celestial ID: {}", indent, cel.id);
    print_body_common(cel, indent2);
    std::println("{}Gravity Model: {}", indent2, gravity_model_str(cel.gravity_model));
    std::println("{}mu: {}", indent2, cel.mu);
    std::println("{}Degree: {}", indent2, cel.degree);
    std::println("{}Order: {}", indent2, cel.order);
    std::println("{}J: {}", indent2, vec_string(cel.J));
    std::println("{}C Dim: ({}, {})", indent2, cel.C.rows(), cel.C.cols());
    std::println("{}S Dim: ({}, {})", indent2, cel.S.rows(), cel.S.cols());
    std::println(
        "{}Attitude Model: {}",
        indent2,
        celestial_attitude_model_str(cel.attitude_model)
    );
    std::println("{}Radiation Model: {}", indent2, radiation_model_str(cel.radiation_model));
    std::println("{}Reference Radius: {}", indent2, cel.ref_radius);
    std::println("{}Mean Radius: {}", indent2, cel.mean_radius);
    std::println("{}Semimajor Axis: {}", indent2, cel.semimajor_axis);
    std::println("{}Semiminor Axis: {}", indent2, cel.semiminor_axis);
    std::println("{}Eccentricity: {}", indent2, cel.eccentricity);
    std::println("{}Flattening: {}", indent2, cel.flattening);
    std::println("{}Spin Rate: {}", indent2, cel.spin_rate());
}

inline void print_satellite(const Satellite& sat, const string& indent = "") {
    const string indent2 = indent + "    ";

    std::println("{}--- Satellite ID: {}", indent, sat.id);
    print_body_common(sat, indent2);
    print_mass_properties(sat.mass_properties, indent2);
}

inline void print_station_instrument(
    const StationInstrument& instrument,
    const string& indent = ""
) {
    std::println("{}--- Instrument ID: {}", indent, instrument.id);
    std::println("{}Name: {}", indent, instrument.name);
    std::println("{}Instrument Type: {}", indent, observation_type_str(instrument.type));
    std::println("{}Enabled: {}", indent, instrument.enabled);
    std::println(
        "{}Measurement Covariance Dim: ({}, {})",
        indent,
        instrument.R.rows(),
        instrument.R.cols()
    );
    std::println("{}Measurement Covariance: R = {}", indent, vec_string(instrument.R));
}

inline void print_station(const Station& stat, const string& indent = "") {
    const string indent2 = indent + "    ";
    const string indent3 = indent2 + "    ";

    std::println("{}--- Station ID: {}", indent, stat.id);
    print_body_common(stat, indent2);
    std::println("{}Anchored: {}", indent2, stat.anchored);
    std::println("{}Anchor ID: {}", indent2, stat.anchor_id);
    std::println("{}r_body_BCBF: {}", indent2, vec_string(stat.r_body_BCBF));
    std::println("{}llh_BCBF: {}", indent2, vec_string(stat.llh_BCBF));

    if (!stat.anchored) {
        print_mass_properties(stat.mass_properties, indent2);
    }

    std::println("{}Instruments Count: {}", indent2, stat.instruments.size());
    std::println(
        "{}Enabled Instrument IDs: {}",
        indent2,
        vec_string(stat.enabled_instrument_ids)
    );
    for (const auto& [id, instrument] : stat.instruments) {
        print_station_instrument(instrument, indent3);
    }
}
