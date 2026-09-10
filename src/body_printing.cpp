// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/body_printing.hpp"

#include "core/body.hpp"
#include "core/instrument.hpp"
#include "core/observation_type.hpp"
#include "util/printing.hpp"

#include <print>

void print_instruments(const InstrumentSuite& suite) {
    for (const auto& [id, instrument] : suite.instruments) {
        std::string enabled_str = instrument.enabled ? "enabled" : "disabled";
        std::println("{} ({})", instrument.name, enabled_str);
    }
}

void print_mass_properties(const MassProperties& mp, const std::string& indent) {
    std::println("{}Mass: {}", indent, mp.mass);
    std::println("{}Inertia Tensor: I = {}", indent, vec_string(mp.I));
    std::println("{}Inverse Inertia Tensor: I_inv = {}", indent, vec_string(mp.I_inv));
    std::println("{}Principle Axes: {}", indent, mp.principal_axes);
    std::println("{}Active: {}", indent, mp.active);
}

void print_body_common(const Body& body, const std::string& indent) {
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

void print_celestial(const Celestial& cel, const std::string& indent) {
    const std::string indent2 = indent + "    ";

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
    std::println(
        "{}Radiation Model: {}",
        indent2,
        radiation_model_str(cel.radiation_model)
    );
    std::println("{}Reference Radius: {}", indent2, cel.ref_radius);
    std::println("{}Mean Radius: {}", indent2, cel.mean_radius);
    std::println("{}Semimajor Axis: {}", indent2, cel.semimajor_axis);
    std::println("{}Semiminor Axis: {}", indent2, cel.semiminor_axis);
    std::println("{}Eccentricity: {}", indent2, cel.eccentricity);
    std::println("{}Flattening: {}", indent2, cel.flattening);
    std::println("{}Spin Rate: {}", indent2, cel.spin_rate());
}

void print_satellite(const Satellite& sat, const std::string& indent) {
    const std::string indent2 = indent + "    ";

    std::println("{}--- Satellite ID: {}", indent, sat.id);
    print_body_common(sat, indent2);
    print_mass_properties(sat.mass_properties, indent2);
}

void print_station_instrument(
    const PlatformInstrument& instrument,
    const std::string& indent
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

void print_station(const Station& stat, const std::string& indent) {
    const std::string indent2 = indent + "    ";
    const std::string indent3 = indent2 + "    ";

    std::println("{}--- Station ID: {}", indent, stat.id);
    print_body_common(stat, indent2);
    std::println("{}Anchored: {}", indent2, stat.anchored);
    std::println("{}Anchor ID: {}", indent2, stat.anchor_id);
    std::println("{}r_body_BCBF: {}", indent2, vec_string(stat.r_body_BCBF));
    std::println("{}llh_BCBF: {}", indent2, vec_string(stat.llh_BCBF));

    if (!stat.anchored) print_mass_properties(stat.mass_properties, indent2);

    std::println(
        "{}Instruments Count: {}",
        indent2,
        stat.instrument_suite.instruments.size()
    );
    std::println(
        "{}Enabled Instrument IDs: {}",
        indent2,
        vec_string(stat.instrument_suite.enabled_ids)
    );
    for (const auto& [id, instrument] : stat.instrument_suite.instruments) {
        print_station_instrument(instrument, indent3);
    }
}
