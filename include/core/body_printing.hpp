// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

struct Body;
struct Celestial;
struct InstrumentSuite;
struct MassProperties;
struct PlatformInstrument;
struct Satellite;
struct Station;

void print_instruments(const InstrumentSuite& suite);
void print_mass_properties(const MassProperties& mp, const std::string& indent = "");
void print_body_common(const Body& body, const std::string& indent = "");
void print_celestial(const Celestial& cel, const std::string& indent = "");
void print_satellite(const Satellite& sat, const std::string& indent = "");
void print_station_instrument(
    const PlatformInstrument& instrument,
    const std::string& indent = ""
);
void print_station(const Station& stat, const std::string& indent = "");
