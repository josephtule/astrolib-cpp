// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/ephemeris.hpp"
#include "util/typedefs.hpp"

struct EphemerisCSVLayout {
    string filepath;
    string delimiter = ",";
    bool has_header = true;
    i32 header_lines = 1;
};

struct EphemerisWriteOptions {
    bool overwrite = false;
    i32 precision = 12;
    EphemerisCSVLayout csv{};
};

StatusCode load_cartesian_ephemeris(
    const string& manifest_filepath,
    CartesianEphemerisTable& out
);

StatusCode load_orientation_ephemeris(
    const string& manifest_filepath,
    OrientationEphemerisTable& out,
    f64 tol = tol12
);

StatusCode save_cartesian_ephemeris(
    const string& manifest_filepath,
    const CartesianEphemerisTable& in,
    const EphemerisWriteOptions& opts
);

StatusCode save_orientation_ephemeris(
    const string& manifest_filepath,
    const OrientationEphemerisTable& in,
    const EphemerisWriteOptions& opts
);
