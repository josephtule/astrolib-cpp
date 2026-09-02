// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/ephemeris.hpp"
#include "core/state.hpp"

struct EphemerisProvider {
    CartesianEphemerisTable table{};
    CartesianSampleOptions options{};
};

struct OrientationProvider {
    OrientationEphemerisTable table{};
    OrientationSampleOptions options{};
};

struct BodyEphemerisProviders {
    sptr<const EphemerisProvider> translation;
    sptr<const OrientationProvider> orientation;
};

StatusCode validate_ephemeris_provider(const EphemerisProvider& provider);

StatusCode validate_orientation_provider(const OrientationProvider& provider);

StatusCode query_ephemeris_provider(
    const EphemerisProvider& provider,
    const JulianDate& query_epoch,
    TimeScale query_scale,
    const TimeOffsets& offsets,
    StateTr& out
);

StatusCode query_orientation_provider(
    const OrientationProvider& provider,
    const JulianDate& query_epoch,
    TimeScale query_scale,
    const TimeOffsets& offsets,
    StateAtt& out
);
