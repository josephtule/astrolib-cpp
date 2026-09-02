// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/ephemeris.hpp"
#include "core/state.hpp"

enum struct ProviderCoverageAction {
    reject_step,
    extrapolate,
    hold_state,
    handoff_to_dynamics,
    stop_world
};

struct ProviderCoveragePolicy {
    // when query out of provider bounds
    ProviderCoverageAction before_start = ProviderCoverageAction::reject_step;
    ProviderCoverageAction after_end = ProviderCoverageAction::reject_step;
};

struct EphemerisProvider {
    CartesianEphemerisTable table{};
    CartesianSampleOptions options{};
    ProviderCoveragePolicy coverage{};
};

struct OrientationProvider {
    OrientationEphemerisTable table{};
    OrientationSampleOptions options{};
    ProviderCoveragePolicy coverage{};
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

StatusCode ephemeris_query_dt(
    const EphemerisEpochMetadata& metadata,
    const JulianDate& query_epoch,
    TimeScale query_scale,
    const TimeOffsets& offsets,
    f64& dt
);
