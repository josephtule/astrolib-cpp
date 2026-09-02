// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/ephemeris_provider.hpp"
#include "core/ephemeris.hpp"
#include "core/interpolation.hpp"
#include "core/state.hpp"
#include "core/status.hpp"
#include "core/time.hpp"
#include "util/math.hpp"
#include <cmath>

static StatusCode validate_sampling_options(const CartesianSampleOptions& opts) {
    if (!finite_pos(opts.tol)) return StatusCode::invalid_input;
    return StatusCode::ok;
}

static StatusCode validate_sampling_options(const OrientationSampleOptions& opts) {
    if (!finite_pos(opts.tol)) return StatusCode::invalid_input;
    return StatusCode::ok;
}

StatusCode validate_ephemeris_provider(const EphemerisProvider& provider) {
    StatusCode status;

    status = validate_sampling_options(provider.options);
    if (status != StatusCode::ok) return status;
    if (!provider.table.has_velocity) {
        if (provider.options.interpolation
            == CartesianInterpolationMethod::cubic_hermite) {
            return StatusCode::invalid_input;
        }
        if (provider.options.extrapolation == ExtrapolationMethod::constant_velocity) {
            return StatusCode::invalid_input;
        }
    }

    status = validate_cartesian_ephemeris_table(provider.table);
    if (status != StatusCode::ok) return status;

    return StatusCode::ok;
}

StatusCode validate_orientation_provider(const OrientationProvider& provider) {
    StatusCode status;

    status = validate_sampling_options(provider.options);
    if (status != StatusCode::ok) return status;
    if (!provider.table.has_angular_velocity
        && provider.options.extrapolation == ExtrapolationMethod::constant_velocity)
        return StatusCode::invalid_input;

    status = validate_orientation_ephemeris_table(provider.table);
    if (status != StatusCode::ok) return status;

    return StatusCode::ok;
}

static StatusCode ephemeris_query_dt(
    const EphemerisEpochMetadata& metadata,
    const JulianDate& query_epoch,
    TimeScale query_scale,
    const TimeOffsets& offsets,
    f64& dt
) {
    if (!std::isfinite(query_epoch.day) || !std::isfinite(query_epoch.frac))
        return StatusCode::invalid_input;

    JulianDate epoch = query_epoch;
    if (epoch.frac >= 1.0 || epoch.frac < 0.0) {
        epoch = normalize_jd(epoch);
    }

    epoch = jd_scale_convert(epoch, query_scale, metadata.time_scale, offsets);
    if (epoch.frac >= 1.0 || epoch.frac < 0.0) {
        epoch = normalize_jd(epoch);
    }

    const JulianDate& ref_epoch = metadata.ref_epoch;
    f64 d_day = epoch.day - ref_epoch.day;
    f64 d_frac = epoch.frac - ref_epoch.frac;
    f64 dt_temp = (d_day + d_frac) * 86400.0;
    if (!std::isfinite(dt_temp)) return StatusCode::non_finite_result;
    dt = dt_temp;

    return StatusCode::ok;
}

StatusCode query_ephemeris_provider(
    const EphemerisProvider& provider,
    const JulianDate& query_epoch,
    TimeScale query_scale,
    const TimeOffsets& offsets,
    StateTr& out
) {
    StatusCode status;

    f64 dt;
    status = ephemeris_query_dt(
        provider.table.metadata.epoch,
        query_epoch,
        query_scale,
        offsets,
        dt
    );
    if (status != StatusCode::ok) return status;

    StateTr temp{};
    status = sample_cartesian_ephemeris(provider.table, dt, temp, provider.options);
    if (status != StatusCode::ok) return status;

    out = temp;
    return StatusCode::ok;
}

StatusCode query_orientation_provider(
    const OrientationProvider& provider,
    const JulianDate& query_epoch,
    TimeScale query_scale,
    const TimeOffsets& offsets,
    StateAtt& out
) {
    f64 dt;
    StatusCode status = ephemeris_query_dt(
        provider.table.metadata.epoch,
        query_epoch,
        query_scale,
        offsets,
        dt
    );
    if (status != StatusCode::ok) return status;

    StateAtt temp{};
    status = sample_orientation_ephemeris(provider.table, dt, temp, provider.options);
    if (status != StatusCode::ok) return status;

    out = temp;
    return StatusCode::ok;
}
