// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/measurement_types.hpp"
#include "core/status.hpp"

StatusCode apply_measurement_noise_diagonal(
    Measurement& meas,
    const MeasurementNoiseOptions& opts,
    ObservationType type
);

StatusCode apply_measurement_noise_cholesky(
    Measurement& meas,
    const MeasurementNoiseOptions& opts,
    ObservationType type
);
