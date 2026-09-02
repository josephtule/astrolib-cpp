// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/measurement_realism.hpp"
#include "core/measurement.hpp"
#include <cmath>

StatusCode evaluate_measurement_realism(
    const MeasurementRealismQuery& query,
    const MeasurementRealismPolicy& policy,
    MeasurementRealismResult& out
) {
    if (!std::isfinite(query.t) || measurement_dim(query.type) <= 0
        || query.observer_id == kInvalidEntityId || query.target_id == kInvalidEntityId
        || query.observer_id == query.target_id) {
        return StatusCode::invalid_input;
    }
    // TODO: implement each requested check; disabled does not mean verified visible
    if (policy.check_station_horizon || policy.check_body_occultation
        || policy.check_target_illumination || policy.check_field_of_view) {
        return StatusCode::unsupported_type;
    }
    out = MeasurementRealismResult{};
    return StatusCode::ok;
}
