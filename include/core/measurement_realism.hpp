// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/entity.hpp"
#include "core/observation_type.hpp"
#include "core/status.hpp"

struct MeasurementRealismPolicy {
    bool check_station_horizon = false;
    bool check_body_occultation = false;
    bool check_target_illumination = false;
    bool check_field_of_view = false;
};

struct MeasurementRealismQuery {
    f64 t = 0.0;
    ObservationType type = ObservationType::radec;
    EntityId observer_id = kInvalidEntityId;
    EntityId target_id = kInvalidEntityId;
};

enum struct MeasurementAvailability { not_evaluated, available, unavailable };

struct MeasurementRealismResult {
    MeasurementAvailability availability = MeasurementAvailability::not_evaluated;
};

StatusCode evaluate_measurement_realism(
    const MeasurementRealismQuery& query,
    const MeasurementRealismPolicy& policy,
    MeasurementRealismResult& out
);
