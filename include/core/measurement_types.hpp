// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/entity.hpp"
#include "core/observation_type.hpp"
#include "core/state.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include <random>

struct Measurement {
    f64 t = 0.0;
    ObservationType type = ObservationType::radec;
    EntityId observer_id = kInvalidEntityId;
    EntityId target_id = kInvalidEntityId;
    vecXd z; // measured value
    matXd R; // measurement covariance
};

enum struct AzelInputFrame : i32 { enu, bcbf };

struct MeasurementContext {
    StateTr x_tr_observer;
    StateTr x_tr_target;

    // local station geometry, used only for azel
    bool has_station_local = false;
    AzelInputFrame azel_frame = AzelInputFrame::enu;
    vec3d station_llh = vec3d0; // [lat, lon, h], angle units from angle_in
};

struct MeasurementNoiseOptions {
    bool diagonal = true; // false = use cholesky
    std::mt19937_64& rng;
    bool enabled = false;
    UAngle u_angle = UAngle::radian;
};

MeasurementContext make_measurement_context(
    const StateTr& x_tr_target,
    const StateTr& x_tr_observer
);

i32 measurement_dim(ObservationType type);
