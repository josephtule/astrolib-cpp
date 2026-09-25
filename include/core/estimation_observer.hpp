// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/entity.hpp"
#include "core/state.hpp"
#include "core/status.hpp"
#include "core/time.hpp"
#include "core/world.hpp"
#include "core/world_history.hpp"
#include "util/constants.hpp"
#include "util/typedefs.hpp"
#include "util/vecdefs.hpp"

enum struct ObserverStateSource { world, provider, estimate };

string observer_state_source_string(ObserverStateSource source);

enum struct ODObserverMissingEpochPolicy { reject, propagate };

string od_observer_missing_epoch_policy_string(ODObserverMissingEpochPolicy policy);

struct ODObserverBinding {
    EntityId observer_id = kInvalidEntityId;
    ObserverStateSource source = ObserverStateSource::world; // runtime state backend
    string source_id; // provider/estimate registry key; empty for world
    ODObserverMissingEpochPolicy missing_epoch = ODObserverMissingEpochPolicy::reject;
    bool known_state = false;
};

StatusCode validate_od_observer_binding(const ODObserverBinding& binding);

StatusCode find_od_observer_binding(
    const svec<ODObserverBinding>& bindings,
    EntityId observer_id,
    const ODObserverBinding*& out
);

StatusCode validate_od_observer_bindings(const svec<ODObserverBinding>& bindings);

struct ODObserverSample {
    EntityId observer_id = kInvalidEntityId;
    f64 t = 0.0;

    StateTr x;
    mat6d P = mat6d0;

    string center;     // origin the observer state is relative to
    string frame_name; // coordinate axes used to express the observer state

    ObserverStateSource source = ObserverStateSource::world; // backend used for this sample
    string source_id; // backend registry key or empty for world
    bool known_state = false;
};

struct ODObserverEpochContext {
    f64 t_ref = 0.0;
    JulianDate jd_ref{};
    TimeScale time_scale = TimeScale::utc;
    TimeOffsets offsets{};
    bool has_calendar_epoch = false;
};

// forward declare
struct EphemerisProvider;
struct ODEstimateHistory;

struct ODObserverQueryContext {
    const World* world = nullptr;
    const WorldHistoryProvider* world_history = nullptr;

    const umap<string, const EphemerisProvider*>* providers = nullptr;
    const umap<string, const ODEstimateHistory*>* estimates = nullptr;

    ODObserverEpochContext epoch;

    string center;     // required origin for the returned observer state
    string frame_name; // required axes for the returned observer state
    f64 tol_time = tol12;
};

StatusCode validate_od_observer_sample(
    const ODObserverSample& sample,
    f64 tol_time = tol12
);

StatusCode query_od_observer(
    const ODObserverBinding& binding,
    const ODObserverQueryContext& ctx,
    f64 t,
    ODObserverSample& out
);
