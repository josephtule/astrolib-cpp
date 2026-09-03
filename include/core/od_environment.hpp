// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/entity.hpp"
#include "core/od_dynamics.hpp"
#include "util/typedefs.hpp"
#include <string>

enum struct ODPropagationMode : i32 {
    lightweight, // no staging
    provider, // staged providers
    full_world, // full staging
};

inline std::string od_propagation_mode_str(ODPropagationMode mode) {
    switch (mode) {
    case ODPropagationMode::lightweight: return "lightweight";
    case ODPropagationMode::provider: return "provider";
    case ODPropagationMode::full_world: return "full_world";
    }
    return "unknown";
}

// configuration for the new context-based propagation path, not the legacy helpers
// context boundary: target/source r and v are absolute in simulation inertial axes
// and simulation units (km, km/s); covariance/STM state order is [r, v]
// force input: r_target_source_I = r_target_I - r_source_I (same epoch)
// source attitude q maps simulation inertial -> source body-fixed
// full_world initially estimates target translation only; sources must not depend
// on the target state. Coupled source sensitivities require an augmented STM later.
struct ODPropagationConfig {
    ODPropagationMode mode = ODPropagationMode::lightweight;
    EntityId target_id = kInvalidEntityId;
    svec<EntityId> source_ids; // force sources, excludes target; runtime world IDs
    f64 t0 = 0.0;              // reference epoch in simulation seconds, not JD
};

struct ODSourceConfig {
    EntityId id = kInvalidEntityId;
    string object;
    string body_frame;
    StateTr x_tr0;
    StateAtt x_att0;
    ODAnchorTrModel tr_model = ODAnchorTrModel::fixed;
    ODAnchorAttModel att_model = ODAnchorAttModel::fixed;
    BodyEphemerisProviders providers;
    GravityModel gravity_model = GravityModel::pointmass;
    f64 mu = 0.0;
    f64 ref_radius = 0.0;
    i32 degree = 0;
    vec7d J = vec7d0;
};

struct ODDynamicsContext {
    ODPropagationConfig propagation;
    svec<ODSourceConfig> sources;
    string inertial_frame;
    string center; // origin label for absolute simulation-inertial states
    bool has_epoch = false;
    JulianDate epoch; // date corresponding to propagation.t0
    TimeScale time_scale = TimeScale::tt;
    TimeOffsets offsets;
};

struct ODSourceState {
    StateTr x_tr;
    StateAtt x_att;
};

StatusCode validate_od_dynamics_context(const ODDynamicsContext& ctx);
StatusCode od_source_state_at_time(
    const ODDynamicsContext& ctx,
    EntityId source_id,
    f64 t,
    ODSourceState& out
);

StatusCode derivtr_od_context(
    const ODDynamicsContext& ctx, f64 t, const StateTr& x_target_I, DerivTr& out
);

// stationary body-fixed anchor; same source query as force evaluation
StatusCode od_anchored_observer_state_at_time(
    const ODDynamicsContext& ctx, EntityId anchor_id, const vec3d& r_observer_body_B,
    f64 t, StateTr& out
);
