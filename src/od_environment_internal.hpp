// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "core/od_environment.hpp"

// caller validates the context once and keeps it unchanged for the propagation call
namespace od_environment_detail {
StatusCode source_state_at_time(const ODDynamicsContext&, EntityId, f64, ODSourceState&);
StatusCode derivative(const ODDynamicsContext&, f64, const StateTr&, DerivTr&);
}
