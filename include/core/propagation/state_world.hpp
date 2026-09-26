#pragma once

#include "core/propagation/state.hpp"
#include "core/status.hpp"

class World;

StatusCode pack_world_propagation_state(
    const World& world,
    const PropagationStateLayout& layout,
    PropagationState& out
);

StatusCode commit_world_propagation_state(
    World& world,
    const PropagationStateLayout& layout,
    const PropagationState& state
);

StatusCode finalize_world_propagation_state_layout(
    const World& world,
    PropagationStateLayout& layout
);

StatusCode validate_world_propagation_state_layout(
    const World& world,
    const PropagationStateLayout& layout
);
