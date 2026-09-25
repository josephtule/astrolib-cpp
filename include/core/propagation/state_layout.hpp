#pragma once

#include "core/entity.hpp"
#include "core/status.hpp"
#include "util/typedefs.hpp"

enum struct PropagationStateDomain {
    integrated, // numerically integrated
    prescribed, // provider or simple model
    algebraic,  // reconstructed from state/geometry
    inactive,
};
inline string propagation_state_domain_string(PropagationStateDomain domain) {
    switch (domain) {
    case PropagationStateDomain::integrated: return "integrated";
    case PropagationStateDomain::prescribed: return "prescribed";
    case PropagationStateDomain::algebraic: return "algebraic";
    case PropagationStateDomain::inactive: return "inactive";
    }
    return "unknown";
}

enum struct PropagationStateKind {
    translation,
    attitude,
    mass,
    resource,
    parameter,
};
inline string propagation_state_kind_string(PropagationStateKind kind) {
    switch (kind) {
    case PropagationStateKind::translation: return "translation";
    case PropagationStateKind::attitude: return "attitude";
    case PropagationStateKind::mass: return "mass";
    case PropagationStateKind::resource: return "resource";
    case PropagationStateKind::parameter: return "parameter";
    }
    return "unknown";
}

struct PropagationStateBlock {
    PropagationStateKind kind = PropagationStateKind::translation;
    string key = "primary";

    PropagationStateDomain domain = PropagationStateDomain::inactive;
    i32 offset = -1; // index in contiguous propagation state
    i32 size = 0;    // translation: 6, attitude: 7, mass: 1, resource/parameter: variable
};

struct PropagationBodyLayout {
    EntityId entity_id = kInvalidEntityId;
    i32 body_index = -1;

    svec<PropagationStateBlock> blocks;
};

struct PropagationStateLayout {
    // layout for dense/contiguous vector of propagation states
    svec<PropagationBodyLayout> bodies;
    i32 state_size = 0;
};

StatusCode validate_propagation_state_layout(const PropagationStateLayout& layout);

const PropagationStateBlock* find_state_block(
    const PropagationBodyLayout& body,
    PropagationStateKind kind,
    const string& key
);

svec<const PropagationStateBlock*> find_state_blocks(
    const PropagationBodyLayout& body,
    PropagationStateKind kind
);