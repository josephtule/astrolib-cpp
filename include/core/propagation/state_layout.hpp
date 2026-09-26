#pragma once

#include "core/entity.hpp"
#include "core/status.hpp"
#include "core/world/world_revision.hpp"
#include "util/typedefs.hpp"

using PropagationBodyIndex = i32;
inline constexpr PropagationBodyIndex kInvalidPropagationBodyIndex = -1;

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

inline constexpr string default_block_key = "primary";
struct PropagationStateBlock {
    PropagationStateKind kind = PropagationStateKind::translation;
    string key = default_block_key;

    PropagationStateDomain domain = PropagationStateDomain::inactive;
    i32 offset = -1; // index in contiguous propagation state
    i32 size = 0;    // translation: 6, attitude: 7, mass: 1, resource/parameter: variable
};

struct PropagationBodyLayout {
    EntityId entity_id = kInvalidEntityId;
    PropagationBodyIndex body_index = kInvalidPropagationBodyIndex;

    svec<PropagationStateBlock> blocks;
};

struct PropagationLayoutRevisions {
    WorldRevision topology = 0;
    WorldRevision dynamics = 0;
    WorldRevision providers = 0;
    bool bound = false;
};

struct PropagationStateLayout {
    // layout for dense/contiguous vector of propagation states
    svec<PropagationBodyLayout> bodies;
    umap<EntityId, PropagationBodyIndex> body_indices;

    PropagationLayoutRevisions revisions{};

    i32 state_size = 0;
};

StatusCode initialize_propagation_state_layout(
    const svec<EntityId>& entity_ids,
    PropagationStateLayout& out
);

StatusCode validate_propagation_state_block(const PropagationStateBlock& block);
StatusCode validate_propagation_body_layout(const PropagationBodyLayout& body);
StatusCode validate_propagation_state_layout(const PropagationStateLayout& layout);

i32 count_layout_domain(
    const PropagationStateLayout& layout,
    PropagationStateDomain domain
);
i32 count_layout_domain_size(
    const PropagationStateLayout& layout,
    PropagationStateDomain domain
);
i32 count_layout_integrated(const PropagationStateLayout& layout);
i32 count_layout_nonintegrated(const PropagationStateLayout& layout);
i32 count_layout_integrated_size(const PropagationStateLayout& layout);
i32 count_layout_nonintegrated_size(const PropagationStateLayout& layout);

PropagationBodyIndex find_propagation_body_index(
    const PropagationStateLayout& layout,
    EntityId entity_id
);

const PropagationBodyLayout* find_propagation_body_layout(
    const PropagationStateLayout& layout,
    EntityId entity_id
);

PropagationBodyLayout* find_propagation_body_layout(
    PropagationStateLayout& layout,
    EntityId entity_id
);

StatusCode add_propagation_state_block(
    PropagationStateLayout& layout,
    EntityId entity_id,
    PropagationStateKind kind,
    const string& key,
    PropagationStateDomain domain,
    i32 size
);

const PropagationStateBlock* find_state_block(
    const PropagationBodyLayout& body,
    PropagationStateKind kind,
    const string& key
);

svec<const PropagationStateBlock*> find_state_blocks(
    const PropagationBodyLayout& body,
    PropagationStateKind kind
);
