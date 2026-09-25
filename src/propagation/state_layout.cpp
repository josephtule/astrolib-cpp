#include "core/propagation/state_layout.hpp"
#include "core/entity.hpp"
#include "core/status.hpp"
#include "util/typedefs.hpp"

#include <algorithm>
#include <limits>
#include <set>

StatusCode validate_propagation_state_block(const PropagationStateBlock& block) {
    if (block.key.empty()) return StatusCode::invalid_input;

    if (block.domain == PropagationStateDomain::integrated) {
        if (block.size <= 0) return StatusCode::invalid_input;
        if (block.offset < 0) return StatusCode::invalid_input;
        if (block.offset > std::numeric_limits<i32>::max() - block.size) {
            return StatusCode::invalid_input;
        }
    } else {
        if (block.size != 0) return StatusCode::invalid_input;
        if (block.offset != -1) return StatusCode::invalid_input;
    }

    return StatusCode::ok;
}

StatusCode validate_propagation_body_layout(const PropagationBodyLayout& body) {
    if (body.entity_id == kInvalidEntityId) return StatusCode::body_not_found;
    if (body.body_index < 0) return StatusCode::invalid_input;

    StatusCode status;

    std::set<std::pair<PropagationStateKind, string>> seen_blocks;

    for (const auto& block : body.blocks) {
        if (!seen_blocks.insert({block.kind, block.key}).second) {
            return StatusCode::duplicate_id;
        }

        status = validate_propagation_state_block(block);
        if (status != StatusCode::ok) return status;
    }

    return StatusCode::ok;
}

StatusCode initialize_propagation_state_layout(
    const svec<EntityId>& entity_ids,
    PropagationStateLayout& out
) {
    if (entity_ids.empty()) return StatusCode::invalid_input;
    PropagationStateLayout temp;

    svec<EntityId> sorted_ids;
    for (i32 i = 0; i < entity_ids.size(); ++i) {
        EntityId id = entity_ids[i];
        if (id == kInvalidEntityId) return StatusCode::invalid_input;
        for (i32 j = 0; j < sorted_ids.size(); ++j) {
            if (id == sorted_ids[j]) return StatusCode::invalid_input;
        }
        sorted_ids.push_back(id);
    }
    std::sort(sorted_ids.begin(), sorted_ids.end());

    for (i32 i = 0; i < sorted_ids.size(); ++i) {
        const auto body_index = static_cast<PropagationBodyIndex>(i);

        temp.bodies.push_back({
            .entity_id = sorted_ids[i],
            .body_index = body_index,
        });
        temp.body_indices.emplace(sorted_ids[i], body_index);
    }

    StatusCode status = validate_propagation_state_layout(temp);
    if (status != StatusCode::ok) return status;

    out = std::move(temp);
    return StatusCode::ok;
}

StatusCode validate_propagation_state_layout(const PropagationStateLayout& layout) {
    StatusCode status;

    if (layout.state_size < 0) return StatusCode::invalid_input;

    if (layout.body_indices.size() != layout.bodies.size())
        return StatusCode::size_mismatch;

    std::set<EntityId> seen_entity_ids;
    std::set<i32> seen_body_indices;
    svec<std::pair<i32, i32>> ranges;
    for (size_t i = 0; i < layout.bodies.size(); ++i) {
        const PropagationBodyLayout& body = layout.bodies[i];

        if (body.body_index != static_cast<PropagationBodyIndex>(i))
            return StatusCode::invalid_input;
        if (find_propagation_body_index(layout, body.entity_id) != body.body_index)
            return StatusCode::invalid_input;

        if (!seen_entity_ids.insert(body.entity_id).second)
            return StatusCode::invalid_input;
        if (!seen_body_indices.insert(body.body_index).second)
            return StatusCode::invalid_input;
        if (body.body_index >= layout.bodies.size()) return StatusCode::invalid_input;

        status = validate_propagation_body_layout(body);
        if (status != StatusCode::ok) return status;

        for (const auto& block : body.blocks) {
            if (block.domain != PropagationStateDomain::integrated) continue;

            i32 end = block.offset + block.size;
            ranges.push_back({block.offset, end});
        }
    }

    for (const auto& [id, idx] : layout.body_indices) {
        if (id == kInvalidEntityId || idx == kInvalidPropagationBodyIndex)
            return StatusCode::invalid_input;
        if (idx >= layout.bodies.size()) return StatusCode::invalid_input;
        if (layout.bodies[idx].entity_id != id) return StatusCode::invalid_input;
    }

    std::sort(ranges.begin(), ranges.end());
    i32 expected_offset = 0;
    for (const auto& [begin, end] : ranges) {
        if (begin != expected_offset) { // overlap or gap
            return StatusCode::invalid_input;
        }

        expected_offset = end;
    }

    if (expected_offset != layout.state_size) return StatusCode::invalid_input;

    return StatusCode::ok;
}

i32 count_layout_domain(
    const PropagationStateLayout& layout,
    PropagationStateDomain domain
) {
    i32 count = 0;
    for (const auto& body : layout.bodies) {
        for (const auto& block : body.blocks) {
            if (block.domain == domain) ++count;
        }
    }

    return count;
}

i32 count_layout_domain_size(
    const PropagationStateLayout& layout,
    PropagationStateDomain domain
) {
    i32 size = 0;
    for (const auto& body : layout.bodies) {
        for (const auto& block : body.blocks) {
            if (block.domain == domain) size += block.size;
        }
    }

    return size;
}

i32 count_layout_integrated(const PropagationStateLayout& layout) {
    i32 count = 0;
    for (const auto& body : layout.bodies) {
        for (const auto& block : body.blocks) {
            if (block.domain == PropagationStateDomain::integrated) ++count;
        }
    }

    return count;
}

i32 count_layout_nonintegrated(const PropagationStateLayout& layout) {
    i32 count = 0;
    for (const auto& body : layout.bodies) {
        for (const auto& block : body.blocks) {
            if (block.domain != PropagationStateDomain::integrated) ++count;
        }
    }

    return count;
}

i32 count_layout_integrated_size(const PropagationStateLayout& layout) {
    i32 size = 0;
    for (const auto& body : layout.bodies) {
        for (const auto& block : body.blocks) {
            if (block.domain == PropagationStateDomain::integrated) size += block.size;
        }
    }

    return size;
}

i32 count_layout_nonintegrated_size(const PropagationStateLayout& layout) {
    i32 size = 0;
    for (const auto& body : layout.bodies) {
        for (const auto& block : body.blocks) {
            if (block.domain != PropagationStateDomain::integrated) size += block.size;
        }
    }

    return size;
}

PropagationBodyIndex find_propagation_body_index(
    const PropagationStateLayout& layout,
    EntityId entity_id
) {
    auto it = layout.body_indices.find(entity_id);
    if (it == layout.body_indices.end()) return kInvalidPropagationBodyIndex;

    return it->second;
}

const PropagationBodyLayout* find_propagation_body_layout(
    const PropagationStateLayout& layout,
    EntityId entity_id
) {
    PropagationBodyIndex body_idx = find_propagation_body_index(layout, entity_id);
    if (body_idx >= layout.bodies.size()) return nullptr;

    const PropagationBodyLayout* body = &layout.bodies[body_idx];
    if (body->entity_id != entity_id) return nullptr;

    return body;
}

PropagationBodyLayout* find_propagation_body_layout(
    PropagationStateLayout& layout,
    EntityId entity_id
) {
    PropagationBodyIndex body_idx = find_propagation_body_index(layout, entity_id);
    if (body_idx >= layout.bodies.size()) return nullptr;

    PropagationBodyLayout* body = &layout.bodies[body_idx];
    if (body->entity_id != entity_id) return nullptr;

    return body;
}

StatusCode add_propagation_state_block(
    PropagationStateLayout& layout,
    EntityId entity_id,
    PropagationStateKind kind,
    const string& key,
    PropagationStateDomain domain,
    i32 size
) {
    PropagationBodyLayout* body = find_propagation_body_layout(layout, entity_id);
    if (body == nullptr) return StatusCode::body_not_found;

    const bool integrated = domain == PropagationStateDomain::integrated;
    if (integrated && size <= 0) return StatusCode::invalid_input;
    if (!integrated && size != 0) return StatusCode::invalid_input;

    PropagationStateBlock block{
        .kind = kind,
        .key = key,
        .domain = domain,
        .offset = integrated ? layout.state_size : -1,
        .size = integrated ? size : 0
    };
    StatusCode status = validate_propagation_state_block(block);
    if (status != StatusCode::ok) return status;

    if (find_state_block(*body, kind, key) != nullptr) return StatusCode::duplicate_id;

    body->blocks.push_back(block);
    layout.state_size += block.size;

    return StatusCode::ok;
}

const PropagationStateBlock* find_state_block(
    const PropagationBodyLayout& body,
    PropagationStateKind kind,
    const string& key
) {
    for (const PropagationStateBlock& block : body.blocks) {
        if (kind == block.kind) {
            if (key == block.key) return &block;
        }
    }

    return nullptr;
}

svec<const PropagationStateBlock*> find_state_blocks(
    const PropagationBodyLayout& body,
    PropagationStateKind kind
) {
    svec<const PropagationStateBlock*> blocks;

    for (const PropagationStateBlock& block : body.blocks) {
        if (kind == block.kind) {
            blocks.push_back(&block);
        }
    }

    return blocks;
}
