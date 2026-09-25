#include "core/propagation/state_layout.hpp"
#include "core/entity.hpp"
#include "core/status.hpp"

#include <algorithm>
#include <set>
#include <limits>

static StatusCode validate_propagation_state_block(const PropagationStateBlock& block) {
    if (block.key.empty()) return StatusCode::invalid_input;

    if (block.domain == PropagationStateDomain::integrated) {
        if (block.size <= 0) return StatusCode::invalid_input;
        if (block.offset < 0) return StatusCode::invalid_input;
    } else {
        if (block.size != 0) return StatusCode::invalid_input;
        if (block.offset != -1) return StatusCode::invalid_input;
    }

    return StatusCode::ok;
}

static StatusCode validate_propagation_body_layout(const PropagationBodyLayout& body) {
    if (body.entity_id == kInvalidEntityId) return StatusCode::body_not_found;
    if (body.body_index < 0) return StatusCode::invalid_input;

    StatusCode status;

    std::set<std::pair<PropagationStateKind, string>> seen_blocks;

    for (const auto& block : body.blocks) {
        if (!seen_blocks.insert({block.kind, block.key}).second) { // duplicate
            return StatusCode::invalid_input;
        }

        status = validate_propagation_state_block(block);
        if (status != StatusCode::ok) return status;
    }

    return StatusCode::ok;
}

StatusCode validate_propagation_state_layout(const PropagationStateLayout& layout) {
    StatusCode status;

    std::set<EntityId> seen_entity_ids;
    std::set<i32> seen_body_indices;
    svec<std::pair<i32, i32>> ranges;
    for (const auto& body : layout.bodies) {
        if (!seen_entity_ids.insert(body.entity_id).second)
            return StatusCode::invalid_input;
        if (!seen_body_indices.insert(body.body_index).second)
            return StatusCode::invalid_input;
        if (body.body_index >= layout.bodies.size()) return StatusCode::invalid_input;

        status = validate_propagation_body_layout(body);
        if (status != StatusCode::ok) return status;

        for (const auto& block : body.blocks) {
            if (block.domain != PropagationStateDomain::integrated) continue;

            if (block.offset > std::numeric_limits<i32>::max() - block.size) {
                return StatusCode::invalid_input;
            }
            i32 end = block.offset + block.size;
            ranges.push_back({block.offset, end});
        }
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
