#include "core/propagation/state_world.hpp"

#include "core/body.hpp"
#include "core/propagation/state.hpp"
#include "core/propagation/state_layout.hpp"
#include "core/state.hpp"
#include "core/status.hpp"
#include "core/world.hpp"
#include "core/world/world_revision.hpp"
#include "util/math.hpp"

#include <cassert>
#include <utility>

static StatusCode validate_world_state_block_binding(const PropagationStateBlock& block) {
    switch (block.kind) {
    case PropagationStateKind::translation: [[fallthrough]];
    case PropagationStateKind::attitude: [[fallthrough]];
    case PropagationStateKind::mass: {
        if (block.key != default_block_key) return StatusCode::unsupported_method;
    } break;
    case PropagationStateKind::resource: [[fallthrough]];
    case PropagationStateKind::parameter: break;
    default: return StatusCode::unsupported_type;
    }

    return StatusCode::ok;
}

// Packing Helpers --------------------------------------------------------------

static StatusCode pack_world_translation_state(
    const Body& body,
    const PropagationStateBlock& block,
    PropagationState& state
) {
    PropagationTranslationView view{};
    StatusCode status = make_propagation_translation_view(state, block, view);
    if (status != StatusCode::ok) return status;

    view.r() = body.x_tr.r;
    view.v() = body.x_tr.v;

    return StatusCode::ok;
}

static StatusCode pack_world_attitude_state(
    const Body& body,
    const PropagationStateBlock& block,
    PropagationState& state
) {
    PropagationAttitudeView view{};
    StatusCode status = make_propagation_attitude_view(state, block, view);
    if (status != StatusCode::ok) return status;

    view.q() = body.x_att.q;
    view.w() = body.x_att.w;

    return StatusCode::ok;
}

static StatusCode pack_world_mass_state(
    const Body& body,
    const PropagationStateBlock& block,
    PropagationState& state
) {
    const MassProperties* mass_properties = nullptr;
    StatusCode status = mass_properties_from_body(body, mass_properties);
    if (status != StatusCode::ok) return status;
    if (mass_properties == nullptr || !mass_properties->active) {
        return StatusCode::invalid_mass_properties;
    }
    if (!finite_pos(mass_properties->mass)) {
        return StatusCode::invalid_mass_properties;
    }

    PropagationScalarView view{};
    status = make_propagation_scalar_view(state, block, view);
    if (status != StatusCode::ok) return status;

    view.value() = mass_properties->mass;
    return StatusCode::ok;
}

static StatusCode pack_world_state_block(
    const Body& body,
    const PropagationStateBlock& block,
    PropagationState& state
) {
    if (block.domain != PropagationStateDomain::integrated) return StatusCode::ok;

    StatusCode status = validate_world_state_block_binding(block);
    if (status != StatusCode::ok) return status;

    switch (block.kind) {
    case PropagationStateKind::translation:
        return pack_world_translation_state(body, block, state);
    case PropagationStateKind::attitude:
        return pack_world_attitude_state(body, block, state);
    case PropagationStateKind::mass: return pack_world_mass_state(body, block, state);
    case PropagationStateKind::resource: [[fallthrough]];
    case PropagationStateKind::parameter: return StatusCode::unsupported_method;
    default: return StatusCode::unsupported_type;
    }
}

static StatusCode pack_world_body_state(
    const Body& body,
    const PropagationBodyLayout& body_layout,
    PropagationState& state
) {
    if (body.id != body_layout.entity_id) return StatusCode::invalid_input;

    for (const auto& block : body_layout.blocks) {
        StatusCode status = pack_world_state_block(body, block, state);
        if (status != StatusCode::ok) return status;
    }

    return StatusCode::ok;
}

// Commit Validation Helpers ----------------------------------------------------

static StatusCode validate_world_translation_commit(
    const PropagationStateBlock& block,
    const PropagationState& state
) {
    PropagationConstTranslationView view{};
    StatusCode status = make_propagation_translation_view(state, block, view);
    if (status != StatusCode::ok) return status;

    StateTr temp{.r = view.r(), .v = view.v()};
    if (!finite_state_tr(temp)) return StatusCode::invalid_state;

    return StatusCode::ok;
}

static StatusCode validate_world_attitude_commit(
    const PropagationStateBlock& block,
    const PropagationState& state
) {
    PropagationConstAttitudeView view{};
    StatusCode status = make_propagation_attitude_view(state, block, view);
    if (status != StatusCode::ok) return status;

    StateAtt temp{.q = view.q(), .w = view.w()};
    if (!finite_state_att(temp)) return StatusCode::invalid_att_state;

    return StatusCode::ok;
}

static StatusCode validate_world_mass_commit(
    const Body& body,
    const PropagationStateBlock& block,
    const PropagationState& state
) {
    const MassProperties* mass_properties = nullptr;
    StatusCode status = mass_properties_from_body(body, mass_properties);
    if (status != StatusCode::ok) return status;
    if (mass_properties == nullptr || !mass_properties->active) {
        return StatusCode::invalid_mass_properties;
    }

    PropagationConstScalarView view{};
    status = make_propagation_scalar_view(state, block, view);
    if (status != StatusCode::ok) return status;
    if (!finite_pos(view.value())) return StatusCode::invalid_mass_properties;

    return StatusCode::ok;
}

static StatusCode validate_world_state_commit_block(
    const Body& body,
    const PropagationStateBlock& block,
    const PropagationState& state
) {
    if (block.domain != PropagationStateDomain::integrated) return StatusCode::ok;

    StatusCode status = validate_world_state_block_binding(block);
    if (status != StatusCode::ok) return status;

    switch (block.kind) {
    case PropagationStateKind::translation:
        return validate_world_translation_commit(block, state);
    case PropagationStateKind::attitude:
        return validate_world_attitude_commit(block, state);
    case PropagationStateKind::mass:
        return validate_world_mass_commit(body, block, state);
    case PropagationStateKind::resource: [[fallthrough]];
    case PropagationStateKind::parameter: return StatusCode::unsupported_method;
    default: return StatusCode::unsupported_type;
    }
}

static StatusCode validate_world_body_commit(
    const Body& body,
    const PropagationBodyLayout& body_layout,
    const PropagationState& state
) {
    if (body.id != body_layout.entity_id) return StatusCode::invalid_input;

    for (const auto& block : body_layout.blocks) {
        StatusCode status = validate_world_state_commit_block(body, block, state);
        if (status != StatusCode::ok) return status;
    }

    return StatusCode::ok;
}

// Commit Application Helpers ---------------------------------------------------

static void apply_world_translation_commit(
    Body& body,
    const PropagationStateBlock& block,
    const PropagationState& state
) {
    PropagationConstTranslationView view{state.values.data() + block.offset};
    body.x_tr.r = view.r();
    body.x_tr.v = view.v();
}

static void apply_world_attitude_commit(
    Body& body,
    const PropagationStateBlock& block,
    const PropagationState& state
) {
    PropagationConstAttitudeView view{state.values.data() + block.offset};
    body.x_att.q = view.q();
    normalize_quaternion_inplace(body.x_att.q);
    body.x_att.w = view.w();
}

static void apply_world_mass_commit(
    Body& body,
    const PropagationStateBlock& block,
    const PropagationState& state
) {
    MassProperties* mass_properties = nullptr;
    StatusCode status = mass_properties_from_body(body, mass_properties);
    assert(status == StatusCode::ok);
    assert(mass_properties != nullptr);

    if (status != StatusCode::ok || mass_properties == nullptr) return;

    PropagationConstScalarView view{state.values.data() + block.offset};
    mass_properties->mass = view.value();
}

static void apply_world_state_commit_block(
    Body& body,
    const PropagationStateBlock& block,
    const PropagationState& state
) {
    if (block.domain != PropagationStateDomain::integrated) return;

    switch (block.kind) {
    case PropagationStateKind::translation:
        apply_world_translation_commit(body, block, state);
        break;
    case PropagationStateKind::attitude:
        apply_world_attitude_commit(body, block, state);
        break;
    case PropagationStateKind::mass: apply_world_mass_commit(body, block, state); break;
    case PropagationStateKind::resource: [[fallthrough]];
    case PropagationStateKind::parameter: [[fallthrough]];
    default: return;
    }
}

static void apply_world_body_commit(
    Body& body,
    const PropagationBodyLayout& body_layout,
    const PropagationState& state
) {
    for (const auto& block : body_layout.blocks) {
        apply_world_state_commit_block(body, block, state);
    }
}

// Public Conversion Functions --------------------------------------------------

StatusCode pack_world_propagation_state(
    const World& world,
    const PropagationStateLayout& layout,
    PropagationState& out
) {
    StatusCode status;

    status = validate_world_propagation_state_layout(world, layout);
    if (status != StatusCode::ok) return status;

    PropagationState temp;
    status = initialize_propagation_state(layout, temp);
    if (status != StatusCode::ok) return status;

    for (const auto& body_layout : layout.bodies) {
        const Body* body = world.body(body_layout.entity_id);
        if (body == nullptr) return StatusCode::body_not_found;

        status = pack_world_body_state(*body, body_layout, temp);
        if (status != StatusCode::ok) return status;
    }

    status = validate_propagation_state(layout, temp);
    if (status != StatusCode::ok) return status;

    out = std::move(temp);
    return StatusCode::ok;
}

StatusCode commit_world_propagation_state(
    World& world,
    const PropagationStateLayout& layout,
    const PropagationState& state
) {
    StatusCode status;

    status = validate_world_propagation_state_layout(world, layout);
    if (status != StatusCode::ok) return status;

    status = validate_propagation_state(layout, state);
    if (status != StatusCode::ok) return status;

    for (const auto& body_layout : layout.bodies) {
        const Body* body = world.body(body_layout.entity_id);
        if (body == nullptr) return StatusCode::body_not_found;

        status = validate_world_body_commit(*body, body_layout, state);
        if (status != StatusCode::ok) return status;
    }

    for (const auto& body_layout : layout.bodies) {
        Body* body = world.body(body_layout.entity_id);
        assert(body != nullptr);
        if (body == nullptr) return StatusCode::body_not_found;

        apply_world_body_commit(*body, body_layout, state);
    }
    world.mark_state_changed();

    return StatusCode::ok;
}

StatusCode finalize_world_propagation_state_layout(
    const World& world,
    PropagationStateLayout& layout
) {
    StatusCode status;
    status = validate_propagation_state_layout(layout);
    if (status != StatusCode::ok) return status;

    const WorldRevisions& revisions = world.get_revisions();
    layout.revisions.topology = revisions.topology;
    layout.revisions.dynamics = revisions.dynamics;
    layout.revisions.providers = revisions.providers;
    layout.revisions.bound = true;

    return StatusCode::ok;
}

StatusCode validate_world_propagation_state_layout(
    const World& world,
    const PropagationStateLayout& layout
) {
    if (!layout.revisions.bound) return StatusCode::stale_revision;

    StatusCode status;
    status = validate_propagation_state_layout(layout);
    if (status != StatusCode::ok) return status;

    const WorldRevisions& revisions = world.get_revisions();
    if (layout.revisions.topology != revisions.topology
        || layout.revisions.dynamics != revisions.dynamics
        || layout.revisions.providers != revisions.providers) {
        return StatusCode::stale_revision;
    }

    return StatusCode::ok;
}
