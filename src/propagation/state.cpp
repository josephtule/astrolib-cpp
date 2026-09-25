#include "core/propagation/state.hpp"
#include "core/propagation/state_layout.hpp"
#include "core/status.hpp"

StatusCode validate_propagation_state(
    const PropagationStateLayout& layout,
    const PropagationState& state
) {
    if (state.values.size() != layout.state_size) return StatusCode::size_mismatch;
    if (!state.values.allFinite()) return StatusCode::non_finite_result;

    return StatusCode::ok;
}

StatusCode initialize_propagation_state(
    const PropagationStateLayout& layout,
    PropagationState& out
) {
    StatusCode status;

    status = validate_propagation_state_layout(layout);
    if (status != StatusCode::ok) return status;

    PropagationState temp;
    temp.values = vecXd::Zero(layout.state_size);

    status = validate_propagation_state(layout, temp);
    if (status != StatusCode::ok) return status;

    // TODO: complete this

    out = std::move(temp);
    return StatusCode::ok;
}

static StatusCode validate_propagation_block_view(
    const PropagationState& state,
    const PropagationStateBlock& block
) {
    if (block.domain != PropagationStateDomain::integrated)
        return StatusCode::invalid_input;
    if (block.offset < 0 || block.size <= 0) return StatusCode::invalid_input;
    if (block.offset > state.values.size() - block.size) return StatusCode::size_mismatch;

    return StatusCode::ok;
}
StatusCode make_propagation_block_view(
    PropagationState& state,
    const PropagationStateBlock& block,
    PropagationBlockView& out
) {
    StatusCode status;
    out = {};

    status = validate_propagation_block_view(state, block);
    if (status != StatusCode::ok) return status;

    PropagationBlockView temp;
    temp.data = state.values.data() + block.offset;
    temp.size = block.size;

    out = std::move(temp);
    return StatusCode::ok;
}

StatusCode make_propagation_block_view(
    const PropagationState& state,
    const PropagationStateBlock& block,
    PropagationConstBlockView& out
) {
    StatusCode status;
    out = {};

    status = validate_propagation_block_view(state, block);
    if (status != StatusCode::ok) return status;

    PropagationConstBlockView temp;
    temp.data = state.values.data() + block.offset;
    temp.size = block.size;

    out = std::move(temp);
    return StatusCode::ok;
}

StatusCode make_propagation_translation_view(
    PropagationState& state,
    const PropagationStateBlock& block,
    PropagationTranslationView& out
) {
    out = {};

    if (block.kind != PropagationStateKind::translation) return StatusCode::invalid_input;
    if (block.domain != PropagationStateDomain::integrated)
        return StatusCode::invalid_input;
    if (block.size != 6) return StatusCode::invalid_input;

    StatusCode status;
    status = validate_propagation_block_view(state, block);
    if (status != StatusCode::ok) return status;

    PropagationTranslationView temp;
    temp.data = state.values.data() + block.offset;

    out = std::move(temp);
    return StatusCode::ok;
}

StatusCode make_propagation_translation_view(
    const PropagationState& state,
    const PropagationStateBlock& block,
    PropagationConstTranslationView& out
) {
    out = {};

    if (block.kind != PropagationStateKind::translation) return StatusCode::invalid_input;
    if (block.domain != PropagationStateDomain::integrated)
        return StatusCode::invalid_input;
    if (block.size != 6) return StatusCode::invalid_input;

    StatusCode status;
    status = validate_propagation_block_view(state, block);
    if (status != StatusCode::ok) return status;

    PropagationConstTranslationView temp;
    temp.data = state.values.data() + block.offset;

    out = std::move(temp);
    return StatusCode::ok;
}

StatusCode make_propagation_attitude_view(
    PropagationState& state,
    const PropagationStateBlock& block,
    PropagationAttitudeView& out
) {
    out = {};

    if (block.kind != PropagationStateKind::attitude) return StatusCode::invalid_input;
    if (block.domain != PropagationStateDomain::integrated)
        return StatusCode::invalid_input;
    if (block.size != 7) return StatusCode::invalid_input;

    StatusCode status;
    status = validate_propagation_block_view(state, block);
    if (status != StatusCode::ok) return status;

    PropagationAttitudeView temp;
    temp.data = state.values.data() + block.offset;

    out = std::move(temp);
    return StatusCode::ok;
}

StatusCode make_propagation_attitude_view(
    const PropagationState& state,
    const PropagationStateBlock& block,
    PropagationConstAttitudeView& out
) {
    out = {};

    if (block.kind != PropagationStateKind::attitude) return StatusCode::invalid_input;
    if (block.domain != PropagationStateDomain::integrated)
        return StatusCode::invalid_input;
    if (block.size != 7) return StatusCode::invalid_input;

    StatusCode status;
    status = validate_propagation_block_view(state, block);
    if (status != StatusCode::ok) return status;

    PropagationConstAttitudeView temp;
    temp.data = state.values.data() + block.offset;

    out = std::move(temp);
    return StatusCode::ok;
}

StatusCode make_propagation_scalar_view(
    PropagationState& state,
    const PropagationStateBlock& block,
    PropagationScalarView& out
) {
    out = {};

    if (block.domain != PropagationStateDomain::integrated)
        return StatusCode::invalid_input;
    if (block.size != 1) return StatusCode::invalid_input;
    if (block.kind != PropagationStateKind::mass
        && block.kind != PropagationStateKind::resource
        && block.kind != PropagationStateKind::parameter)
        return StatusCode::invalid_input;

    StatusCode status;
    status = validate_propagation_block_view(state, block);
    if (status != StatusCode::ok) return status;

    PropagationScalarView temp;
    temp.data = state.values.data() + block.offset;

    out = std::move(temp);
    return StatusCode::ok;
}

StatusCode make_propagation_scalar_view(
    const PropagationState& state,
    const PropagationStateBlock& block,
    PropagationConstScalarView& out
) {
    out = {};

    if (block.domain != PropagationStateDomain::integrated)
        return StatusCode::invalid_input;
    if (block.size != 1) return StatusCode::invalid_input;
    if (block.kind != PropagationStateKind::mass
        && block.kind != PropagationStateKind::resource
        && block.kind != PropagationStateKind::parameter)
        return StatusCode::invalid_input;

    StatusCode status;
    status = validate_propagation_block_view(state, block);
    if (status != StatusCode::ok) return status;

    PropagationConstScalarView temp;
    temp.data = state.values.data() + block.offset;

    out = std::move(temp);
    return StatusCode::ok;
}

emap<vec3d> PropagationTranslationView::r() const { return emap<vec3d>(data); }
emap<vec3d> PropagationTranslationView::v() const { return emap<vec3d>(data + 3); }
emap<const vec3d> PropagationConstTranslationView::r() const {
    return emap<const vec3d>(data);
}
emap<const vec3d> PropagationConstTranslationView::v() const {
    return emap<const vec3d>(data + 3);
}

emap<vec4d> PropagationAttitudeView::q() const { return emap<vec4d>(data); }
emap<vec3d> PropagationAttitudeView::w() const { return emap<vec3d>(data + 4); }
emap<const vec4d> PropagationConstAttitudeView::q() const {
    return emap<const vec4d>(data);
}
emap<const vec3d> PropagationConstAttitudeView::w() const {
    return emap<const vec3d>(data + 4);
}

f64& PropagationScalarView::value() const { return *data; }
const f64& PropagationConstScalarView::value() const { return *data; }
