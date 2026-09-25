#pragma once

#include "core/propagation/state_layout.hpp"
#include "core/status.hpp"
#include "util/vecdefs.hpp"

struct PropagationState {
    vecXd values;
};

struct PropagationBlockView {
    f64* data = nullptr;
    i32 size = 0;

    emap<vecXd> vector() const;
};

struct PropagationConstBlockView {
    const f64* data = nullptr;
    i32 size = 0;

    emap<const vecXd> vector() const;
};

inline emap<vecXd> PropagationBlockView::vector() const {
    return emap<vecXd>{data, size};
}

inline emap<const vecXd> PropagationConstBlockView::vector() const {
    return emap<const vecXd>{data, size};
}

struct PropagationTranslationView {
    f64* data = nullptr;

    emap<vec3d> r() const;
    emap<vec3d> v() const;
};

struct PropagationConstTranslationView {
    const f64* data = nullptr;

    emap<const vec3d> r() const;
    emap<const vec3d> v() const;
};

struct PropagationAttitudeView {
    f64* data = nullptr;

    emap<vec4d> q() const;
    emap<vec3d> w() const;
};

struct PropagationConstAttitudeView {
    const f64* data = nullptr;

    emap<const vec4d> q() const;
    emap<const vec3d> w() const;
};

struct PropagationScalarView {
    f64* data = nullptr;

    f64& value() const;
};

struct PropagationConstScalarView {
    const f64* data = nullptr;

    const f64& value() const;
};

StatusCode validate_propagation_state(
    const PropagationStateLayout& layout,
    const PropagationState& state
);

StatusCode initialize_propagation_state(
    const PropagationStateLayout& layout,
    PropagationState& out
);

StatusCode make_propagation_block_view(
    PropagationState& state,
    const PropagationStateBlock& block,
    PropagationBlockView& out
);

StatusCode make_propagation_block_view(
    const PropagationState& state,
    const PropagationStateBlock& block,
    PropagationConstBlockView& out
);

StatusCode make_propagation_translation_view(
    PropagationState& state,
    const PropagationStateBlock& block,
    PropagationTranslationView& out
);

StatusCode make_propagation_translation_view(
    const PropagationState& state,
    const PropagationStateBlock& block,
    PropagationConstTranslationView& out
);

StatusCode make_propagation_attitude_view(
    PropagationState& state,
    const PropagationStateBlock& block,
    PropagationAttitudeView& out
);

StatusCode make_propagation_attitude_view(
    const PropagationState& state,
    const PropagationStateBlock& block,
    PropagationConstAttitudeView& out
);

StatusCode make_propagation_scalar_view(
    PropagationState& state,
    const PropagationStateBlock& block,
    PropagationScalarView& out
);

StatusCode make_propagation_scalar_view(
    const PropagationState& state,
    const PropagationStateBlock& block,
    PropagationConstScalarView& out
);