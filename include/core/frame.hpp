// Copyright 2025-2026 Joseph Tu Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/state.hpp"
#include "util/vecdefs.hpp"

enum struct FrameType {
    inertial,
    body_fixed,
    two_point_rotating,
    topocentric,
    vector_aligned
};

struct FrameTransform {
    vec4d q = q_identity;
    vec3d r = vec3d0;
    // optional:
    vec3d v = vec3d0;
    vec3d w = vec3d0;
};
