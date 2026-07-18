// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

class World;
struct RenderLoopConfig;
struct RenderLoopState;

void render_loop_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state);
