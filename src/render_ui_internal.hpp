// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/body.hpp"
#include "core/estimation_common.hpp"
#include "core/scenario_io.hpp"
#include "core/world.hpp"
#include "graphics/render_loop.hpp"
#include "imgui.h"
#include "raylib.h"

namespace render_ui_detail {

ImVec4 rl_to_im(const Vector4& v);
ImVec4 rl_to_im(const Color& c);
void status_text(const char* label, StatusCode code);

void render_simulation_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state);
void render_renderer_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state);
void render_camera_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state);
void render_performance_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state);
void render_scenario_file_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state);

void render_body_stats_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state);
void render_add_body(World& world, RenderLoopConfig& cfg, RenderLoopState& state);
void render_body_lists(World& world, RenderLoopConfig& cfg, RenderLoopState& state);

string make_unique_scenario_body_id(
    const ScenarioSession& session,
    BodyType type,
    EntityId runtime_id
);
StatusCode register_scenario_body_mapping(
    ScenarioBuildResult& result,
    BodyType type,
    const string& config_id,
    EntityId runtime_id
);
StatusCode sync_scenario_body_active(
    ScenarioSession& scenario,
    EntityId id,
    bool active,
    const World& world
);
StatusCode sync_scenario_body(
    ScenarioSession& scenario,
    const Body& body,
    const World& world
);
StatusCode make_scenario_celestial_config(
    ScenarioConfig& scenario,
    const Celestial& cel,
    bool active,
    ScenarioCelestialConfig& out
);
StatusCode make_scenario_satellite_config(
    const Satellite& sat,
    bool active,
    ScenarioSatelliteConfig& out
);
StatusCode make_scenario_station_config(
    const Station& stat,
    bool active,
    const ScenarioBuildResult& mappings,
    ScenarioStationConfig& out
);

} // namespace render_ui_detail
