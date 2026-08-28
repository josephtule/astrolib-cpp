// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/status.hpp"
#include "core/observation_type.hpp"
#include "core/scenario_io.hpp"
#include "core/world.hpp"
#include "core/world_stepper.hpp"

#include "graphics/camera.hpp"
#include "graphics/render_assets.hpp"
#include "graphics/renderer.hpp"

#include "util/typedefs.hpp"

struct RenderLoopConfig {
    i32 screen_width = 1600;
    i32 screen_height = 900;
    string window_title = "astrolib-cpp";

    bool set_target_fps = true;
    i32 target_fps = 60;

    bool realtime = false;
    WorldStepperConfig stepper_cfg{};
    bool step_single = false;
    bool display_ui = true;

    bool display_body_stats = true;
    bool edit_body_stats = false;
    EntityId body_stats_id = kInvalidEntityId; // TODO: figure out where to put this

    RenderCameraConfig camera{};
    RenderDrawOptions draw{};
    RenderAssetConfig assets{};
};

enum struct SimpleSpinEditMode : i32 { rate = 0, axis = 1 };

struct BodyEditDraft {
    BodyType edit_body_type = BodyType::unknown;
    EntityId edit_body_id = kInvalidEntityId;
    Celestial edit_celestial = Celestial{};
    Satellite edit_satellite = Satellite{};
    Station edit_station = Station{};
    StatusCode edit_body_status = StatusCode::invalid_input;
    SimpleSpinEditMode simple_spin_edit_mode = SimpleSpinEditMode::rate;
    f64 simple_spin_rate = 0.0;
    vec3d simple_spin_axis = vec3d{0.0, 0.0, 1.0};
};

enum struct ScenarioFileOperation { none, load, save_as };

struct ScenarioFileUIState {
    string path_text;
    string resolved_path;
    StatusCode status = StatusCode::file_not_found;
    ScenarioFileOperation operation = ScenarioFileOperation::none;
    bool relative_path = true;
    bool show_status = false;
    bool overwrite_pending = false;
};

struct RenderLoopState {
    WorldStepperStats stats;
    WorldStepperWorkspace wksp;
    Camera3D camera;
    RenderAssets assets;

    bool relative_step = true;
    f64 step_to_time = 0;
    f64 step_by_delta = 0;

    f64 dt = 0.0; // dt used (either dt0 or frametime)

    f32 frame_time = 0.0f;
    f32 fps = 0.0f;
    svec<f32> frame_time_ms;
    svec<f32> fps_history;
    i32 frame_history_max = 240;

    bool add_body = false;
    BodyType add_body_type = BodyType::unknown;
    Celestial temp_celestial;
    Satellite temp_satellite;
    Station temp_station;
    StatusCode add_body_status = StatusCode::ok;

    BodyEditDraft draft;

    i32 add_instrument_type = static_cast<i32>(ObservationType::radec);
    string add_instrument_name = "New Instrument";
    vecXd add_instrument_R_diag = vecXd::Ones(2);
    StatusCode add_instrument_status = StatusCode::ok;

    ScenarioFileUIState scenario_file;
    ScenarioSession scenario;

    BodyFilterMode list_filter = BodyFilterMode::all;
};

void run_world_render_loop(
    World& world,
    RenderLoopConfig& cfg,
    f64 dt0,
    ScenarioSession scenario = {}
);
inline vec3f camera_pivot_from_mode(const RenderCameraConfig& cfg, const World& world) {
    switch (cfg.mode) {
    case RenderCameraMode::locked: return cfg.target;
    case RenderCameraMode::target: {
        if (cfg.target_id == kInvalidEntityId) return vec3f0;
        const Body* body = world.body(cfg.target_id);
        if (body == nullptr) {
            return vec3f0;
        }
        switch (body->body_type) {
        case BodyType::unknown: return vec3f0;
        case BodyType::celestial: return body->x_tr.r.cast<f32>();
        case BodyType::satellite: return body->x_tr.r.cast<f32>();
        case BodyType::station: return world.stat_r_inertial(cfg.target_id).cast<f32>();
        }
    }
    case RenderCameraMode::origin: return vec3f0;
    case RenderCameraMode::free: return cfg.position;
    }

    return vec3f0;
}
inline void sync_camera_tracking(RenderCameraConfig& cfg, const World& world) {
    switch (cfg.mode) {
    case RenderCameraMode::locked:
    case RenderCameraMode::free: break;
    case RenderCameraMode::target:
    case RenderCameraMode::origin: cfg.target = camera_pivot_from_mode(cfg, world); break;
    }
}
inline void cycle_id(EntityId& id, const svec<EntityId>& ids, i32 step) {
    i32 n = static_cast<i32>(ids.size());
    if (n == 0) {
        id = kInvalidEntityId;
        return;
    }

    i32 idx = -1;
    for (i32 i = 0; i < n; ++i) {
        if (id == ids[i]) {
            idx = i;
            break;
        }
    }

    if (idx < 0) {
        id = ids[0];
        return;
    }

    i32 next_idx = (idx + step) % n;
    if (next_idx < 0) next_idx += n;
    id = ids[next_idx];
}
inline void cycle_active_id(EntityId& id, const World& world, i32 step) {
    svec<EntityId> ids = world.active_entity_ids();
    cycle_id(id, ids, step);
}
inline EntityId first_celestial_id(const World& world) {
    svec<EntityId> ids = world.active_celestial_ids();
    if (ids.empty()) return kInvalidEntityId;
    return ids[0];
}
