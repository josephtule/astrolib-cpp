#pragma once

#include "core/world.hpp"
#include "core/world_stepper.hpp"

#include "graphics/camera.hpp"
#include "graphics/renderer.hpp"

#include "util/typedefs.hpp"

struct RenderLoopConfig {
    i32 screen_width = 1280;
    i32 screen_height = 720;
    string window_title = "astrolib-cpp";

    bool set_target_fps = false;
    i32 target_fps = 60;

    bool realtime = false;
    WorldStepperConfig stepper_cfg{};

    RenderCameraConfig camera{};
    RenderDrawOptions draw{};
    RenderAssetConfig assets{};
};

struct RenderLoopState {
    WorldStepperStats stats;
    WorldStepperWorkspace wksp;
    Camera3D camera;
    RenderAssets assets;

    f64 dt = 0.0; // dt used (either dt0 or frametime)

    f32 frame_time = 0.0f;
    f32 fps = 0.0f;
    svec<f32> frame_time_ms;
    svec<f32> fps_history;
    i32 frame_history_max = 240;
};

void run_world_render_loop(World& world, RenderLoopConfig& cfg, f64 dt0);