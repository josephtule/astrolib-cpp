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
    bool step_world = true;
    bool paused = false;
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
};

void run_world_render_loop(World& world, RenderLoopConfig& cfg, f64 dt0);