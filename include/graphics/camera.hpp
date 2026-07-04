#pragma once

#include "core/entity.hpp"

#include "graphics/raygen.hpp"

#include "util/vecdefs.hpp"

enum struct RenderCameraMode {
    locked, // camera+target locked
    target, // track target
    origin, // track origin
    free,   // track desired direction
    // track_axis
    // track_attitude
};

struct RenderCameraConfig {
    RenderCameraMode mode = RenderCameraMode::origin;

    vec3f position = {0.0f, -20.0f, 10.0f};
    vec3f target = {0.0f, 0.0f, 0.0f};
    EntityId target_id = kInvalidEntityId; // TODO: should this live here?
    vec3f up = {0.0f, 0.0f, 1.0f};
    CameraProjection projection = CAMERA_PERSPECTIVE;
    f32 fovy = 45.0f;

    bool invert_mousewheel = true; // true for "natural scrolling", false otherwise
    f32 zoom_rate = 5.0f;

    f32 orbit_speed = 1.0f;
    f32 pan_speed = 1.0f;
    f32 fly_speed = 10000.0f;
};

inline Camera3D make_render_camera(const RenderCameraConfig& cfg) {
    Camera3D camera{};
    camera.position = eig_to_rl(cfg.position);
    camera.target = eig_to_rl(cfg.target);
    camera.up = eig_to_rl(cfg.up);
    camera.fovy = cfg.fovy;
    camera.projection = cfg.projection;
    return camera;
}
