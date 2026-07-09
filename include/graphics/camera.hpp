#pragma once

#include "core/entity.hpp"

#include "graphics/raygen.hpp"

#include "util/vecdefs.hpp"

enum struct RenderCameraMode {
    locked = 0, // camera+target locked
    target = 1, // track target
    origin = 2, // track origin
    free = 3,   // track desired direction
    // track_axis
    // track_attitude
};
inline bool is_orbit(RenderCameraMode mode) {
    switch (mode) {
    case RenderCameraMode::locked: return false;
    case RenderCameraMode::target: return true;
    case RenderCameraMode::origin: return true;
    case RenderCameraMode::free: return false;
    }
}
inline string camera_mode_str(const RenderCameraMode& mode) {
    switch (mode) {
    case RenderCameraMode::locked: return "locked";
    case RenderCameraMode::target: return "target";
    case RenderCameraMode::origin: return "origin";
    case RenderCameraMode::free: return "free";
    }
    return "unknoown";
}
struct RenderCameraConfig {
    RenderCameraMode mode = RenderCameraMode::target;

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
    f32 mouse_sensitivity = 0.001f;
    f32 scroll_sensitivity = 0.5f;
    f32 settings_scroll_sensitivity = 2.5f;
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
