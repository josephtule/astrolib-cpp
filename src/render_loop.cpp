#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/transform.hpp"
#include "raylib.h"

#include "graphics/camera.hpp"
#include "graphics/raygen.hpp"
#include "graphics/rdraw.hpp"
#include "graphics/render_assets.hpp"
#include "graphics/render_loop.hpp"
#include "graphics/renderer.hpp"

#include "raymath.h"
#include "util/lightweight_tools.hpp"
#include "util/vecdefs.hpp"

#include <algorithm>
#include <print>

static void render_grids(const RenderDrawOptions& cfg) {
    if (cfg.draw_inertial_axes) {
        DrawLine3D(rlvec30, rlaxis_x * cfg.inertial_axes_scale, RED);
        DrawLine3D(rlvec30, rlaxis_y * cfg.inertial_axes_scale, GREEN);
        DrawLine3D(rlvec30, rlaxis_z * cfg.inertial_axes_scale, BLUE);
        DrawLine3D(rlvec30, rlaxis_x * -cfg.inertial_axes_scale, MAGENTA);
        DrawLine3D(rlvec30, rlaxis_y * -cfg.inertial_axes_scale, YELLOW);
        DrawLine3D(rlvec30, rlaxis_z * -cfg.inertial_axes_scale, CYAN);
    }

    if (cfg.draw_grid_xy) {
        vec2f min = vec2f{cfg.grid_x_min, cfg.grid_y_min};
        vec2f max = vec2f{cfg.grid_x_max, cfg.grid_y_max};
        vec2f inc = vec2f{cfg.u_inc, cfg.v_inc};
        draw_grid(min, max, inc, GridPlane::XY, cfg.color_axes, cfg.grid_color);
    }

    if (cfg.draw_grid_xz) {
        vec2f min = vec2f{cfg.grid_x_min, cfg.grid_z_min};
        vec2f max = vec2f{cfg.grid_x_max, cfg.grid_z_max};
        vec2f inc = vec2f{cfg.u_inc, cfg.v_inc};
        draw_grid(min, max, inc, GridPlane::XZ, cfg.color_axes, cfg.grid_color);
    }

    if (cfg.draw_grid_zy) {
        vec2f min = vec2f{cfg.grid_z_min, cfg.grid_y_min};
        vec2f max = vec2f{cfg.grid_z_max, cfg.grid_y_max};
        vec2f inc = vec2f{cfg.u_inc, cfg.v_inc};
        draw_grid(min, max, inc, GridPlane::ZY, cfg.color_axes, cfg.grid_color);
    }
}

static void draw_selected_body_marker(
    const World& world,
    EntityId selected_id,
    const RenderDrawOptions& draw
) {
    if (!draw.draw_selected_body) return;

    const Body* body = world.body(selected_id);
    if (body == nullptr) return;

    f32 radius = 1.0f;
    switch (body->body_type) {
    case BodyType::unknown: return;
    case BodyType::celestial: {
        const Celestial* cel = world.celestial(selected_id);
        if (cel == nullptr) return;
        radius = std::max(cel->semimajor_axis, cel->semiminor_axis)
                 * draw.selected_marker_scale;
        DrawSphereWires(
            eig_to_rl(cel->x_tr.r),
            radius,
            draw.selected_segements,
            draw.selected_segements,
            draw.selected_color
        );
    } break;
    case BodyType::satellite: {
        radius = 400.0f * draw.selected_marker_scale;
        DrawSphereWires(
            eig_to_rl(body->x_tr.r),
            radius,
            draw.selected_segements,
            draw.selected_segements,
            draw.selected_color
        );
    } break;
    case BodyType::station: {
        radius = 250.0f * draw.selected_marker_scale;
        DrawSphereWires(
            eig_to_rl(world.stat_r_inertial(selected_id)),
            radius,
            draw.selected_segements,
            draw.selected_segements,
            draw.selected_color
        );
        // TODO: center of cylinder model is at the bottom, which is negligible for
        // anchored stations but is wrong for free stations, fix later
    } break;
    }
}

static void render_single_frame(
    const World& world,
    RenderAssets& assets,
    const RenderLoopConfig& cfg,
    Camera3D& camera
) {
    BeginDrawing();
    ClearBackground(cfg.draw.background);
    BeginMode3D(camera);

    render_grids(cfg.draw);
    if (cfg.draw.draw_body_axes) {
        svec<EntityId> ids = world.active_entity_ids();
        for (const EntityId id : ids) {
            const Body* body = world.body(id);
            if (body != nullptr) {
                f32 scale = 1.0f;
                switch (body->body_type) {
                case BodyType::unknown: break;
                case BodyType::station: {
                    scale = 1000.0f;
                    draw_axes(
                        world.stat_x_tr_inertial(id),
                        world.stat_x_att_inertial(id),
                        scale
                    );
                } break;
                case BodyType::celestial: {
                    const Celestial* cel = world.celestial(id);
                    scale = 1.5f
                            * static_cast<f32>(
                                std::max(cel->semimajor_axis, cel->semiminor_axis)
                            );
                    draw_axes(body->x_tr.r, body->x_att.q, scale);
                } break;
                case BodyType::satellite: {
                    scale = 1000.0f;
                    draw_axes(body->x_tr.r, body->x_att.q, scale);
                } break;
                }
            }
        }
    }
    RenderSceneSnapshot scene = build_render_scene_snapshot(world, assets.builtin);
    draw_render_scene(scene, assets);

    draw_selected_body_marker(world, cfg.camera.target_id, cfg.draw);

    EndMode3D();
    if (cfg.draw.draw_fps) DrawFPS(0, 0);
    EndDrawing();
}

static void update_camera(RenderCameraConfig& cfg, Camera3D& camera) {
    camera.fovy = cfg.fovy;
    camera.position = eig_to_rl(cfg.position);
    camera.target = eig_to_rl(cfg.target);
    camera.up = eig_to_rl(cfg.up);
    camera.projection = cfg.projection;
}

static vec3f camera_pivot_from_mode(const RenderCameraConfig& cfg, const World& world) {
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

static void sync_camera_tracking(RenderCameraConfig& cfg, const World& world) {
    switch (cfg.mode) {
    case RenderCameraMode::locked:
    case RenderCameraMode::free: break;
    case RenderCameraMode::target:
    case RenderCameraMode::origin: cfg.target = camera_pivot_from_mode(cfg, world); break;
    }
}

static void cycle_camera_target(
    RenderCameraConfig& camera,
    const World& world,
    i32 step
) {
    svec<EntityId> ids = world.active_entity_ids();
    i32 n = static_cast<i32>(ids.size());
    if (n == 0) {
        camera.target_id = kInvalidEntityId;
        return;
    }

    i32 idx = -1;
    for (i32 i = 0; i < n; ++i) {
        if (camera.target_id == ids[i]) {
            idx = i;
            break;
        }
    }

    if (idx < 0) {
        camera.target_id = ids[0];
        return;
    }

    i32 next_idx = (idx + step) % n;
    if (next_idx < 0) next_idx += n;
    camera.target_id = ids[next_idx];
}

static void init_camera(RenderCameraConfig& camera, const World& world) {
    if (camera.mode != RenderCameraMode::target) return;
    if (camera.target_id == kInvalidEntityId) {
        cycle_camera_target(camera, world, 1);
    }
    sync_camera_tracking(camera, world);
}

static void orbit_camera_about_pivot(
    RenderCameraConfig& cfg,
    const World& world,
    f32 d_az,
    f32 d_el
) {
    vec3f pivot = camera_pivot_from_mode(cfg, world);
    vec3f offset = cfg.position - pivot;

    f32 radius = offset.norm();
    if (radius <= tol6) {
        return;
    }

    vec3f azelr = cart_to_sph<f32>(offset, tol6);
    if (!azelr.allFinite()) {
        return;
    }

    f32 bound = pio2 - tol6;
    azelr(0) += d_az * cfg.orbit_speed;
    azelr(1) = std::clamp(azelr(1) + d_el * cfg.orbit_speed, -bound, bound);

    vec3f new_offset = sph_to_cart(azelr);
    cfg.position = pivot + new_offset;
    cfg.target = pivot;
}

static void orbit_camera_zoom(RenderCameraConfig& cfg, const World& world, f32 wheel) {
    vec3f pivot = camera_pivot_from_mode(cfg, world);
    vec3f offset = cfg.position - pivot;

    f32 radius = offset.norm();
    if (radius < tol6) {
        return;
    }

    radius *= std::exp(-wheel * cfg.zoom_rate);
    cfg.position = pivot + offset.normalized() * radius;
    cfg.target = pivot;
}

static void rotate_free_camera_direction(RenderCameraConfig& cfg, f32 d_az, f32 d_el) {
    vec3f look = cfg.target - cfg.position;
    f32 focus_dist = look.norm();
    if (focus_dist <= tol6) {
        return;
    }

    vec3f azelr = cart_to_sph<f32>(look, tol6);

    f32 bound = pio2 - tol6;
    azelr(0) += d_az * cfg.pan_speed;
    azelr(1) = std::clamp(azelr(1) + d_el * cfg.pan_speed, -bound, bound);
    azelr(2) = 1.0f;

    vec3f new_dir = sph_to_cart(azelr);
    cfg.target = cfg.position + new_dir * focus_dist;
}

static void move_free_camera(
    RenderCameraConfig& cfg,
    f32 dt,
    f32 forward_input,
    f32 right_input,
    f32 up_input
) {
    vec3f look = cfg.target - cfg.position;
    f32 focus_dist = look.norm();
    if (focus_dist <= tol6) {
        return;
    }

    vec3f forward = look.normalized();
    vec3f world_up = cfg.up.normalized();
    vec3f right = forward.cross(world_up);
    f32 right_norm = right.norm();
    if (right_norm <= tol6) {
        return;
    }
    right /= right_norm;

    vec3f delta = cfg.fly_speed * dt
                  * (forward_input * forward + right_input * right + up_input * world_up);

    cfg.position += delta;
    cfg.target += delta;
}

static void handle_sim_input(RenderLoopConfig& cfg) {
    if (IsKeyPressed(KEY_SPACE)) {
        toggle(cfg.paused);
    }
    if (IsKeyPressed(KEY_UP)) {
        cfg.stepper_cfg.ticks += 1;
    } else if (IsKeyPressed(KEY_DOWN)) {
        cfg.stepper_cfg.ticks -= 1;
        if (cfg.stepper_cfg.ticks < 1) cfg.stepper_cfg.ticks = 1;
    } else if (IsKeyPressed(KEY_RIGHT)) {
        cfg.stepper_cfg.substeps += 1;
    } else if (IsKeyPressed(KEY_LEFT)) {
        cfg.stepper_cfg.substeps -= 1;
        if (cfg.stepper_cfg.substeps < 1) cfg.stepper_cfg.substeps = 1;
    }
}

static void handle_render_toggle_input(RenderDrawOptions& draw) {
    if (IsKeyPressed(KEY_COMMA)) {
        toggle(draw.draw_body_axes);
    } else if (IsKeyPressed(KEY_PERIOD)) {
        toggle(draw.draw_inertial_axes);
    } else if (IsKeyPressed(KEY_SLASH)) {
        toggle(draw.color_axes);
    } else if (IsKeyPressed(KEY_APOSTROPHE)) {
        toggle(draw.draw_selected_body);
    }

    if (IsKeyPressed(KEY_ONE)) {
        toggle(draw.draw_grid_xy);
    } else if (IsKeyPressed(KEY_TWO)) {
        toggle(draw.draw_grid_xz);
    } else if (IsKeyPressed(KEY_THREE)) {
        toggle(draw.draw_grid_zy);
    }

    if (IsKeyPressed(KEY_GRAVE)) {
        toggle(draw.draw_fps);
    }

    if (IsKeyPressed(KEY_BACKSPACE)) {
        RenderDrawOptions default_cfg{};
        draw.draw_body_axes = default_cfg.draw_body_axes;
        draw.draw_inertial_axes = default_cfg.draw_inertial_axes;
        draw.draw_fps = default_cfg.draw_fps;
        draw.draw_grid_xy = default_cfg.draw_grid_xy;
        draw.draw_grid_xz = default_cfg.draw_grid_xz;
        draw.draw_grid_zy = default_cfg.draw_grid_zy;
        draw.draw_labels = default_cfg.draw_labels;
        draw.color_axes = default_cfg.color_axes;
        draw.draw_selected_body = default_cfg.draw_selected_body;
    }
}

static void handle_camera_settings_input(
    RenderCameraConfig& camera,
    const World& world,
    f32 dt
) {
    if (IsKeyPressed(KEY_ONE)) {
        camera.mode = RenderCameraMode::locked;
    } else if (IsKeyPressed(KEY_TWO)) {
        camera.mode = RenderCameraMode::target;
        if (camera.target_id == kInvalidEntityId) cycle_camera_target(camera, world, 1);
    } else if (IsKeyPressed(KEY_THREE)) {
        camera.mode = RenderCameraMode::origin;
    } else if (IsKeyPressed(KEY_FOUR)) {
        camera.mode = RenderCameraMode::free;
    }

    f32 wheel = GetMouseWheelMove();
    if (wheel != 0.0f) {
        f32 dir = camera.invert_mousewheel ? -1.0f : 1.0f;
        camera.zoom_rate += camera.settings_scroll_sensitivity * dt * dir * wheel;
        camera.zoom_rate = std::clamp(camera.zoom_rate, 0.05f, 15.0f);
    }

    if (camera.mode == RenderCameraMode::free) {
        if (IsKeyPressed(KEY_UP)) {
            camera.fly_speed += 1000.0f;
        } else if (IsKeyPressed(KEY_RIGHT)) {
            camera.fly_speed += 100.0f;
        } else if (IsKeyPressed(KEY_DOWN)) {
            camera.fly_speed = std::max(100.0f, camera.fly_speed - 1000.0f);
        } else if (IsKeyPressed(KEY_LEFT)) {
            camera.fly_speed = std::max(100.0f, camera.fly_speed - 100.0f);
        }
    } else {
        if (IsKeyPressed(KEY_UP)) {
            camera.orbit_speed += 1.0f;
        } else if (IsKeyPressed(KEY_RIGHT)) {
            camera.orbit_speed += 0.1f;
        } else if (IsKeyPressed(KEY_DOWN)) {
            camera.orbit_speed = std::max(0.1f, camera.orbit_speed - 1.0f);
        } else if (IsKeyPressed(KEY_LEFT)) {
            camera.orbit_speed = std::max(0.1f, camera.orbit_speed - 0.1f);
        }
    }

    if (camera.mode == RenderCameraMode::target) {
        if (IsKeyPressed(KEY_RIGHT_BRACKET)) {
            cycle_camera_target(camera, world, 1);
        } else if (IsKeyPressed(KEY_LEFT_BRACKET)) {
            cycle_camera_target(camera, world, -1);
        }
    }

    if (IsKeyPressed(KEY_GRAVE)) {
        toggle(camera.invert_mousewheel);
    }

    if (IsKeyPressed(KEY_MINUS)) {
        camera.fovy = std::clamp<f32>(camera.fovy - camera.zoom_rate, 1, 179);
    } else if (IsKeyPressed(KEY_EQUAL)) {
        camera.fovy = std::clamp<f32>(camera.fovy + camera.zoom_rate, 1, 179);
    }

    if (IsKeyPressed(KEY_BACKSLASH)) {
        RenderCameraConfig default_cfg{};
        camera.fovy = default_cfg.fovy;
        camera.zoom_rate = default_cfg.zoom_rate;
    }

    if (IsKeyPressed(KEY_BACKSPACE)) {
        camera = RenderCameraConfig{};
    }
}

static void handle_camera_movement_input(
    RenderCameraConfig& camera,
    const World& world,
    f32 dt
) {
    f32 d_az = 0.0f;
    f32 d_el = 0.0f;
    f32 d_forward = 0.0f;
    f32 d_right = 0.0f;
    f32 d_up = 0.0f;

    if (camera.mode == RenderCameraMode::free) {
        if (IsKeyDown(KEY_W)) d_forward += 1.0f;
        if (IsKeyDown(KEY_S)) d_forward -= 1.0f;
        if (IsKeyDown(KEY_D)) d_right += 1.0f;
        if (IsKeyDown(KEY_A)) d_right -= 1.0f;
        if (IsKeyDown(KEY_E)) d_up += 1.0f;
        if (IsKeyDown(KEY_Q)) d_up -= 1.0f;
    } else {
        if (IsKeyDown(KEY_W)) d_el += dt;
        if (IsKeyDown(KEY_S)) d_el -= dt;
        if (IsKeyDown(KEY_D)) d_az += dt;
        if (IsKeyDown(KEY_A)) d_az -= dt;
    }

    f32 wheel = GetMouseWheelMove();
    if (wheel != 0.0f) {
        f32 dir = camera.invert_mousewheel ? -1.0f : 1.0f;
        if (camera.mode == RenderCameraMode::free) {
            f32 wheel_step = wheel * camera.scroll_sensitivity * dir;
            d_forward += wheel_step;
        } else if (
            camera.mode == RenderCameraMode::target
            || camera.mode == RenderCameraMode::origin
        ) {
            f32 wheel_step = wheel * dt * camera.scroll_sensitivity * dir;
            orbit_camera_zoom(camera, world, wheel_step);
        }
    }

    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        Vector2 mouse_delta = GetMouseDelta();
        d_az -= mouse_delta.x * camera.mouse_sensitivity;
        d_el += mouse_delta.y * camera.mouse_sensitivity;
    }

    switch (camera.mode) {
    case RenderCameraMode::locked: break;
    case RenderCameraMode::target:
        orbit_camera_about_pivot(camera, world, d_az, d_el);
        break;
    case RenderCameraMode::origin:
        orbit_camera_about_pivot(camera, world, d_az, d_el);
        break;
    case RenderCameraMode::free: {
        move_free_camera(camera, dt, d_forward, d_right, d_up);
        rotate_free_camera_direction(camera, -d_az, d_el);
    } break;
    }
}

static void handle_input(RenderLoopConfig& cfg, const World& world, f32 dt) {
    // TODO: create a key map and key map settings page where users can change map

    bool render_modifier = IsKeyDown(KEY_LEFT_SHIFT);
    bool camera_modifier = IsKeyDown(KEY_RIGHT_SHIFT);
    bool no_modifier = !render_modifier && !camera_modifier;
    bool control_modifier = no_modifier && IsKeyDown(KEY_LEFT_CONTROL);

    if (camera_modifier) {
        handle_camera_settings_input(cfg.camera, world, dt);
    } else if (render_modifier) {
        handle_render_toggle_input(cfg.draw);
    } else if (control_modifier) {
        // reserved for future non-render shortcuts
    } else {
        handle_sim_input(cfg);
        handle_camera_movement_input(cfg.camera, world, dt);
    }

    sync_camera_tracking(cfg.camera, world);
}

struct RenderLoopState {
    WorldStepperStats stats;
    WorldStepperWorkspace wksp;
    Camera3D camera;
    RenderAssets assets;
};

static void init_render_loop_state(
    const World& world,
    RenderLoopConfig& cfg,
    RenderLoopState& state
) {
    state.stats.success = true;

    state.camera = make_render_camera(cfg.camera);
    state.assets = load_render_assets(cfg.assets);
    init_camera(cfg.camera, world);
    update_camera(cfg.camera, state.camera);
}

static void shutdown_render_loop_state(RenderLoopState& state) {
    unload_render_assets(state.assets);
}

void run_world_render_loop(World& world, RenderLoopConfig& cfg, f64 dt0) {
    InitWindow(cfg.screen_width, cfg.screen_height, cfg.window_title.c_str());

    RenderLoopState state;
    init_render_loop_state(world, cfg, state);

    f64 dt = dt0;
    f32 dt_window;

    render_single_frame(world, state.assets, cfg, state.camera);
    while (!WindowShouldClose()) {
        if (cfg.set_target_fps) SetTargetFPS(cfg.target_fps);
        dt_window = GetFrameTime();

        handle_input(cfg, world, dt_window);
        update_camera(cfg.camera, state.camera);

        dt = cfg.realtime ? dt_window : dt0;
        if (!cfg.paused) {
            state.stats += step_world(world, dt, cfg.stepper_cfg, state.wksp);
            if (!state.stats.success) {
                std::println("Simulation Failure");
                break;
            }
        }
        render_single_frame(world, state.assets, cfg, state.camera);
    }

    shutdown_render_loop_state(state);
    CloseWindow();
}

static void render_ui_placeholder(
    World& world,
    RenderLoopConfig& cfg,
    RenderLoopState& state
) {}