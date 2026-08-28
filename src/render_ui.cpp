// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "imgui_internal.h"
#include "render_ui_internal.hpp"

#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/status.hpp"
#include "core/time.hpp"
#include "core/transform.hpp"
#include "core/world.hpp"
#include "core/world_stepper.hpp"
#include "graphics/camera.hpp"
#include "graphics/render_loop.hpp"
#include "graphics/ui.hpp"

#include "imgui.h"
#include "implot.h"

#include "raylib.h"
#include "util/constants.hpp"
#include "util/lightweight_tools.hpp"
#include "util/math.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

#include <algorithm>
#include <string>

namespace render_ui_detail {

namespace im = ImGui;
namespace imp = ImPlot;

static Color status_color(const StatusCode code) {
    switch (code) {
    // Success
    case StatusCode::ok: return DARKGREEN;

    // Informational non-failure states
    case StatusCode::file_overwritten:
    case StatusCode::prediction_only: return SKYBLUE;

    // Warning states that may still produce a usable result
    case StatusCode::max_iters_reached:
    case StatusCode::correction_rejected: return GOLD;

    // Input/configuration errors
    case StatusCode::invalid_input:
    case StatusCode::validation_failed:
    case StatusCode::unsupported_type:
    case StatusCode::missing_reference:
    case StatusCode::duplicate_id:
    case StatusCode::inactive_entity:
    case StatusCode::invalid_state:
    case StatusCode::invalid_attitude_state:
    case StatusCode::invalid_mass_properties:
    case StatusCode::invalid_shape:
    case StatusCode::anchor_not_found:
    case StatusCode::size_mismatch:
    case StatusCode::time_mismatch:
    case StatusCode::invalid_covariance:
    case StatusCode::unsupported_method:
    case StatusCode::parse_failed: return ORANGE;

    // Missing data/reference lookups
    case StatusCode::observer_not_found:
    case StatusCode::target_not_found:
    case StatusCode::instrument_not_found:
    case StatusCode::body_not_found:
    case StatusCode::gravity_model_not_found:
    case StatusCode::attitude_type_not_found:
    case StatusCode::celestial_model_not_found:
    case StatusCode::sample_not_found:
    case StatusCode::empty_measurements:
    case StatusCode::empty_history:
    case StatusCode::empty_events: return YELLOW;

    // Runtime/solver/IO failures
    case StatusCode::propagation_failed:
    case StatusCode::singular_normal_matrix:
    case StatusCode::singular_innovation:
    case StatusCode::interp_failed:
    case StatusCode::file_not_found:
    case StatusCode::file_write_failed:
    case StatusCode::file_close_failed:
    case StatusCode::file_open_failed:
    case StatusCode::file_already_exists:
    case StatusCode::matrix_invert_failed:
    case StatusCode::step_size_underflow:
    case StatusCode::max_steps_reached:
    case StatusCode::max_rejections_reached:
    case StatusCode::non_finite_derivative:
    case StatusCode::non_finite_result: return RED;
    }

    return RAYWHITE;
}
ImVec4 rl_to_im(const Vector4& v) { return ImVec4{v.x, v.y, v.z, v.w}; }
ImVec4 rl_to_im(const Color& c) {
    f32 r = static_cast<f32>(c.r) / 255.0f;
    f32 g = static_cast<f32>(c.g) / 255.0f;
    f32 b = static_cast<f32>(c.b) / 255.0f;
    f32 a = static_cast<f32>(c.a) / 255.0f;
    return ImVec4{r, g, b, a};
}

void status_text(const char* label, const StatusCode code) {
    im::TextColored(
        rl_to_im(status_color(code)),
        "%s: %s",
        label,
        status_string(code).c_str()
    );
}
void render_simulation_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state) {
    im::Begin("Simulation");
    WorldStepperConfig& stepper = cfg.stepper_cfg;

    if (im::Button("Run/Pause")) {
        toggle(stepper.paused);
    }
    // TODO: hide single step while sim running
    if (!stepper.paused && cfg.step_single) {
        cfg.step_single = false;
        stepper.paused = true;
    } else {
        im::SameLine();
        if (im::Button("Step")) {
            cfg.step_single = true;
            stepper.paused = false;
        }
    }

    im::Checkbox("Realtime", &cfg.realtime);

    im::Text("Time = %10.4f", world.t_sim());
    DHMStime dhms = sec_to_dhms(world.t_sim());
    im::Text("T0+%04d-%02d:%02d:%05.4f", dhms.day, dhms.hour, dhms.minute, dhms.second);
    if (world.is_date_active()) {
        im::Text("JD: %.6lf", jd_to_scalar(world.get_date_jd()));
        im::Text("Date: %s", cal_str(world.get_date_cal()).c_str());
    }

    // TODO: allow set date, enum on date type and copy "now"
    im::Text("dt = %.3f", state.dt);
    im::Text(
        "Effective dt = %.3f (%d steps)",
        state.dt * stepper.ticks * stepper.dt_scale,
        stepper.ticks * stepper.substeps
    );

    // im::Checkbox("Paused", &stepper.paused);

    im::InputInt("Ticks", &stepper.ticks);
    stepper.ticks = std::max(1, stepper.ticks);

    im::InputInt("Substeps", &stepper.substeps);
    stepper.substeps = std::max(1, stepper.substeps);

    im::InputDouble("Time Scale", &stepper.dt_scale);
    if (!finite_pos(stepper.dt_scale)) {
        stepper.dt_scale = 1.0;
    }

    // TODO: add history/snapshot setting, capturing, setting/resetting
    // only allow while paused
    // if using history make vector and dropdown based on time or name?, resetting clears
    // history

    if (!cfg.realtime && stepper.paused) {
        im::Separator();
        // step to a certain world sim time or by delta amount of time using the current
        // dt, ticks, substeps, and dt_scale settings
        im::Text("Offline step"); // TODO: rename this
        im::Checkbox("Relative", &state.relative_step);

        // TODO: allow other date/time types later
        if (state.relative_step) {
            im::TextSL("Step by:");
            im::InputDouble("sec", &state.step_by_delta);
            state.step_to_time = world.t_sim() + state.step_by_delta;
        } else {
            im::TextSL("Step to:");
            im::InputDouble("sec", &state.step_to_time);
            state.step_by_delta = state.step_to_time - world.t_sim();
        }
        if (im::Button("Run")) {
            // TODO: make ui still responsive while this runs

            const f64 advance_per_call = state.dt * stepper.dt_scale * stepper.ticks;
            if (finite_pos(state.step_by_delta) && finite_pos(advance_per_call)) {
                f64 remaining = state.step_by_delta;
                while (remaining > tol12) {
                    const f64 dt_call = std::min(
                        state.dt,
                        remaining / (stepper.dt_scale * stepper.ticks)
                    );
                    WorldStepperStats step_stats
                        = step_world(world, dt_call, stepper, state.wksp);
                    state.stats += step_stats;
                    if (!step_stats.success) break;

                    remaining -= step_stats.dt_sim_advanced;
                    if (remaining < tol12) remaining = 0.0;
                }
            }
        }
    }

    render_scenario_file_ui(world, cfg, state);
    im::End();
}

void render_renderer_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state) {
    // render
    // grids, axes, selected marker, FPS toggle

    im::Begin("Renderer");

    im::Checkbox("Show FPS", &cfg.draw.draw_fps);

    bool changed_lock = im::Checkbox("Lock FPS", &cfg.set_target_fps);
    bool changed_target_fps = false;
    if (cfg.set_target_fps) {
        changed_target_fps = im::InputInt("Target FPS", &cfg.target_fps);
        if (!finite_pos(cfg.target_fps)) {
            cfg.target_fps = 60;
        }
    }
    if (changed_lock || changed_target_fps) {
        SetTargetFPS(cfg.set_target_fps ? cfg.target_fps : 0);
    }

    im::Checkbox("Draw Grids", &cfg.draw.draw_grids);
    if (cfg.draw.draw_grids) {
        if (im::CollapsingHeader("Grids")) {
            im::Checkbox("Show XY Grid", &cfg.draw.draw_grid_xy);
            im::Checkbox("Show ZY Grid", &cfg.draw.draw_grid_zy);
            im::Checkbox("Show XZ Grid", &cfg.draw.draw_grid_xz);
        };
        im::Separator();
    }

    if (im::CollapsingHeader("Axes")) {
        im::Checkbox("Show Inertial Axes", &cfg.draw.draw_inertial_axes);
        im::Checkbox("Show Body Axes", &cfg.draw.draw_body_axes);
        im::Checkbox("Color Axes", &cfg.draw.color_axes);
        im::Separator();
    }

    im::Checkbox("Highlight Selected Body", &cfg.draw.draw_selected_body);
    // im::Checkbox("Draw Labels", &cfg.draw.draw_labels);

    im::End();
}

void render_camera_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state) {
    // camera
    im::Begin("Camera");

    RenderCameraConfig& camera = cfg.camera;
    const char* mode_names[] = {"Locked", "Target", "Origin", "Free"};
    i32 mode_idx = static_cast<i32>(camera.mode);
    if (im::Combo("Camera Mode", &mode_idx, mode_names, 4)) {
        // TODO: update as more camera modes are added
        camera.mode = static_cast<RenderCameraMode>(mode_idx);

        if (camera.mode == RenderCameraMode::target) {
            if (camera.target_id == kInvalidEntityId) {
                cycle_active_id(camera.target_id, world, 1);
            }
            sync_camera_tracking(camera, world);
        }
    }

    if (camera.mode == RenderCameraMode::target) {
        im::Indent();
        const Body* body = world.body(camera.target_id);
        string name = "";
        if (body != nullptr) name = body->name + " ";

        im::Text("Target: %s(ID: %llu)", name.c_str(), camera.target_id);
        if (im::Button("Prev Target")) {
            cycle_active_id(camera.target_id, world, -1);
            sync_camera_tracking(camera, world);
        }
        im::SameLine();
        if (im::Button("Next Target")) {
            cycle_active_id(camera.target_id, world, 1);
            sync_camera_tracking(camera, world);
        }
        im::Unindent();
    }

    im::Text(
        "Position: [%.3f, %.3f, %.3f]",
        camera.position(0),
        camera.position(1),
        camera.position(2)
    );
    im::Text(
        "Target Position: [%.3f, %.3f, %.3f]",
        camera.target(0),
        camera.target(1),
        camera.target(2)
    );

    im::Checkbox("Invert Mouse Wheel", &camera.invert_mousewheel);

    im::SliderFloat("FOV", &camera.fovy, 1.0f, 179.0f);
    camera.fovy = std::clamp(camera.fovy, 1.0f, 179.0f); // TODO: might not be needed

    im::InputFloat("Zoom Rate", &camera.zoom_rate);
    if (!finite_pos(camera.zoom_rate)) {
        RenderCameraConfig default_cfg{};
        camera.zoom_rate = default_cfg.zoom_rate;
    }
    im::InputFloat("Fly Speed", &camera.fly_speed);
    if (!finite_pos(camera.fly_speed)) {
        RenderCameraConfig default_cfg{};
        camera.fly_speed = default_cfg.fly_speed;
    }
    im::InputFloat("Orbit Speed", &camera.orbit_speed);
    if (!finite_pos(camera.orbit_speed)) {
        RenderCameraConfig default_cfg{};
        camera.orbit_speed = default_cfg.orbit_speed;
    }
    im::InputFloat("Pan Speed", &camera.pan_speed);
    if (!finite_pos(camera.pan_speed)) {
        RenderCameraConfig default_cfg{};
        camera.pan_speed = default_cfg.pan_speed;
    }
    if (im::Button("Reset Camera")) {
        RenderCameraConfig default_cfg{};
        camera.invert_mousewheel = default_cfg.invert_mousewheel;
        camera.fovy = default_cfg.fovy;
        camera.zoom_rate = default_cfg.zoom_rate;
        camera.fly_speed = default_cfg.fly_speed;
        camera.orbit_speed = default_cfg.orbit_speed;
        camera.pan_speed = default_cfg.pan_speed;
    }

    if (is_orbit(camera.mode)) {
        vec3f azelr = cart_to_sph<f32>(camera.position, tol9, UAngle::degree);
        vec2f azel = azelr.segment<2>(0);
        if (im::SliderFloat2(
                "[azimuth, elevation]",
                azel,
                vec2f{-179.9f, -89.9f},
                vec2f{179.9f, 89.9f},
                "%.1f"
            )) {
            azelr(0) = azel(0);
            azelr(1) = azel(1);
            camera.position = sph_to_cart(azelr, UAngle::degree);
        }
    }

    im::End();
}

void render_performance_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state) {
    im::Begin("Performance");

    im::Checkbox("Plot Performance", &cfg.draw.plot_performance);
    if (cfg.draw.plot_performance) {
        if (state.frame_time_ms.size() >= state.frame_history_max) {
            state.frame_time_ms.erase(state.frame_time_ms.begin());
        }
        state.frame_time_ms.push_back(state.frame_time * 1000.0f);

        if (state.fps_history.size() >= state.frame_history_max) {
            state.fps_history.erase(state.fps_history.begin());
        }
        state.fps_history.push_back(state.fps);

        // TODO: maybe add input field for plot limits
        if (imp::BeginPlot("Frame Time")) {
            imp::SetupAxes("Sample", "ms");
            imp::SetupAxisLimits(
                ImAxis_X1,
                0.0,
                static_cast<f64>(state.frame_time_ms.size()),
                ImGuiCond_Always
            );
            imp::SetupAxisLimits(ImAxis_Y1, 0.0, 30.0, ImGuiCond_Once);
            // TODO: decide between bars and line
            imp::PlotLine(
                "Frame Time",
                state.frame_time_ms.data(),
                state.frame_time_ms.size()
            );

            imp::EndPlot();
        }

        if (imp::BeginPlot("FPS")) {
            imp::SetupAxes("Sample", "FPS");
            imp::SetupAxisLimits(
                ImAxis_X1,
                0.0,
                static_cast<f64>(state.fps_history.size()),
                ImGuiCond_Always
            );
            imp::SetupAxisLimits(ImAxis_Y1, 0.0, 500.0, ImGuiCond_Once);

            imp::PlotLine("FPS", state.fps_history.data(), state.fps_history.size());
            imp::EndPlot();
        }
    } else {
        if (!state.frame_time_ms.empty()) state.frame_time_ms.clear();
        if (!state.fps_history.empty()) state.fps_history.clear();
    }

    im::End();
}

constexpr ImGuiWindowFlags dockspace_flags
    = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse
      | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove
      | ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus
      | ImGuiWindowFlags_NoBackground | ImGuiWindowFlags_NoDocking;

static void build_default_dock_layout(
    ImGuiID dockspace_id,
    const ImVec2& dockspace_size
) {
    ImGui::DockBuilderRemoveNode(dockspace_id);
    ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
    ImGui::DockBuilderSetNodeSize(dockspace_id, dockspace_size);

    ImGuiID dock_center = dockspace_id;

    ImGuiID dock_left = ImGui::DockBuilderSplitNode(
        dock_center,
        ImGuiDir_Left,
        0.20f, // percentage of main viewport
        nullptr,
        &dock_center
    );

    ImGuiID dock_right = ImGui::DockBuilderSplitNode(
        dock_center,
        ImGuiDir_Right,
        0.33f,
        nullptr,
        &dock_center
    );

    ImGuiID dock_bottom = ImGui::DockBuilderSplitNode(
        dock_center,
        ImGuiDir_Down,
        0.33f,
        nullptr,
        &dock_center
    );
    ImGui::DockBuilderDockWindow("Bodies", dock_left);
    ImGui::DockBuilderDockWindow("Body Statistics", dock_right);
    ImGui::DockBuilderDockWindow("Add Body", dock_right);
    ImGui::DockBuilderDockWindow("Simulation", dock_bottom);
    ImGui::DockBuilderDockWindow("Camera", dock_bottom);
    ImGui::DockBuilderDockWindow("Renderer", dock_bottom);
    ImGui::DockBuilderDockWindow("Performance", dock_bottom);

    ImGui::DockBuilderFinish(dockspace_id);
}

static void render_main_dockspace() {
    const ImGuiViewport* viewport = ImGui::GetMainViewport();

    ImGui::SetNextWindowPos(viewport->WorkPos);
    ImGui::SetNextWindowSize(viewport->WorkSize);
    ImGui::SetNextWindowViewport(viewport->ID);

    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2{0.0f, 0.0f});
    ImGui::Begin("MainDockspaceHost", nullptr, dockspace_flags);
    ImGui::PopStyleVar();

    ImGuiID dockspace_id = ImGui::GetID("MainDockspace");
    if (ImGui::DockBuilderGetNode(dockspace_id) == nullptr) {
        build_default_dock_layout(dockspace_id, viewport->WorkSize);
    }

    ImGui::DockSpace(
        dockspace_id,
        ImVec2{0.0f, 0.0f},
        ImGuiDockNodeFlags_PassthruCentralNode
    );

    im::End();
}
} // namespace render_ui_detail

void render_loop_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state) {
    begin_render_ui_frame();

    render_ui_detail::render_main_dockspace();

    render_ui_detail::render_performance_ui(world, cfg, state);
    render_ui_detail::render_simulation_ui(world, cfg, state);
    render_ui_detail::render_camera_ui(world, cfg, state);
    render_ui_detail::render_renderer_ui(world, cfg, state);
    render_ui_detail::render_body_stats_ui(world, cfg, state);
    render_ui_detail::render_add_body(world, cfg, state);
    render_ui_detail::render_body_lists(world, cfg, state);
    // render_world_history_ui(world, cfg, state);

    end_render_ui_frame();
}
