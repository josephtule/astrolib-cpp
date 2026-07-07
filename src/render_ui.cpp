#include "graphics/render_ui.hpp"
#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/estimation_common.hpp"
#include "core/measurement.hpp"
#include "graphics/render_loop.hpp"
#include "graphics/ui.hpp"
#include "imgui.h"
#include "implot.h"
#include "misc/cpp/imgui_stdlib.h"
#include "util/lightweight_tools.hpp"

namespace im = ImGui;
namespace imp = ImPlot;

static void init_add_body_draft_defaults(
    RenderLoopState& state,
    const World& world,
    BodyType type
);
static void reset_add_body_draft(RenderLoopState& state, BodyType type);
static StatusCode validate_add_body_draft(
    const RenderLoopState& state,
    const World& world
);
static StatusCode validate_mass_properties(const MassProperties& mp, f64 tol = tol12);
static StatusCode validate_body(const Body& body, const World& world);
static StatusCode load_body_edit_draft(
    EntityId id,
    BodyEditDraft& draft,
    const World& world
);
static StatusCode validate_body_edit_draft(
    const BodyEditDraft& draft,
    const World& world
);
static StatusCode apply_body_edit_draft(BodyEditDraft& draft, World& world);
static void cancel_body_edit_draft(BodyEditDraft& draft);

static void render_simulation_ui(
    World& world,
    RenderLoopConfig& cfg,
    RenderLoopState& state
);
static void render_renderer_ui(
    World& world,
    RenderLoopConfig& cfg,
    RenderLoopState& state
);
static void render_camera_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state);
static void render_body_stats_ui(
    World& world,
    RenderLoopConfig& cfg,
    RenderLoopState& state
);
static void render_add_body(World& world, RenderLoopConfig& cfg, RenderLoopState& state);
static void render_performance_ui(
    World& world,
    RenderLoopConfig& cfg,
    RenderLoopState& state
);

static bool render_state_tr_ui(Body& body);
static bool render_state_att_ui(Body& body);
static bool render_state_ui(Body& body);
static bool render_body_propagation_ui(Body& body);
static bool render_stat_state_ui(Station& stat, World& world);
static void render_celestial_stats_ui(Celestial& cel);
static void render_satellite_stats_ui(Satellite& sat);
static void render_mass_properties_ui(MassProperties& mp);
static void render_station_stats_ui(
    Station& stat,
    RenderLoopConfig& cfg,
    RenderLoopState& state,
    World& world
);
static void render_station_draft_ui(Station& stat, RenderLoopState& state, World& world);
static void sync_add_instrument_diag(RenderLoopState& state, i32 dim);
static void render_station_instruments_ui(Station& stat, RenderLoopState& state);

static Body* draft_body_ptr(BodyEditDraft& draft) {
    switch (draft.edit_body_type) {
    case BodyType::unknown: return nullptr;
    case BodyType::celestial: return &draft.edit_celestial;
    case BodyType::satellite: return &draft.edit_satellite;
    case BodyType::station: return &draft.edit_station;
    }

    return nullptr;
}

static bool render_state_tr_ui(Body& body) {
    return im::InputDouble3("r", body.x_tr.r) || im::InputDouble3("v", body.x_tr.v);
}

static bool render_state_att_ui(Body& body) {
    return im::InputDouble4("q", body.x_att.q) || im::InputDouble3("w", body.x_att.w);
}

static bool render_state_ui(Body& body) {
    return render_state_tr_ui(body) || render_state_att_ui(body);
}

static void init_add_body_draft_defaults(
    RenderLoopState& state,
    const World& world,
    BodyType type
) {
    switch (type) {
    case BodyType::unknown: break;
    case BodyType::celestial: {
        Celestial& cel = state.temp_celestial;
        if (cel.name.empty()) cel.name = "New Celestial";
        if (cel.semimajor_axis == 0.0 && cel.semiminor_axis == 0.0
            && cel.mean_radius == 0.0) {
            cel.semimajor_axis = 1000.0;
            cel.semiminor_axis = 1000.0;
            cel.mean_radius = mean_from_semiaxes(cel.semimajor_axis, cel.semiminor_axis);
            cel.eccentricity = ecc_from_semiaxes(cel.semimajor_axis, cel.semiminor_axis);
            cel.flattening = flat_from_semiaxes(cel.semimajor_axis, cel.semiminor_axis);
            if (cel.ref_radius == 0.0) cel.ref_radius = cel.semimajor_axis;
        }
    } break;
    case BodyType::satellite: {
        Satellite& sat = state.temp_satellite;
        if (sat.name.empty()) sat.name = "New Satellite";
        if (sat.propagate_att && !sat.mass_properties.active) {
            sat.propagate_att = false;
        }
    } break;
    case BodyType::station: {
        Station& stat = state.temp_station;
        if (stat.name.empty()) stat.name = "New Station";
        if (stat.anchored && world.celestial(stat.anchor_id) == nullptr) {
            EntityId anchor_id = first_celestial_id(world);
            if (anchor_id == kInvalidEntityId) {
                stat.anchored = false;
            } else {
                stat.anchor_id = anchor_id;
            }
        }
    } break;
    }
}

static void render_simulation_ui(
    World& world,
    RenderLoopConfig& cfg,
    RenderLoopState& state
) {
    im::Begin("Simulation");
    WorldStepperConfig& stepper = cfg.stepper_cfg;

    if (im::Button("Play/Pause")) toggle(stepper.paused);

    im::Checkbox("Realtime", &cfg.realtime);

    im::Text("Time = %.3f", world.t_sim());
    im::Text("dt = %.3f", state.dt);
    im::Text("Effective dt = %.3f", state.dt * stepper.ticks * stepper.dt_scale);

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

    if (im::CollapsingHeader("Scenario")) {
        im::Indent();
        im::Text("Filepath:");
        im::SameLine();
        string path;
        im::InputText("##filepath", &path);
        if (im::Button("Save")) {
            // TODO:
            // TODO: if file already exists, pop up warning
        }
        im::SameLine();
        if (im::Button("Load")) {
            // TEMP: use raylib file explorer or imgui-filebrowser
            // TODO: if file doesn't exist, pop up warning
        }
        im::Unindent();
    }

    if (im::Button("Add Body")) {
        state.add_body = true;
        if (state.add_body_type == BodyType::unknown) {
            state.add_body_type = BodyType::celestial;
        }
        init_add_body_draft_defaults(state, world, state.add_body_type);
        state.add_body_status = StatusCode::ok;
    }

    im::End();
}

static void render_renderer_ui(
    World& world,
    RenderLoopConfig& cfg,
    RenderLoopState& state
) {
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
        im::Indent();
        if (im::CollapsingHeader("Grids")) {
            im::Checkbox("Show XY Grid", &cfg.draw.draw_grid_xy);
            im::Checkbox("Show ZY Grid", &cfg.draw.draw_grid_zy);
            im::Checkbox("Show XZ Grid", &cfg.draw.draw_grid_xz);
        };
        im::Unindent();
    }

    if (im::CollapsingHeader("Axes")) {
        im::Indent();
        im::Checkbox("Show Inertial Axes", &cfg.draw.draw_inertial_axes);
        im::Checkbox("Show Body Axes", &cfg.draw.draw_body_axes);
        im::Checkbox("Color Axes", &cfg.draw.color_axes);
        im::Unindent();
    }

    im::Checkbox("Highlight Selected Body", &cfg.draw.draw_selected_body);
    // im::Checkbox("Draw Labels", &cfg.draw.draw_labels);

    im::End();
}

static void render_camera_ui(
    World& world,
    RenderLoopConfig& cfg,
    RenderLoopState& state
) {
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

    // TODO: add dropdown or +- or separate submenu
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

    im::End();
}

static void render_body_stats_ui(
    World& world,
    RenderLoopConfig& cfg,
    RenderLoopState& state
) {
    im::Begin("Body Statistics");
    im::Checkbox("Display Statistics", &cfg.display_body_stats);
    if (cfg.display_body_stats) {
        if (cfg.body_stats_id == kInvalidEntityId)
            cycle_active_id(cfg.body_stats_id, world, 1);

        Body* body = world.body(cfg.body_stats_id);
        if (body == nullptr) {
            im::End();
            return;
        }

        bool active = world.is_active(body->id);

        im::Text(
            "%s: %s (id: %llu, active: %s)",
            body_type_str(body->body_type).c_str(),
            body->name.c_str(),
            body->id,
            bool_str(active).c_str()
        );
        im::Indent();
        if (im::Button("Prev Body")) {
            svec<EntityId> ids = world.all_entity_ids();
            cycle_id(cfg.body_stats_id, ids, -1);
            body = world.body(cfg.body_stats_id);
            if (body == nullptr) {
                im::End();
                return;
            }
        }
        im::SameLine();
        if (im::Button("Next Body")) {
            svec<EntityId> ids = world.all_entity_ids();
            cycle_id(cfg.body_stats_id, ids, 1);
            body = world.body(cfg.body_stats_id);
            if (body == nullptr) {
                im::End();
                return;
            }
        }
        // TODO: add dropdown + body type filtering?
        im::Unindent();

        bool anchored = false;
        Station* stat = nullptr;
        if (body->body_type == BodyType::station) {
            stat = world.station(cfg.body_stats_id);
            if (stat == nullptr) {
                im::End();
                return;
            }
            anchored = stat->anchored;
        }

        // TODO: add ability to select/copy text but not edit
        bool editable = cfg.edit_body_stats && cfg.stepper_cfg.paused;
        // im::Checkbox("Edit Body Stats", &cfg.edit_body_stats);
        if (!editable && im::Button("Edit")) {
            state.draft.edit_body_status
                = load_body_edit_draft(body->id, state.draft, world);
            if (state.draft.edit_body_status == StatusCode::ok) {
                cfg.edit_body_stats = true;
                editable = cfg.stepper_cfg.paused;
            }
        }
        if (editable) {
            if (im::Button("Save")) {
                state.draft.edit_body_status = apply_body_edit_draft(state.draft, world);
                if (state.draft.edit_body_status == StatusCode::ok) {
                    cfg.edit_body_stats = false;
                    cancel_body_edit_draft(state.draft);
                }
            }
            im::SameLine();
            if (im::Button("Cancel")) {
                cancel_body_edit_draft(state.draft);
                cfg.edit_body_stats = false;
            }
        }
        if (state.draft.edit_body_id != kInvalidEntityId
            && state.draft.edit_body_status != StatusCode::ok) {
            im::Text(
                "Edit Body Status: %s",
                status_string(state.draft.edit_body_status).c_str()
            );
        }

        if (!editable) im::BeginDisabled();
        im::TextSL("Name:");
        Body* edit_body = editable ? draft_body_ptr(state.draft) : body;
        if (edit_body != nullptr) {
            im::InputText("##Name", &edit_body->name);
        }

        switch (body->body_type) {
        case BodyType::unknown: break;
        case BodyType::celestial: {
            Celestial* cel = editable ? &state.draft.edit_celestial
                                      : world.celestial(cfg.body_stats_id);
            if (cel == nullptr) {
                im::End();
                return;
            }
            render_celestial_stats_ui(*cel);
        } break;
        case BodyType::satellite: {
            Satellite* sat = editable ? &state.draft.edit_satellite
                                      : world.satellite(cfg.body_stats_id);
            if (sat == nullptr) {
                im::End();
                return;
            }
            render_satellite_stats_ui(*sat);
        } break;
        case BodyType::station: {
            Station* stat
                = editable ? &state.draft.edit_station : world.station(cfg.body_stats_id);
            if (stat == nullptr) {
                im::End();
                return;
            }
            render_station_stats_ui(*stat, cfg, state, world);
        } break;
        }

        // im::Checkbox("Emits Gravity", &body->emits_gravity);
        // im::Checkbox("Emits Radiation", &body->emits_radiation);
        // im::Checkbox("Has Atmosphere", &body->has_atmosphere);

        if (!editable) im::EndDisabled();

        if (active && im::Button("Set as Camera Target")) {
            cfg.camera.target_id = cfg.body_stats_id;
            sync_camera_tracking(cfg.camera, world);
        }
        if (active && im::Button("Make Inactive")) {
            world.make_inactive(cfg.body_stats_id);
            state.wksp.dirty = true;
            if (cfg.camera.target_id == body->id) {
                cycle_active_id(cfg.camera.target_id, world, 1);
                sync_camera_tracking(cfg.camera, world);
            }
        }
        if (!active && im::Button("Make Active")) {
            world.make_active(cfg.body_stats_id);
            state.wksp.dirty = true;
        }
    }

    im::End();
}

static void render_performance_ui(
    World& world,
    RenderLoopConfig& cfg,
    RenderLoopState& state
) {
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

static void render_add_body(World& world, RenderLoopConfig& cfg, RenderLoopState& state) {
    im::Begin("Add Body");
    if (cfg.stepper_cfg.paused && state.add_body) {
        ImGui::SetWindowCollapsed(false);

        const char* type_names[] = {"Unknown", "Celestial", "Satellite", "Station"};
        i32 type_idx = static_cast<i32>(state.add_body_type);
        if (im::Combo("Body Type", &type_idx, type_names, 4)) {
            state.add_body_type = static_cast<BodyType>(type_idx);
            init_add_body_draft_defaults(state, world, state.add_body_type);
            state.add_body_status = StatusCode::ok;
        }

        switch (state.add_body_type) {
        case BodyType::unknown: break;
        case BodyType::celestial: render_celestial_stats_ui(state.temp_celestial); break;
        case BodyType::satellite: render_satellite_stats_ui(state.temp_satellite); break;
        case BodyType::station:
            render_station_draft_ui(state.temp_station, state, world);
            break;
        }

        // TODO: add copy from existing bodies, maybe allow cross types for states, etc.

        if (state.add_body_type != BodyType::unknown) {
            if (im::Button("Save")) {
                EntityId id = kInvalidEntityId;
                state.add_body_status = validate_add_body_draft(state, world);

                if (state.add_body_status == StatusCode::ok) {
                    switch (state.add_body_type) {
                    case BodyType::unknown: break;
                    case BodyType::celestial: {
                        auto cel = std::make_unique<Celestial>(state.temp_celestial);
                        id = world.insert_celestial(std::move(cel));
                        state.temp_celestial = Celestial{};
                    } break;
                    case BodyType::satellite: {
                        auto sat = std::make_unique<Satellite>(state.temp_satellite);
                        id = world.insert_satellite(std::move(sat));
                        state.temp_satellite = Satellite{};
                    } break;
                    case BodyType::station: {
                        auto stat = std::make_unique<Station>(state.temp_station);
                        id = world.insert_station(std::move(stat));
                        state.temp_station = Station{};
                    } break;
                    }

                    if (id != kInvalidEntityId) {
                        state.wksp.dirty = true;
                        cfg.body_stats_id = id;
                        cfg.camera.target_id = id;
                    } else {
                        state.add_body_status = StatusCode::body_not_found;
                    }

                    if (state.add_body_status == StatusCode::ok) {
                        state.add_body = false;
                        im::SetWindowCollapsed(true);
                    }
                }
            }

            im::SameLine();
            if (im::Button("Reset Draft")) {
                reset_add_body_draft(state, state.add_body_type);
                init_add_body_draft_defaults(state, world, state.add_body_type);
            }

            im::SameLine();
            if (im::Button("Cancel")) {
                state.add_body = false;
                state.add_body_status = StatusCode::ok;
                im::SetWindowCollapsed(true);
            }

            if (state.add_body_status != StatusCode::ok) {
                im::Text(
                    "Add Body Status: %s",
                    status_string(state.add_body_status).c_str()
                );
            }
        }
    } else {
        im::SetWindowCollapsed(true);
    }
    im::End();
}

static bool render_body_propagation_ui(Body& body) {
    return im::Checkbox("Propagate Translation", &body.propagate_tr)
           || im::Checkbox("Propagate Attitude", &body.propagate_att);
}

static bool render_stat_state_ui(Station& stat, World& world) {
    vec3d r = world.stat_r_inertial(stat.id);
    vec4d q = world.stat_q_inertial(stat.id);
    return im::InputDouble3("r", r) || im::InputDouble4("q", q);
}

static void render_celestial_stats_ui(Celestial& cel) {
    render_state_ui(cel);

    if (im::CollapsingHeader("Celestial Model", ImGuiTreeNodeFlags_DefaultOpen)) {
        bool edit_a = im::InputDouble("Semimajor Axis", &cel.semimajor_axis);
        bool edit_b = im::InputDouble("Semiminor Axis", &cel.semiminor_axis);
        if (edit_a || edit_b) {
            cel.mean_radius = mean_from_semiaxes(cel.semimajor_axis, cel.semiminor_axis);
            cel.eccentricity = ecc_from_semiaxes(cel.semimajor_axis, cel.semiminor_axis);
            cel.flattening = flat_from_semiaxes(cel.semimajor_axis, cel.semiminor_axis);
        }
        // TODO: add recomputation of the below
        im::BeginDisabled();
        im::InputDouble("Mean Radius", &cel.mean_radius);
        im::InputDouble("Eccentricity", &cel.eccentricity);
        im::InputDouble("Flattening", &cel.flattening);
        im::EndDisabled();
    }

    if (im::CollapsingHeader("Gravity Model", ImGuiTreeNodeFlags_DefaultOpen)) {
        im::Indent();
        i32 model_idx = static_cast<i32>(cel.gravity_model);
        const char* model_names[] = {"Pointmass", "Zonal", "Spherical Harmonics"};
        if (im::Combo("##Gravity Model", &model_idx, model_names, 3)) {
            // TODO: update as new models are added
            cel.gravity_model = static_cast<GravityModel>(model_idx);
            // TODO: set to fall back to zonal or pointmass if coefs not available
            // (depending on J, C, S active bools?)
        }
        im::InputDouble("Mu", &cel.mu);
        switch (cel.gravity_model) {
        case GravityModel::pointmass: break;
        case GravityModel::zonal: {
            im::InputDouble("Reference Radius", &cel.ref_radius);
            im::InputInt("Degree", &cel.degree);
            im::InputDouble7("Zonal Coefficients (J)", cel.J);
        } break;
        case GravityModel::spherical_harmonics: {
            im::InputDouble("Reference Radius", &cel.ref_radius);
            im::InputInt("Degree", &cel.degree);
            im::InputInt("Order", &cel.order);
            bool csactive = cel.C.cols() > 0 && cel.C.rows() > 0 && cel.S.cols() > 0
                            && cel.S.rows() > 0;
            im::Text("C/S Matrices Active: %s", bool_str(csactive).c_str());
        } break;
        }
        im::Unindent();
    }

    if (im::CollapsingHeader("Attitude Model")) {
        i32 att_idx = static_cast<i32>(cel.attitude_model);
        const char* att_name[] = {"Fixed", "Simple Spin", "Provider"};
        if (im::Combo("##Attitude Model", &att_idx, att_name, 3)) {
            if (static_cast<CelestialAttitudeModel>(att_idx)
                == CelestialAttitudeModel::provider) {
                att_idx = static_cast<i32>(CelestialAttitudeModel::fixed);
            } // TEMP: no attitude providers yet, fallback to fixed
            cel.attitude_model = static_cast<CelestialAttitudeModel>(att_idx);
        }
        if (cel.attitude_model == CelestialAttitudeModel::simple_spin) {
            im::BeginDisabled();
            f64 w_mag = cel.x_att.w.norm();
            im::InputDouble("Spin Rate", &w_mag);
            im::InputDouble3("Spin Axis", cel.x_att.w);
            im::EndDisabled();
        }
    }

    render_body_propagation_ui(cel);
    // TODO: allow change of coefs (dropdown of providers)
}

static void render_mass_properties_ui(MassProperties& mp) {
    if (im::CollapsingHeader("Mass Properties", ImGuiTreeNodeFlags_DefaultOpen)) {
        im::InputDouble("Mass", &mp.mass);

        if (im::Checkbox("Principal Axes", &mp.principal_axes)) {
            if (mp.principal_axes) {
                vec3d I_diag = mp.I.diagonal();
                mp.I = I_diag.asDiagonal();
                mp.I_inv = mp.I.inverse();
            }
        }
        if (mp.principal_axes) {
            vec3d I_diag = mp.I.diagonal();
            if (im::InputDouble3("Inertia Tensor", I_diag)) {
                mp.I = I_diag.asDiagonal();
                mp.I_inv = mp.I.inverse();
            }
        } else {
            if (im::InputDouble3x3("Inertia Tensor", mp.I)) {
                mp.I_inv = mp.I.inverse();
            }
        }

        im::Indent();
        im::InputDouble3("Offset", mp.offset_body);
        if (im::Button("Recompute Inertia Tensor about COM")) {
            mp.I = inertia_PAT(mp.I, mp.mass, mp.offset_body);
            if (!mp.I.isDiagonal()) mp.principal_axes = false;
            mp.offset_body = vec3d0;
            mp.I_inv = mp.I.inverse();
        }
        im::Unindent();

        im::Text("Inertia Matrix Active: %s", bool_str(mp.active).c_str());
        im::Indent();
        if (!finite_mat(mp.I_inv)) {
            im::Text("Inertia Tensor Cannot Be Inverted");
            mp.active = false;
        } else {
            mp.active = true;
        }
        im::Unindent();
    }
}

static void render_satellite_stats_ui(Satellite& sat) {
    render_state_ui(sat);
    render_mass_properties_ui(sat.mass_properties);
    render_body_propagation_ui(sat);
}

static void sync_add_instrument_diag(RenderLoopState& state, i32 dim) {
    if (dim <= 0) {
        state.add_instrument_R_diag.resize(0);
        return;
    }

    if (state.add_instrument_R_diag.size() == dim) return;

    vecXd diag = vecXd::Ones(dim);
    i32 n_copy = std::min<i32>(dim, static_cast<i32>(state.add_instrument_R_diag.size()));
    for (i32 i = 0; i < n_copy; ++i) {
        diag(i) = state.add_instrument_R_diag(i);
    }
    state.add_instrument_R_diag = diag;
}

static void render_station_instruments_ui(Station& stat, RenderLoopState& state) {
    if (im::CollapsingHeader("Instruments")) {
        im::Indent();
        bool enabled_changed = false;
        for (auto& [id, instr] : stat.instruments) {
            im::PushID(&instr.name);
            im::TextSL("Name:");
            im::InputText("##name", &instr.name);
            im::Text("ID: %u", instr.id);
            im::Text("Observation Type: %s", observation_type_str(instr.type).c_str());
            im::InputMatXd("Covariance", instr.R, "%.3E");
            if (im::Checkbox("Enabled", &instr.enabled)) {
                enabled_changed = true;
            }
            im::PopID();
            im::Separator();
        }

        if (enabled_changed) {
            stat.enabled_instrument_ids = enabled_station_instrument_ids(stat);
        }

        if (im::CollapsingHeader("Add Instrument")) {
            const char* type_names[]
                = {"Right-Ascension + Declination",
                   "Azimuth + Elevation",
                   "Range",
                   "Range-Rate",
                   "Position",
                   "Position + Velocity",
                   "Relative Position",
                   "Relative Position + Velocity"};

            im::Combo("Type", &state.add_instrument_type, type_names, 8);
            im::TextSL("Name:");
            im::InputText("##new instrument name", &state.add_instrument_name);

            ObservationType type
                = static_cast<ObservationType>(state.add_instrument_type);
            i32 dim = measurement_dim(type);
            sync_add_instrument_diag(state, dim);

            for (i32 i = 0; i < state.add_instrument_R_diag.size(); ++i) {
                im::PushID(i);
                f64 variance = state.add_instrument_R_diag(i);
                string label = "R(" + std::to_string(i) + "," + std::to_string(i) + ")";
                im::InputDouble(label.c_str(), &variance, 0.0, 0.0, "%.3E");
                if (!finite_pos(variance)) variance = 1.0;
                state.add_instrument_R_diag(i) = variance;
                im::PopID();
            }

            if (im::Button("Add")) {
                if (dim <= 0) {
                    state.add_instrument_status = StatusCode::invalid_input;
                } else {
                    StationInstrument instrument;
                    instrument.name = state.add_instrument_name;
                    instrument.type = type;
                    instrument.enabled = true;
                    instrument.R = state.add_instrument_R_diag.asDiagonal();
                    state.add_instrument_status
                        = add_station_instrument(stat, instrument);
                    if (state.add_instrument_status == StatusCode::ok) {
                        state.add_instrument_name = "New Instrument";
                    }
                }
            }
            // im::SameLine();
            // if (im::Button("Reset")) {
            // }
            if (state.add_instrument_status != StatusCode::ok) {
                im::Text(
                    "Add Instrument Status: %s",
                    status_string(state.add_instrument_status).c_str()
                );
            }
        }
        im::Unindent();
    }
}

static void render_station_stats_ui(
    Station& stat,
    RenderLoopConfig& cfg,
    RenderLoopState& state,
    World& world
) {
    im::Checkbox("Anchored", &stat.anchored);

    if (stat.anchored) {
        Celestial* cel = world.celestial(stat.anchor_id);
        string anchor_name;
        if (cel == nullptr) {
            anchor_name = "unknown";
        } else {
            anchor_name = cel->name;
        }
        render_stat_state_ui(stat, world);
        im::Text("Anchor: %s (ID: %llu)", anchor_name.c_str(), stat.anchor_id);
        im::Indent();
        bool prev_anchor = im::Button("Prev Anchor");
        bool next_anchor = im::Button("Next Anchor");
        if (prev_anchor || next_anchor) {
            svec<EntityId> cel_ids = world.active_celestial_ids();
            i32 step = next_anchor ? 1 : -1;
            EntityId trial_id = stat.anchor_id;
            cycle_id(trial_id, cel_ids, step);
            Celestial* trial_cel = world.celestial(trial_id);
            if (trial_cel != nullptr) {
                stat.anchor_id = trial_id;
                stat.anchored
                    = world.set_stat_anchor_detic(stat.id, stat.anchor_id, stat.llh_BCBF);
            }
        }
        im::Unindent();
        if (im::Button("Switch to Anchor")) {
            cfg.body_stats_id = stat.anchor_id;
        }
        if (cel == nullptr) {
            im::Text("Invalid Anchor, Cannot Convert Station Geometry");
        } else {
            if (im::InputDouble3("r_body (BCBF)", stat.r_body_BCBF)) {
                stat.llh_BCBF = bcbf_to_detic(stat.r_body_BCBF, *cel);
            }
            if (im::InputDouble3("Planetodetic LLH", stat.llh_BCBF)) {
                stat.r_body_BCBF = detic_to_bcbf(stat.llh_BCBF, *cel);
            }
        }

    } else {
        render_state_ui(stat);
        render_mass_properties_ui(stat.mass_properties);
    }

    render_station_instruments_ui(stat, state);
}

static void render_station_draft_ui(Station& stat, RenderLoopState& state, World& world) {
    im::Checkbox("Anchored", &stat.anchored);

    if (stat.anchored) {
        Celestial* cel = world.celestial(stat.anchor_id);
        string anchor_name = cel == nullptr ? "unknown" : cel->name;

        im::Text("Anchor: %s (ID: %llu)", anchor_name.c_str(), stat.anchor_id);
        im::Indent();
        bool prev_anchor = im::Button("Prev Anchor");
        bool next_anchor = im::Button("Next Anchor");
        if (prev_anchor || next_anchor) {
            svec<EntityId> cel_ids = world.active_celestial_ids();
            i32 step = next_anchor ? 1 : -1;
            EntityId trial_id = stat.anchor_id;
            cycle_id(trial_id, cel_ids, step);
            Celestial* trial_cel = world.celestial(trial_id);
            if (trial_cel != nullptr) {
                stat.anchor_id = trial_id;
                cel = trial_cel;
            }
        }
        im::Unindent();

        if (cel == nullptr) {
            im::Text("Select a Valid Anchor to Edit LLH/r_body");
        } else {
            if (im::InputDouble3("r_body (BCBF)", stat.r_body_BCBF)) {
                stat.llh_BCBF = bcbf_to_detic(stat.r_body_BCBF, *cel);
            }
            if (im::InputDouble3("Planetodetic LLH", stat.llh_BCBF)) {
                stat.r_body_BCBF = detic_to_bcbf(stat.llh_BCBF, *cel);
            }
        }
    } else {
        render_state_ui(stat);
        render_mass_properties_ui(stat.mass_properties);
    }

    render_body_propagation_ui(stat);
    render_station_instruments_ui(stat, state);
}

static void reset_add_body_draft(RenderLoopState& state, BodyType type) {
    switch (type) {
    case BodyType::unknown: break;
    case BodyType::celestial: state.temp_celestial = Celestial{}; break;
    case BodyType::satellite: state.temp_satellite = Satellite{}; break;
    case BodyType::station: state.temp_station = Station{}; break;
    }
    state.add_body_status = StatusCode::ok;
}

static StatusCode validate_add_body_draft(
    const RenderLoopState& state,
    const World& world
) {
    switch (state.add_body_type) {
    case BodyType::unknown: return StatusCode::invalid_input;
    case BodyType::celestial: {
        const Celestial& cel = state.temp_celestial;
        if (!finite_state_tr(cel.x_tr) || !finite_state_att(cel.x_att)) {
            return StatusCode::invalid_input;
        }
        if (!finite_nonneg(cel.mu)) return StatusCode::invalid_input;
        if (!finite_nonneg(cel.semimajor_axis) || !finite_nonneg(cel.semiminor_axis)) {
            return StatusCode::invalid_input;
        }
        if (cel.semimajor_axis == 0.0 && cel.semiminor_axis == 0.0) {
            return StatusCode::invalid_input;
        }
        if (cel.semimajor_axis > 0.0 && cel.semiminor_axis > cel.semimajor_axis) {
            return StatusCode::invalid_input;
        }
        if (cel.ref_radius < 0.0 || !std::isfinite(cel.ref_radius)) {
            return StatusCode::invalid_input;
        }
    } break;
    case BodyType::satellite: {
        const Satellite& sat = state.temp_satellite;
        if (!finite_state_tr(sat.x_tr) || !finite_state_att(sat.x_att)) {
            return StatusCode::invalid_input;
        }
        if (sat.propagate_att) {
            if (!finite_pos(sat.mass_properties.mass)) return StatusCode::invalid_input;
            if (!finite_mat(sat.mass_properties.I)
                || !finite_mat(sat.mass_properties.I_inv)) {
                return StatusCode::invalid_input;
            }
            if (!sat.mass_properties.active) return StatusCode::invalid_input;
        }
    } break;
    case BodyType::station: {
        const Station& stat = state.temp_station;
        if (stat.anchored) {
            if (world.celestial(stat.anchor_id) == nullptr)
                return StatusCode::body_not_found;
            if (!finite_vec(stat.r_body_BCBF) || !finite_vec(stat.llh_BCBF)) {
                return StatusCode::invalid_input;
            }
        } else {
            if (!finite_state_tr(stat.x_tr) || !finite_state_att(stat.x_att)) {
                return StatusCode::invalid_input;
            }
            if (stat.propagate_att && !stat.mass_properties.active) {
                return StatusCode::invalid_input;
            }
        }
    } break;
    }

    return StatusCode::ok;
}

void render_body_lists(World& world, RenderLoopConfig& cfg, RenderLoopState& state) {
    im::Begin("Bodies");
    // TODO: default this to closed

    if (ImGui::CollapsingHeader("Celestials")) {
        svec<EntityId> ids = world.all_celestial_ids();
        for (const EntityId id : ids) {
            Celestial* cel = world.celestial(id);
            if (cel == nullptr) continue;
            im::PushID(&cel->name);
            im::Text(
                "Name: %s (id: %llu, active: %s)",
                cel->name.c_str(),
                cel->id,
                active_str(world.is_active(id)).c_str()
            );
            im::PopID();
        }
    }

    if (ImGui::CollapsingHeader("Satellites")) {
        svec<EntityId> ids = world.all_satellite_ids();
        for (const EntityId id : ids) {
            Satellite* sat = world.satellite(id);
            if (sat == nullptr) continue;
            im::PushID(&sat->name);
            im::Text(
                "Name: %s (id: %llu, active: %s)",
                sat->name.c_str(),
                sat->id,
                active_str(world.is_active(id)).c_str()
            );
            im::PopID();
        }
    }

    if (ImGui::CollapsingHeader("Stations")) {
        svec<EntityId> ids = world.all_station_ids();
        for (const EntityId id : ids) {
            Station* stat = world.station(id);
            if (stat == nullptr) continue;
            im::PushID(&stat->name);
            im::Text(
                "Name: %s (id: %llu, active: %s)",
                stat->name.c_str(),
                stat->id,
                active_str(world.is_active(id)).c_str()
            );
            im::PopID();
        }
    }

    im::End();
}

void render_loop_ui(World& world, RenderLoopConfig& cfg, RenderLoopState& state) {
    begin_render_ui_frame();

    // TODO: make these panels dockable
    render_performance_ui(world, cfg, state);
    render_simulation_ui(world, cfg, state);
    render_camera_ui(world, cfg, state);
    render_renderer_ui(world, cfg, state);
    render_body_stats_ui(world, cfg, state);
    render_add_body(world, cfg, state);
    render_body_lists(world, cfg, state);
    // render_world_history_ui(world, cfg, state);

    end_render_ui_frame();
}

static StatusCode validate_mass_properties(const MassProperties& mp, f64 tol) {
    if (mp.mass == 0.0) return StatusCode::ok;
    if (!finite_pos(mp.mass)) return StatusCode::invalid_input;
    if (!finite_mat(mp.I)) return StatusCode::invalid_input;

    for (i32 i = 0; i < 3; ++i) {
        if (!finite_pos(mp.I(i, i))) return StatusCode::invalid_input;
    }
    if (!finite_nonzero(mp.I.determinant(), tol)) return StatusCode::invalid_input;

    return StatusCode::ok;
}

static StatusCode validate_body(const Body& body, const World& world) {
    switch (body.body_type) {
    case BodyType::unknown: return StatusCode::invalid_input;
    case BodyType::celestial: {
        const Celestial* cel = dynamic_cast<const Celestial*>(&body);
        if (cel == nullptr) return StatusCode::invalid_input;
        if (!finite_state_tr(cel->x_tr) || !finite_state_att(cel->x_att)) {
            return StatusCode::invalid_input;
        }
        if (!finite_nonneg(cel->mu)) return StatusCode::invalid_input;
        if (!finite_nonneg(cel->semimajor_axis) || !finite_nonneg(cel->semiminor_axis)) {
            return StatusCode::invalid_input;
        }
        if (cel->ref_radius < 0.0 || !std::isfinite(cel->ref_radius)) {
            return StatusCode::invalid_input;
        }
    } break;
    case BodyType::satellite: {
        const Satellite* sat = dynamic_cast<const Satellite*>(&body);
        if (sat == nullptr) return StatusCode::invalid_input;
        if (!finite_state_tr(sat->x_tr) || !finite_state_att(sat->x_att)) {
            return StatusCode::invalid_input;
        }
        if (sat->propagate_att) {
            // TODO: use validate mass properties here
            if (!finite_pos(sat->mass_properties.mass)) return StatusCode::invalid_input;
            if (!finite_mat(sat->mass_properties.I)
                || !finite_mat(sat->mass_properties.I_inv)) {
                return StatusCode::invalid_input;
            }
            if (!sat->mass_properties.active) return StatusCode::invalid_input;
        }
    } break;
    case BodyType::station: {
        const Station* stat = dynamic_cast<const Station*>(&body);
        if (stat == nullptr) return StatusCode::invalid_input;
        if (stat->anchored) {
            if (world.celestial(stat->anchor_id) == nullptr)
                return StatusCode::body_not_found;
            if (!finite_vec(stat->r_body_BCBF) || !finite_vec(stat->llh_BCBF)) {
                return StatusCode::invalid_input;
            }
        } else {
            // TODO: use validate mass properties here
            if (!finite_state_tr(stat->x_tr) || !finite_state_att(stat->x_att)) {
                return StatusCode::invalid_input;
            }
            if (stat->propagate_att && !stat->mass_properties.active) {
                return StatusCode::invalid_input;
            }
        }
    } break;
    }

    return StatusCode::ok;
}

static StatusCode load_body_edit_draft(
    EntityId id,
    BodyEditDraft& draft,
    const World& world
) {
    draft.edit_body_id = id;
    const Body* body = world.body(id);
    if (body == nullptr) return StatusCode::body_not_found;
    draft.edit_body_type = body->body_type;

    switch (draft.edit_body_type) {
    case BodyType::unknown: return StatusCode::body_not_found;
    case BodyType::celestial: {
        const Celestial* cel = world.celestial(draft.edit_body_id);
        if (cel == nullptr) return StatusCode::body_not_found;
        draft.edit_celestial = *cel;
    } break;
    case BodyType::satellite: {
        const Satellite* sat = world.satellite(draft.edit_body_id);
        if (sat == nullptr) return StatusCode::body_not_found;
        draft.edit_satellite = *sat;
    } break;
    case BodyType::station: {
        const Station* stat = world.station(draft.edit_body_id);
        if (stat == nullptr) return StatusCode::body_not_found;
        draft.edit_station = *stat;
    } break;
    }

    return StatusCode::ok;
}
static StatusCode validate_body_edit_draft(
    const BodyEditDraft& draft,
    const World& world
) {
    StatusCode status = StatusCode::ok;
    switch (draft.edit_body_type) {
    case BodyType::unknown: return StatusCode::invalid_input;
    case BodyType::celestial: {
        status = validate_body(draft.edit_celestial, world);
    } break;
    case BodyType::satellite: {
        status = validate_body(draft.edit_satellite, world);
    } break;
    case BodyType::station: {
        status = validate_body(draft.edit_station, world);
    } break;
    }
    return status;
}
static StatusCode apply_body_edit_draft(BodyEditDraft& draft, World& world) {
    StatusCode status = validate_body_edit_draft(draft, world);
    if (status != StatusCode::ok) return status;

    switch (draft.edit_body_type) {
    case BodyType::unknown: return StatusCode::invalid_input;
    case BodyType::celestial: {
        Celestial* cel = world.celestial(draft.edit_body_id);
        if (cel == nullptr) return StatusCode::body_not_found;
        *cel = draft.edit_celestial;
        return StatusCode::ok;
    } break;
    case BodyType::satellite: {
        Satellite* sat = world.satellite(draft.edit_body_id);
        if (sat == nullptr) return StatusCode::body_not_found;
        *sat = draft.edit_satellite;
    } break;
    case BodyType::station: {
        Station* stat = world.station(draft.edit_body_id);
        if (stat == nullptr) return StatusCode::body_not_found;
        *stat = draft.edit_station;
    } break;
    }

    return StatusCode::ok;
}
static void cancel_body_edit_draft(BodyEditDraft& draft) { draft = BodyEditDraft{}; }
