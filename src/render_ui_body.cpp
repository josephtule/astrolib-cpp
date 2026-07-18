// Copyright 2025-2026 Joseph Tu Le
// SPDX-License-Identifier: Apache-2.0

#include "render_ui_internal.hpp"

#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/estimation_common.hpp"
#include "core/measurement.hpp"
#include "core/scenario_io.hpp"
#include "core/time.hpp"
#include "core/transform.hpp"
#include "core/world.hpp"
#include "graphics/camera.hpp"
#include "graphics/render_loop.hpp"
#include "graphics/ui.hpp"

#include "ImGuiFD.h"
#include "imgui.h"
#include "implot.h"
#include "misc/cpp/imgui_stdlib.h"

#include "raylib.h"
#include "util/lightweight_tools.hpp"
#include "util/math.hpp"
#include "util/units.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <memory>
#include <string>

namespace im = ImGui;
namespace imp = ImPlot;
namespace imfd = ImGuiFD;

namespace render_ui_detail {
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
static bool render_state_tr_ui(Body& body, ImGuiInputTextFlags flags = 0);
static bool render_state_att_ui(Body& body, ImGuiInputTextFlags flags = 0);
static bool render_state_ui(Body& body, ImGuiInputTextFlags flags = 0);
static bool render_body_propagation_ui(Body& body, bool editable = true);
static bool render_stat_state_ui(
    Station& stat,
    World& world,
    ImGuiInputTextFlags flags = 0
);
static void render_celestial_stats_ui(
    Celestial& cel,
    ImGuiInputTextFlags flags = 0,
    BodyEditDraft* draft = nullptr
);
static void render_satellite_stats_ui(
    Satellite& sat,
    ImGuiInputTextFlags flags = 0
);
static void render_mass_properties_ui(
    MassProperties& mp,
    ImGuiInputTextFlags flags = 0
);
static void render_station_stats_ui(
    Station& stat,
    RenderLoopConfig& cfg,
    RenderLoopState& state,
    World& world,
    ImGuiInputTextFlags flags = 0
);
static void render_station_draft_ui(
    Station& stat,
    RenderLoopState& state,
    World& world
);
static void render_station_instruments_ui(
    Station& stat,
    RenderLoopState& state,
    ImGuiInputTextFlags flags = 0
);
static bool body_edit_draft_changed(const BodyEditDraft& draft, const World& world);
static void render_body_list_row(
    const Body& body,
    World& world,
    RenderLoopConfig& cfg
);
static Body* draft_body_ptr(BodyEditDraft& draft) {
    switch (draft.edit_body_type) {
    case BodyType::unknown: return nullptr;
    case BodyType::celestial: return &draft.edit_celestial;
    case BodyType::satellite: return &draft.edit_satellite;
    case BodyType::station: return &draft.edit_station;
    }

    return nullptr;
}

static bool fields_readonly(ImGuiInputTextFlags flags) {
    return (flags & ImGuiInputTextFlags_ReadOnly) != 0;
}

static void push_edit_field_style() {
    im::PushStyleColor(ImGuiCol_FrameBg, ImVec4(0.31f, 0.25f, 0.36f, 1.0f));
    im::PushStyleColor(ImGuiCol_FrameBgHovered, ImVec4(0.38f, 0.30f, 0.46f, 1.0f));
    im::PushStyleColor(ImGuiCol_FrameBgActive, ImVec4(0.44f, 0.34f, 0.54f, 1.0f));
    im::PushStyleColor(ImGuiCol_Border, ImVec4(0.72f, 0.58f, 0.90f, 1.0f));
    im::PushStyleColor(ImGuiCol_CheckMark, ImVec4(0.84f, 0.68f, 0.96f, 1.0f));
    im::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
}

static void pop_edit_field_style() {
    im::PopStyleVar();
    im::PopStyleColor(5);
}

static void set_draft_simple_spin_from_w(BodyEditDraft& draft, ecref<vec3d> w) {
    draft.simple_spin_rate = w.norm();
    if (draft.simple_spin_rate > tol12) {
        draft.simple_spin_axis = w / draft.simple_spin_rate;
    } else {
        draft.simple_spin_axis = axis_z;
    }
    draft.simple_spin_edit_mode = SimpleSpinEditMode::rate;
}

static StatusCode apply_draft_simple_spin(BodyEditDraft& draft) {
    if (draft.edit_body_type != BodyType::celestial) return StatusCode::ok;
    if (draft.edit_celestial.attitude_model != CelestialAttitudeModel::simple_spin) {
        return StatusCode::ok;
    }
    if (!std::isfinite(draft.simple_spin_rate) || draft.simple_spin_rate < 0.0) {
        return StatusCode::invalid_input;
    }
    if (!finite_vec(draft.simple_spin_axis)) return StatusCode::invalid_input;

    f64 axis_norm = draft.simple_spin_axis.norm();
    if (draft.simple_spin_rate <= tol12) {
        draft.edit_celestial.x_att.w = vec3d0;
        if (axis_norm <= tol12) draft.simple_spin_axis = axis_z;
        return StatusCode::ok;
    }
    if (axis_norm <= tol12) return StatusCode::invalid_input;

    draft.simple_spin_axis /= axis_norm;
    draft.edit_celestial.x_att.w = draft.simple_spin_axis * draft.simple_spin_rate;

    return StatusCode::ok;
}

static bool render_state_tr_ui(Body& body, ImGuiInputTextFlags flags) {
    return im::InputDouble3("r", body.x_tr.r, "%.3f", flags)
           || im::InputDouble3("v", body.x_tr.v, "%.3f", flags);
}

static bool render_state_att_ui(Body& body, ImGuiInputTextFlags flags) {
    return im::InputDouble4("q", body.x_att.q, "%.3f", flags)
           || im::InputDouble3("w", body.x_att.w, "%.6E", flags);
}

static bool render_state_ui(Body& body, ImGuiInputTextFlags flags) {
    return render_state_tr_ui(body, flags) || render_state_att_ui(body, flags);
}

static bool scalar_changed(f64 a, f64 b, f64 tol = tol12) {
    return std::abs(a - b) > tol;
}

static bool vec3_changed(ecref<vec3d> a, ecref<vec3d> b, f64 tol = tol12) {
    return !a.isApprox(b, tol);
}

static bool vec4_changed(ecref<vec4d> a, ecref<vec4d> b, f64 tol = tol12) {
    return !a.isApprox(b, tol);
}

static bool vec7_changed(ecref<vec7d> a, ecref<vec7d> b, f64 tol = tol12) {
    return !a.isApprox(b, tol);
}

static bool simple_spin_cache_changed(
    const BodyEditDraft& draft,
    const Celestial& cel,
    f64 tol = tol12
) {
    if (draft.edit_celestial.attitude_model != CelestialAttitudeModel::simple_spin) {
        return false;
    }

    f64 rate = cel.x_att.w.norm();
    if (scalar_changed(draft.simple_spin_rate, rate, tol)) return true;
    if (rate <= tol && draft.simple_spin_rate <= tol) return false;

    vec3d axis = axis_z;
    if (rate > tol) axis = cel.x_att.w / rate;
    return vec3_changed(draft.simple_spin_axis, axis, tol);
}

static bool mat3_changed(ecref<mat3d> a, ecref<mat3d> b, f64 tol = tol12) {
    return !a.isApprox(b, tol);
}

static bool matx_changed(ecref<matXd> a, ecref<matXd> b, f64 tol = tol12) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return true;
    return !a.isApprox(b, tol);
}

static bool state_tr_changed(const StateTr& a, const StateTr& b, f64 tol = tol12) {
    return vec3_changed(a.r, b.r, tol) || vec3_changed(a.v, b.v, tol);
}

static bool state_att_changed(const StateAtt& a, const StateAtt& b, f64 tol = tol12) {
    return vec4_changed(a.q, b.q, tol) || vec3_changed(a.w, b.w, tol);
}

static bool body_common_changed(const Body& a, const Body& b, f64 tol = tol12) {
    return a.name != b.name || a.propagate_tr != b.propagate_tr
           || a.propagate_att != b.propagate_att || state_tr_changed(a.x_tr, b.x_tr, tol)
           || state_att_changed(a.x_att, b.x_att, tol);
}

static bool mass_properties_changed(
    const MassProperties& a,
    const MassProperties& b,
    f64 tol = tol12
) {
    return scalar_changed(a.mass, b.mass, tol) || mat3_changed(a.I, b.I, tol)
           || mat3_changed(a.I_inv, b.I_inv, tol) || a.principal_axes != b.principal_axes
           || vec3_changed(a.offset_body, b.offset_body, tol) || a.active != b.active;
}

static bool station_instrument_changed(
    const StationInstrument& a,
    const StationInstrument& b,
    f64 tol = tol12
) {
    return a.id != b.id || a.name != b.name || a.type != b.type || a.enabled != b.enabled
           || matx_changed(a.R, b.R, tol);
}

static bool station_instruments_changed(
    const umap<u32, StationInstrument>& a,
    const umap<u32, StationInstrument>& b,
    f64 tol = tol12
) {
    if (a.size() != b.size()) return true;

    for (const auto& [id, instr_a] : a) {
        auto it = b.find(id);
        if (it == b.end()) return true;
        if (station_instrument_changed(instr_a, it->second, tol)) return true;
    }

    return false;
}

static bool celestial_changed(const Celestial& a, const Celestial& b, f64 tol = tol12) {
    return body_common_changed(a, b, tol) || a.gravity_model != b.gravity_model
           || scalar_changed(a.mu, b.mu, tol) || a.degree != b.degree
           || a.order != b.order || vec7_changed(a.J, b.J, tol)
           || matx_changed(a.C, b.C, tol) || matx_changed(a.S, b.S, tol)
           || a.attitude_model != b.attitude_model
           || a.radiation_model != b.radiation_model
           || scalar_changed(a.ref_radius, b.ref_radius, tol)
           || scalar_changed(a.semimajor_axis, b.semimajor_axis, tol)
           || scalar_changed(a.semiminor_axis, b.semiminor_axis, tol)
           || scalar_changed(a.mean_radius, b.mean_radius, tol)
           || scalar_changed(a.eccentricity, b.eccentricity, tol)
           || scalar_changed(a.flattening, b.flattening, tol);
}

static bool satellite_changed(const Satellite& a, const Satellite& b, f64 tol = tol12) {
    return body_common_changed(a, b, tol)
           || mass_properties_changed(a.mass_properties, b.mass_properties, tol);
}

static bool station_changed(const Station& a, const Station& b, f64 tol = tol12) {
    return body_common_changed(a, b, tol) || a.anchored != b.anchored
           || a.anchor_id != b.anchor_id
           || vec3_changed(a.r_body_BCBF, b.r_body_BCBF, tol)
           || vec3_changed(a.llh_BCBF, b.llh_BCBF, tol)
           || a.next_instrument_id != b.next_instrument_id
           || a.enabled_instrument_ids != b.enabled_instrument_ids
           || station_instruments_changed(a.instruments, b.instruments, tol)
           || mass_properties_changed(a.mass_properties, b.mass_properties, tol);
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

void render_body_stats_ui(
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

        bool edit_locked = !cfg.stepper_cfg.paused;
        bool editable = cfg.edit_body_stats && !edit_locked;
        ImGuiInputTextFlags field_flags = editable ? 0 : ImGuiInputTextFlags_ReadOnly;
        bool draft_loaded = state.draft.edit_body_id == body->id
                            && state.draft.edit_body_type == body->body_type;
        bool draft_dirty = draft_loaded && body_edit_draft_changed(state.draft, world);

        if (!cfg.edit_body_stats && edit_locked) {
            im::BeginDisabled();
            im::Button("Edit");
            im::EndDisabled();
            im::SameLine();
            im::Text("Pause simulation to edit");
        } else if (!cfg.edit_body_stats && im::Button("Edit")) {
            state.draft.edit_body_status
                = load_body_edit_draft(body->id, state.draft, world);
            if (state.draft.edit_body_status == StatusCode::ok) {
                cfg.edit_body_stats = true;
                editable = !edit_locked;
                field_flags = editable ? 0 : ImGuiInputTextFlags_ReadOnly;
                draft_loaded = true;
                draft_dirty = false;
            }
        }

        if (cfg.edit_body_stats) {
            if (edit_locked) im::BeginDisabled();
            im::TextSL("Edit:");
            if (im::Button("Confirm")) {
                EntityId edit_body_id = state.draft.edit_body_id;
                state.draft.edit_body_status = apply_body_edit_draft(state.draft, world);
                if (state.draft.edit_body_status == StatusCode::ok) {
                    const Body* edited_body = world.body(edit_body_id);
                    if (edited_body == nullptr) {
                        state.draft.edit_body_status = StatusCode::body_not_found;
                    } else {
                        state.draft.edit_body_status = sync_scenario_body(
                            state.scenario,
                            *edited_body,
                            world
                        );
                    }
                }
                if (state.draft.edit_body_status == StatusCode::ok) {
                    state.scenario.dirty = true;
                    state.wksp.dirty = true;
                    if (cfg.camera.target_id == edit_body_id) {
                        sync_camera_tracking(cfg.camera, world);
                    }
                    cfg.edit_body_stats = false;
                    cancel_body_edit_draft(state.draft);
                }
            }
            if (edit_locked) im::EndDisabled();
            im::SameLine();
            if (edit_locked) im::BeginDisabled();
            if (im::Button("Reset")) {
                state.draft.edit_body_status
                    = load_body_edit_draft(body->id, state.draft, world);
                if (state.draft.edit_body_status == StatusCode::ok) {
                    cfg.edit_body_stats = true;
                    editable = !edit_locked;
                    field_flags = editable ? 0 : ImGuiInputTextFlags_ReadOnly;
                    draft_loaded = true;
                    draft_dirty = false;
                }
            }
            if (edit_locked) im::EndDisabled();
            im::SameLine();
            if (im::Button("Cancel")) {
                cancel_body_edit_draft(state.draft);
                cfg.edit_body_stats = false;
            }
        }
        if (state.draft.edit_body_id != kInvalidEntityId
            && state.draft.edit_body_status != StatusCode::ok) {
            status_text("Edit Body Status", state.draft.edit_body_status);
        }

        if (editable) {
            im::TextColored(ImVec4(0.84f, 0.68f, 0.96f, 1.0f), "Mode: Editing Draft");
            im::SameLine();
            if (draft_dirty) {
                im::TextColored(ImVec4(1.0f, 0.82f, 0.34f, 1.0f), "Draft: Modified");
            } else {
                im::Text("Draft: Clean");
            }
            push_edit_field_style();
        } else {
            im::Text("Mode: Read Only");
            if (cfg.edit_body_stats && edit_locked) {
                im::Text("Edit draft locked until simulation is paused");
            }
        }

        im::TextSL("Name:");
        Body* edit_body = editable ? draft_body_ptr(state.draft) : body;
        if (edit_body != nullptr) {
            im::InputText("##Name", &edit_body->name, field_flags);
        }

        switch (body->body_type) {
        case BodyType::unknown: break;
        case BodyType::celestial: {
            Celestial* cel = editable ? &state.draft.edit_celestial
                                      : world.celestial(cfg.body_stats_id);
            if (cel == nullptr) {
                im::Text("Invalid Celestial");
            } else {
                render_celestial_stats_ui(
                    *cel,
                    field_flags,
                    editable ? &state.draft : nullptr
                );
            }
        } break;
        case BodyType::satellite: {
            Satellite* sat = editable ? &state.draft.edit_satellite
                                      : world.satellite(cfg.body_stats_id);
            if (sat == nullptr) {
                im::Text("Invalid Satellite");
            } else {
                render_satellite_stats_ui(*sat, field_flags);
            }
        } break;
        case BodyType::station: {
            Station* stat
                = editable ? &state.draft.edit_station : world.station(cfg.body_stats_id);
            if (stat == nullptr) {
                im::Text("Invalid Station");
            } else {
                render_station_stats_ui(*stat, cfg, state, world, field_flags);
            }
        } break;
        }

        // im::Checkbox("Emits Gravity", &body->emits_gravity);
        // im::Checkbox("Emits Radiation", &body->emits_radiation);
        // im::Checkbox("Has Atmosphere", &body->has_atmosphere);

        if (editable) pop_edit_field_style();

        if (active && im::Button("Set as Camera Target")) {
            cfg.camera.target_id = cfg.body_stats_id;
            sync_camera_tracking(cfg.camera, world);
        }
        if (cfg.stepper_cfg.paused) {
            if (active && im::Button("Make Inactive")) {
                if (world.make_inactive(cfg.body_stats_id)) {
                    StatusCode status = sync_scenario_body_active(
                        state.scenario,
                        cfg.body_stats_id,
                        false,
                        world
                    );
                    if (status == StatusCode::ok) {
                        state.scenario.dirty = true;
                        state.wksp.dirty = true;
                        if (cfg.camera.target_id == body->id) {
                            cycle_active_id(cfg.camera.target_id, world, 1);
                            sync_camera_tracking(cfg.camera, world);
                        }
                    } else {
                        world.make_active(cfg.body_stats_id);
                    }
                }
            }
            if (!active && im::Button("Make Active")) {
                if (world.make_active(cfg.body_stats_id)) {
                    StatusCode status = sync_scenario_body_active(
                        state.scenario,
                        cfg.body_stats_id,
                        true,
                        world
                    );
                    if (status == StatusCode::ok) {
                        state.scenario.dirty = true;
                        state.wksp.dirty = true;
                    } else {
                        world.make_inactive(cfg.body_stats_id);
                    }
                }
            }
        }
    }

    im::End();
}

void render_add_body(World& world, RenderLoopConfig& cfg, RenderLoopState& state) {
    im::Begin("Add Body");
    if (cfg.stepper_cfg.paused && !state.add_body) {
        im::Separator();
        if (im::Button("Add Body")) {
            state.add_body = true;
            if (state.add_body_type == BodyType::unknown) {
                state.add_body_type = BodyType::celestial;
            }
            init_add_body_draft_defaults(state, world, state.add_body_type);
            state.add_body_status = StatusCode::ok;
        }
    }

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
            if (im::Button("Create")) {
                EntityId id = kInvalidEntityId;
                state.add_body_status = validate_add_body_draft(state, world);

                if (state.add_body_status == StatusCode::ok) {
                    switch (state.add_body_type) {
                    case BodyType::unknown: break;
                    case BodyType::celestial: {
                        ScenarioCelestialConfig cel_cfg;
                        state.add_body_status = make_scenario_celestial_config(
                            state.scenario.config,
                            state.temp_celestial,
                            true,
                            cel_cfg
                        );
                        if (state.add_body_status != StatusCode::ok) break;

                        auto cel = std::make_unique<Celestial>(state.temp_celestial);
                        id = world.insert_celestial(std::move(cel));
                        if (id == kInvalidEntityId) {
                            state.add_body_status = StatusCode::body_not_found;
                            break;
                        }

                        const string config_id = make_unique_scenario_body_id(
                            state.scenario,
                            BodyType::celestial,
                            id
                        );
                        if (config_id.empty()) {
                            state.add_body_status = StatusCode::invalid_input;
                            break;
                        }

                        cel_cfg.id = config_id;
                        state.add_body_status = register_scenario_body_mapping(
                            state.scenario.build_result,
                            BodyType::celestial,
                            config_id,
                            id
                        );
                        if (state.add_body_status != StatusCode::ok) break;

                        state.scenario.config.celestials.push_back(std::move(cel_cfg));
                    } break;
                    case BodyType::satellite: {
                        ScenarioSatelliteConfig sat_cfg;
                        state.add_body_status = make_scenario_satellite_config(
                            state.temp_satellite,
                            true,
                            sat_cfg
                        );
                        if (state.add_body_status != StatusCode::ok) break;

                        auto sat = std::make_unique<Satellite>(state.temp_satellite);
                        id = world.insert_satellite(std::move(sat));
                        if (id == kInvalidEntityId) {
                            state.add_body_status = StatusCode::body_not_found;
                            break;
                        }

                        const string config_id = make_unique_scenario_body_id(
                            state.scenario,
                            BodyType::satellite,
                            id
                        );
                        if (config_id.empty()) {
                            state.add_body_status = StatusCode::invalid_input;
                            break;
                        }

                        sat_cfg.id = config_id;
                        state.add_body_status = register_scenario_body_mapping(
                            state.scenario.build_result,
                            BodyType::satellite,
                            config_id,
                            id
                        );
                        if (state.add_body_status != StatusCode::ok) break;

                        state.scenario.config.satellites.push_back(std::move(sat_cfg));
                    } break;
                    case BodyType::station: {
                        ScenarioStationConfig stat_cfg;
                        state.add_body_status = make_scenario_station_config(
                            state.temp_station,
                            true,
                            state.scenario.build_result,
                            stat_cfg
                        );
                        if (state.add_body_status != StatusCode::ok) break;

                        auto stat = std::make_unique<Station>(state.temp_station);
                        id = world.insert_station(std::move(stat));
                        if (id == kInvalidEntityId) {
                            state.add_body_status = StatusCode::body_not_found;
                            break;
                        }

                        const string config_id = make_unique_scenario_body_id(
                            state.scenario,
                            BodyType::station,
                            id
                        );
                        if (config_id.empty()) {
                            state.add_body_status = StatusCode::invalid_input;
                            break;
                        }

                        stat_cfg.id = config_id;
                        state.add_body_status = register_scenario_body_mapping(
                            state.scenario.build_result,
                            BodyType::station,
                            config_id,
                            id
                        );
                        if (state.add_body_status != StatusCode::ok) break;

                        state.scenario.config.stations.push_back(std::move(stat_cfg));
                    } break;
                    }

                    if (state.add_body_status == StatusCode::ok
                        && id != kInvalidEntityId) {
                        state.scenario.dirty = true;
                        state.wksp.dirty = true;
                        cfg.body_stats_id = id;
                        cfg.camera.target_id = id;
                        switch (state.add_body_type) {
                        case BodyType::unknown: break;
                        case BodyType::celestial:
                            state.temp_celestial = Celestial{};
                            break;
                        case BodyType::satellite:
                            state.temp_satellite = Satellite{};
                            break;
                        case BodyType::station: state.temp_station = Station{}; break;
                        }
                    } else if (state.add_body_status == StatusCode::ok) {
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
                status_text("Add Body Status", state.add_body_status);
            }
        }
    } else {
        im::SetWindowCollapsed(true);
    }
    im::End();
}

static bool render_body_propagation_ui(Body& body, bool editable) {
    bool changed = false;

    if (im::CollapsingHeader("Propagation")) {
        if (!editable) im::BeginDisabled();
        changed = im::Checkbox("Propagate Translation", &body.propagate_tr)
                  || im::Checkbox("Propagate Attitude", &body.propagate_att);
        if (!editable) im::EndDisabled();
        im::Separator();
    }

    return changed;
}

static bool render_stat_state_ui(Station& stat, World& world, ImGuiInputTextFlags flags) {
    flags |= ImGuiInputTextFlags_ReadOnly;
    vec3d r = world.stat_r_inertial(stat.id);
    vec4d q = world.stat_q_inertial(stat.id);
    return im::InputDouble3("r", r, "%.3f", flags)
           || im::InputDouble4("q", q, "%.3f", flags);
}

static void render_celestial_stats_ui(
    Celestial& cel,
    ImGuiInputTextFlags flags,
    BodyEditDraft* draft
) {
    bool editable = !fields_readonly(flags);

    render_state_tr_ui(cel, flags);
    im::InputDouble4("q", cel.x_att.q, "%.3f", flags);
    if (cel.attitude_model != CelestialAttitudeModel::simple_spin) {
        im::InputDouble3("w", cel.x_att.w, "%.6E", flags);
    }

    if (im::CollapsingHeader("Attitude Model")) {
        i32 att_idx = static_cast<i32>(cel.attitude_model);
        const char* att_name[] = {"Fixed", "Simple Spin", "Provider"};
        if (!editable) im::BeginDisabled();
        if (im::Combo("##Attitude Model", &att_idx, att_name, 3)) {
            if (static_cast<CelestialAttitudeModel>(att_idx)
                == CelestialAttitudeModel::provider) {
                att_idx = static_cast<i32>(CelestialAttitudeModel::fixed);
            } // TEMP: no attitude providers yet, fallback to fixed
            cel.attitude_model = static_cast<CelestialAttitudeModel>(att_idx);
        }
        if (!editable) im::EndDisabled();
        if (cel.attitude_model == CelestialAttitudeModel::simple_spin) {
            if (editable && draft != nullptr) {
                i32 spin_mode = static_cast<i32>(draft->simple_spin_edit_mode);
                if (im::RadioButton("Edit Spin Rate", &spin_mode, 0)) {
                    draft->simple_spin_edit_mode = SimpleSpinEditMode::rate;
                }
                im::SameLine();
                if (im::RadioButton("Edit Spin Axis", &spin_mode, 1)) {
                    draft->simple_spin_edit_mode = SimpleSpinEditMode::axis;
                }

                ImGuiInputTextFlags rate_flags = flags;
                ImGuiInputTextFlags axis_flags = flags;
                if (draft->simple_spin_edit_mode != SimpleSpinEditMode::rate) {
                    rate_flags |= ImGuiInputTextFlags_ReadOnly;
                }
                if (draft->simple_spin_edit_mode != SimpleSpinEditMode::axis) {
                    axis_flags |= ImGuiInputTextFlags_ReadOnly;
                }

                im::InputDouble(
                    "Spin Rate",
                    &draft->simple_spin_rate,
                    0.0,
                    0.0,
                    "%.6E",
                    rate_flags
                );
                im::InputDouble3(
                    "Spin Axis",
                    draft->simple_spin_axis,
                    "%.6E",
                    axis_flags
                );
                im::Text("Spin angular velocity is applied when Confirm is pressed.");
            } else {
                f64 w_mag = cel.x_att.w.norm();
                vec3d spin_axis = axis_z;
                if (w_mag > tol12) {
                    spin_axis = cel.x_att.w / w_mag;
                }
                ImGuiInputTextFlags spin_flags = flags | ImGuiInputTextFlags_ReadOnly;
                im::InputDouble("Spin Rate", &w_mag, 0.0, 0.0, "%.6E", spin_flags);
                im::InputDouble3("Spin Axis", spin_axis, "%.6E", spin_flags);
            }
        }
        im::Separator();
    }

    if (im::CollapsingHeader("Celestial Model", ImGuiTreeNodeFlags_DefaultOpen)) {
        bool edit_a = im::InputDouble(
            "Semimajor Axis",
            &cel.semimajor_axis,
            0.0,
            0.0,
            "%.3f",
            flags
        );
        bool edit_b = im::InputDouble(
            "Semiminor Axis",
            &cel.semiminor_axis,
            0.0,
            0.0,
            "%.3f",
            flags
        );
        if (edit_a || edit_b) {
            cel.mean_radius = mean_from_semiaxes(cel.semimajor_axis, cel.semiminor_axis);
            cel.eccentricity = ecc_from_semiaxes(cel.semimajor_axis, cel.semiminor_axis);
            cel.flattening = flat_from_semiaxes(cel.semimajor_axis, cel.semiminor_axis);
        }

        ImGuiInputTextFlags computed_flags = flags | ImGuiInputTextFlags_ReadOnly;
        im::InputDouble(
            "Mean Radius",
            &cel.mean_radius,
            0.0,
            0.0,
            "%.3f",
            computed_flags
        );
        im::InputDouble(
            "Eccentricity",
            &cel.eccentricity,
            0.0,
            0.0,
            "%.3f",
            computed_flags
        );
        im::InputDouble("Flattening", &cel.flattening, 0.0, 0.0, "%.3f", computed_flags);
        im::Separator();
    }

    if (im::CollapsingHeader("Gravity Model", ImGuiTreeNodeFlags_DefaultOpen)) {
        i32 model_idx = static_cast<i32>(cel.gravity_model);
        const char* model_names[] = {"Pointmass", "Zonal", "Spherical Harmonics"};
        if (!editable) im::BeginDisabled();
        if (im::Combo("##Gravity Model", &model_idx, model_names, 3)) {
            // TODO: update as new models are added
            cel.gravity_model = static_cast<GravityModel>(model_idx);
            // TODO: set to fall back to zonal or pointmass if coefs not available
            // (depending on J, C, S active bools?)
        }
        if (!editable) im::EndDisabled();

        im::InputDouble("Mu", &cel.mu, 0.0, 0.0, "%.3f", flags);
        switch (cel.gravity_model) {
        case GravityModel::pointmass: break;
        case GravityModel::zonal: {
            im::InputDouble("Reference Radius", &cel.ref_radius, 0.0, 0.0, "%.3f", flags);
            im::InputInt("Degree", &cel.degree, 1, 100, flags);
            im::InputDouble7("Zonal Coefficients (J)", cel.J, "%.3f", flags);
        } break;
        case GravityModel::spherical_harmonics: {
            im::InputDouble("Reference Radius", &cel.ref_radius, 0.0, 0.0, "%.3f", flags);
            im::InputInt("Degree", &cel.degree, 1, 100, flags);
            im::InputInt("Order", &cel.order, 1, 100, flags);
            bool csactive = cel.C.cols() > 0 && cel.C.rows() > 0 && cel.S.cols() > 0
                            && cel.S.rows() > 0;
            im::Text("C/S Matrices Active: %s", bool_str(csactive).c_str());
        } break;
        }
        im::Separator();
    }

    render_body_propagation_ui(cel, editable);
    // TODO: allow change of coefs (dropdown of providers)
}

static void render_mass_properties_ui(MassProperties& mp, ImGuiInputTextFlags flags) {
    bool editable = !fields_readonly(flags);

    if (im::CollapsingHeader("Mass Properties", ImGuiTreeNodeFlags_DefaultOpen)) {
        im::InputDouble("Mass", &mp.mass, 0.0, 0.0, "%.3f", flags);

        if (!editable) im::BeginDisabled();
        if (im::Checkbox("Principal Axes", &mp.principal_axes)) {
            if (mp.principal_axes) {
                vec3d I_diag = mp.I.diagonal();
                mp.I = I_diag.asDiagonal();
                mp.I_inv = mp.I.inverse();
            }
        }
        if (!editable) im::EndDisabled();
        if (mp.principal_axes) {
            vec3d I_diag = mp.I.diagonal();
            if (im::InputDouble3("Inertia Tensor", I_diag, "%.3f", flags)) {
                mp.I = I_diag.asDiagonal();
                mp.I_inv = mp.I.inverse();
            }
        } else {
            if (im::InputDouble3x3("Inertia Tensor", mp.I, "%.3f", flags)) {
                mp.I_inv = mp.I.inverse();
            }
        }

        im::Indent();
        im::InputDouble3("Offset", mp.offset_body, "%.3f", flags);
        if (!editable) im::BeginDisabled();
        if (im::Button("Recompute Inertia Tensor about COM")) {
            mp.I = inertia_PAT(mp.I, mp.mass, mp.offset_body);
            if (!mp.I.isDiagonal()) mp.principal_axes = false;
            mp.offset_body = vec3d0;
            mp.I_inv = mp.I.inverse();
        }
        if (!editable) im::EndDisabled();
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

static void render_satellite_stats_ui(Satellite& sat, ImGuiInputTextFlags flags) {
    bool editable = !fields_readonly(flags);

    render_state_ui(sat, flags);
    render_mass_properties_ui(sat.mass_properties, flags);
    render_body_propagation_ui(sat, editable);
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

static void render_station_instruments_ui(
    Station& stat,
    RenderLoopState& state,
    ImGuiInputTextFlags flags
) {
    bool editable = !fields_readonly(flags);

    if (im::CollapsingHeader("Instruments")) {
        im::Indent();
        bool enabled_changed = false;
        for (auto& [id, instr] : stat.instruments) {
            im::PushID(&instr.name);
            im::TextSL("Name:");
            im::InputText("##name", &instr.name, flags);
            im::Text("ID: %u", instr.id);
            im::Text("Observation Type: %s", observation_type_str(instr.type).c_str());
            im::InputMatXd("Covariance", instr.R, "%.3E", flags);
            if (!editable) im::BeginDisabled();
            if (im::Checkbox("Enabled", &instr.enabled)) {
                enabled_changed = true;
            }
            if (!editable) im::EndDisabled();
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

            if (!editable) im::BeginDisabled();
            im::Combo("Type", &state.add_instrument_type, type_names, 8);
            im::TextSL("Name:");
            im::InputText("##new instrument name", &state.add_instrument_name, flags);

            ObservationType type
                = static_cast<ObservationType>(state.add_instrument_type);
            i32 dim = measurement_dim(type);
            sync_add_instrument_diag(state, dim);

            for (i32 i = 0; i < state.add_instrument_R_diag.size(); ++i) {
                im::PushID(i);
                f64 variance = state.add_instrument_R_diag(i);
                string label = "R(" + std::to_string(i) + "," + std::to_string(i) + ")";
                im::InputDouble(label.c_str(), &variance, 0.0, 0.0, "%.3E", flags);
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
            if (!editable) im::EndDisabled();
            // im::SameLine();
            // if (im::Button("Reset")) {
            // }
            if (state.add_instrument_status != StatusCode::ok) {
                status_text("Add Instrument Status", state.add_instrument_status);
            }
        }
        im::Unindent();
    }
}

static void render_station_stats_ui(
    Station& stat,
    RenderLoopConfig& cfg,
    RenderLoopState& state,
    World& world,
    ImGuiInputTextFlags flags
) {
    bool editable = !fields_readonly(flags);

    if (!editable) im::BeginDisabled();
    im::Checkbox("Anchored", &stat.anchored);
    if (!editable) im::EndDisabled();

    if (stat.anchored) {
        Celestial* cel = world.celestial(stat.anchor_id);
        string anchor_name;
        if (cel == nullptr) {
            anchor_name = "unknown";
        } else {
            anchor_name = cel->name;
        }
        render_stat_state_ui(stat, world, flags);
        im::Text("Anchor: %s (ID: %llu)", anchor_name.c_str(), stat.anchor_id);
        im::Indent();
        if (!editable) im::BeginDisabled();
        bool prev_anchor = im::Button("Prev Anchor");
        bool next_anchor = im::Button("Next Anchor");
        if (!editable) im::EndDisabled();
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
            if (im::InputDouble3("r_body (BCBF)", stat.r_body_BCBF, "%.3f", flags)) {
                stat.llh_BCBF = bcbf_to_detic(stat.r_body_BCBF, *cel);
            }
            if (im::InputDouble3("Planetodetic LLH", stat.llh_BCBF, "%.3f", flags)) {
                stat.r_body_BCBF = detic_to_bcbf(stat.llh_BCBF, *cel);
            }
        }

    } else {
        render_state_ui(stat, flags);
        render_mass_properties_ui(stat.mass_properties, flags);
    }

    render_body_propagation_ui(stat, editable);
    render_station_instruments_ui(stat, state, flags);
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
    case BodyType::unknown: return StatusCode::unsupported_type;
    case BodyType::celestial: {
        const Celestial& cel = state.temp_celestial;
        if (!finite_state_tr(cel.x_tr) || !finite_state_att(cel.x_att)) {
            return StatusCode::invalid_state;
        }
        if (!finite_nonneg(cel.mu)) return StatusCode::invalid_input;
        if (!finite_nonneg(cel.semimajor_axis) || !finite_nonneg(cel.semiminor_axis)) {
            return StatusCode::invalid_shape;
        }
        if (cel.semimajor_axis == 0.0 && cel.semiminor_axis == 0.0) {
            return StatusCode::invalid_shape;
        }
        if (cel.semimajor_axis > 0.0 && cel.semiminor_axis > cel.semimajor_axis) {
            return StatusCode::invalid_shape;
        }
        if (cel.ref_radius < 0.0 || !std::isfinite(cel.ref_radius)) {
            return StatusCode::invalid_shape;
        }
    } break;
    case BodyType::satellite: {
        const Satellite& sat = state.temp_satellite;
        if (!finite_state_tr(sat.x_tr) || !finite_state_att(sat.x_att)) {
            return StatusCode::invalid_state;
        }
        if (sat.propagate_att) {
            if (!finite_pos(sat.mass_properties.mass))
                return StatusCode::invalid_mass_properties;
            if (!finite_mat(sat.mass_properties.I)
                || !finite_mat(sat.mass_properties.I_inv)) {
                return StatusCode::invalid_mass_properties;
            }
            if (!sat.mass_properties.active) return StatusCode::invalid_mass_properties;
        }
    } break;
    case BodyType::station: {
        const Station& stat = state.temp_station;
        if (stat.anchored) {
            if (world.celestial(stat.anchor_id) == nullptr)
                return StatusCode::anchor_not_found;
            if (!finite_vec(stat.r_body_BCBF) || !finite_vec(stat.llh_BCBF)) {
                return StatusCode::anchor_not_found;
            }
        } else {
            if (!finite_state_tr(stat.x_tr) || !finite_state_att(stat.x_att)) {
                return StatusCode::invalid_state;
            }
            if (stat.propagate_att && !stat.mass_properties.active) {
                return StatusCode::invalid_mass_properties;
            }
        }
    } break;
    }

    return StatusCode::ok;
}

void render_body_lists(World& world, RenderLoopConfig& cfg, RenderLoopState& state) {
    im::Begin("Bodies");
    // TODO: default this to closed

    const char* filter_names[] = {"all", "active", "inactive"};
    i32 filter_idx = static_cast<i32>(state.list_filter);
    if (im::Combo("Filter", &filter_idx, filter_names, 3)) {
        state.list_filter = static_cast<BodyFilterMode>(filter_idx);
    }

    if (ImGui::CollapsingHeader("Celestials")) {
        svec<EntityId> ids = world.celestial_ids(state.list_filter);
        for (const EntityId id : ids) {
            Celestial* cel = world.celestial(id);
            if (cel == nullptr) continue;
            render_body_list_row(*cel, world, cfg);
        }
    }

    if (ImGui::CollapsingHeader("Satellites")) {
        svec<EntityId> ids = world.satellite_ids(state.list_filter);
        for (const EntityId id : ids) {
            Satellite* sat = world.satellite(id);
            if (sat == nullptr) continue;
            render_body_list_row(*sat, world, cfg);
        }
    }

    if (ImGui::CollapsingHeader("Stations")) {
        svec<EntityId> ids = world.station_ids(state.list_filter);
        for (const EntityId id : ids) {
            Station* stat = world.station(id);
            if (stat == nullptr) continue;
            render_body_list_row(*stat, world, cfg);
        }
    }

    im::End();
}

static StatusCode validate_mass_properties(const MassProperties& mp, f64 tol) {
    if (mp.mass == 0.0) return StatusCode::ok;
    if (!finite_pos(mp.mass)) return StatusCode::invalid_mass_properties;
    if (!finite_mat(mp.I)) return StatusCode::invalid_mass_properties;

    for (i32 i = 0; i < 3; ++i) {
        if (!finite_pos(mp.I(i, i))) return StatusCode::invalid_mass_properties;
    }
    if (!finite_nonzero(mp.I.determinant(), tol))
        return StatusCode::invalid_mass_properties;

    return StatusCode::ok;
}

static StatusCode validate_body(const Body& body, const World& world) {
    switch (body.body_type) {
    case BodyType::unknown: return StatusCode::unsupported_type;
    case BodyType::celestial: {
        const Celestial* cel = dynamic_cast<const Celestial*>(&body);
        if (cel == nullptr) return StatusCode::unsupported_type;
        if (!finite_state_tr(cel->x_tr) || !finite_state_att(cel->x_att)) {
            return StatusCode::invalid_state;
        }
        if (!finite_nonneg(cel->mu)) return StatusCode::invalid_input;
        if (!finite_nonneg(cel->semimajor_axis) || !finite_nonneg(cel->semiminor_axis)) {
            return StatusCode::invalid_shape;
        }
        if (cel->ref_radius < 0.0 || !std::isfinite(cel->ref_radius)) {
            return StatusCode::invalid_shape;
        }
    } break;
    case BodyType::satellite: {
        const Satellite* sat = dynamic_cast<const Satellite*>(&body);
        if (sat == nullptr) return StatusCode::unsupported_type;
        if (!finite_state_tr(sat->x_tr) || !finite_state_att(sat->x_att)) {
            return StatusCode::invalid_state;
        }
        if (sat->propagate_att) {
            if (!finite_pos(sat->mass_properties.mass))
                return StatusCode::invalid_mass_properties;
            if (!finite_mat(sat->mass_properties.I)
                || !finite_mat(sat->mass_properties.I_inv)) {
                return StatusCode::invalid_mass_properties;
            }
            if (!sat->mass_properties.active) return StatusCode::invalid_mass_properties;
        }
    } break;
    case BodyType::station: {
        const Station* stat = dynamic_cast<const Station*>(&body);
        if (stat == nullptr) return StatusCode::unsupported_type;
        if (stat->anchored) {
            if (world.celestial(stat->anchor_id) == nullptr)
                return StatusCode::anchor_not_found;
            if (!finite_vec(stat->r_body_BCBF) || !finite_vec(stat->llh_BCBF)) {
                return StatusCode::anchor_not_found;
            }
        } else {
            if (!finite_state_tr(stat->x_tr) || !finite_state_att(stat->x_att)) {
                return StatusCode::invalid_state;
            }
            if (stat->propagate_att && !stat->mass_properties.active) {
                return StatusCode::invalid_mass_properties;
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
    case BodyType::unknown: return StatusCode::unsupported_type;
    case BodyType::celestial: {
        const Celestial* cel = world.celestial(draft.edit_body_id);
        if (cel == nullptr) return StatusCode::body_not_found;
        draft.edit_celestial = *cel;
        set_draft_simple_spin_from_w(draft, cel->x_att.w);
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

static bool body_edit_draft_changed(const BodyEditDraft& draft, const World& world) {
    switch (draft.edit_body_type) {
    case BodyType::unknown: return false;
    case BodyType::celestial: {
        const Celestial* cel = world.celestial(draft.edit_body_id);
        if (cel == nullptr) return false;
        return celestial_changed(*cel, draft.edit_celestial)
               || simple_spin_cache_changed(draft, *cel);
    }
    case BodyType::satellite: {
        const Satellite* sat = world.satellite(draft.edit_body_id);
        if (sat == nullptr) return false;
        return satellite_changed(*sat, draft.edit_satellite);
    }
    case BodyType::station: {
        const Station* stat = world.station(draft.edit_body_id);
        if (stat == nullptr) return false;
        return station_changed(*stat, draft.edit_station);
    }
    }

    return false;
}

static void render_body_list_row(const Body& body, World& world, RenderLoopConfig& cfg) {
    const bool active = world.is_active(body.id);
    const bool selected = cfg.body_stats_id == body.id;

    im::PushID(static_cast<int>(body.id));
    float row_height = ImGui::GetTextLineHeightWithSpacing();
    if (ImGui::Selectable(
            "##row",
            selected,
            ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap,
            ImVec2(0.0f, row_height)

        )) {
        cfg.body_stats_id = body.id;
    }
    if (im::IsItemHovered()) {
        if (im::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
            cfg.camera.target_id = body.id;
            sync_camera_tracking(cfg.camera, world);
        }
    }

    im::SameLine();

    if (active) {
        im::TextColored(rl_to_im(DARKGREEN), "active");
    } else {
        im::TextColored(rl_to_im(RED), "inactive");
    }
    im::SameLine();

    string camera_target_marker = "";
    if (cfg.camera.target_id == body.id) camera_target_marker = " <";
    im::Text("%s (id: %llu)%s", body.name.c_str(), body.id, camera_target_marker.c_str());

    im::PopID();
}

static StatusCode validate_body_edit_draft(
    const BodyEditDraft& draft,
    const World& world
) {
    StatusCode status = StatusCode::ok;
    switch (draft.edit_body_type) {
    case BodyType::unknown: return StatusCode::unsupported_type;
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
    StatusCode status = apply_draft_simple_spin(draft);
    if (status != StatusCode::ok) return status;

    status = validate_body_edit_draft(draft, world);
    if (status != StatusCode::ok) return status;

    switch (draft.edit_body_type) {
    case BodyType::unknown: return StatusCode::unsupported_type;
    case BodyType::celestial: {
        Celestial* cel = world.celestial(draft.edit_body_id);
        if (cel == nullptr) return StatusCode::body_not_found;
        *cel = draft.edit_celestial;
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

} // namespace render_ui_detail
