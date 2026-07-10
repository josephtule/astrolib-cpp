#include "graphics/render_ui.hpp"
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
#include <string>

namespace im = ImGui;
namespace imp = ImPlot;
namespace imfd = ImGuiFD;

// imfd::settings.ascii

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
    case StatusCode::matrix_invert_failed: return RED;
    }

    return RAYWHITE;
}

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
static void open_scenario_file_dialog(RenderLoopState& state);
static void render_open_scenario_file_dialog(RenderLoopState& state);

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

static bool render_state_tr_ui(Body& body, ImGuiInputTextFlags flags = 0);
static bool render_state_att_ui(Body& body, ImGuiInputTextFlags flags = 0);
static bool render_state_ui(Body& body, ImGuiInputTextFlags flags = 0);
static bool render_body_propagation_ui(Body& body, bool editable = true);
static bool render_stat_state_ui(
    Station& stat,
    World& world,
    ImGuiInputTextFlags flags = 0
);
static void render_celestial_stats_ui(Celestial& cel, ImGuiInputTextFlags flags = 0);
static void render_satellite_stats_ui(Satellite& sat, ImGuiInputTextFlags flags = 0);
static void render_mass_properties_ui(MassProperties& mp, ImGuiInputTextFlags flags = 0);
static void render_station_stats_ui(
    Station& stat,
    RenderLoopConfig& cfg,
    RenderLoopState& state,
    World& world,
    ImGuiInputTextFlags flags = 0
);
static void render_station_draft_ui(Station& stat, RenderLoopState& state, World& world);
static void sync_add_instrument_diag(RenderLoopState& state, i32 dim);
static void render_station_instruments_ui(
    Station& stat,
    RenderLoopState& state,
    ImGuiInputTextFlags flags = 0
);
static bool body_edit_draft_changed(const BodyEditDraft& draft, const World& world);
static void render_body_list_row(const Body& body, World& world, RenderLoopConfig& cfg);
static void request_scenario_save(
    RenderLoopState& state,
    const RenderLoopConfig& cfg,
    const World& world,
    bool overwrite = false
);

static const ScenarioCelestialConfig* find_celestial_config(
    const ScenarioConfig& cfg,
    const string& id
) {
    for (const auto& cel : cfg.celestials) {
        if (cel.id == id) return &cel;
    }
    return nullptr;
}
static ScenarioCelestialConfig* find_celestial_config_mut(
    ScenarioConfig& cfg,
    const string& id
) {
    for (auto& cel : cfg.celestials) {
        if (cel.id == id) return &cel;
    }
    return nullptr;
}
static const ScenarioCelestialConfig* find_celestial_config(
    const ScenarioSession& scenario,
    const EntityId id
) {
    auto it = scenario.build_result.celestial_config_ids.find(id);
    if (it == scenario.build_result.celestial_config_ids.end()) return nullptr;
    string id_str = it->second;
    return find_celestial_config(scenario.config, id_str);
}
static ScenarioCelestialConfig* find_celestial_config_mut(
    ScenarioSession& scenario,
    const EntityId id
) {
    auto it = scenario.build_result.celestial_config_ids.find(id);
    if (it == scenario.build_result.celestial_config_ids.end()) return nullptr;
    string id_str = it->second;
    return find_celestial_config_mut(scenario.config, id_str);
}
static const ScenarioSatelliteConfig* find_satellite_config(
    const ScenarioConfig& cfg,
    const string& id
) {
    for (const auto& sat : cfg.satellites) {
        if (sat.id == id) return &sat;
    }
    return nullptr;
}
static ScenarioSatelliteConfig* find_satellite_config_mut(
    ScenarioConfig& cfg,
    const string& id
) {
    for (auto& sat : cfg.satellites) {
        if (sat.id == id) return &sat;
    }
    return nullptr;
}
static const ScenarioSatelliteConfig* find_satellite_config(
    const ScenarioSession& scenario,
    const EntityId id
) {
    auto it = scenario.build_result.satellite_config_ids.find(id);
    if (it == scenario.build_result.satellite_config_ids.end()) return nullptr;
    string id_str = it->second;
    return find_satellite_config(scenario.config, id_str);
}
static ScenarioSatelliteConfig* find_satellite_config_mut(
    ScenarioSession& scenario,
    const EntityId id
) {
    auto it = scenario.build_result.satellite_config_ids.find(id);
    if (it == scenario.build_result.satellite_config_ids.end()) return nullptr;
    string id_str = it->second;
    return find_satellite_config_mut(scenario.config, id_str);
}
static const ScenarioStationConfig* find_station_config(
    const ScenarioConfig& cfg,
    const string& id
) {
    for (const auto& stat : cfg.stations) {
        if (stat.id == id) return &stat;
    }
    return nullptr;
}
static ScenarioStationConfig* find_station_config_mut(
    ScenarioConfig& cfg,
    const string& id
) {
    for (auto& stat : cfg.stations) {
        if (stat.id == id) return &stat;
    }
    return nullptr;
}
static const ScenarioStationConfig* find_station_config(
    const ScenarioSession& scenario,
    const EntityId id
) {
    auto it = scenario.build_result.station_config_ids.find(id);
    if (it == scenario.build_result.station_config_ids.end()) return nullptr;
    string id_str = it->second;
    return find_station_config(scenario.config, id_str);
}
static ScenarioStationConfig* find_station_config_mut(
    ScenarioSession& scenario,
    const EntityId id
) {
    auto it = scenario.build_result.station_config_ids.find(id);
    if (it == scenario.build_result.station_config_ids.end()) return nullptr;
    string id_str = it->second;
    return find_station_config_mut(scenario.config, id_str);
}

static StatusCode sync_scenario_body_active(
    ScenarioSession& scenario,
    const EntityId id,
    const bool active,
    const World& world
) {
    const auto* body = world.body(id);
    switch (body->body_type) {
    case BodyType::celestial: {
        ScenarioCelestialConfig* cel_cfg = find_celestial_config_mut(scenario, id);
        if (cel_cfg == nullptr) return StatusCode::body_not_found;
        cel_cfg->active = active;
        return StatusCode::ok;
    } break;
    case BodyType::satellite: {
        ScenarioSatelliteConfig* sat_cfg = find_satellite_config_mut(scenario, id);
        if (sat_cfg == nullptr) return StatusCode::body_not_found;
        sat_cfg->active = active;
        return StatusCode::ok;
    } break;
    case BodyType::station: {
        ScenarioStationConfig* stat_cfg = find_station_config_mut(scenario, id);
        if (stat_cfg == nullptr) return StatusCode::body_not_found;
        stat_cfg->active = active;
        return StatusCode::ok;
    } break;
    case BodyType::unknown: return StatusCode::body_not_found;
    }

    return StatusCode::ok;
}

static StatusCode sync_scenario_x_tr(ScenarioStateTrConfig& cfg, const StateTr& x_tr) {
    cfg.r = x_tr.r;
    cfg.v = x_tr.v;
    return StatusCode::ok;
}
static StatusCode sync_scenario_x_att(
    ScenarioStateAttConfig& cfg,
    const StateAtt& x_att
) {
    cfg.q = x_att.q;
    cfg.w = x_att.w;
    return StatusCode::ok;
}
static StatusCode sync_scenario_body_state(
    ScenarioSession& scenario,
    const Body& body,
    const World& world
) {
    ScenarioConfig& cfg = scenario.config;
    ScenarioBuildResult& build_result = scenario.build_result;

    switch (body.body_type) {
    case BodyType::unknown: return StatusCode::body_not_found;
    case BodyType::celestial: {
        auto* cfg = find_celestial_config_mut(scenario, body.id);
        if (cfg == nullptr) return StatusCode::body_not_found;
        cfg->x_tr.r = body.x_tr.r;
        cfg->x_tr.v = body.x_tr.v;
        cfg->x_att.q = body.x_att.q;
        cfg->x_att.w = body.x_att.w;
    } break;
    case BodyType::satellite: {
        auto* body_cfg = find_satellite_config_mut(scenario, body.id);
        if (body_cfg == nullptr) return StatusCode::body_not_found;
        body_cfg->x_tr.r = body.x_tr.r;
        body_cfg->x_tr.v = body.x_tr.v;
        body_cfg->x_att.q = body.x_att.q;
        body_cfg->x_att.w = body.x_att.w;
    } break;
    case BodyType::station: {
        auto* stat = world.station(body.id);
        if (stat == nullptr) return StatusCode::body_not_found;
        auto* body_cfg = find_station_config_mut(scenario, body.id);
        if (body_cfg == nullptr) return StatusCode::body_not_found;
        if (stat->anchored) {
            auto it = scenario.build_result.celestial_config_ids.find(stat->anchor_id);
            if (it == scenario.build_result.celestial_config_ids.end())
                return StatusCode::anchor_not_found;
            body_cfg->anchor = it->second;
            body_cfg->llh = stat->llh_BCBF;
            body_cfg->r_body = stat->r_body_BCBF;
        } else {
            body_cfg->x_tr.r = body.x_tr.r;
            body_cfg->x_tr.v = body.x_tr.v;
            body_cfg->x_att.q = body.x_att.q;
            body_cfg->x_att.w = body.x_att.w;
        }
    } break;
    }
    return StatusCode::ok;
}

static StatusCode sync_scenario_mass_properties(
    ScenarioMassPropertiesConfig& cfg,
    const MassProperties& mp
) {
    cfg.mass = mp.mass;
    cfg.inertia = mp.I;
    cfg.offset_body = mp.offset_body;
    cfg.principle_axes = mp.principal_axes;
    return StatusCode::ok;
}

static StatusCode sync_scenario_gravity_model(
    ScenarioGravityConfig& cfg,
    const Celestial& cel
) {
    cfg.model = cel.gravity_model;
    cfg.mu = cel.mu;
    cfg.radius = cel.ref_radius;
    cfg.degree = cel.degree;
    cfg.order = cel.order;
    cfg.J = cel.J;
    cfg.C = cel.C;
    cfg.S = cel.S;
    return StatusCode::ok;
}
static StatusCode sync_scenario_celestial_model(
    ScenarioCelestialModelConfig& cfg,
    const Celestial& cel
) {
    cfg.semimajor_axis = cel.semimajor_axis;
    cfg.semiminor_axis = cel.semiminor_axis;
    cfg.mean_radius = cel.mean_radius;
    cfg.eccentricity = cel.eccentricity;
    return StatusCode::ok;
}

static StatusCode sync_scenario_celestial(
    ScenarioSession& scenario,
    const Celestial& cel
) {
    ScenarioConfig& cfg = scenario.config;
    ScenarioBuildResult& build_result = scenario.build_result;

    auto* cel_cfg = find_celestial_config_mut(scenario, cel.id);

    StatusCode status;

    status = sync_scenario_x_tr(cel_cfg->x_tr, cel.x_tr);
    if (status != StatusCode::ok) return status;

    cel_cfg->attitude_model = cel.attitude_model;
    status = sync_scenario_x_att(cel_cfg->x_att, cel.x_att);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_celestial_model(cel_cfg->model, cel);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_gravity_model(cel_cfg->model.gravity_model, cel);
    if (status != StatusCode::ok) return status;

    cel_cfg->propagation.translation = cel.propagate_tr;
    cel_cfg->propagation.attitude = cel.propagate_att;

    return StatusCode::ok;
}

static StatusCode sync_scenario_add_celestial(
    ScenarioSession& scenario,
    const Celestial& cel
) {
    ScenarioCelestialConfig cfg;
    StatusCode status;

    status = sync_scenario_celestial_model(cfg.model, cel);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_gravity_model(cfg.model.gravity_model, cel);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_x_tr(cfg.x_tr, cel.x_tr);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_x_att(cfg.x_att, cel.x_att);
    if (status != StatusCode::ok) return status;

    cfg.propagation.translation = cel.propagate_tr;
    cfg.propagation.attitude = cel.propagate_att;

    scenario.config.celestials.push_back(cfg);
    cfg.name = cel.name;
    string id = cel.name.empty() ? "celestial " + std::to_string(cel.id) : cel.name;
    scenario.build_result.celestial_config_ids.at(cel.id) = id;
    scenario.build_result.celestial_ids.at(id) = cel.id;
    scenario.build_result.body_config_ids.at(cel.id) = id;
    scenario.build_result.body_ids.at(id) = cel.id;

    return StatusCode::ok;
}

static StatusCode sync_scenario_satellite(
    ScenarioSession& scenario,
    const Satellite& sat
) {
    ScenarioConfig& cfg = scenario.config;
    ScenarioBuildResult& build_result = scenario.build_result;

    auto* sat_cfg = find_satellite_config_mut(scenario, sat.id);

    StatusCode status;

    status = sync_scenario_x_tr(sat_cfg->x_tr, sat.x_tr);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_x_att(sat_cfg->x_att, sat.x_att);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_mass_properties(sat_cfg->mass_properties, sat.mass_properties);
    if (status != StatusCode::ok) return status;

    sat_cfg->propagation.translation = sat.propagate_tr;
    sat_cfg->propagation.attitude = sat.propagate_att;

    return StatusCode::ok;
}

static StatusCode sync_scenario_add_satellite(
    ScenarioSession& scenario,
    const Satellite& sat
) {
    ScenarioSatelliteConfig cfg;
    StatusCode status;

    status = sync_scenario_x_tr(cfg.x_tr, sat.x_tr);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_x_att(cfg.x_att, sat.x_att);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_mass_properties(cfg.mass_properties, sat.mass_properties);
    if (status != StatusCode::ok) return status;

    cfg.propagation.translation = sat.propagate_tr;
    cfg.propagation.attitude = sat.propagate_att;

    scenario.config.satellites.push_back(cfg);
    cfg.name = sat.name;
    string id = sat.name.empty() ? "satellite " + std::to_string(sat.id) : sat.name;
    scenario.build_result.satellite_config_ids.at(sat.id) = id;
    scenario.build_result.satellite_ids.at(id) = sat.id;
    scenario.build_result.body_config_ids.at(sat.id) = id;
    scenario.build_result.body_ids.at(id) = sat.id;

    return StatusCode::ok;
}

static StatusCode sync_scenario_station_anchor(
    ScenarioStationConfig& cfg,
    const Station& stat,
    ScenarioBuildResult& build_result
) {
    auto it = build_result.celestial_config_ids.find(stat.anchor_id);
    if (it == build_result.celestial_config_ids.end())
        return StatusCode::anchor_not_found;

    cfg.anchor = it->second;
    cfg.llh = stat.llh_BCBF;
    cfg.r_body = stat.r_body_BCBF;
    cfg.local_frame = "ENU"; // TODO: make enum for this?
    cfg.units_angle = UAngle::radian;
    cfg.units_length = ULength::kilometer;

    return StatusCode::ok;
}

static StatusCode sync_scenario_instrument(
    ScenarioInstrumentConfig& cfg,
    const StationInstrument& instr
) {
    cfg.covariance_cfg.covariance = instr.R;
    cfg.enabled = instr.enabled;
    cfg.type = instr.type;

    return StatusCode::ok;
}
static StatusCode sync_scenario_add_station(
    ScenarioSession& scenario,
    const Station& stat
) {
    ScenarioStationConfig cfg;
    StatusCode status;

    cfg.anchored = stat.anchored;
    if (stat.anchored) {
        status = sync_scenario_station_anchor(cfg, stat, scenario.build_result);
        if (status != StatusCode::ok) return status;
    } else {
        status = sync_scenario_x_tr(cfg.x_tr, stat.x_tr);
        if (status != StatusCode::ok) return status;

        status = sync_scenario_x_att(cfg.x_att, stat.x_att);
        if (status != StatusCode::ok) return status;

        status = sync_scenario_mass_properties(cfg.mass_properties, stat.mass_properties);
        if (status != StatusCode::ok) return status;
    }
    cfg.propagation.translation = stat.propagate_tr;
    cfg.propagation.attitude = stat.propagate_att;

    for (const auto& [id, instr] : stat.instruments) {
        ScenarioInstrumentConfig instrument;
        sync_scenario_instrument(instrument, instr);
        instrument.id = "Instrument " + std::to_string(instr.id) + " ("
                        + observation_type_str_simple(instr.type) + ")";
        cfg.instruments.push_back(instrument);
    }

    scenario.config.stations.push_back(cfg);
    cfg.name = stat.name;
    string id = stat.name.empty() ? "station" + std::to_string(stat.id) : stat.name;
    scenario.build_result.station_config_ids.at(stat.id) = id;
    scenario.build_result.station_ids.at(id) = stat.id;
    scenario.build_result.body_config_ids.at(stat.id) = id;
    scenario.build_result.body_ids.at(id) = stat.id;

    return StatusCode::ok;
}

static void sync_scenario_runtime_state_from_world(
    RenderLoopState& state,
    const World& world
) {
    if (state.scenario.dirty) {}
    state.scenario.dirty = false;
}

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

static ImVec4 rl_to_im(const Vector4& v) { return ImVec4{v.x, v.y, v.z, v.w}; }
static ImVec4 rl_to_im(const Color& c) {
    f32 r = static_cast<f32>(c.r) / 255.0f;
    f32 g = static_cast<f32>(c.g) / 255.0f;
    f32 b = static_cast<f32>(c.b) / 255.0f;
    f32 a = static_cast<f32>(c.a) / 255.0f;
    return ImVec4{r, g, b, a};
}

static void status_text(const char* label, const StatusCode code) {
    im::TextColored(
        rl_to_im(status_color(code)),
        "%s: %s",
        label,
        status_string(code).c_str()
    );
}

static string normalize_file_dialog_path(const char* path_in) {
    string path = path_in;
    std::replace(path.begin(), path.end(), '\\', '/');
    return path;
}

static std::filesystem::path resolve_textbox_path(
    const string& filepath,
    bool relative_path
) {
    std::filesystem::path path = filepath;

    if (relative_path) {
        return std::filesystem::path(pwd) / path.relative_path();
    }

    if (path.is_absolute()) return path;
    return std::filesystem::path("/") / path.relative_path();
}

static string scenario_dialog_start_path(const RenderLoopState& state) {
    std::filesystem::path path;

    if (state.filepath.empty()) {
        path = std::filesystem::path(PROJECT_ROOT) / "scenarios";
    } else {
        path = resolve_textbox_path(state.filepath, state.relative_path);
    }

    if (std::filesystem::is_regular_file(path)) {
        path = path.parent_path();
    } else if (!std::filesystem::exists(path) && !path.parent_path().empty()) {
        path = path.parent_path();
    }

    if (path.empty()) {
        path = std::filesystem::current_path();
    }

    return normalize_file_dialog_path(path.string().c_str());
}

static string scenario_file_path(const RenderLoopState& state) {
    std::filesystem::path path
        = resolve_textbox_path(state.filepath, state.relative_path);
    return normalize_file_dialog_path(path.string().c_str());
}

static bool path_exists(const string& filepath) {
    std::error_code ec;
    return std::filesystem::exists(filepath, ec);
}

static bool is_existing_regular_file(const string& filepath) {
    std::error_code ec;
    return std::filesystem::exists(filepath, ec)
           && std::filesystem::is_regular_file(filepath, ec);
}

static bool is_existing_non_regular_file(const string& filepath) {
    std::error_code ec;
    return std::filesystem::exists(filepath, ec)
           && !std::filesystem::is_regular_file(filepath, ec);
}

static void set_scenario_dialog_style() {
    imfd::settings.displayMode = imfd::GlobalSettings::DisplayMode_List;
    imfd::settings.asciiArtIcons = false;
}

static void open_scenario_file_dialog(RenderLoopState& state) {
    const string path = scenario_dialog_start_path(state);
    set_scenario_dialog_style();
    imfd::OpenDialog("Open Scenario", ImGuiFDMode_LoadFile, path.c_str(), "*.json", 0, 1);
}

static void open_scenario_save_dialog(RenderLoopState& state) {
    const string path = scenario_dialog_start_path(state);
    set_scenario_dialog_style();
    imfd::OpenDialog("Save Scenario", ImGuiFDMode_SaveFile, path.c_str(), "*.json", 0, 1);
}

static void open_folder_dialog(RenderLoopState& state) {
    const string path = scenario_dialog_start_path(state);
    set_scenario_dialog_style();
    imfd::OpenDialog("Open Folder", ImGuiFDMode_OpenDir, path.c_str(), NULL, 0, 1);
}

static void render_open_scenario_file_dialog(RenderLoopState& state) {
    set_scenario_dialog_style();
    if (imfd::BeginDialog("Open Scenario")) {
        if (imfd::ActionDone()) {
            if (imfd::SelectionMade()) {
                state.filepath
                    = normalize_file_dialog_path(imfd::GetSelectionPathString(0));
                state.relative_path = false;
                state.file_attempt = false;
            }
            imfd::CloseCurrentDialog();
        }
        imfd::EndDialog();
    }
}

static void render_save_scenario_file_dialog(
    RenderLoopState& state,
    const RenderLoopConfig& cfg,
    const World& world
) {
    set_scenario_dialog_style();
    if (imfd::BeginDialog("Save Scenario")) {
        if (imfd::ActionDone()) {
            if (imfd::SelectionMade()) {
                state.filepath
                    = normalize_file_dialog_path(imfd::GetSelectionPathString(0));
                state.relative_path = false;
                state.save_filepath = state.filepath;
                state.file_attempt = true;

                if (is_existing_non_regular_file(state.save_filepath)) {
                    state.file_status = StatusCode::invalid_input;
                } else if (is_existing_regular_file(state.save_filepath)) {
                    state.file_status = StatusCode::file_already_exists;
                    im::OpenPopup("File Warning");
                } else {
                    request_scenario_save(state, cfg, world);
                }
            }
            imfd::CloseCurrentDialog();
        }
        imfd::EndDialog();
    }
}

static void render_open_folder_dialog(RenderLoopState& state) {
    set_scenario_dialog_style();
    if (imfd::BeginDialog("Open Folder")) {
        if (imfd::ActionDone()) {
            if (imfd::SelectionMade()) {
                state.filepath
                    = normalize_file_dialog_path(imfd::GetSelectionPathString(0));
                state.relative_path = false;
                state.file_attempt = false;
            }
            imfd::CloseCurrentDialog();
        }
        imfd::EndDialog();
    }
}

static void request_scenario_save(
    RenderLoopState& state,
    const RenderLoopConfig& cfg,
    const World& world,
    bool overwrite
) {
    state.file_attempt = true;
    state.save_filepath = scenario_file_path(state);

    if (state.save_filepath.empty()) {
        state.file_status = StatusCode::invalid_input;
        return;
    }

    if (is_existing_non_regular_file(state.save_filepath)) {
        state.file_status = StatusCode::invalid_input;
        return;
    }

    if (!overwrite && is_existing_regular_file(state.save_filepath)) {
        state.file_status = StatusCode::file_already_exists;
        im::OpenPopup("File Warning");
        return;
    }

    // TODO: use as fallback?
    // state.file_status
    //     = build_scenario_config_from_world(state.scenario.config, world,
    //     cfg.stepper_cfg);
    // if (state.file_status != StatusCode::ok) return;

    state.file_status = save_scenario_json(state.save_filepath, state.scenario.config);
    if (state.file_status != StatusCode::ok) return;
}

static void render_scenario_overwrite_popup(
    RenderLoopState& state,
    const RenderLoopConfig& cfg,
    const World& world
) {
    if (im::BeginPopupModal("File Warning")) {
        im::Text("File already exists, overwrite?");
        if (im::Button("Yes")) {
            state.file_status = StatusCode::file_overwritten;
            im::CloseCurrentPopup();
            request_scenario_save(state, cfg, world, true);
        }

        im::SameLine();

        if (im::Button("Cancel")) {
            state.file_status = StatusCode::file_already_exists;
            im::CloseCurrentPopup();
        }

        im::EndPopup();
    }
}
// static void push_selected_field_style() {
//     im::PushStyleColor(ImGuiCol_Border, rl_to_im(YELLOW));
//     im::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
// }

// static void pop_selected_field_style() {
//     im::PopStyleVar();
//     im::PopStyleColor(1);
// }

static bool render_state_tr_ui(Body& body, ImGuiInputTextFlags flags) {
    return im::InputDouble3("r", body.x_tr.r, "%.3f", flags)
           || im::InputDouble3("v", body.x_tr.v, "%.3f", flags);
}

static bool render_state_att_ui(Body& body, ImGuiInputTextFlags flags) {
    return im::InputDouble4("q", body.x_att.q, "%.3f", flags)
           || im::InputDouble3("w", body.x_att.w, "%.3f", flags);
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

static void render_simulation_ui(
    World& world,
    RenderLoopConfig& cfg,
    RenderLoopState& state
) {
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

            bool finite_pos_dt = finite_nonneg(state.step_by_delta);
            if (finite_pos_dt) {
                // TODO: step here
                im::SameLine();
                im::Text("Not yet implemented");
            } else {
                // TODO: use statuscode
                im::SameLine();
                im::Text("The time change must be positive");
            }
        }
    }

    if (im::CollapsingHeader("Scenario")) {
        im::Indent();
        im::Text("Filepath:");
        im::SameLine();
        im::InputText("##filepath", &state.filepath);
        im::SameLine();
        if (im::Button("/")) {
            open_folder_dialog(state);
        }
        im::SameLine();
        if (im::Button("...")) {
            open_scenario_file_dialog(state);
        }
        render_open_scenario_file_dialog(state);
        render_save_scenario_file_dialog(state, cfg, world);
        render_open_folder_dialog(state);

        const std::filesystem::path path(state.filepath);
        const string filename = path.filename().string();
        im::Text("File: %s", filename.c_str());

        im::Checkbox("Relative filepath", &state.relative_path);

        if (im::Button("Save As")) {
            open_scenario_save_dialog(state);
        }
        render_scenario_overwrite_popup(state, cfg, world);

        im::SameLine();
        if (im::Button("Load")) {
            state.file_attempt = true;
            string path = scenario_file_path(state);
            state.load_filepath = path;

            if (!path_exists(path)) {
                state.file_status = StatusCode::file_not_found;
            } else if (!is_existing_regular_file(path)) {
                state.file_status = StatusCode::invalid_input;
            } else {
                state.file_status = load_scenario_json(path, state.scenario.config);
            }

            if (state.file_status == StatusCode::ok) {
                world = World{};
                state.file_status = build_world_from_scenario_config(
                    state.scenario.config,
                    world,
                    state.scenario.build_result,
                    cfg.stepper_cfg
                );
                if (state.file_status == StatusCode::ok) {
                    state.scenario.dirty = false;
                }
            }
        }
        im::SameLine();
        if (im::Button("Cancel")) {
            state.file_attempt = false;
            state.filepath = "";
            state.load_filepath = "";
            state.save_filepath = "";
        }
        if (state.file_attempt) {
            status_text("Scenario Loader/Saver Status", state.file_status);
        }
        im::Unindent();
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
            if (im::Button("Save")) {
                EntityId edit_body_id = state.draft.edit_body_id;
                state.draft.edit_body_status = apply_body_edit_draft(state.draft, world);
                if (state.draft.edit_body_status == StatusCode::ok) {
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
                render_celestial_stats_ui(*cel, field_flags);
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
                world.make_inactive(cfg.body_stats_id);
                sync_scenario_body_active(
                    state.scenario,
                    cfg.body_stats_id,
                    false,
                    world
                );
                state.wksp.dirty = true;
                if (cfg.camera.target_id == body->id) {
                    cycle_active_id(cfg.camera.target_id, world, 1);
                    sync_camera_tracking(cfg.camera, world);
                }
            }
            if (!active && im::Button("Make Active")) {
                world.make_active(cfg.body_stats_id);
                sync_scenario_body_active(state.scenario, cfg.body_stats_id, true, world);
                state.wksp.dirty = true;
            }
        }
    }

    if (cfg.stepper_cfg.paused) {
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
                        state.add_body_status = sync_scenario_add_celestial(
                            state.scenario,
                            state.temp_celestial
                        );
                        if (state.add_body_status != StatusCode::ok) break;
                        id = world.insert_celestial(std::move(cel));
                        state.temp_celestial = Celestial{};
                    } break;
                    case BodyType::satellite: {
                        auto sat = std::make_unique<Satellite>(state.temp_satellite);
                        state.add_body_status = sync_scenario_add_satellite(
                            state.scenario,
                            state.temp_satellite
                        );
                        if (state.add_body_status != StatusCode::ok) break;
                        id = world.insert_satellite(std::move(sat));
                        state.temp_satellite = Satellite{};
                    } break;
                    case BodyType::station: {
                        auto stat = std::make_unique<Station>(state.temp_station);
                        state.add_body_status = sync_scenario_add_station(
                            state.scenario,
                            state.temp_station
                        );
                        if (state.add_body_status != StatusCode::ok) break;
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
                status_text("Add Body Status", state.add_body_status);
            }
        }
    } else {
        im::SetWindowCollapsed(true);
    }
    im::End();
}

static bool render_body_propagation_ui(Body& body, bool editable) {
    if (!editable) im::BeginDisabled();
    bool changed = im::Checkbox("Propagate Translation", &body.propagate_tr)
                   || im::Checkbox("Propagate Attitude", &body.propagate_att);
    if (!editable) im::EndDisabled();

    return changed;
}

static bool render_stat_state_ui(Station& stat, World& world, ImGuiInputTextFlags flags) {
    flags |= ImGuiInputTextFlags_ReadOnly;
    vec3d r = world.stat_r_inertial(stat.id);
    vec4d q = world.stat_q_inertial(stat.id);
    return im::InputDouble3("r", r, "%.3f", flags)
           || im::InputDouble4("q", q, "%.3f", flags);
}

static void render_celestial_stats_ui(Celestial& cel, ImGuiInputTextFlags flags) {
    bool editable = !fields_readonly(flags);

    render_state_ui(cel, flags);

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
    }

    if (im::CollapsingHeader("Gravity Model", ImGuiTreeNodeFlags_DefaultOpen)) {
        im::Indent();
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
        im::Unindent();
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
            ImGuiInputTextFlags spin_flags = flags | ImGuiInputTextFlags_ReadOnly;
            f64 w_mag = cel.x_att.w.norm();
            im::InputDouble("Spin Rate", &w_mag, 0.0, 0.0, "%.3f", spin_flags);
            im::InputDouble3("Spin Axis", cel.x_att.w, "%.3f", spin_flags);
        }
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
        return celestial_changed(*cel, draft.edit_celestial);
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
    StatusCode status = validate_body_edit_draft(draft, world);
    if (status != StatusCode::ok) return status;

    switch (draft.edit_body_type) {
    case BodyType::unknown: return StatusCode::unsupported_type;
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
