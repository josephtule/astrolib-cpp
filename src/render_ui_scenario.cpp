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
    const string& id_str = it->second;
    return find_celestial_config(scenario.config, id_str);
}
static ScenarioCelestialConfig* find_celestial_config_mut(
    ScenarioSession& scenario,
    const EntityId id
) {
    auto it = scenario.build_result.celestial_config_ids.find(id);
    if (it == scenario.build_result.celestial_config_ids.end()) return nullptr;
    const string& id_str = it->second;
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
    const string& id_str = it->second;
    return find_satellite_config(scenario.config, id_str);
}
static ScenarioSatelliteConfig* find_satellite_config_mut(
    ScenarioSession& scenario,
    const EntityId id
) {
    auto it = scenario.build_result.satellite_config_ids.find(id);
    if (it == scenario.build_result.satellite_config_ids.end()) return nullptr;
    const string& id_str = it->second;
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
    const string& id_str = it->second;
    return find_station_config(scenario.config, id_str);
}
static ScenarioStationConfig* find_station_config_mut(
    ScenarioSession& scenario,
    const EntityId id
) {
    auto it = scenario.build_result.station_config_ids.find(id);
    if (it == scenario.build_result.station_config_ids.end()) return nullptr;
    const string& id_str = it->second;
    return find_station_config_mut(scenario.config, id_str);
}

static bool scenario_body_config_id_exists(const ScenarioConfig& cfg, const string& id) {
    return find_celestial_config(cfg, id) != nullptr
           || find_satellite_config(cfg, id) != nullptr
           || find_station_config(cfg, id) != nullptr;
}

static const ScenarioGravityProviderConfig* find_gravity_provider_config(
    const ScenarioConfig& cfg,
    const string& id
) {
    for (const auto& provider : cfg.gravity_providers) {
        if (provider.id == id) return &provider;
    }
    return nullptr;
}

string make_unique_scenario_body_id(
    const ScenarioSession& session,
    BodyType type,
    EntityId runtime_id
) {
    string prefix;

    switch (type) {
    case BodyType::celestial: prefix = "celestial"; break;
    case BodyType::satellite: prefix = "satellite"; break;
    case BodyType::station: prefix = "station"; break;
    case BodyType::unknown: return "";
    }

    const string base = prefix + "_" + std::to_string(runtime_id);
    string candidate = base;
    i32 suffix = 2;

    while (session.build_result.body_ids.contains(candidate)
           || scenario_body_config_id_exists(session.config, candidate)) {
        candidate = base + "_" + std::to_string(suffix);
        ++suffix;
    }

    return candidate;
}

StatusCode register_scenario_body_mapping(
    ScenarioBuildResult& result,
    BodyType type,
    const string& config_id,
    EntityId runtime_id
) {
    if (type == BodyType::unknown) return StatusCode::unsupported_type;
    if (runtime_id == kInvalidEntityId) return StatusCode::body_not_found;
    if (config_id.empty()) return StatusCode::body_not_found;

    const auto it1 = result.body_ids.find(config_id);
    if (it1 != result.body_ids.end()) return StatusCode::body_not_found;
    const auto it2 = result.body_config_ids.find(runtime_id);
    if (it2 != result.body_config_ids.end()) return StatusCode::body_not_found;

    switch (type) {

    case BodyType::unknown: return StatusCode::unsupported_type;
    case BodyType::celestial: {
        const auto it1 = result.celestial_ids.find(config_id);
        if (it1 != result.celestial_ids.end()) return StatusCode::body_not_found;
        const auto it2 = result.celestial_config_ids.find(runtime_id);
        if (it2 != result.celestial_config_ids.end()) return StatusCode::body_not_found;
    } break;
    case BodyType::satellite: {
        const auto it1 = result.satellite_ids.find(config_id);
        if (it1 != result.satellite_ids.end()) return StatusCode::body_not_found;
        const auto it2 = result.satellite_config_ids.find(runtime_id);
        if (it2 != result.satellite_config_ids.end()) return StatusCode::body_not_found;
    } break;
    case BodyType::station: {
        const auto it1 = result.station_ids.find(config_id);
        if (it1 != result.station_ids.end()) return StatusCode::body_not_found;
        const auto it2 = result.station_config_ids.find(runtime_id);
        if (it2 != result.station_config_ids.end()) return StatusCode::body_not_found;
    } break;
    }

    result.body_ids[config_id] = runtime_id;
    result.body_config_ids[runtime_id] = config_id;

    switch (type) {
    case BodyType::unknown: return StatusCode::body_not_found;
    case BodyType::celestial: {
        result.celestial_ids[config_id] = runtime_id;
        result.celestial_config_ids[runtime_id] = config_id;
    } break;
    case BodyType::satellite: {
        result.satellite_ids[config_id] = runtime_id;
        result.satellite_config_ids[runtime_id] = config_id;
    } break;
    case BodyType::station: {
        result.station_ids[config_id] = runtime_id;
        result.station_config_ids[runtime_id] = config_id;
    } break;
    }

    return StatusCode::ok;
}

StatusCode sync_scenario_body_active(
    ScenarioSession& scenario,
    const EntityId id,
    const bool active,
    const World& world
) {
    const auto* body = world.body(id);
    if (body == nullptr) return StatusCode::body_not_found;
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
    case BodyType::unknown: return StatusCode::unsupported_type;
    }

    return StatusCode::ok;
}

static StatusCode sync_scenario_x_tr(ScenarioStateTrConfig& cfg, const StateTr& x_tr) {
    cfg = ScenarioStateTrConfig{}; // clear out unnecessary fields
    cfg.input_type = StateTrInputType::pos_vel;
    cfg.units_length = ULength::kilometer;
    cfg.r = x_tr.r;
    cfg.v = x_tr.v;
    return StatusCode::ok;
}
static StatusCode sync_scenario_x_att(
    ScenarioStateAttConfig& cfg,
    const StateAtt& x_att
) {
    cfg = ScenarioStateAttConfig{}; // clear out unnecessary fields
    cfg.input_type = AttitudeType::quaternion;
    cfg.units_angle = UAngle::radian;
    cfg.q = x_att.q;
    cfg.w = x_att.w;
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
    cfg.offset = !mp.offset_body.isZero(tol12);
    return StatusCode::ok;
}

static StatusCode sync_scenario_gravity_model(
    ScenarioConfig& scenario,
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

    cfg.coefficients.clear();
    switch (cel.gravity_model) {
    case GravityModel::pointmass: {
        cfg.coefficient_source = GravityCoefficientSource::none;
    } break;
    case GravityModel::zonal:
    case GravityModel::spherical_harmonics: {
        if (cel.gravity_provider.empty()) {
            cfg.coefficient_source
                = cel.gravity_model == GravityModel::zonal
                      ? GravityCoefficientSource::direct_zonal
                      : GravityCoefficientSource::direct_spherical_harmonics;
            break;
        }

        if (find_gravity_provider_config(scenario, cel.gravity_provider) == nullptr) {
            if (cel.gravity_provider_format.empty()
                || cel.gravity_provider_filepath.empty()) {
                return StatusCode::gravity_model_not_found;
            }
            if (scenario_body_config_id_exists(scenario, cel.gravity_provider)) {
                return StatusCode::duplicate_id;
            }

            ScenarioGravityProviderConfig provider;
            provider.id = cel.gravity_provider;
            provider.format = cel.gravity_provider_format;
            provider.filepath = cel.gravity_provider_filepath;
            provider.lineskips = cel.gravity_provider_lineskips;
            provider.normalized = cel.gravity_provider_normalized;
            scenario.gravity_providers.push_back(std::move(provider));
        }

        cfg.coefficient_source = GravityCoefficientSource::provider;
        cfg.coefficients = cel.gravity_provider;
    } break;
    }

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
    cfg.flattening = cel.flattening;
    cfg.units_length = ULength::kilometer;
    return StatusCode::ok;
}

static StatusCode sync_scenario_celestial(
    ScenarioSession& scenario,
    const Celestial& cel,
    bool active
) {
    auto* cel_cfg = find_celestial_config_mut(scenario, cel.id);
    if (cel_cfg == nullptr) return StatusCode::body_not_found;

    StatusCode status;
    cel_cfg->name = cel.name;

    status = sync_scenario_celestial_model(cel_cfg->model, cel);
    if (status != StatusCode::ok) return status;

    status
        = sync_scenario_gravity_model(scenario.config, cel_cfg->model.gravity_model, cel);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_x_tr(cel_cfg->x_tr, cel.x_tr);
    if (status != StatusCode::ok) return status;

    cel_cfg->attitude_model = cel.attitude_model;
    status = sync_scenario_x_att(cel_cfg->x_att, cel.x_att);
    if (status != StatusCode::ok) return status;

    cel_cfg->propagation.translation = cel.propagate_tr;
    cel_cfg->propagation.attitude = cel.propagate_att;

    cel_cfg->active = active;

    return StatusCode::ok;
}

static StatusCode sync_scenario_satellite(
    ScenarioSession& scenario,
    const Satellite& sat,
    bool active
) {
    auto* sat_cfg = find_satellite_config_mut(scenario, sat.id);
    if (sat_cfg == nullptr) return StatusCode::body_not_found;

    StatusCode status;
    sat_cfg->name = sat.name;

    status = sync_scenario_x_tr(sat_cfg->x_tr, sat.x_tr);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_x_att(sat_cfg->x_att, sat.x_att);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_mass_properties(sat_cfg->mass_properties, sat.mass_properties);
    if (status != StatusCode::ok) return status;

    sat_cfg->propagation.translation = sat.propagate_tr;
    sat_cfg->propagation.attitude = sat.propagate_att;
    sat_cfg->active = active;

    return StatusCode::ok;
}

StatusCode make_scenario_celestial_config(
    ScenarioConfig& scenario,
    const Celestial& cel,
    bool active,
    ScenarioCelestialConfig& out
) {
    out = ScenarioCelestialConfig{};
    StatusCode status;

    out.name = cel.name;

    status = sync_scenario_x_tr(out.x_tr, cel.x_tr);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_x_att(out.x_att, cel.x_att);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_celestial_model(out.model, cel);
    if (status != StatusCode::ok) return status;

    out.model.id = "custom";

    status = sync_scenario_gravity_model(scenario, out.model.gravity_model, cel);
    if (status != StatusCode::ok) return status;

    out.attitude_model = cel.attitude_model;
    out.has_attitude_model = true;

    out.propagation.translation = cel.propagate_tr;
    out.propagation.attitude = cel.propagate_att;

    out.active = active;

    return StatusCode::ok;
}

StatusCode make_scenario_satellite_config(
    const Satellite& sat,
    bool active,
    ScenarioSatelliteConfig& out
) {
    out = ScenarioSatelliteConfig{};
    StatusCode status;

    status = sync_scenario_x_tr(out.x_tr, sat.x_tr);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_x_att(out.x_att, sat.x_att);
    if (status != StatusCode::ok) return status;

    status = sync_scenario_mass_properties(out.mass_properties, sat.mass_properties);
    if (status != StatusCode::ok) return status;

    out.name = sat.name;

    out.propagation.translation = sat.propagate_tr;
    out.propagation.attitude = sat.propagate_att;

    out.active = active;

    return StatusCode::ok;
}

static StatusCode sync_scenario_station_anchor(
    ScenarioStationConfig& cfg,
    const Station& stat,
    const ScenarioBuildResult& build_result
) {
    auto it = build_result.celestial_config_ids.find(stat.anchor_id);
    if (it == build_result.celestial_config_ids.end())
        return StatusCode::anchor_not_found;

    cfg.anchor = it->second;
    cfg.llh = stat.llh_BCBF;
    cfg.r_body = stat.r_body_BCBF;
    cfg.coordinate_type = "detic_llh";
    cfg.local_frame = "ENU"; // TODO: make enum for this?
    cfg.units_angle = UAngle::radian;
    cfg.units_length = ULength::kilometer;

    return StatusCode::ok;
}

static StatusCode sync_scenario_instrument(
    ScenarioInstrumentConfig& cfg,
    const StationInstrument& instr
) {
    cfg.id = instr.name.empty() ? "instrument_" + std::to_string(instr.id) : instr.name;

    cfg.covariance_cfg.type = "matrix";

    cfg.covariance_cfg.covariance = instr.R;
    cfg.enabled = instr.enabled;
    cfg.type = instr.type;

    return StatusCode::ok;
}

static StatusCode sync_scenario_station_instruments(
    ScenarioStationConfig& cfg,
    const Station& stat
) {
    svec<InstrumentId> instrument_ids;
    instrument_ids.reserve(stat.instruments.size());
    for (const auto& entry : stat.instruments) {
        instrument_ids.push_back(entry.first);
    }
    std::sort(instrument_ids.begin(), instrument_ids.end());

    svec<ScenarioInstrumentConfig> instruments;
    instruments.reserve(instrument_ids.size());
    uset<string> config_ids;

    for (InstrumentId id : instrument_ids) {
        auto it = stat.instruments.find(id);
        if (it == stat.instruments.end()) return StatusCode::instrument_not_found;

        ScenarioInstrumentConfig instr_cfg;
        StatusCode status = sync_scenario_instrument(instr_cfg, it->second);
        if (status != StatusCode::ok) return status;
        if (!config_ids.insert(instr_cfg.id).second) return StatusCode::duplicate_id;
        instruments.push_back(std::move(instr_cfg));
    }

    cfg.instruments = std::move(instruments);
    return StatusCode::ok;
}

static StatusCode sync_scenario_station(
    ScenarioSession& session,
    const Station& stat,
    bool active
) {
    auto* stat_cfg = find_station_config_mut(session, stat.id);
    if (stat_cfg == nullptr) return StatusCode::body_not_found;

    StatusCode status;
    stat_cfg->name = stat.name;

    stat_cfg->anchored = stat.anchored;
    if (stat.anchored) {
        status = sync_scenario_station_anchor(*stat_cfg, stat, session.build_result);
        if (status != StatusCode::ok) return status;
    } else {
        status = sync_scenario_x_tr(stat_cfg->x_tr, stat.x_tr);
        if (status != StatusCode::ok) return status;

        status = sync_scenario_x_att(stat_cfg->x_att, stat.x_att);
        if (status != StatusCode::ok) return status;

        status = sync_scenario_mass_properties(
            stat_cfg->mass_properties,
            stat.mass_properties
        );
        if (status != StatusCode::ok) return status;
    }

    status = sync_scenario_station_instruments(*stat_cfg, stat);
    if (status != StatusCode::ok) return status;

    stat_cfg->propagation.translation = stat.propagate_tr;
    stat_cfg->propagation.attitude = stat.propagate_att;
    stat_cfg->active = active;

    return StatusCode::ok;
}

StatusCode sync_scenario_body(
    ScenarioSession& session,
    const Body& body,
    const World& world
) {
    const bool active = world.is_active(body.id);

    switch (body.body_type) {
    case BodyType::unknown: return StatusCode::unsupported_type;
    case BodyType::celestial: {
        const Celestial* cel = world.celestial(body.id);
        if (cel == nullptr) return StatusCode::body_not_found;
        return sync_scenario_celestial(session, *cel, active);
    }
    case BodyType::satellite: {
        const Satellite* sat = world.satellite(body.id);
        if (sat == nullptr) return StatusCode::body_not_found;
        return sync_scenario_satellite(session, *sat, active);
    }
    case BodyType::station: {
        const Station* stat = world.station(body.id);
        if (stat == nullptr) return StatusCode::body_not_found;
        return sync_scenario_station(session, *stat, active);
    }
    }

    return StatusCode::unsupported_type;
}
StatusCode make_scenario_station_config(
    const Station& stat,
    bool active,
    const ScenarioBuildResult& mappings,
    ScenarioStationConfig& out
) {
    out = ScenarioStationConfig{};

    StatusCode status;
    out.name = stat.name;

    out.anchored = stat.anchored;
    if (stat.anchored) {
        status = sync_scenario_station_anchor(out, stat, mappings);
        if (status != StatusCode::ok) return status;
    } else {
        status = sync_scenario_x_tr(out.x_tr, stat.x_tr);
        if (status != StatusCode::ok) return status;

        status = sync_scenario_x_att(out.x_att, stat.x_att);
        if (status != StatusCode::ok) return status;

        status = sync_scenario_mass_properties(out.mass_properties, stat.mass_properties);
        if (status != StatusCode::ok) return status;
    }

    status = sync_scenario_station_instruments(out, stat);
    if (status != StatusCode::ok) return status;

    out.propagation.translation = stat.propagate_tr;
    out.propagation.attitude = stat.propagate_att;
    out.active = active;

    return StatusCode::ok;
}

static StatusCode sync_scenario_runtime_state_from_world(
    ScenarioSession& session,
    const World& world,
    const WorldStepperConfig& stepper
) {
    StatusCode status;
    ScenarioConfig& cfg = session.config;
    ScenarioBuildResult& mappings = session.build_result;

    cfg.time.t0 = world.t_sim();

    cfg.world_stepper.dt_scale = stepper.dt_scale;
    cfg.world_stepper.substeps = stepper.substeps;
    cfg.world_stepper.ticks = stepper.ticks;
    cfg.world_stepper.paused = stepper.paused;
    cfg.world_stepper.integrator_tr = stepper.integrator_tr;
    cfg.world_stepper.integrator_att = stepper.integrator_att;

    for (const EntityId id : world.all_celestial_ids()) {
        const Celestial* cel = world.celestial(id);
        if (cel == nullptr) return StatusCode::body_not_found;

        const auto mapping = mappings.celestial_config_ids.find(id);
        if (mapping != mappings.celestial_config_ids.end()) {
            // sync already existing
            status = sync_scenario_celestial(session, *cel, world.is_active(id));
            if (status != StatusCode::ok) return status;
            continue;
        }

        // add new
        ScenarioCelestialConfig cel_cfg;
        status = make_scenario_celestial_config(cfg, *cel, world.is_active(id), cel_cfg);
        if (status != StatusCode::ok) return status;

        cel_cfg.id = make_unique_scenario_body_id(session, BodyType::celestial, id);
        if (cel_cfg.id.empty()) return StatusCode::invalid_input;

        status = register_scenario_body_mapping(
            mappings,
            BodyType::celestial,
            cel_cfg.id,
            id
        );
        if (status != StatusCode::ok) return status;

        cfg.celestials.push_back(std::move(cel_cfg));
    }

    for (const EntityId id : world.all_satellite_ids()) {
        const Satellite* sat = world.satellite(id);
        if (sat == nullptr) return StatusCode::body_not_found;

        const auto mapping = mappings.satellite_config_ids.find(id);
        if (mapping != mappings.satellite_config_ids.end()) {
            status = sync_scenario_satellite(session, *sat, world.is_active(id));
            if (status != StatusCode::ok) return status;
            continue;
        }

        ScenarioSatelliteConfig sat_cfg;
        status = make_scenario_satellite_config(*sat, world.is_active(id), sat_cfg);
        if (status != StatusCode::ok) return status;

        sat_cfg.id = make_unique_scenario_body_id(session, BodyType::satellite, id);
        if (sat_cfg.id.empty()) return StatusCode::invalid_input;

        status = register_scenario_body_mapping(
            mappings,
            BodyType::satellite,
            sat_cfg.id,
            id
        );
        if (status != StatusCode::ok) return status;

        cfg.satellites.push_back(std::move(sat_cfg));
    }

    for (const EntityId id : world.all_station_ids()) {
        const Station* stat = world.station(id);
        if (stat == nullptr) return StatusCode::body_not_found;

        const auto mapping = mappings.station_config_ids.find(id);
        if (mapping != mappings.station_config_ids.end()) {
            status = sync_scenario_station(session, *stat, world.is_active(id));
            if (status != StatusCode::ok) return status;
            continue;
        }

        ScenarioStationConfig stat_cfg;
        status = make_scenario_station_config(
            *stat,
            world.is_active(id),
            mappings,
            stat_cfg
        );
        if (status != StatusCode::ok) return status;

        stat_cfg.id = make_unique_scenario_body_id(session, BodyType::station, id);
        if (stat_cfg.id.empty()) return StatusCode::invalid_input;

        status = register_scenario_body_mapping(
            mappings,
            BodyType::station,
            stat_cfg.id,
            id
        );
        if (status != StatusCode::ok) return status;

        cfg.stations.push_back(std::move(stat_cfg));
    }

    return StatusCode::ok;
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

    state.file_status
        = sync_scenario_runtime_state_from_world(state.scenario, world, cfg.stepper_cfg);
    if (state.file_status != StatusCode::ok) return;

    state.file_status = save_scenario_json(state.save_filepath, state.scenario.config);
    if (state.file_status != StatusCode::ok) return;

    state.scenario.filepath = state.save_filepath;
    state.scenario.has_filepath = true;
    state.scenario.dirty = false;
    if (overwrite) state.file_status = StatusCode::file_overwritten;
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

void render_scenario_file_ui(
    World& world,
    RenderLoopConfig& cfg,
    RenderLoopState& state
) {
    if (im::CollapsingHeader("Scenario")) {
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
                ScenarioConfig loaded_config;
                state.file_status = load_scenario_json(path, loaded_config);
                if (state.file_status == StatusCode::ok) {
                    World loaded_world;
                    ScenarioBuildResult loaded_mappings;
                    WorldStepperConfig loaded_stepper = cfg.stepper_cfg;
                    state.file_status = build_world_from_scenario_config(
                        loaded_config,
                        loaded_world,
                        loaded_mappings,
                        loaded_stepper
                    );
                    if (state.file_status == StatusCode::ok) {
                        world = std::move(loaded_world);
                        cfg.stepper_cfg = loaded_stepper;
                        state.scenario.config = std::move(loaded_config);
                        state.scenario.build_result = std::move(loaded_mappings);
                        state.scenario.filepath = path;
                        state.scenario.has_filepath = true;
                        state.scenario.dirty = false;
                        state.wksp.dirty = true;
                    }
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
        im::Separator();
    }
}

} // namespace render_ui_detail
