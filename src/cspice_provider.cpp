// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/cspice_provider.hpp"
#include "core/ephemeris.hpp"
#include "core/state.hpp"
#include "core/status.hpp"
#include "core/time.hpp"
#include "core/transform.hpp"
#include "util/math.hpp"
#include "util/typedefs.hpp"

#include <algorithm>
#include <cmath>
#include <mutex>
#include <sstream>

#if ASTROLIB_HAS_CSPICE
#include "SpiceUsr.h"
#include "SpiceZdf.h"
#endif

namespace {

i32 cspice_body_id(CSpiceBodyPreset body) {
    switch (body) {
    case CSpiceBodyPreset::solar_system_barycenter: return 0;
    case CSpiceBodyPreset::mercury_barycenter: return 1;
    case CSpiceBodyPreset::venus_barycenter: return 2;
    case CSpiceBodyPreset::earth_moon_barycenter: return 3;
    case CSpiceBodyPreset::mars_barycenter: return 4;
    case CSpiceBodyPreset::jupiter_barycenter: return 5;
    case CSpiceBodyPreset::saturn_barycenter: return 6;
    case CSpiceBodyPreset::uranus_barycenter: return 7;
    case CSpiceBodyPreset::neptune_barycenter: return 8;
    case CSpiceBodyPreset::pluto_barycenter: return 9;
    case CSpiceBodyPreset::sun: return 10;
    case CSpiceBodyPreset::mercury: return 199;
    case CSpiceBodyPreset::venus: return 299;
    case CSpiceBodyPreset::earth: return 399;
    case CSpiceBodyPreset::moon: return 301;
    case CSpiceBodyPreset::mars: return 499;
    case CSpiceBodyPreset::jupiter: return 599;
    case CSpiceBodyPreset::saturn: return 699;
    case CSpiceBodyPreset::uranus: return 799;
    case CSpiceBodyPreset::neptune: return 899;
    case CSpiceBodyPreset::pluto: return 999;
    }
    return 0;
}

} // namespace

CSpiceBodyRef::CSpiceBodyRef(const char* value)
    : CSpiceBodyRef::CSpiceBodyRef(string{value == nullptr ? "" : value}) {}

CSpiceBodyRef::CSpiceBodyRef(string value)
    : type(Type::name), name_value(std::move(value)) {}

CSpiceBodyRef CSpiceBodyRef::preset(CSpiceBodyPreset value) {
    CSpiceBodyRef out;
    out.type = Type::preset;
    out.preset_value = value;
    return out;
}

CSpiceBodyRef CSpiceBodyRef::id(i32 value) {
    CSpiceBodyRef out;
    out.type = Type::naif_id;
    out.id_value = value;
    return out;
}

CSpiceBodyRef CSpiceBodyRef::name(string value) {
    return CSpiceBodyRef{std::move(value)};
}

CSpiceFrameRef::CSpiceFrameRef(const char* value)
    : CSpiceFrameRef::CSpiceFrameRef(string{value == nullptr ? "" : value}) {}

CSpiceFrameRef::CSpiceFrameRef(string value)
    : type(Type::name), name_value(std::move(value)) {}

CSpiceFrameRef CSpiceFrameRef::preset(CSpiceFramePreset value) {
    CSpiceFrameRef out;
    out.type = Type::preset;
    out.preset_value = value;
    return out;
}

CSpiceFrameRef CSpiceFrameRef::name(string value) {
    return CSpiceFrameRef{std::move(value)};
}

string cspice_body_string(const CSpiceBodyRef& body) {
    switch (body.type) {
    case CSpiceBodyRef::Type::preset:
        return std::to_string(cspice_body_id(body.preset_value));
    case CSpiceBodyRef::Type::naif_id: return std::to_string(body.id_value);
    case CSpiceBodyRef::Type::name: return body.name_value;
    }
    return {};
}

string cspice_frame_string(const CSpiceFrameRef& frame) {
    if (frame.type == CSpiceFrameRef::Type::name) return frame.name_value;

    switch (frame.preset_value) {
    case CSpiceFramePreset::j2000: return "J2000";
    case CSpiceFramePreset::eclipj2000: return "ECLIPJ2000";
    case CSpiceFramePreset::iau_sun: return "IAU_SUN";
    case CSpiceFramePreset::iau_earth: return "IAU_EARTH";
    case CSpiceFramePreset::iau_moon: return "IAU_MOON";
    case CSpiceFramePreset::iau_mars: return "IAU_MARS";
    }
    return {};
}

string cspice_aberration_string(CSpiceAberrationCorrection correction) {
    switch (correction) {
    case CSpiceAberrationCorrection::none: return "NONE";
    case CSpiceAberrationCorrection::lt: return "LT";
    case CSpiceAberrationCorrection::lt_stellar: return "LT+S";
    case CSpiceAberrationCorrection::converged_newtonian: return "CN";
    case CSpiceAberrationCorrection::converged_newtonian_stellar: return "CN+S";
    case CSpiceAberrationCorrection::transmission_lt: return "XLT";
    case CSpiceAberrationCorrection::transmission_lt_stellar: return "XLT+S";
    case CSpiceAberrationCorrection::transmission_converged_newtonian: return "XCN";
    case CSpiceAberrationCorrection::transmission_converged_newtonian_stellar:
        return "XCN+S";
    }
    return "NONE";
}

StatusCode validate_earth_orientation_source_selection(
    bool use_cspice_orientation,
    bool use_internal_eop
) {
    if (use_cspice_orientation && use_internal_eop) return StatusCode::invalid_input;
    return StatusCode::ok;
}

StatusCode cspice_et_from_jd(
    const JulianDate& epoch,
    TimeScale scale,
    const TimeOffsets& offsets,
    f64& et
) {
    if (!std::isfinite(epoch.day) || !std::isfinite(epoch.frac)) {
        return StatusCode::invalid_input;
    }

    JulianDate epoch_tdb
        = jd_scale_convert(normalize_jd(epoch), scale, TimeScale::tdb, offsets);
    f64 et_temp = ((epoch_tdb.day - 2451545.0) + epoch_tdb.frac) * 86400.0;
    if (!std::isfinite(et_temp)) return StatusCode::non_finite_result;

    et = et_temp;
    return StatusCode::ok;
}

#if ASTROLIB_HAS_CSPICE
StatusCode cspice_et_from_utccal_str(const string& utccal_str, f64& et) {
    // TODO: add utc cal string validation

    SpiceDouble et_temp;
    str2et_c(utccal_str.c_str(), &et_temp);
    if (!std::isfinite(et_temp)) return StatusCode::non_finite_result;

    et = f64(et_temp);
    return StatusCode::ok;
}
#endif

#if ASTROLIB_HAS_CSPICE

namespace {

constexpr SpiceInt cspice_message_length = 1841;
constexpr SpiceInt cspice_filepath_length = 4096;
constexpr SpiceInt cspice_type_length = 64;

std::mutex cspice_mutex;

void prepare_cspice_call() {
    if (failed_c()) reset_c();

    SpiceChar action[] = "RETURN";
    SpiceChar print_list[] = "NONE";
    erract_c("SET", 0, action);
    errprt_c("SET", 0, print_list);
}

StatusCode finish_cspice_call(CSpiceError* error) {
    if (!failed_c()) {
        if (error != nullptr) *error = CSpiceError{};
        return StatusCode::ok;
    }

    SpiceChar short_message[cspice_message_length]{};
    SpiceChar long_message[cspice_message_length]{};
    getmsg_c("SHORT", cspice_message_length, short_message);
    getmsg_c("LONG", cspice_message_length, long_message);

    if (error != nullptr) {
        error->short_message = short_message;
        error->long_message = long_message;
    }

    reset_c();
    return StatusCode::external_library_error;
}

bool valid_state_query(const CSpiceStateQuery& query) {
    return !cspice_body_string(query.target).empty()
           && !cspice_body_string(query.observer).empty()
           && !cspice_frame_string(query.frame).empty();
}

bool valid_orientation_query(const CSpiceOrientationQuery& query) {
    return !cspice_frame_string(query.source_frame).empty()
           && !cspice_frame_string(query.target_frame).empty();
}

JulianDate sample_epoch(const CSpiceSampleGrid& grid, f64 dt) {
    JulianDate epoch = grid.epoch.ref_epoch;
    epoch.frac += dt / 86400.0;
    return normalize_jd(epoch);
}

StatusCode validate_sample_grid(const CSpiceSampleGrid& grid) {
    if (grid.dt.empty()) return StatusCode::empty_ephemeris;
    if (!std::isfinite(grid.epoch.ref_epoch.day)
        || !std::isfinite(grid.epoch.ref_epoch.frac)
        || grid.epoch.offset_unit != UTime::second) {
        return StatusCode::invalid_input;
    }

    for (i32 i = 0; i < grid.dt.size(); ++i) {
        if (!std::isfinite(grid.dt[i])) return StatusCode::invalid_input;
        if (i > 0 && grid.dt[i] <= grid.dt[i - 1]) {
            return StatusCode::non_monotonic_time;
        }
    }

    return StatusCode::ok;
}

string loaded_kernel_paths() {
    SpiceInt count = 0;
    ktotal_c("ALL", &count);
    if (failed_c()) return {};

    std::ostringstream paths;
    for (SpiceInt i = 0; i < count; ++i) {
        SpiceChar filepath[cspice_filepath_length]{};
        SpiceChar type[cspice_type_length]{};
        SpiceChar source[cspice_filepath_length]{};
        SpiceInt handle = 0;
        SpiceBoolean found = SPICEFALSE;
        kdata_c(
            i,
            "ALL",
            cspice_filepath_length,
            cspice_type_length,
            cspice_filepath_length,
            filepath,
            type,
            source,
            &handle,
            &found
        );
        if (failed_c() || !found) return {};
        if (i > 0) paths << ';';
        paths << filepath;
    }
    return paths.str();
}

} // namespace

bool cspice_available() { return true; }

string cspice_toolkit_version() {
    std::scoped_lock lock(cspice_mutex);
    prepare_cspice_call();
    const SpiceChar* version = tkvrsn_c("TOOLKIT");
    if (failed_c() || version == nullptr) {
        reset_c();
        return {};
    }
    return version;
}

StatusCode CSpiceKernelPool::load(const string& filepath, CSpiceError* error) const {
    if (filepath.empty()) return StatusCode::invalid_input;

    std::scoped_lock lock(cspice_mutex);
    prepare_cspice_call();
    furnsh_c(filepath.c_str());
    return finish_cspice_call(error);
}

StatusCode CSpiceKernelPool::unload(const string& filepath, CSpiceError* error) const {
    if (filepath.empty()) return StatusCode::invalid_input;

    std::scoped_lock lock(cspice_mutex);
    prepare_cspice_call();
    unload_c(filepath.c_str());
    return finish_cspice_call(error);
}

StatusCode CSpiceKernelPool::clear(CSpiceError* error) const {
    std::scoped_lock lock(cspice_mutex);
    prepare_cspice_call();
    kclear_c();
    return finish_cspice_call(error);
}

StatusCode CSpiceKernelPool::loaded(
    svec<CSpiceKernelInfo>& out,
    CSpiceError* error
) const {
    std::scoped_lock lock(cspice_mutex);
    prepare_cspice_call();

    SpiceInt count = 0;
    ktotal_c("ALL", &count);
    StatusCode status = finish_cspice_call(error);
    if (status != StatusCode::ok) return status;

    svec<CSpiceKernelInfo> temp;
    temp.reserve(count);
    for (SpiceInt i = 0; i < count; ++i) {
        SpiceChar filepath[cspice_filepath_length]{};
        SpiceChar type[cspice_type_length]{};
        SpiceChar source[cspice_filepath_length]{};
        SpiceInt handle = 0;
        SpiceBoolean found = SPICEFALSE;
        kdata_c(
            i,
            "ALL",
            cspice_filepath_length,
            cspice_type_length,
            cspice_filepath_length,
            filepath,
            type,
            source,
            &handle,
            &found
        );
        status = finish_cspice_call(error);
        if (status != StatusCode::ok) return status;
        if (!found) return StatusCode::kernel_not_loaded;

        temp.push_back(
            CSpiceKernelInfo{
                .filepath = filepath,
                .type = type,
                .source = source,
                .handle = handle
            }
        );
    }

    out = std::move(temp);
    return StatusCode::ok;
}

StatusCode query_cspice_state_et(
    const CSpiceStateQuery& query,
    f64 et,
    StateTr& out,
    f64* light_time,
    CSpiceError* error
) {
    if (!valid_state_query(query) || !std::isfinite(et)) {
        return StatusCode::invalid_input;
    }

    std::scoped_lock lock(cspice_mutex);
    prepare_cspice_call();

    SpiceDouble state[6]{};
    SpiceDouble lt = 0.0;
    string target = cspice_body_string(query.target);
    string observer = cspice_body_string(query.observer);
    string frame = cspice_frame_string(query.frame);
    string aberration = cspice_aberration_string(query.aberration);
    spkezr_c(
        target.c_str(),
        et,
        frame.c_str(),
        aberration.c_str(),
        observer.c_str(),
        state,
        &lt
    );
    StatusCode status = finish_cspice_call(error);
    if (status != StatusCode::ok) return status;

    StateTr temp{
        .r = vec3d{state[0], state[1], state[2]},
        .v = vec3d{state[3], state[4], state[5]}
    };
    if (!finite_state(temp) || !std::isfinite(lt)) return StatusCode::non_finite_result;

    out = temp;
    if (light_time != nullptr) *light_time = lt;
    return StatusCode::ok;
}

StatusCode query_cspice_state(
    const CSpiceStateQuery& query,
    const JulianDate& epoch,
    TimeScale scale,
    const TimeOffsets& offsets,
    StateTr& out,
    f64* light_time,
    CSpiceError* error
) {
    f64 et = 0.0;
    StatusCode status = cspice_et_from_jd(epoch, scale, offsets, et);
    if (status != StatusCode::ok) return status;
    return query_cspice_state_et(query, et, out, light_time, error);
}

StatusCode query_cspice_orientation_et(
    const CSpiceOrientationQuery& query,
    f64 et,
    StateAtt& out,
    CSpiceError* error
) {
    if (!valid_orientation_query(query) || !std::isfinite(et)) {
        return StatusCode::invalid_input;
    }

    std::scoped_lock lock(cspice_mutex);
    prepare_cspice_call();

    string source_frame = cspice_frame_string(query.source_frame);
    string target_frame = cspice_frame_string(query.target_frame);
    SpiceDouble transform[6][6]{};
    sxform_c(source_frame.c_str(), target_frame.c_str(), et, transform);
    StatusCode status = finish_cspice_call(error);
    if (status != StatusCode::ok) return status;

    SpiceDouble rotation[3][3]{};
    SpiceDouble angular_velocity[3]{};
    xf2rav_c(transform, rotation, angular_velocity);

    mat3d R_target_source;
    vec3d w_target_source_source;
    for (i32 i = 0; i < 3; ++i) {
        w_target_source_source(i) = angular_velocity[i];
        for (i32 j = 0; j < 3; ++j) {
            R_target_source(i, j) = rotation[i][j];
        }
    }

    StateAtt temp{
        .q = dcm_to_ep(R_target_source),
        .w = R_target_source * w_target_source_source
    };
    if (!finite_state_att(temp) || temp.q.norm() <= tol12) {
        return StatusCode::non_finite_result;
    }
    temp.q.normalize();

    out = temp;
    return StatusCode::ok;
}

StatusCode query_cspice_orientation(
    const CSpiceOrientationQuery& query,
    const JulianDate& epoch,
    TimeScale scale,
    const TimeOffsets& offsets,
    StateAtt& out,
    CSpiceError* error
) {
    f64 et = 0.0;
    StatusCode status = cspice_et_from_jd(epoch, scale, offsets, et);
    if (status != StatusCode::ok) return status;
    return query_cspice_orientation_et(query, et, out, error);
}

StatusCode sample_cspice_ephemeris(
    const CSpiceStateQuery& query,
    const CSpiceSampleGrid& grid,
    const TimeOffsets& offsets,
    CartesianEphemerisTable& out,
    CSpiceError* error
) {
    if (!valid_state_query(query)) return StatusCode::invalid_input;
    StatusCode status = validate_sample_grid(grid);
    if (status != StatusCode::ok) return status;

    CartesianEphemerisTable temp;
    string target = cspice_body_string(query.target);
    string observer = cspice_body_string(query.observer);
    string frame = cspice_frame_string(query.frame);
    temp.metadata.epoch = grid.epoch;
    temp.metadata.frame = {.object = target, .center = observer, .frame = frame};
    temp.metadata.units = {.length = ULength::kilometer, .time = UTime::second};
    {
        std::scoped_lock lock(cspice_mutex);
        prepare_cspice_call();
        string paths = loaded_kernel_paths();
        status = finish_cspice_call(error);
        if (status != StatusCode::ok) return status;
        temp.metadata.source
            = {.source_type = "cspice_spk",
               .source_name = target + " relative to " + observer,
               .source_path = std::move(paths),
               .description = "Materialized from a CSPICE state query"};
    }
    temp.dt = grid.dt;
    temp.states.reserve(grid.dt.size());

    for (f64 dt : grid.dt) {
        StateTr state;
        status = query_cspice_state(
            query,
            sample_epoch(grid, dt),
            grid.epoch.time_scale,
            offsets,
            state,
            nullptr,
            error
        );
        if (status != StatusCode::ok) return status;
        temp.states.push_back(state);
    }

    status = validate_cartesian_ephemeris_table(temp);
    if (status != StatusCode::ok) return status;
    out = std::move(temp);
    return StatusCode::ok;
}

StatusCode sample_cspice_orientation(
    const CSpiceOrientationQuery& query,
    const CSpiceSampleGrid& grid,
    const TimeOffsets& offsets,
    OrientationEphemerisTable& out,
    CSpiceError* error
) {
    if (!valid_orientation_query(query)) return StatusCode::invalid_input;
    StatusCode status = validate_sample_grid(grid);
    if (status != StatusCode::ok) return status;

    OrientationEphemerisTable temp;
    string object = cspice_body_string(query.object);
    string source_frame = cspice_frame_string(query.source_frame);
    string target_frame = cspice_frame_string(query.target_frame);
    temp.metadata.epoch = grid.epoch;
    temp.metadata.frame
        = {.object = object.empty() ? target_frame : object,
           .source_frame = source_frame,
           .target_frame = target_frame};
    temp.metadata.source
        = {.source_type = "cspice_orientation",
           .source_name = source_frame + " to " + target_frame,
           .description = "Materialized from a CSPICE frame query"};
    {
        std::scoped_lock lock(cspice_mutex);
        prepare_cspice_call();
        temp.metadata.source.source_path = loaded_kernel_paths();
        status = finish_cspice_call(error);
        if (status != StatusCode::ok) return status;
    }
    temp.dt = grid.dt;
    temp.states.reserve(grid.dt.size());

    for (f64 dt : grid.dt) {
        StateAtt state;
        status = query_cspice_orientation(
            query,
            sample_epoch(grid, dt),
            grid.epoch.time_scale,
            offsets,
            state,
            error
        );
        if (status != StatusCode::ok) return status;
        temp.states.push_back(state);
    }

    status = canonicalize_orientation_ephemeris_samples(temp);
    if (status != StatusCode::ok) return status;
    out = std::move(temp);
    return StatusCode::ok;
}

#else

namespace {

StatusCode cspice_unavailable(CSpiceError* error) {
    if (error != nullptr) {
        error->short_message = "CSPICE support is disabled";
        error->long_message = "Configure with -DASTROLIB_ENABLE_CSPICE=ON";
    }
    return StatusCode::external_library_unavailable;
}

} // namespace

bool cspice_available() { return false; }

string cspice_toolkit_version() { return {}; }

StatusCode CSpiceKernelPool::load(const string&, CSpiceError* error) const {
    return cspice_unavailable(error);
}

StatusCode CSpiceKernelPool::unload(const string&, CSpiceError* error) const {
    return cspice_unavailable(error);
}

StatusCode CSpiceKernelPool::clear(CSpiceError* error) const {
    return cspice_unavailable(error);
}

StatusCode CSpiceKernelPool::loaded(svec<CSpiceKernelInfo>&, CSpiceError* error) const {
    return cspice_unavailable(error);
}

StatusCode query_cspice_state_et(
    const CSpiceStateQuery&,
    f64,
    StateTr&,
    f64*,
    CSpiceError* error
) {
    return cspice_unavailable(error);
}

StatusCode query_cspice_state(
    const CSpiceStateQuery&,
    const JulianDate&,
    TimeScale,
    const TimeOffsets&,
    StateTr&,
    f64*,
    CSpiceError* error
) {
    return cspice_unavailable(error);
}

StatusCode query_cspice_orientation_et(
    const CSpiceOrientationQuery&,
    f64,
    StateAtt&,
    CSpiceError* error
) {
    return cspice_unavailable(error);
}

StatusCode query_cspice_orientation(
    const CSpiceOrientationQuery&,
    const JulianDate&,
    TimeScale,
    const TimeOffsets&,
    StateAtt&,
    CSpiceError* error
) {
    return cspice_unavailable(error);
}

StatusCode sample_cspice_ephemeris(
    const CSpiceStateQuery&,
    const CSpiceSampleGrid&,
    const TimeOffsets&,
    CartesianEphemerisTable&,
    CSpiceError* error
) {
    return cspice_unavailable(error);
}

StatusCode sample_cspice_orientation(
    const CSpiceOrientationQuery&,
    const CSpiceSampleGrid&,
    const TimeOffsets&,
    OrientationEphemerisTable&,
    CSpiceError* error
) {
    return cspice_unavailable(error);
}

#endif
