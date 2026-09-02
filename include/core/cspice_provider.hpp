// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/ephemeris.hpp"
#include "core/status.hpp"
#include "core/time.hpp"
#include "util/typedefs.hpp"

struct CSpiceError {
    string short_message;
    string long_message;
};

struct CSpiceKernelInfo {
    string filepath;
    string type;
    string source;
    i32 handle = 0;
};

// CSPICE owns one process-wide kernel pool. Destruction does not clear it.
class CSpiceKernelPool {
  public:
    StatusCode load(const string& filepath, CSpiceError* error = nullptr) const;
    StatusCode unload(const string& filepath, CSpiceError* error = nullptr) const;
    StatusCode clear(CSpiceError* error = nullptr) const;
    StatusCode loaded(svec<CSpiceKernelInfo>& out, CSpiceError* error = nullptr) const;
};

enum struct CSpiceAberrationCorrection {
    none,
    lt,
    lt_stellar,
    converged_newtonian,
    converged_newtonian_stellar,
    transmission_lt,
    transmission_lt_stellar,
    transmission_converged_newtonian,
    transmission_converged_newtonian_stellar
};

string cspice_aberration_string(CSpiceAberrationCorrection correction);

enum struct CSpiceBodyPreset {
    solar_system_barycenter,
    mercury_barycenter,
    venus_barycenter,
    earth_moon_barycenter,
    mars_barycenter,
    jupiter_barycenter,
    saturn_barycenter,
    uranus_barycenter,
    neptune_barycenter,
    pluto_barycenter,
    sun,
    mercury,
    venus,
    earth,
    moon,
    mars,
    jupiter,
    saturn,
    uranus,
    neptune,
    pluto
};

enum struct CSpiceFramePreset {
    j2000,
    eclipj2000,
    iau_sun,
    iau_earth,
    iau_moon,
    iau_mars

    // itrf93, // requires earth_latest_high_prec.bpc and earth_assoc_itrf93.tf
    // moon_pa, // requires moon_pa_de440_200625.bpc and a lunar FK
    // moon_me // requires moon_pa_de440_200625.bpc and a lunar FK
};

struct CSpiceBodyRef {
    enum struct Type { preset, naif_id, name };

    Type type = Type::name;
    CSpiceBodyPreset preset_value = CSpiceBodyPreset::earth;
    i32 id_value = 0;
    string name_value;

    CSpiceBodyRef() = default;
    CSpiceBodyRef(const char* value);
    CSpiceBodyRef(string value);

    static CSpiceBodyRef preset(CSpiceBodyPreset value);
    static CSpiceBodyRef id(i32 value);
    static CSpiceBodyRef name(string value);
};

struct CSpiceFrameRef {
    enum struct Type { preset, name };

    Type type = Type::name;
    CSpiceFramePreset preset_value = CSpiceFramePreset::j2000;
    string name_value;

    CSpiceFrameRef() = default;
    CSpiceFrameRef(const char* value);
    CSpiceFrameRef(string value);

    static CSpiceFrameRef preset(CSpiceFramePreset value);
    static CSpiceFrameRef name(string value);
};

string cspice_body_string(const CSpiceBodyRef& body);
string cspice_frame_string(const CSpiceFrameRef& frame);

struct CSpiceStateQuery {
    CSpiceBodyRef target;
    CSpiceBodyRef observer;
    CSpiceFrameRef frame = CSpiceFrameRef::preset(CSpiceFramePreset::j2000);
    CSpiceAberrationCorrection aberration = CSpiceAberrationCorrection::none;
};

struct CSpiceOrientationQuery {
    CSpiceBodyRef object;
    CSpiceFrameRef source_frame = CSpiceFrameRef::preset(CSpiceFramePreset::j2000);
    CSpiceFrameRef target_frame;
};

struct CSpiceSampleGrid {
    EphemerisEpochMetadata epoch{};
    svec<f64> dt;
};

bool cspice_available();
string cspice_toolkit_version();

StatusCode cspice_et_from_jd(
    const JulianDate& epoch,
    TimeScale scale,
    const TimeOffsets& offsets,
    f64& et
);

StatusCode query_cspice_state_et(
    const CSpiceStateQuery& query,
    f64 et,
    StateTr& out,
    f64* light_time = nullptr,
    CSpiceError* error = nullptr
);

StatusCode query_cspice_state(
    const CSpiceStateQuery& query,
    const JulianDate& epoch,
    TimeScale scale,
    const TimeOffsets& offsets,
    StateTr& out,
    f64* light_time = nullptr,
    CSpiceError* error = nullptr
);

StatusCode query_cspice_orientation_et(
    const CSpiceOrientationQuery& query,
    f64 et,
    StateAtt& out,
    CSpiceError* error = nullptr
);

StatusCode query_cspice_orientation(
    const CSpiceOrientationQuery& query,
    const JulianDate& epoch,
    TimeScale scale,
    const TimeOffsets& offsets,
    StateAtt& out,
    CSpiceError* error = nullptr
);

StatusCode sample_cspice_ephemeris(
    const CSpiceStateQuery& query,
    const CSpiceSampleGrid& grid,
    const TimeOffsets& offsets,
    CartesianEphemerisTable& out,
    CSpiceError* error = nullptr
);

StatusCode sample_cspice_orientation(
    const CSpiceOrientationQuery& query,
    const CSpiceSampleGrid& grid,
    const TimeOffsets& offsets,
    OrientationEphemerisTable& out,
    CSpiceError* error = nullptr
);

StatusCode validate_earth_orientation_source_selection(
    bool use_cspice_orientation,
    bool use_internal_eop
);
