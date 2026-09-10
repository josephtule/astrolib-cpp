// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/observation_type.hpp"
#include "core/status.hpp"
#include "util/typedefs.hpp"
#include "util/vecdefs.hpp"

#include <string>

using InstrumentId = u32;
inline constexpr InstrumentId kInvalidInstrumentId = 0;

struct PlatformInstrument {
    InstrumentId id = kInvalidInstrumentId;
    std::string name;
    ObservationType type = ObservationType::radec;
    matXd R;
    bool enabled = true;
};

struct InstrumentSuite {
    InstrumentId next_id = 1;
    umap<InstrumentId, PlatformInstrument> instruments;
    svec<InstrumentId> enabled_ids;
};

StatusCode measurement_covariance(
    const InstrumentSuite& suite,
    InstrumentId instrument_id,
    matXd& R
);

StatusCode set_instrument(InstrumentSuite& suite, const PlatformInstrument& instrument);

StatusCode add_instrument(
    InstrumentSuite& suite,
    const PlatformInstrument& instrument,
    InstrumentId& out_id
);

StatusCode add_instrument(InstrumentSuite& suite, const PlatformInstrument& instrument);

StatusCode get_instrument(
    const InstrumentSuite& suite,
    InstrumentId id,
    PlatformInstrument& out
);

svec<InstrumentId> enabled_instrument_ids(const InstrumentSuite& suite);

StatusCode enable_instrument(InstrumentSuite& suite, InstrumentId id);

StatusCode disable_instrument(InstrumentSuite& suite, InstrumentId id);

StatusCode add_radec_instrument(
    InstrumentSuite& suite,
    const mat2d& R,
    std::string name = "Ra/Dec Instrument"
);

StatusCode add_radec_instrument(
    InstrumentSuite& suite,
    const mat2d& R,
    InstrumentId& out_id,
    std::string name = "Ra/Dec Instrument"
);

StatusCode add_azel_instrument(
    InstrumentSuite& suite,
    const mat2d& R,
    std::string name = "Az/El Instrument"
);

StatusCode add_azel_instrument(
    InstrumentSuite& suite,
    const mat2d& R,
    InstrumentId& out_id,
    std::string name = "Az/El Instrument"
);

StatusCode add_range_instrument(
    InstrumentSuite& suite,
    const matXd& R,
    std::string name = "Range Instrument"
);

StatusCode add_range_instrument(
    InstrumentSuite& suite,
    const matXd& R,
    InstrumentId& out_id,
    std::string name = "Range Instrument"
);

StatusCode add_range_rate_instrument(
    InstrumentSuite& suite,
    const matXd& R,
    std::string name = "Range-Rate Instrument"
);

StatusCode add_range_rate_instrument(
    InstrumentSuite& suite,
    const matXd& R,
    InstrumentId& out_id,
    std::string name = "Range-Rate Instrument"
);

StatusCode add_pos_instrument(
    InstrumentSuite& suite,
    const mat3d& R,
    std::string name = "Simulation Position Instrument"
);

StatusCode add_pos_instrument(
    InstrumentSuite& suite,
    const mat3d& R,
    InstrumentId& out_id,
    std::string name = "Inertial Position Instrument"
);

StatusCode add_posvel_instrument(
    InstrumentSuite& suite,
    const mat6d& R,
    std::string name = "Simulation State Instrument"
);

StatusCode add_posvel_instrument(
    InstrumentSuite& suite,
    const mat6d& R,
    InstrumentId& out_id,
    std::string name = "Inertial State Instrument"
);

StatusCode add_rel_pos_instrument(
    InstrumentSuite& suite,
    const mat3d& R,
    std::string name = "Relative Position Instrument"
);

StatusCode add_rel_pos_instrument(
    InstrumentSuite& suite,
    const mat3d& R,
    InstrumentId& out_id,
    std::string name = "Relative Position Instrument"
);

StatusCode add_rel_posvel_instrument(
    InstrumentSuite& suite,
    const mat6d& R,
    std::string name = "Relative State Instrument"
);

StatusCode add_rel_posvel_instrument(
    InstrumentSuite& suite,
    const mat6d& R,
    InstrumentId& out_id,
    std::string name = "Relative State Instrument"
);
