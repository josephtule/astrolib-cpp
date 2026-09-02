// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/body.hpp"
#include "core/ephemeris_provider.hpp"
#include "core/measurement.hpp"
#include "core/observation_type.hpp"
#include "core/status.hpp"

#include <cmath>
#include <limits>
#include <utility>

static StatusCode validate_instrument(const PlatformInstrument& instrument) {
    i32 dim = measurement_dim(instrument.type);
    if (dim <= 0) return StatusCode::unsupported_type;
    if (instrument.R.rows() != dim || instrument.R.cols() != dim) {
        return StatusCode::size_mismatch;
    }
    if (!instrument.R.allFinite()) return StatusCode::invalid_covariance;

    for (i32 i = 0; i < dim; ++i) {
        if (instrument.R(i, i) < 0.0) return StatusCode::invalid_covariance;
        for (i32 j = i + 1; j < dim; ++j) {
            if (std::abs(instrument.R(i, j) - instrument.R(j, i)) >= tol12) {
                return StatusCode::invalid_covariance;
            }
        }
    }

    return StatusCode::ok;
}

static StatusCode add_typed_instrument(
    InstrumentSuite& suite,
    ObservationType type,
    const matXd& R,
    InstrumentId& out_id,
    std::string name
) {
    i32 dim = measurement_dim(type);
    if (dim <= 0) return StatusCode::unsupported_type;

    PlatformInstrument instrument;
    instrument.type = type;
    instrument.enabled = true;
    instrument.R = R;
    instrument.name = name;

    return add_instrument(suite, instrument, out_id);
}

StatusCode measurement_covariance(
    const InstrumentSuite& suite,
    InstrumentId instrument_id,
    matXd& R
) {
    auto it = suite.instruments.find(instrument_id);
    if (it == suite.instruments.end()) {
        return StatusCode::instrument_not_found;
    }

    R = it->second.R;

    return StatusCode::ok;
}

StatusCode set_instrument(InstrumentSuite& suite, const PlatformInstrument& instrument) {
    StatusCode status = validate_instrument(instrument);
    if (status != StatusCode::ok) return status;
    if (instrument.id == kInvalidInstrumentId) {
        return StatusCode::instrument_not_found;
    }
    auto it = suite.instruments.find(instrument.id);
    if (it == suite.instruments.end()) {
        return StatusCode::instrument_not_found;
    }

    it->second = instrument;

    suite.enabled_ids = enabled_instrument_ids(suite);

    return StatusCode::ok;
}

StatusCode add_instrument(
    InstrumentSuite& suite,
    const PlatformInstrument& instrument,
    InstrumentId& out_id
) {
    PlatformInstrument copy = instrument;
    StatusCode status = validate_instrument(copy);
    if (status != StatusCode::ok) return status;

    if (copy.id == kInvalidInstrumentId) {
        if (suite.next_id == kInvalidInstrumentId
            || suite.next_id == std::numeric_limits<InstrumentId>::max()) {
            return StatusCode::invalid_input;
        }
        copy.id = suite.next_id++;
    } else if (copy.id >= suite.next_id) {
        if (copy.id == std::numeric_limits<InstrumentId>::max()) {
            return StatusCode::invalid_input;
        }
        suite.next_id = copy.id + 1;
    }

    auto it = suite.instruments.find(copy.id);
    if (it != suite.instruments.end()) {
        // id already exists
        out_id = copy.id;
        return set_instrument(suite, copy);
    }

    suite.instruments.emplace(copy.id, copy);
    out_id = copy.id;

    suite.enabled_ids = enabled_instrument_ids(suite);

    return StatusCode::ok;
}

StatusCode add_instrument(InstrumentSuite& suite, const PlatformInstrument& instrument) {
    InstrumentId _;
    return add_instrument(suite, instrument, _);
}

StatusCode get_instrument(
    const InstrumentSuite& suite,
    InstrumentId id,
    PlatformInstrument& out
) {
    auto it = suite.instruments.find(id);
    if (it == suite.instruments.end()) {
        return StatusCode::instrument_not_found;
    }
    out = it->second;

    return StatusCode::ok;
}

StatusCode add_radec_instrument(
    InstrumentSuite& suite,
    const mat2d& R,
    std::string name
) {
    InstrumentId _;
    return add_radec_instrument(suite, R, _, name);
}

StatusCode add_radec_instrument(
    InstrumentSuite& suite,
    const mat2d& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_instrument(suite, ObservationType::radec, R, out_id, name);
}

StatusCode add_azel_instrument(InstrumentSuite& suite, const mat2d& R, std::string name) {
    InstrumentId _;
    return add_azel_instrument(suite, R, _, name);
}

StatusCode add_azel_instrument(
    InstrumentSuite& suite,
    const mat2d& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_instrument(suite, ObservationType::azel, R, out_id, name);
}

StatusCode add_range_instrument(
    InstrumentSuite& suite,
    const matXd& R,
    std::string name
) {
    InstrumentId _;
    return add_range_instrument(suite, R, _, name);
}

StatusCode add_range_instrument(
    InstrumentSuite& suite,
    const matXd& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_instrument(suite, ObservationType::range, R, out_id, name);
}

StatusCode add_range_rate_instrument(
    InstrumentSuite& suite,
    const matXd& R,
    std::string name
) {
    InstrumentId _;
    return add_range_rate_instrument(suite, R, _, name);
}

StatusCode add_range_rate_instrument(
    InstrumentSuite& suite,
    const matXd& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_instrument(suite, ObservationType::range_rate, R, out_id, name);
}

StatusCode add_pos_instrument(InstrumentSuite& suite, const mat3d& R, std::string name) {
    InstrumentId _;
    return add_pos_instrument(suite, R, _, name);
}

StatusCode add_pos_instrument(
    InstrumentSuite& suite,
    const mat3d& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_instrument(suite, ObservationType::pos, R, out_id, name);
}

StatusCode add_posvel_instrument(
    InstrumentSuite& suite,
    const mat6d& R,
    std::string name
) {
    InstrumentId _;
    return add_posvel_instrument(suite, R, _, name);
}

StatusCode add_posvel_instrument(
    InstrumentSuite& suite,
    const mat6d& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_instrument(suite, ObservationType::pos_vel, R, out_id, name);
}

StatusCode add_rel_pos_instrument(
    InstrumentSuite& suite,
    const mat3d& R,
    std::string name
) {
    InstrumentId _;
    return add_rel_pos_instrument(suite, R, _, name);
}

StatusCode add_rel_pos_instrument(
    InstrumentSuite& suite,
    const mat3d& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_instrument(suite, ObservationType::rel_pos, R, out_id, name);
}

StatusCode add_rel_posvel_instrument(
    InstrumentSuite& suite,
    const mat6d& R,
    std::string name
) {
    InstrumentId _;
    return add_rel_posvel_instrument(suite, R, _, name);
}

StatusCode add_rel_posvel_instrument(
    InstrumentSuite& suite,
    const mat6d& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_instrument(suite, ObservationType::rel_pos_vel, R, out_id, name);
}

svec<InstrumentId> enabled_instrument_ids(const InstrumentSuite& suite) {
    svec<InstrumentId> ids;
    ids.reserve(suite.instruments.size());

    for (auto& [instr_id, instrument] : suite.instruments) {
        if (instrument.enabled) {
            ids.push_back(instr_id);
        }
    }
    std::sort(ids.begin(), ids.end());

    return ids;
}

StatusCode enable_instrument(InstrumentSuite& suite, InstrumentId instrument_id) {
    auto it = suite.instruments.find(instrument_id);
    if (it == suite.instruments.end()) {
        return StatusCode::instrument_not_found;
    }
    it->second.enabled = true;

    suite.enabled_ids = enabled_instrument_ids(suite);

    return StatusCode::ok;
}

StatusCode disable_instrument(InstrumentSuite& suite, InstrumentId instrument_id) {
    auto it = suite.instruments.find(instrument_id);
    if (it == suite.instruments.end()) {
        return StatusCode::instrument_not_found;
    }
    it->second.enabled = false;

    suite.enabled_ids = enabled_instrument_ids(suite);

    return StatusCode::ok;
}

void print_instruments(const InstrumentSuite& suite) {
    for (const auto& [id, instrument] : suite.instruments) {
        std::string enabled_str = instrument.enabled ? "enabled" : "disabled";
        std::println("{} ({})", instrument.name, enabled_str);
    }
}

StatusCode set_celestial_ephemeris_providers(
    Celestial& celestial,
    BodyEphemerisProviders providers
) {
    StatusCode status;
    if (providers.translation != nullptr) {
        status = validate_ephemeris_provider(*providers.translation);
        if (status != StatusCode::ok) return status;
    }

    if (providers.orientation != nullptr) {
        status = validate_orientation_provider(*providers.orientation);
        if (status != StatusCode::ok) return status;
    }
    if (providers.orientation == nullptr
        && celestial.attitude_model == CelestialAttitudeModel::provider) {
        return StatusCode::invalid_input;
    }

    if (providers.translation != nullptr) celestial.propagate_tr = false;
    if (providers.orientation != nullptr) {
        celestial.attitude_model = CelestialAttitudeModel::provider;
        celestial.propagate_att = true;
    }
    celestial.ephemeris_providers = std::move(providers);

    return StatusCode::ok;
}

StatusCode instrument_suite_from_body(Body& body, InstrumentSuite*& out) {
    out = nullptr;
    switch (body.body_type) {
    case BodyType::station: {
        out = &static_cast<Station&>(body).instrument_suite;
    } break;
    case BodyType::satellite: {
        out = &static_cast<Satellite&>(body).instrument_suite;
    } break;
    case BodyType::celestial: [[fallthrough]];
    case BodyType::unknown: {
        return StatusCode::unsupported_type;
    } break;
    }

    return StatusCode::ok;
}

StatusCode instrument_suite_from_body(const Body& body, const InstrumentSuite*& out) {
    out = nullptr;
    switch (body.body_type) {
    case BodyType::station: {
        out = &static_cast<const Station&>(body).instrument_suite;
    } break;
    case BodyType::satellite: {
        out = &static_cast<const Satellite&>(body).instrument_suite;
    } break;
    case BodyType::celestial: [[fallthrough]];
    case BodyType::unknown: {
        return StatusCode::unsupported_type;
    } break;
    }

    return StatusCode::ok;
}
