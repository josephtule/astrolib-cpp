// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/instrument.hpp"

#include "core/measurement_types.hpp"
#include "util/constants.hpp"

#include <algorithm>
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
    instrument.name = std::move(name);

    return add_instrument(suite, instrument, out_id);
}

StatusCode measurement_covariance(
    const InstrumentSuite& suite,
    InstrumentId instrument_id,
    matXd& R
) {
    auto it = suite.instruments.find(instrument_id);
    if (it == suite.instruments.end()) return StatusCode::instrument_not_found;

    R = it->second.R;
    return StatusCode::ok;
}

StatusCode set_instrument(InstrumentSuite& suite, const PlatformInstrument& instrument) {
    StatusCode status = validate_instrument(instrument);
    if (status != StatusCode::ok) return status;
    if (instrument.id == kInvalidInstrumentId) return StatusCode::instrument_not_found;

    auto it = suite.instruments.find(instrument.id);
    if (it == suite.instruments.end()) return StatusCode::instrument_not_found;

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
        out_id = copy.id;
        return set_instrument(suite, copy);
    }

    suite.instruments.emplace(copy.id, copy);
    out_id = copy.id;
    suite.enabled_ids = enabled_instrument_ids(suite);
    return StatusCode::ok;
}

StatusCode add_instrument(InstrumentSuite& suite, const PlatformInstrument& instrument) {
    InstrumentId out_id;
    return add_instrument(suite, instrument, out_id);
}

StatusCode get_instrument(
    const InstrumentSuite& suite,
    InstrumentId id,
    PlatformInstrument& out
) {
    auto it = suite.instruments.find(id);
    if (it == suite.instruments.end()) return StatusCode::instrument_not_found;

    out = it->second;
    return StatusCode::ok;
}

#define DEFINE_TYPED_INSTRUMENT(NAME, TYPE, MATRIX_TYPE)                                 \
    StatusCode add_##NAME##_instrument(                                                  \
        InstrumentSuite& suite,                                                          \
        const MATRIX_TYPE& R,                                                            \
        std::string name                                                                 \
    ) {                                                                                  \
        InstrumentId out_id;                                                             \
        return add_##NAME##_instrument(suite, R, out_id, std::move(name));               \
    }                                                                                    \
    StatusCode add_##NAME##_instrument(                                                  \
        InstrumentSuite& suite,                                                          \
        const MATRIX_TYPE& R,                                                            \
        InstrumentId& out_id,                                                            \
        std::string name                                                                 \
    ) {                                                                                  \
        return add_typed_instrument(                                                     \
            suite,                                                                       \
            ObservationType::TYPE,                                                       \
            R,                                                                           \
            out_id,                                                                      \
            std::move(name)                                                              \
        );                                                                               \
    }

DEFINE_TYPED_INSTRUMENT(radec, radec, mat2d)
DEFINE_TYPED_INSTRUMENT(azel, azel, mat2d)
DEFINE_TYPED_INSTRUMENT(range, range, matXd)
DEFINE_TYPED_INSTRUMENT(range_rate, range_rate, matXd)
DEFINE_TYPED_INSTRUMENT(pos, pos, mat3d)
DEFINE_TYPED_INSTRUMENT(posvel, pos_vel, mat6d)
DEFINE_TYPED_INSTRUMENT(rel_pos, rel_pos, mat3d)
DEFINE_TYPED_INSTRUMENT(rel_posvel, rel_pos_vel, mat6d)

#undef DEFINE_TYPED_INSTRUMENT

svec<InstrumentId> enabled_instrument_ids(const InstrumentSuite& suite) {
    svec<InstrumentId> ids;
    ids.reserve(suite.instruments.size());

    for (const auto& [instrument_id, instrument] : suite.instruments) {
        if (instrument.enabled) ids.push_back(instrument_id);
    }
    std::sort(ids.begin(), ids.end());
    return ids;
}

StatusCode enable_instrument(InstrumentSuite& suite, InstrumentId instrument_id) {
    auto it = suite.instruments.find(instrument_id);
    if (it == suite.instruments.end()) return StatusCode::instrument_not_found;

    it->second.enabled = true;
    suite.enabled_ids = enabled_instrument_ids(suite);
    return StatusCode::ok;
}

StatusCode disable_instrument(InstrumentSuite& suite, InstrumentId instrument_id) {
    auto it = suite.instruments.find(instrument_id);
    if (it == suite.instruments.end()) return StatusCode::instrument_not_found;

    it->second.enabled = false;
    suite.enabled_ids = enabled_instrument_ids(suite);
    return StatusCode::ok;
}
