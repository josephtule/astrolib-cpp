// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/body.hpp"
#include "core/estimation_common.hpp"
#include "core/measurement.hpp"
#include "core/observation_type.hpp"

static StatusCode add_typed_station_instrument(
    Station& station,
    ObservationType type,
    const matXd& R,
    InstrumentId& out_id,
    std::string name
) {
    i32 dim = measurement_dim(type);
    if (R.rows() != dim || R.cols() != dim) {
        return StatusCode::size_mismatch;
    }

    StationInstrument instrument;
    instrument.type = type;
    instrument.enabled = true;
    instrument.R = R;
    instrument.name = name;

    return add_station_instrument(station, instrument, out_id);
}

StatusCode station_measurement_covariance(
    const Station& station,
    InstrumentId instrument_id,
    matXd& R
) {
    auto it = station.instruments.find(instrument_id);
    if (it == station.instruments.end()) {
        return StatusCode::instrument_not_found;
    }

    R = it->second.R;

    return StatusCode::ok;
}

StatusCode station_measurement_covariance(
    const Station& station,
    ObservationType type,
    matXd& R
) {
    // only works if there is only one measurement of each type tied to the station
    i32 found_count = 0;
    matXd R_match;

    for (const auto& [id, instrument] : station.instruments) {
        if (!instrument.enabled) continue;
        if (instrument.type != type) continue;

        ++found_count;
        R_match = instrument.R;

        if (found_count > 1) {
            return StatusCode::duplicate_id;
        }
    }

    if (found_count == 0) {
        return StatusCode::instrument_not_found;
    }

    R = R_match;
    return StatusCode::ok;
}

StatusCode set_station_instrument(Station& station, const StationInstrument& instrument) {
    if (instrument.id == kInvalidInstrumentId) {
        return StatusCode::instrument_not_found;
    }
    auto it = station.instruments.find(instrument.id);
    if (it == station.instruments.end()) {
        return StatusCode::instrument_not_found;
    }

    it->second = instrument;

    station.enabled_instrument_ids = enabled_station_instrument_ids(station);

    return StatusCode::ok;
}

StatusCode add_station_instrument(
    Station& station,
    const StationInstrument& instrument,
    InstrumentId& out_id
) {
    StationInstrument copy = instrument;
    if (copy.id == kInvalidInstrumentId) {
        copy.id = station.next_instrument_id++;
    }

    auto it = station.instruments.find(copy.id);
    if (it != station.instruments.end()) {
        // id already exists
        out_id = copy.id;
        return set_station_instrument(station, copy);
    }

    station.instruments.emplace(copy.id, copy);
    out_id = copy.id;

    if (copy.enabled) {
        station.enabled_instrument_ids = enabled_station_instrument_ids(station);
    }

    return StatusCode::ok;
}

StatusCode add_station_instrument(Station& station, const StationInstrument& instrument) {
    InstrumentId _;
    return add_station_instrument(station, instrument, _);
}

StatusCode get_station_instrument(
    const Station& station,
    StationInstrument& instrument,
    InstrumentId id
) {

    auto it = station.instruments.find(id);
    if (it == station.instruments.end()) {
        return StatusCode::instrument_not_found;
    }
    instrument = it->second;

    return StatusCode::ok;
}

StatusCode add_radec_instrument(Station& station, const mat2d& R, std::string name) {
    InstrumentId _;
    return add_radec_instrument(station, R, _, name);
}

StatusCode add_radec_instrument(
    Station& station,
    const mat2d& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_station_instrument(station, ObservationType::radec, R, out_id, name);
}

StatusCode add_azel_instrument(Station& station, const mat2d& R, std::string name) {
    InstrumentId _;
    return add_azel_instrument(station, R, _, name);
}

StatusCode add_azel_instrument(
    Station& station,
    const mat2d& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_station_instrument(station, ObservationType::azel, R, out_id, name);
}

StatusCode add_range_instrument(Station& station, const matXd& R, std::string name) {
    InstrumentId _;
    return add_range_instrument(station, R, _, name);
}

StatusCode add_range_instrument(
    Station& station,
    const matXd& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_station_instrument(station, ObservationType::range, R, out_id, name);
}

StatusCode add_range_rate_instrument(Station& station, const matXd& R, std::string name) {
    InstrumentId _;
    return add_range_rate_instrument(station, R, _, name);
}

StatusCode add_range_rate_instrument(
    Station& station,
    const matXd& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_station_instrument(
        station,
        ObservationType::range_rate,
        R,
        out_id,
        name
    );
}

StatusCode add_pos_instrument(Station& station, const mat3d& R, std::string name) {
    InstrumentId _;
    return add_pos_instrument(station, R, _, name);
}

StatusCode add_pos_instrument(
    Station& station,
    const mat3d& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_station_instrument(station, ObservationType::pos, R, out_id, name);
}

StatusCode add_posvel_instrument(Station& station, const mat6d& R, std::string name) {
    InstrumentId _;
    return add_posvel_instrument(station, R, _, name);
}

StatusCode add_posvel_instrument(
    Station& station,
    const mat6d& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_station_instrument(
        station,
        ObservationType::pos_vel,
        R,
        out_id,
        name
    );
}

StatusCode add_rel_pos_instrument(Station& station, const mat3d& R, std::string name) {
    InstrumentId _;
    return add_rel_pos_instrument(station, R, _, name);
}

StatusCode add_rel_pos_instrument(
    Station& station,
    const mat3d& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_station_instrument(
        station,
        ObservationType::rel_pos,
        R,
        out_id,
        name
    );
}

StatusCode add_rel_posvel_instrument(Station& station, const mat6d& R, std::string name) {
    InstrumentId _;
    return add_rel_posvel_instrument(station, R, _, name);
}

StatusCode add_rel_posvel_instrument(
    Station& station,
    const mat6d& R,
    InstrumentId& out_id,
    std::string name
) {
    return add_typed_station_instrument(
        station,
        ObservationType::rel_pos_vel,
        R,
        out_id,
        name
    );
}

svec<InstrumentId> enabled_station_instrument_ids(const Station& station) {
    svec<InstrumentId> ids;
    ids.reserve(station.instruments.size());

    for (auto& [instr_id, instrument] : station.instruments) {
        if (instrument.enabled) {
            ids.push_back(instr_id);
        }
    }
    std::sort(ids.begin(), ids.end());

    return ids;
}

StatusCode enable_station_instrument(Station& station, InstrumentId instrument_id) {
    auto it = station.instruments.find(instrument_id);
    if (it == station.instruments.end()) {
        return StatusCode::instrument_not_found;
    }
    it->second.enabled = true;

    station.enabled_instrument_ids = enabled_station_instrument_ids(station);

    return StatusCode::ok;
}

StatusCode disable_station_instrument(Station& station, InstrumentId instrument_id) {
    auto it = station.instruments.find(instrument_id);
    if (it == station.instruments.end()) {
        return StatusCode::instrument_not_found;
    }
    it->second.enabled = false;

    station.enabled_instrument_ids = enabled_station_instrument_ids(station);

    return StatusCode::ok;
}

void print_station_instruments(const Station& station) {
    for (const auto& [id, instrument] : station.instruments) {
        std::string enabled_str = instrument.enabled ? "enabled" : "disabled";
        std::println("{} ({})", instrument.name, enabled_str);
    }
}
