#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/estimation_common.hpp"
#include "core/observation_type.hpp"

ODStatus stat_meas_cov(const Station& station, u32 instrument_id, matXd& R) {
    auto it = station.instruments.find(instrument_id);
    if (it == station.instruments.end()) {
        return ODStatus::instrument_not_found;
    }

    R = it->second.R;

    return ODStatus::ok;
}

ODStatus stat_meas_cov(const Station& station, ObservationType type, matXd& R) {
    // only works if there is only one measurement of each type tied to the station
    i32 found_count = 0;
    matXd R_match;

    for (const auto& [id, instrument] : station.instruments) {
        if (!instrument.enabled) continue;
        if (instrument.type != type) continue;

        ++found_count;
        R_match = instrument.R;

        if (found_count > 1) {
            return ODStatus::invalid_input;
        }
    }

    if (found_count == 0) {
        return ODStatus::instrument_not_found;
    }

    R = R_match;
    return ODStatus::ok;
}

ODStatus set_station_instrument(Station& station, const StationInstrument& instrument) {
    if (instrument.id == kInvalidInstrumentId) {
        return ODStatus::instrument_not_found;
    }
    auto it = station.instruments.find(instrument.id);
    if (it == station.instruments.end()) {
        return ODStatus::instrument_not_found;
    }

    it->second = instrument;

    return ODStatus::ok;
}

ODStatus add_station_instrument(
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
        return set_station_instrument(station, copy);
    }

    station.instruments.emplace(copy.id, copy);
    out_id = copy.id;
    return ODStatus::ok;
}

ODStatus add_station_instrument(Station& station, const StationInstrument& instrument) {
    InstrumentId _;
    return add_station_instrument(station, instrument, _);
}
