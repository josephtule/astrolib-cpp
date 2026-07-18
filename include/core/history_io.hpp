// Copyright 2025-2026 Joseph Tu Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/estimation_common.hpp"
#include "core/state.hpp"
#include "core/world_history.hpp"
#include <fstream>
#include <iomanip>
#include <set>
#include <string>

struct HistoryCSVExportOptions {
    bool include_attitude = false;
    bool require_all_samples = true;
    bool write_header = true;
    std::string separator = ",";
    i32 precision = 17;
};

inline StatusCode open_history_csv(std::ofstream& file, const std::string& filepath) {
    file.open(filepath);
    if (!file.is_open()) {
        return StatusCode::file_not_found;
    }

    return StatusCode::ok;
}

inline void write_history_csv_header(
    std::ofstream& file,
    const HistoryCSVExportOptions& opts
) {
    const std::string& sep = opts.separator;

    file << "t" << sep << "entity_id" << sep << "body_type" << sep << "rx" << sep << "ry"
         << sep << "rz" << sep << "vx" << sep << "vy" << sep << "vz";
    if (opts.include_attitude) {
        file << sep << "qx" << sep << "qy" << sep << "qz" << sep << "qw" << sep << "wx"
             << sep << "wy" << sep << "wz";
    }
    file << "\n";
}

inline StatusCode close_history_csv(std::ofstream& file) {
    file.close();
    if (file.fail()) {
        return StatusCode::file_close_failed;
    }

    return StatusCode::ok;
}

inline StatusCode write_history_csv_row(
    std::ofstream& file,
    f64 t,
    EntityId id,
    const StateTr& x_tr,
    const StateAtt* x_att,
    const HistoryCSVExportOptions& opts,
    const BodyType type = BodyType::unknown
) {
    if (opts.include_attitude && x_att == nullptr) {
        return StatusCode::sample_not_found;
    }

    const vec6d x_tr_vec = statetr_to_vec6d(x_tr);

    file << t << opts.separator << id << opts.separator << body_type_str(type);
    for (i32 j = 0; j < x_tr_vec.size(); ++j) {
        file << opts.separator << x_tr_vec(j);
    }

    if (opts.include_attitude) {
        const vec7d x_att_vec = stateatt_to_vec7d(*x_att);
        for (i32 j = 0; j < x_att_vec.size(); ++j) {
            file << opts.separator << x_att_vec(j);
        }
    }

    file << "\n";

    if (file.fail()) {
        return StatusCode::file_write_failed;
    }

    return StatusCode::ok;
}

inline StatusCode validate_history_for_export(const WorldHistory& history) {
    if (history.samples.size() == 0) {
        return StatusCode::empty_history;
    }

    return StatusCode::ok;
}

inline StatusCode write_body_history_csv(
    const WorldHistory& history,
    EntityId body_id,
    const std::string& filepath,
    const HistoryCSVExportOptions& opts = HistoryCSVExportOptions{}
) {
    StatusCode status = validate_history_for_export(history);
    if (status != StatusCode::ok) {
        return status;
    }

    std::ofstream file;
    status = open_history_csv(file, filepath);
    if (status != StatusCode::ok) {
        return status;
    }
    file << std::setprecision(opts.precision);

    if (opts.write_header) {
        write_history_csv_header(file, opts);
    }

    for (i32 i = 0; i < history.samples.size(); ++i) {
        const WorldHistorySample& sample = history.samples[i];
        const auto it_tr = sample.x_tr.find(body_id);
        if (it_tr == sample.x_tr.end()) {
            if (opts.require_all_samples) {
                file.close();
                return StatusCode::sample_not_found;
            }
            continue;
        }

        const StateAtt* x_att = nullptr;
        if (opts.include_attitude) {
            const auto it_att = sample.x_att.find(body_id);
            if (it_att == sample.x_att.end()) {
                if (opts.require_all_samples) {
                    file.close();
                    return StatusCode::sample_not_found;
                }
                continue;
            }
            x_att = &it_att->second;
        }

        status
            = write_history_csv_row(file, sample.t, body_id, it_tr->second, x_att, opts);
        if (status != StatusCode::ok) {
            file.close();
            return status;
        }
    }

    return close_history_csv(file);
}

inline std::set<EntityId> sample_tr_ids(const WorldHistorySample& sample) {
    std::set<EntityId> ids;
    for (const auto& [id, _] : sample.x_tr) {
        ids.insert(id);
    }

    return ids;
}

inline bool sample_has_all_tr_ids(
    const WorldHistorySample& sample,
    const std::set<EntityId>& ids
) {
    for (EntityId id : ids) {
        if (sample.x_tr.find(id) == sample.x_tr.end()) {
            return false;
        }
    }

    return true;
}

inline StatusCode validate_world_history_csv_samples(
    const WorldHistory& history,
    const HistoryCSVExportOptions& opts
) {
    StatusCode status = validate_history_for_export(history);
    if (status != StatusCode::ok) {
        return status;
    }

    if (!opts.require_all_samples) {
        return StatusCode::ok;
    }

    const std::set<EntityId> ids = sample_tr_ids(history.samples.front());
    for (const WorldHistorySample& sample : history.samples) {
        if (!sample_has_all_tr_ids(sample, ids)) {
            return StatusCode::sample_not_found;
        }
        if (opts.include_attitude) {
            for (EntityId id : ids) {
                if (sample.x_att.find(id) == sample.x_att.end()) {
                    return StatusCode::sample_not_found;
                }
            }
        }
    }

    return StatusCode::ok;
}

inline StatusCode write_world_history_csv(
    const WorldHistory& history,
    const std::string& filepath,
    const HistoryCSVExportOptions& opts = HistoryCSVExportOptions{}
) {
    // does not write anchored station states
    StatusCode status = validate_world_history_csv_samples(history, opts);
    if (status != StatusCode::ok) {
        return status;
    }

    std::ofstream file;
    status = open_history_csv(file, filepath);
    if (status != StatusCode::ok) {
        return status;
    }
    file << std::setprecision(opts.precision);

    if (opts.write_header) {
        write_history_csv_header(file, opts);
    }

    for (const WorldHistorySample& sample : history.samples) {
        for (const auto& [id, x_tr] : sample.x_tr) {
            const StateAtt* x_att = nullptr;
            if (opts.include_attitude) {
                const auto it_att = sample.x_att.find(id);
                if (it_att == sample.x_att.end()) {
                    if (opts.require_all_samples) {
                        file.close();
                        return StatusCode::sample_not_found;
                    }
                    continue;
                }
                x_att = &it_att->second;
            }

            status = write_history_csv_row(file, sample.t, id, x_tr, x_att, opts);
            if (status != StatusCode::ok) {
                file.close();
                return status;
            }
        }
    }

    return close_history_csv(file);
}

inline StatusCode write_station_history_csv(
    const World& world,
    const WorldHistory& history,
    EntityId station_id,
    const std::string& filepath,
    const HistoryCSVExportOptions& opts = HistoryCSVExportOptions{}
) {
    StatusCode status = validate_history_for_export(history);
    if (status != StatusCode::ok) {
        return status;
    }

    if (world.station(station_id) == nullptr) {
        return StatusCode::observer_not_found;
    }

    std::ofstream file;
    status = open_history_csv(file, filepath);
    if (status != StatusCode::ok) {
        return status;
    }
    file << std::setprecision(opts.precision);

    if (opts.write_header) {
        write_history_csv_header(file, opts);
    }

    for (const WorldHistorySample& sample : history.samples) {
        StateTr x_tr;
        status
            = sample_station_tr_interp_linear(world, history, station_id, sample.t, x_tr);
        if (status != StatusCode::ok) {
            if (opts.require_all_samples) {
                file.close();
                return status;
            }
            continue;
        }

        StateAtt x_att;
        StateAtt* x_att_ptr = nullptr;
        if (opts.include_attitude) {
            status = sample_station_att_interp_linear(
                world,
                history,
                station_id,
                sample.t,
                x_att
            );
            if (status != StatusCode::ok) {
                if (opts.require_all_samples) {
                    file.close();
                    return status;
                }
                continue;
            }
            x_att_ptr = &x_att;
        }
        status = write_history_csv_row(
            file,
            sample.t,
            station_id,
            x_tr,
            x_att_ptr,
            opts,
            BodyType::station
        );
        if (status != StatusCode::ok) {
            file.close();
            return status;
        }
    }

    return close_history_csv(file);
}

inline StatusCode write_world_history_csv(
    const World& world,
    const WorldHistory& history,
    const std::string& filepath,
    const HistoryCSVExportOptions& opts = HistoryCSVExportOptions{}
) {
    // writes all bodies (linear interp for stations)

    StatusCode status = validate_world_history_csv_samples(history, opts);
    if (status != StatusCode::ok) {
        return status;
    }

    std::ofstream file;
    status = open_history_csv(file, filepath);
    if (status != StatusCode::ok) {
        return status;
    }
    file << std::setprecision(opts.precision);

    if (opts.write_header) {
        write_history_csv_header(file, opts);
    }

    for (const WorldHistorySample& sample : history.samples) {
        for (const auto& [id, x_tr] : sample.x_tr) {
            const Body* body = world.body(id);
            if (body == nullptr) {
                if (opts.require_all_samples) {
                    file.close();
                    return StatusCode::sample_not_found;
                }
                continue;
            }

            if (body->body_type != BodyType::station) {
                // celestials + satellites
                const StateAtt* x_att = nullptr;
                if (opts.include_attitude) {
                    const auto it_att = sample.x_att.find(id);
                    if (it_att == sample.x_att.end()) {
                        if (opts.require_all_samples) {
                            file.close();
                            return StatusCode::sample_not_found;
                        }
                        continue;
                    }
                    x_att = &it_att->second;
                }

                status = write_history_csv_row(
                    file,
                    sample.t,
                    id,
                    x_tr,
                    x_att,
                    opts,
                    body->body_type
                );
                if (status != StatusCode::ok) {
                    file.close();
                    return status;
                }
            } else {
                // stations
                StateTr x_tr;
                status
                    = sample_station_tr_interp_linear(world, history, id, sample.t, x_tr);
                if (status != StatusCode::ok) {
                    if (opts.require_all_samples) {
                        file.close();
                        return status;
                    }
                    continue;
                }

                StateAtt x_att;
                StateAtt* x_att_ptr = nullptr;
                if (opts.include_attitude) {
                    status = sample_station_att_interp_linear(
                        world,
                        history,
                        id,
                        sample.t,
                        x_att
                    );
                    if (status != StatusCode::ok) {
                        if (opts.require_all_samples) {
                            file.close();
                            return status;
                        }
                        continue;
                    }
                    x_att_ptr = &x_att;
                    status = write_history_csv_row(
                        file,
                        sample.t,
                        id,
                        x_tr,
                        x_att_ptr,
                        opts,
                        body->body_type
                    );
                    if (status != StatusCode::ok) {
                        file.close();
                        return status;
                    }
                }
            }
        }
    }

    return close_history_csv(file);
}