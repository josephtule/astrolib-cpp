#pragma once

#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/estimation_common.hpp"
#include "core/state.hpp"
#include "core/world_history.hpp"
#include <fstream>

struct HistoryCSVExportOptions {
    bool include_attitude = false;
    bool require_all_samples = true;
    bool write_header = true;
    std::string separator = ",";
};

inline StatusCode write_body_history_csv(
    const WorldHistory& history,
    EntityId body_id,
    std::string filepath,
    const HistoryCSVExportOptions& opts
) {
    if (history.samples.size() == 0) {
        return ODStatus::empty_history;
    }

    std::ofstream file(filepath);
    if (!file.is_open()) {
        return StatusCode::file_not_found;
    }

    if (opts.write_header) {
        file << "t,entity_id,rx,ry,rz,vx,vy,vz";
        if (opts.include_attitude) {
            file << ",qx,qy,qz,qw,wx,wy,wz";
        }
        file << "\n";
    }

    for (i32 i = 0; i < history.samples.size(); ++i) {
        const WorldHistorySample& sample = history.samples[i];
        const auto it_tr = sample.x_tr.find(body_id);
        if (it_tr == sample.x_tr.end()) {
            if (opts.require_all_samples) {
                file.close();
                return ODStatus::sample_not_found;
            }
            continue;
        }
        const vec6d& x_tr_vec = statetr_to_vec6d(it_tr->second);

        file << sample.t << opts.separator << body_id;
        for (i32 j = 0; j < x_tr_vec.size(); ++j) {
            file << opts.separator << x_tr_vec(j);
        }

        if (opts.include_attitude) {
            const auto it_att = sample.x_att.find(body_id);
            if (it_att == sample.x_att.end()) {
                if (opts.require_all_samples) {
                    file.close();
                    return ODStatus::sample_not_found;
                }
                continue;
            }
            const vec7d& x_att_vec = stateatt_to_vec7d(it_att->second);

            for (i32 j = 0; j < x_att_vec.size(); ++j) {
                file << opts.separator << x_att_vec(j);
            }
        }

        file << "\n";

        if (file.fail()) {
            file.close();
            return ODStatus::file_write_failed;
        }
    }

    file.close();
    if (file.fail()) {
        return ODStatus::file_close_failed;
    }

    return StatusCode::ok;
}
