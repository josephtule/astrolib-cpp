// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/estimation_observer.hpp"
#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/ephemeris.hpp"
#include "core/ephemeris_provider.hpp"
#include "core/state.hpp"
#include "core/status.hpp"
#include "core/time.hpp"

#include "core/world_history.hpp"
#include "util/math.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"
#include <cmath>

string observer_state_source_string(ObserverStateSource source) {
    switch (source) {
    case ObserverStateSource::world: return "world";
    case ObserverStateSource::provider: return "provider";
    case ObserverStateSource::estimate: return "estimate";
    }
    return "unknown";
}

string od_observer_missing_epoch_policy_string(ODObserverMissingEpochPolicy policy) {
    switch (policy) {
    case ODObserverMissingEpochPolicy::reject: return "reject";
    case ODObserverMissingEpochPolicy::propagate: return "propagate";
    }
    return "unknown";
}

StatusCode validate_od_observer_binding(const ODObserverBinding& binding) {
    if (binding.observer_id == kInvalidEntityId) return StatusCode::invalid_input;

    switch (binding.source) {
    case ObserverStateSource::world:
    case ObserverStateSource::provider:
    case ObserverStateSource::estimate: break;
    default: return StatusCode::unsupported_type;
    }

    switch (binding.missing_epoch) {
    case ODObserverMissingEpochPolicy::reject:
    case ODObserverMissingEpochPolicy::propagate: break;
    default: return StatusCode::unsupported_type;
    }

    if (binding.source != ObserverStateSource::world) {
        if (binding.source_id.empty()) return StatusCode::invalid_input;
    }

    if (binding.known_state && binding.source == ObserverStateSource::estimate)
        return StatusCode::invalid_input;

    return StatusCode::ok;
}

StatusCode find_od_observer_binding(
    const svec<ODObserverBinding>& bindings,
    EntityId observer_id,
    const ODObserverBinding*& out
) {
    out = nullptr;
    const ODObserverBinding* temp = nullptr;

    if (observer_id == kInvalidEntityId) return StatusCode::invalid_input;

    StatusCode status;
    bool found = false;
    for (const auto& binding : bindings) {
        if (observer_id == binding.observer_id) {
            if (!found) {
                temp = &binding;
                found = true;
            } else {
                return StatusCode::duplicate_id;
            }
        }
    }

    if (!found) {
        return StatusCode::observer_not_found;
    }

    status = validate_od_observer_binding(*temp);
    if (status != StatusCode::ok) return status;

    out = temp;

    return StatusCode::ok;
}

StatusCode validate_od_observer_bindings(const svec<ODObserverBinding>& bindings) {
    StatusCode status;
    uset<EntityId> observer_ids;

    for (auto binding : bindings) {
        auto it = observer_ids.insert(binding.observer_id);
        if (!it.second) return StatusCode::duplicate_id;

        status = validate_od_observer_binding(binding);
        if (status != StatusCode::ok) return status;
    }

    return StatusCode::ok;
}

static StatusCode observer_query_epoch(
    const ODObserverEpochContext& epoch,
    f64 t,
    JulianDate& out
) {
    if (!epoch.has_calendar_epoch) return StatusCode::invalid_input;
    if (!std::isfinite(t)) return StatusCode::invalid_input;
    if (!std::isfinite(epoch.t_ref)) return StatusCode::invalid_input;

    f64 dt = t - epoch.t_ref;
    JulianDate temp = jd_from_scalar(jd_to_scalar(epoch.jd_ref) + dt / seconds_per_day);
    temp = normalize_jd(temp);

    StatusCode status = validate_julian_date(temp);
    if (status != StatusCode::ok) return status;

    out = temp;

    return StatusCode::ok;
}

StatusCode validate_od_observer_sample(
    const ODObserverSample& sample,
    f64 tol_covariance
) {
    if (sample.observer_id == kInvalidEntityId) return StatusCode::invalid_input;

    if (!std::isfinite(sample.t) || !finite_state(sample.x)
        || !finite_nonneg(tol_covariance)) {
        return StatusCode::invalid_input;
    }

    if (sample.center.empty() || sample.frame_name.empty())
        return StatusCode::invalid_input;

    switch (sample.source) {
    case ObserverStateSource::world:
    case ObserverStateSource::provider:
    case ObserverStateSource::estimate: break;
    default: return StatusCode::unsupported_type;
    }
    if (sample.source != ObserverStateSource::world) {
        if (sample.source_id.empty()) return StatusCode::invalid_input;
    }

    if (!covariance_is_psd(sample.P, tol_covariance))
        return StatusCode::invalid_covariance;

    if (sample.known_state) {
        if (sample.source == ObserverStateSource::estimate)
            return StatusCode::invalid_input;

        if (!sample.P.isZero(tol_covariance)) {
            return StatusCode::invalid_covariance;
        }
    }

    return StatusCode::ok;
}

static StatusCode validate_od_observer_query_context(const ODObserverQueryContext& ctx) {
    if (ctx.center.empty() || ctx.frame_name.empty()) return StatusCode::invalid_input;
    if (!finite_nonneg(ctx.tol_time) || !finite_nonneg(ctx.tol_covariance)
        || !std::isfinite(ctx.epoch.t_ref)) {
        return StatusCode::invalid_input;
    }

    if (ctx.epoch.has_calendar_epoch) {
        StatusCode status = validate_julian_date(ctx.epoch.jd_ref);
        if (status != StatusCode::ok) return status;
    }

    if (ctx.world != nullptr && ctx.world_history != nullptr) {
        if (ctx.world_history->world != ctx.world) {
            return StatusCode::invalid_input;
        }
    }

    return StatusCode::ok;
}

static StatusCode query_world_observer(
    const ODObserverBinding& binding,
    const ODObserverQueryContext& ctx,
    f64 t,
    ODObserverSample& out
) {
    StatusCode status;

    if (binding.source != ObserverStateSource::world) return StatusCode::invalid_input;

    if (ctx.world == nullptr) return StatusCode::invalid_input;
    if (ctx.world_center.empty() || ctx.world_frame_name.empty())
        return StatusCode::invalid_input;
    if (ctx.center != ctx.world_center || ctx.frame_name != ctx.world_frame_name)
        return StatusCode::unsupported_method;

    ODObserverSample temp;
    const Body* body = ctx.world->body(binding.observer_id);
    if (body == nullptr) return StatusCode::observer_not_found;
    if (!ctx.world->is_active(body->id)) return StatusCode::inactive_entity;
    if (body->body_type == BodyType::unknown) return StatusCode::invalid_input;

    bool at_world_epoch = std::abs(t - ctx.world->t_sim()) <= ctx.tol_time;
    if (at_world_epoch) {
        if (body->body_type == BodyType::station) {
            const Station* station = ctx.world->station(binding.observer_id);
            if (station == nullptr) return StatusCode::observer_not_found;
            if (station->anchored) {
                if (station->anchor_id == kInvalidEntityId
                    || ctx.world->body(station->anchor_id) == nullptr) {
                    return StatusCode::anchor_not_found;
                }

                if (!ctx.world->is_active(station->anchor_id)) {
                    return StatusCode::inactive_entity;
                }
            }
            temp.x = ctx.world->stat_x_tr_inertial(binding.observer_id);
        } else {
            temp.x = body->x_tr;
        }
    } else {
        if (ctx.world_history == nullptr) {
            if (binding.missing_epoch == ODObserverMissingEpochPolicy::propagate) {
                return StatusCode::unsupported_method;
            }
            return StatusCode::missing_reference;
        }
        if (body->body_type == BodyType::station) {
            status = provider_station_tr(
                *ctx.world_history,
                binding.observer_id,
                t,
                temp.x,
                ctx.tol_time
            );
        } else {
            status = provider_body_tr(
                *ctx.world_history,
                binding.observer_id,
                t,
                temp.x,
                ctx.tol_time
            );
        }
        if (status != StatusCode::ok) {
            bool missing_epoch = status == StatusCode::time_mismatch
                                 || status == StatusCode::sample_not_found
                                 || status == StatusCode::empty_history;

            if (missing_epoch
                && binding.missing_epoch == ODObserverMissingEpochPolicy::propagate) {
                return StatusCode::unsupported_method;
            }

            return status;
        }
    }

    if (!binding.known_state) return StatusCode::invalid_covariance;

    temp.observer_id = binding.observer_id;
    temp.t = t;
    temp.P = mat6d0;
    temp.center = ctx.world_center;
    temp.frame_name = ctx.world_frame_name;
    temp.source = ObserverStateSource::world;
    temp.source_id = binding.source_id;
    temp.known_state = true;

    out = temp;
    return StatusCode::ok;
}

static StatusCode query_provider_observer(
    const ODObserverBinding& binding,
    const ODObserverQueryContext& ctx,
    f64 t,
    ODObserverSample& out
) {
    StatusCode status;

    if (binding.source != ObserverStateSource::provider) return StatusCode::invalid_input;
    if (ctx.providers == nullptr) return StatusCode::missing_reference;
    auto it = ctx.providers->find(binding.source_id);
    if (it == ctx.providers->end()) return StatusCode::missing_reference;
    if (it->second == nullptr) return StatusCode::missing_reference;
    const EphemerisProvider* provider = it->second;
    if (!ctx.epoch.has_calendar_epoch) return StatusCode::invalid_input;

    if (ctx.center != provider->table.metadata.frame.center)
        return StatusCode::unsupported_method;
    if (ctx.frame_name != provider->table.metadata.frame.frame_name)
        return StatusCode::unsupported_method;

    JulianDate query_epoch;
    status = observer_query_epoch(ctx.epoch, t, query_epoch);
    if (status != StatusCode::ok) return status;

    ODObserverSample temp;

    bool covered = false;
    status = query_in_coverage(
        provider->table.metadata.epoch,
        query_epoch,
        provider->table,
        ctx.epoch.time_scale,
        ctx.epoch.offsets,
        provider->options,
        covered
    );
    if (status != StatusCode::ok) return status;
    if (!covered && binding.missing_epoch == ODObserverMissingEpochPolicy::reject) {
        // prevent hold_state or extrapolate
        return StatusCode::sample_not_found;
    }

    status = query_ephemeris_provider(
        *provider,
        query_epoch,
        ctx.epoch.time_scale,
        ctx.epoch.offsets,
        temp.x
    );
    if (status != StatusCode::ok) return status;

    if (!binding.known_state) {
        return StatusCode::invalid_covariance;
    }

    temp.observer_id = binding.observer_id;
    temp.t = t;
    temp.P = mat6d0;
    temp.center = ctx.center;
    temp.frame_name = ctx.frame_name;
    temp.source = ObserverStateSource::provider;
    temp.source_id = binding.source_id;
    temp.known_state = true;

    out = temp;
    return StatusCode::ok;
}

StatusCode query_od_observer(
    const ODObserverBinding& binding,
    const ODObserverQueryContext& ctx,
    f64 t,
    ODObserverSample& out
) {
    if (!std::isfinite(t)) {
        return StatusCode::invalid_input;
    }

    StatusCode status = validate_od_observer_binding(binding);
    if (status != StatusCode::ok) return status;

    status = validate_od_observer_query_context(ctx);
    if (status != StatusCode::ok) return status;

    ODObserverSample temp;

    switch (binding.source) {
    case ObserverStateSource::world: {
        status = query_world_observer(binding, ctx, t, temp);
    } break;
    case ObserverStateSource::provider: {
        status = query_provider_observer(binding, ctx, t, temp);
    } break;
    case ObserverStateSource::estimate: {
        // status = query_estimate_observer(binding, ctx, t, out);
        return StatusCode::unsupported_method;
    } break;
    }
    if (status != StatusCode::ok) return status;

    status = validate_od_observer_sample(temp, ctx.tol_covariance);
    if (status != StatusCode::ok) return status;

    if (temp.observer_id != binding.observer_id) return StatusCode::observer_not_found;

    if (std::abs(temp.t - t) > ctx.tol_time) return StatusCode::time_mismatch;

    if (temp.source != binding.source) return StatusCode::invalid_input;
    if (temp.source_id != binding.source_id) return StatusCode::missing_reference;

    if (temp.center != ctx.center || temp.frame_name != ctx.frame_name)
        return StatusCode::unsupported_method;

    out = temp;
    return StatusCode::ok;
}
