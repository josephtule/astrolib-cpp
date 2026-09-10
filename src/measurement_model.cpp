// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/measurement_model.hpp"

#include "core/observations.hpp"
#include "core/station_geometry.hpp"
#include "util/math.hpp"
#include "util/units.hpp"

MeasurementContext make_measurement_context(
    const StateTr& x_tr_target,
    const StateTr& x_tr_observer
) {
    // for non-azel measurements only, TODO: create one for azel
    MeasurementContext ctx;
    ctx.x_tr_observer = x_tr_observer;
    ctx.x_tr_target = x_tr_target;
    return ctx;
}

i32 measurement_dim(ObservationType type) {
    switch (type) {
    case ObservationType::radec: return 2;
    case ObservationType::azel: return 2;
    case ObservationType::range: return 1;
    case ObservationType::range_rate: return 1;
    case ObservationType::rel_pos: [[fallthrough]];
    case ObservationType::pos: return 3;
    case ObservationType::rel_pos_vel: [[fallthrough]];
    case ObservationType::pos_vel: return 6;
    default: return 0;
    }
}

vecXd predict_measurement(
    ObservationType type,
    const MeasurementContext& ctx,
    UAngle angle_in,
    UAngle angle_out,
    f64 tol
) {
    const StateTr& x_target = ctx.x_tr_target;
    const StateTr& x_observer = ctx.x_tr_observer;
    i32 dim = measurement_dim(type);
    vecXd measurement(dim);

    switch (type) {
    case ObservationType::radec: {
        vec3d radec = radec_from_pos(x_target.r, x_observer.r, angle_out, tol);
        measurement = radec.segment(0, dim);
    } break;
    case ObservationType::azel: {
        // NOTE: requires local context
        if (!ctx.has_station_local) return vecXd{};
        switch (ctx.azel_frame) {
        case AzelInputFrame::enu: {
            vec3d azel
                = azel_from_enu(ctx.x_tr_target.r - ctx.x_tr_observer.r, angle_out, tol);
            measurement = azel.segment(0, dim);
        } break;
        case AzelInputFrame::bcbf: {
            vec3d azel = azel_from_bcbf(
                ctx.x_tr_target.r,
                ctx.x_tr_observer.r,
                ctx.station_llh(0),
                ctx.station_llh(1),
                angle_in,
                angle_out,
                tol
            );
            measurement = azel.segment(0, dim);
        } break;
        }
    } break;
    case ObservationType::range: {
        measurement << (x_target.r - x_observer.r).norm();
    } break;
    case ObservationType::range_rate: {
        StateTr x_rel = x_target - x_observer;
        f64 rho = x_rel.r.norm();
        if (rho <= tol) return measurement.setZero();
        measurement << x_rel.r.dot(x_rel.v) / rho;
    } break;
    case ObservationType::pos: {
        measurement = x_target.r;
    } break;
    case ObservationType::pos_vel: {
        measurement = statetr_to_vec6d(x_target);
    } break;
    case ObservationType::rel_pos: {
        measurement = x_target.r - x_observer.r;
    } break;
    case ObservationType::rel_pos_vel: {
        measurement = statetr_to_vec6d(x_target - x_observer);
    } break;
    }

    return measurement;
}

vecXd predict_measurement(
    ObservationType type,
    const StateTr& x_target,
    const StateTr& x_observer,
    UAngle angle_out,
    f64 tol
) {
    if (type == ObservationType::azel) return vecXd{}; // NOTE: unsupported here

    MeasurementContext ctx;
    ctx.x_tr_target = x_target;
    ctx.x_tr_observer = x_observer;
    return predict_measurement(type, ctx, UAngle::degree, angle_out, tol);
}

vecXd measurement_residual(
    ObservationType type,
    ecref<vecXd> z_obs,
    ecref<vecXd> z_pred,
    UAngle angle_in
) {
    vecXd residual = z_obs - z_pred;

    if (type == ObservationType::radec || type == ObservationType::azel) {
        f64 wrap_min = convert_angle(-pi, UAngle::radian, angle_in);
        f64 wrap_max = convert_angle(pi, UAngle::radian, angle_in);
        residual(0) = wrap_angle(residual(0), wrap_min, wrap_max, angle_in, angle_in);
    }

    return residual;
}
