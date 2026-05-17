#pragma once

#include "core/entity.hpp"
#include "core/observation_type.hpp"
#include "core/observations.hpp"
#include "core/state.hpp"
#include "core/station_geometry.hpp"
#include "util/constants.hpp"
#include "util/typedefs.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

struct Measurement {
    f64 t = 0.0;
    ObservationType type = ObservationType::radec;
    EntityId observer_id = kInvalidEntityId;
    EntityId target_id = kInvalidEntityId;
    vecXd z; // measured value
    matXd R; // measurement covariance
};

enum struct AzelInputFrame : i32 { enu, bcbf };
struct MeasurementContext {
    StateTr x_tr_observer;
    StateTr x_tr_target;

    // local station geometry, used only for azel
    bool has_station_local = false;
    AzelInputFrame azel_frame = AzelInputFrame::enu;
    vec3d station_llh = vec3d0; // [lat, lon, h], angle units from angle_in
};

inline MeasurementContext make_measurement_context(
    const StateTr& x_tr_target,
    const StateTr& x_tr_observer
) {
    // for non-azel measurements only, TODO: create one for azel

    MeasurementContext ctx;
    ctx.x_tr_observer = x_tr_observer;
    ctx.x_tr_target = x_tr_target;

    return ctx;
}

inline i32 measurement_dim(ObservationType type) {
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

inline vecXd predict_measurement(
    ObservationType type,
    const MeasurementContext& ctx,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
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
        if (!ctx.has_station_local) {
            measurement.resize(0);
            return vecXd{};
        }
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
        measurement = statetr_to_vec6(x_target);
    } break;
    case ObservationType::rel_pos: {
        measurement = x_target.r - x_observer.r;
    } break;
    case ObservationType::rel_pos_vel:
        measurement = statetr_to_vec6(x_target - x_observer);
    }

    return measurement;
}

inline vecXd predict_measurement(
    ObservationType type,
    const StateTr& x_target,
    const StateTr& x_observer,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
) {
    if (type == ObservationType::azel)
        return vecXd{}; // NOTE: azel unsupported in this overload

    MeasurementContext ctx;
    ctx.x_tr_target = x_target;
    ctx.x_tr_observer = x_observer;
    return predict_measurement(type, ctx, UAngle::degree, angle_out, tol);
}

inline vecXd measurement_residual(
    ObservationType type,
    ecref<vecXd> z_obs,
    ecref<vecXd> z_pred,
    UAngle angle_in = UAngle::radian
) {
    vecXd residual = z_obs - z_pred;

    if (type == ObservationType::radec || type == ObservationType::azel) {
        i32 i_max = std::min(2, static_cast<i32>(residual.size()));
        f64 wrap_min = convert_angle(-pi, UAngle::radian, angle_in);
        f64 wrap_max = convert_angle(pi, UAngle::radian, angle_in);
        for (i32 i = 0; i < i_max; ++i) {
            residual(i) = wrap_angle(residual(i), wrap_min, wrap_max, angle_in, angle_in);
        }
    }

    return residual;
}

inline matXd jacobian_fd_measurement(
    ObservationType type,
    const MeasurementContext& ctx,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 eps_pos = 1e-3,
    f64 eps_vel = 1e-6,
    f64 tol = tol12
) {
    i32 rows = measurement_dim(type);
    i32 cols = 6;
    matXd H(rows, cols);

    vecXd z0 = predict_measurement(type, ctx, angle_in, angle_out, tol);
    if (z0.size() != rows) return matXd{};

    MeasurementContext ctx_pert = ctx;
    for (i32 i = 0; i < cols; ++i) {
        // forward differencing
        f64 eps_i = i < 3 ? eps_pos : eps_vel;

        vec6d x_target_pert_vec = statetr_to_vec6(ctx.x_tr_target);
        x_target_pert_vec(i) += eps_i;
        StateTr x_target_pert = vec6_to_statetr(x_target_pert_vec);
        ctx_pert.x_tr_target = x_target_pert;
        vecXd z_pert = predict_measurement(type, ctx_pert, angle_in, angle_out, tol);
        if (z_pert.size() != rows) return matXd{};
        vecXd residual = measurement_residual(type, z_pert, z0, angle_out);
        H.col(i) = residual / eps_i;
    }

    return H;
}

inline matXd jacobian_radec(
    const MeasurementContext& ctx,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
) {
    vec3d r_rel = ctx.x_tr_target.r - ctx.x_tr_observer.r;
    f64 rho2 = r_rel.squaredNorm();
    f64 rho = std::sqrt(rho2);
    f64 rhoxy2 = r_rel.segment<2>(0).squaredNorm();
    f64 rhoxy = std::sqrt(rhoxy2);
    f64 denom = rhoxy * rho2;

    if (rho2 <= tol || rhoxy2 <= tol || denom <= tol) return matXd{};

    matd<2, 6> H;
    H.row(0) = vec6d{-r_rel(1) / rhoxy2, r_rel(0) / rhoxy2, 0.0, 0.0, 0.0, 0.0};
    H.row(1) = vec6d{
        -r_rel(2) * r_rel(0) / denom,
        -r_rel(2) * r_rel(1) / denom,
        rhoxy / rho2,
        0.0,
        0.0,
        0.0
    };
    f64 converter = convert_angle(1.0, UAngle::radian, angle_out);
    H *= converter;

    return H;
}

inline matXd jacobian_azel_enu(
    const MeasurementContext& ctx,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
) {
    vec3d r_rel_enu = ctx.x_tr_target.r - ctx.x_tr_observer.r;
    f64 e = r_rel_enu(0);
    f64 n = r_rel_enu(1);
    f64 u = r_rel_enu(2);

    f64 rho2 = r_rel_enu.squaredNorm();
    f64 rhoen2 = r_rel_enu.segment<2>(0).squaredNorm();
    f64 rhoen = std::sqrt(rhoen2);
    f64 denom = rhoen * rho2;

    if (rho2 <= tol || rhoen2 <= tol || denom <= tol) return matXd{};

    matd<2, 6> H;
    H.row(0) = vec6d{n / rhoen2, -e / rhoen2, 0.0, 0.0, 0.0, 0.0};
    H.row(1) = vec6d{-u * e / denom, -u * n / denom, rhoen / rho2, 0.0, 0.0, 0.0};
    f64 converter = convert_angle(1.0, UAngle::radian, angle_out);
    H *= converter;

    return H;
}

inline matXd jacobian_azel_inertial_from_enu(
    const MeasurementContext& ctx,
    const mat3d& R_ENU_I,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
) {
    matXd H_ENU = jacobian_azel_enu(ctx, angle_out, tol);
    matXd H = H_ENU;
    if (H.rows() != 2 || H.cols() != 6) return H;
    H.block<2, 3>(0, 0) = H.block<2, 3>(0, 0) * R_ENU_I;

    return H;
}

inline matXd jacobian_azel_bcbf(
    const MeasurementContext& ctx,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
) {
    f64 lat = ctx.station_llh(0);
    f64 lon = ctx.station_llh(1);
    mat3d R_ENU_BCBF = rot_enu_from_bcbf(lat, lon, angle_in);
    MeasurementContext enu_ctx = ctx;
    enu_ctx.x_tr_observer = R_ENU_BCBF * enu_ctx.x_tr_observer;
    enu_ctx.x_tr_target = R_ENU_BCBF * enu_ctx.x_tr_target;

    matXd H = jacobian_azel_enu(enu_ctx, angle_out, tol);
    if (H.rows() != 2 || H.cols() != 6) return H;
    H.block<2, 3>(0, 0) = H.block<2, 3>(0, 0) * R_ENU_BCBF;

    return H;
}

inline matXd jacobian_azel_inertial_from_bcbf(
    const MeasurementContext& ctx,
    const mat3d& R_BCBF_I,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
) {
    matXd H_BCBF = jacobian_azel_bcbf(ctx, angle_in, angle_out, tol);
    matXd H = H_BCBF;
    if (H.rows() != 2 || H.cols() != 6) return H;
    H.block<2, 3>(0, 0) = H.block<2, 3>(0, 0) * R_BCBF_I;

    return H;
}

inline matXd measurement_jacobian(
    ObservationType type,
    const MeasurementContext& ctx,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 eps_pos = 1e-3,
    f64 eps_vel = 1e-6,
    f64 tol = tol12
) {
    i32 rows = measurement_dim(type);
    i32 cols = 6;
    if (rows == 0) return matXd{};

    matXd H(rows, cols);

    switch (type) {
    case ObservationType::radec: {
        H = jacobian_radec(ctx, angle_out, tol);
    } break;
    // TODO: add jacobians for the below if possible
    case ObservationType::azel: [[fallthrough]]; // only works in world context
    case ObservationType::range: [[fallthrough]];
    case ObservationType::range_rate: [[fallthrough]];
    case ObservationType::pos: [[fallthrough]];
    case ObservationType::pos_vel: [[fallthrough]];
    case ObservationType::rel_pos: [[fallthrough]];
    case ObservationType::rel_pos_vel: [[fallthrough]];
    default:
        H = jacobian_fd_measurement(
            type,
            ctx,
            angle_in,
            angle_out,
            eps_pos,
            eps_vel,
            tol
        );
    }
    return H;
}
