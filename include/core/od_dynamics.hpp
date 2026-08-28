// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/body.hpp"
#include "core/dynamics_rotational.hpp"
#include "core/dynamics_translational.hpp"
#include "core/integrator_adaptive.hpp"
#include "core/integrator_fixed.hpp"
#include "core/state.hpp"
#include "core/transform.hpp"
#include "util/typedefs.hpp"
#include "util/vecdefs.hpp"

enum struct ODTrDynamicsModel : i32 { two_body, zonal };

inline ODTrDynamicsModel worldtrmodel_to_odtrmodel(GravityModel model) {
    switch (model) {
    case GravityModel::pointmass: return ODTrDynamicsModel::two_body;
    case GravityModel::zonal: return ODTrDynamicsModel::zonal;
    case GravityModel::spherical_harmonics:
        return ODTrDynamicsModel::zonal; // TODO: replace or keep
    }
}

enum struct ODAnchorAttModel : i32 { fixed, simple_spin };
enum struct ODAnchorTrModel : i32 { fixed, constant_velocity };

inline ODAnchorAttModel worldattmodel_to_odattmodel(CelestialAttitudeModel model) {
    switch (model) {
    case CelestialAttitudeModel::fixed: return ODAnchorAttModel::fixed;
    case CelestialAttitudeModel::simple_spin: return ODAnchorAttModel::simple_spin;
    case CelestialAttitudeModel::provider:
        return ODAnchorAttModel::fixed; // TODO: replace this
    }
}

struct ODDynamicsConfig {
    // what the station knows of the dynamics
    ODTrDynamicsModel tr_model = ODTrDynamicsModel::two_body;
    ODAnchorAttModel att_model = ODAnchorAttModel::fixed;
    ODAnchorTrModel anchor_tr_model = ODAnchorTrModel::fixed;
    f64 t0 = 0.0;
    f64 mu = 0.0;
    vec7d J = vec7d0;
    f64 R_cb_ref = 0.0; // dependent on supplied gravity model, usually semimajor-axis
    i32 zonal_degree = 0;
    StateTr x_cb0;
    vec4d q_cb0 = q_identity;
    vec3d w_cb = vec3d0;
    bool update_body_attitude = false;
    IntegratorType integrator = IntegratorType::rk4;
};

inline ODDynamicsConfig make_od_cfg_from_celestial(const Celestial& cel) {
    ODDynamicsConfig cfg;
    cfg.tr_model = worldtrmodel_to_odtrmodel(cel.gravity_model);
    cfg.att_model = worldattmodel_to_odattmodel(cel.attitude_model);
    cfg.mu = cel.mu;
    cfg.J = cel.J;
    cfg.R_cb_ref = cel.semimajor_axis;
    cfg.zonal_degree = cel.degree;
    cfg.x_cb0 = cel.x_tr;
    cfg.q_cb0 = cel.x_att.q;
    cfg.w_cb = cel.x_att.w;

    return cfg;
}

inline StateTr od_anchor_tr_at_time(const ODDynamicsConfig& cfg, f64 dt) {
    StateTr x_cb = cfg.x_cb0;
    if (cfg.anchor_tr_model == ODAnchorTrModel::constant_velocity) {
        x_cb.r += x_cb.v * dt;
    }
    return x_cb;
}

inline DerivTr derivtr_two_body(const StateTr& x_rel, f64 mu) {
    // x is state of target relative to the central body at the OD model origin
    DerivTr dx;
    dx.dr = x_rel.v;
    dx.dv = accel_gravity_pointmass(x_rel.r, mu);
    return dx;
}
inline DerivTr derivtr_zonal(
    const StateTr& x_rel, // in sim inertial
    f64 mu,
    f64 R_cb,
    i32 degree,
    const vec4d& q_cb, // orientation of source
    const vec7d& J
) {
    DerivTr dx;
    dx.dr = x_rel.v;
    vec3d r_bcbf = ep_rotate_fast_passive(q_cb, x_rel.r);
    vec3d a_bcbf = accel_gravity_zonal(r_bcbf, mu, R_cb, degree, J);
    dx.dv = ep_rotate_fast_passive(ep_conj(q_cb), a_bcbf);
    return dx;
}

inline DerivTr derivtr_od(f64 t, const StateTr& x, const ODDynamicsConfig& cfg) {

    switch (cfg.tr_model) {
        // NOTE: OD only considers a two body problem, third bodies are considered
        // perturbations (no staging used)
    case ODTrDynamicsModel::two_body: return derivtr_two_body(x, cfg.mu);
    case ODTrDynamicsModel::zonal: {
        vec4d q_cb = cfg.q_cb0;
        if (cfg.update_body_attitude) {
            if (cfg.att_model == ODAnchorAttModel::simple_spin) {
                q_cb = step_q_simple_spin(StateAtt{.q = q_cb, .w = cfg.w_cb}, t - cfg.t0);
            }
        }
        return derivtr_zonal(x, cfg.mu, cfg.R_cb_ref, cfg.zonal_degree, q_cb, cfg.J);
    }
    default: return DerivTr{};
    }

    // TODO: NOTE: lower fidelity model for nbody would probably be taking the state of
    // third bodies and compute position update assuming constant velocity to avoid
    // staging within estimators or even simpler count them as static bodies, may need
    // third body model to switch between, will need to test

    // same for zonal (need celestial orientation) compute simple spin update on the spot
}

inline StateTr rk4_steptr_od(
    f64 t,
    const StateTr& x,
    f64 dt,
    const ODDynamicsConfig& cfg
) {
    DerivTr k1 = derivtr_od(t, x, cfg);
    DerivTr k2 = derivtr_od(t + dt / 2.0, x + k1 * (dt / 2.0), cfg);
    DerivTr k3 = derivtr_od(t + dt / 2.0, x + k2 * (dt / 2.0), cfg);
    DerivTr k4 = derivtr_od(t + dt, x + k3 * dt, cfg);

    StateTr x_next = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    return x_next;
}

inline StateTr propagate_tr_od(
    f64 t0,
    const StateTr& x0,
    f64 t_interval,
    i32 n_steps,
    const ODDynamicsConfig& cfg
) {
    StateTr x = x0;
    if (n_steps <= 0 || t_interval == 0) return x;

    f64 dt_step = t_interval / static_cast<f64>(n_steps);
    f64 t = t0;
    auto f = [&cfg](f64 t, const StateTr& x) -> DerivTr { return derivtr_od(t, x, cfg); };
    for (i32 i = 0; i < n_steps; ++i) {
        // x = rk4_steptr_od(t, x, dt_step, cfg);
        auto tx = step_integrator<StateTr, DerivTr>(f, t, x, dt_step, cfg.integrator);
        t = tx.first;
        x = tx.second;
    }
    return x;
}

inline AdaptivePropagationResult<StateTr> propagate_tr_od_adaptive(
    f64 t0,
    const StateTr& x0,
    f64 t_interval,
    const ODDynamicsConfig& dyn_cfg,
    const AdaptiveIntegratorConfig& integrator_cfg
) {
    f64 tf = t0 + t_interval;

    AdaptivePropagationResult<StateTr> res;

    res.status = StatusCode::invalid_state;
    res.t = t0;
    res.x = x0;
    if (!finite_state_tr(x0)) return res;

    auto f = [&dyn_cfg](f64 t, const StateTr& x) -> DerivTr {
        return derivtr_od(t, x, dyn_cfg);
    };

    auto error_norm = [&integrator_cfg](
                          const StateTr& x,
                          const StateTr& x_high,
                          const DerivTr& error_delta
                      ) -> f64 {
        return adaptive_error_norm_state_tr(x, x_high, error_delta, integrator_cfg);
    };

    return propagate_dopri54<StateTr, DerivTr>(f, error_norm, t0, x0, tf, integrator_cfg);
}

// State and Derivatives for orbit determination (state + STM)
struct VarStateTr {
    StateTr x;
    mat6d Phi = mat6d1; // STM
};
struct VarDerivTr {
    DerivTr dx;
    mat6d dPhi = mat6d0;
};
inline VarStateTr operator+(const VarStateTr& y, const VarDerivTr& dy) {
    return VarStateTr{.x = y.x + dy.dx, .Phi = y.Phi + dy.dPhi};
}
inline VarStateTr operator+(const VarDerivTr& dy, const VarStateTr& y) { return y + dy; }
inline VarStateTr operator-(const VarStateTr& y, const VarDerivTr& dy) {
    return {.x = y.x - dy.dx, .Phi = y.Phi - dy.dPhi};
}
inline VarStateTr& operator+=(VarStateTr& y, const VarDerivTr& dy) {
    y.x += dy.dx;
    y.Phi += dy.dPhi;
    return y;
}
inline VarStateTr operator-(const VarStateTr& y) {
    return VarStateTr{.x = -y.x, .Phi = -y.Phi};
}
inline VarDerivTr operator+(const VarDerivTr& dy1, const VarDerivTr& dy2) {
    return VarDerivTr{.dx = dy1.dx + dy2.dx, .dPhi = dy1.dPhi + dy2.dPhi};
}
inline VarDerivTr operator-(const VarDerivTr& dy1, const VarDerivTr& dy2) {
    return VarDerivTr{.dx = dy1.dx - dy2.dx, .dPhi = dy1.dPhi - dy2.dPhi};
}
inline VarDerivTr operator*(const VarDerivTr& dy, f64 scalar) {
    return VarDerivTr{.dx = dy.dx * scalar, .dPhi = dy.dPhi * scalar};
}
inline VarDerivTr operator*(f64 scalar, const VarDerivTr& dy) { return dy * scalar; }
inline VarDerivTr operator/(const VarDerivTr& dy, f64 scalar) {
    return dy * (1.0 / scalar);
}
inline VarDerivTr operator-(const VarDerivTr& dy) {
    return VarDerivTr{.dx = -dy.dx, .dPhi = -dy.dPhi};
}

// Jacobians
inline mat6d jacobian_tr_two_body(const StateTr& x_rel, f64 mu, f64 tol = tol12) {
    mat6d J = mat6d0;
    if (mu <= 0.0) return J;

    f64 r_mag2 = x_rel.r.squaredNorm();
    f64 r_mag = std::sqrt(r_mag2);
    if (r_mag <= tol) return J;

    f64 r_mag3 = r_mag2 * r_mag;
    f64 r_mag5 = r_mag2 * r_mag2 * r_mag;

    J.block<3, 3>(0, 0) = mat3d0;
    J.block<3, 3>(0, 3) = mat3d1;

    J.block<3, 3>(3, 0)
        = -mu * (mat3d1 / r_mag3 - 3 * x_rel.r * x_rel.r.transpose() / r_mag5);
    J.block<3, 3>(3, 3) = mat3d0;

    return J;
}

// TODO: add jacobians for other perturbations, finite-diff jacobians as well?

inline mat6d jacobian_tr_zonal(
    const StateTr& x_tr_rel, // inertial
    f64 mu,
    f64 R,
    i32 degree,
    const vec4d& q, // cb orientation
    const vec7d& J,
    f64 tol = tol12
) {
    f64 R2 = R * R;
    f64 R3 = R * R2;
    f64 R4 = R2 * R2;
    f64 R5 = R2 * R3;
    f64 R6 = R3 * R3;
    f64 r_mag2 = x_tr_rel.r.squaredNorm();
    f64 r_mag = std::sqrt(r_mag2);
    if (r_mag <= tol) {
        return mat6d0;
    }
    vec3d r_bcbf = ep_rotate_fast_passive(q, x_tr_rel.r);
    mat3d R_B_N = ep_to_dcm(q);
    mat3d R_N_B = R_B_N.transpose();
    f64 x = r_bcbf(0);
    f64 y = r_bcbf(1);
    f64 z = r_bcbf(2);
    f64 x2 = x * x;
    f64 y2 = y * y;
    f64 z2 = z * z;
    f64 x4 = x2 * x2;
    f64 y4 = y2 * y2;
    f64 z4 = z2 * z2;
    f64 x6 = x2 * x4;
    f64 y6 = y2 * y4;
    f64 z6 = z2 * z4;

    // jacobians higher than j2 took too long to do by hand, see ./python/jacobians.py

    mat6d G = mat6d0;
    G.block<3, 3>(0, 0) = mat3d0;
    G.block<3, 3>(0, 3) = mat3d1;
    G.block<3, 3>(3, 3) = mat3d0;

    switch (degree) {
    case 6: {
        mat3d G_j6;
        f64 x8 = x4 * x4;
        f64 y8 = y4 * y4;
        f64 z8 = z4 * z4;
        f64 r_mag17 = std::pow(r_mag, 17);
        G_j6(0, 0)
            = (7.0 / 16.0) * J(6) * R6 * mu
              * (40 * x8 + 115 * x6 * y2 - 1235 * x6 * z2 + 105 * x4 * y4
                 - 2355 * x4 * y2 * z2 + 3480 * x4 * z4 + 25 * x2 * y6
                 - 1005 * x2 * y4 * z2 + 3360 * x2 * y2 * z4 - 1616 * x2 * z6 - 5 * y8
                 + 115 * y6 * z2 - 120 * y4 * z4 - 176 * y2 * z6 + 64 * z8)
              / r_mag17;
        G_j6(0, 1) = (315.0 / 16.0) * J(6) * R6 * mu * x * y
                     * (x6 + 3 * x4 * y2 - 30 * x4 * z2 + 3 * x2 * y4 - 60 * x2 * y2 * z2
                        + 80 * x2 * z4 + y6 - 30 * y4 * z2 + 80 * y2 * z4 - 32 * z6)
                     / r_mag17;
        G_j6(0, 2) = (63.0 / 16.0) * J(6) * R6 * mu * x * z
                     * (35 * x6 + 105 * x4 * y2 - 280 * x4 * z2 + 105 * x2 * y4
                        - 560 * x2 * y2 * z2 + 336 * x2 * z4 + 35 * y6 - 280 * y4 * z2
                        + 336 * y2 * z4 - 64 * z6)
                     / r_mag17;
        G_j6(1, 0) = (315.0 / 16.0) * J(6) * R6 * mu * x * y
                     * (x6 + 3 * x4 * y2 - 30 * x4 * z2 + 3 * x2 * y4 - 60 * x2 * y2 * z2
                        + 80 * x2 * z4 + y6 - 30 * y4 * z2 + 80 * y2 * z4 - 32 * z6)
                     / r_mag17;
        G_j6(1, 1)
            = (7.0 / 16.0) * J(6) * R6 * mu
              * (-5 * x8 + 25 * x6 * y2 + 115 * x6 * z2 + 105 * x4 * y4
                 - 1005 * x4 * y2 * z2 - 120 * x4 * z4 + 115 * x2 * y6
                 - 2355 * x2 * y4 * z2 + 3360 * x2 * y2 * z4 - 176 * x2 * z6 + 40 * y8
                 - 1235 * y6 * z2 + 3480 * y4 * z4 - 1616 * y2 * z6 + 64 * z8)
              / r_mag17;
        G_j6(1, 2) = (63.0 / 16.0) * J(6) * R6 * mu * y * z
                     * (35 * x6 + 105 * x4 * y2 - 280 * x4 * z2 + 105 * x2 * y4
                        - 560 * x2 * y2 * z2 + 336 * x2 * z4 + 35 * y6 - 280 * y4 * z2
                        + 336 * y2 * z4 - 64 * z6)
                     / r_mag17;
        G_j6(2, 0) = (63.0 / 16.0) * J(6) * R6 * mu * x * z
                     * (35 * x6 + 105 * x4 * y2 - 280 * x4 * z2 + 105 * x2 * y4
                        - 560 * x2 * y2 * z2 + 336 * x2 * z4 + 35 * y6 - 280 * y4 * z2
                        + 336 * y2 * z4 - 64 * z6)
                     / r_mag17;
        G_j6(2, 1) = (63.0 / 16.0) * J(6) * R6 * mu * y * z
                     * (35 * x6 + 105 * x4 * y2 - 280 * x4 * z2 + 105 * x2 * y4
                        - 560 * x2 * y2 * z2 + 336 * x2 * z4 + 35 * y6 - 280 * y4 * z2
                        + 336 * y2 * z4 - 64 * z6)
                     / r_mag17;
        G_j6(2, 2)
            = (7.0 / 16.0) * J(6) * R6 * mu
              * (-35 * x8 - 140 * x6 * y2 + 1120 * x6 * z2 - 210 * x4 * y4
                 + 3360 * x4 * y2 * z2 - 3360 * x4 * z4 - 140 * x2 * y6
                 + 3360 * x2 * y4 * z2 - 6720 * x2 * y2 * z4 + 1792 * x2 * z6 - 35 * y8
                 + 1120 * y6 * z2 - 3360 * y4 * z4 + 1792 * y2 * z6 - 128 * z8)
              / r_mag17;
        G.block<3, 3>(3, 0) += R_N_B * G_j6 * R_B_N;
    }
    case 5: {
        mat3d G_j5;
        f64 r_mag15 = std::pow(r_mag, 15);
        G_j5(0, 0) = (21.0 / 8.0) * J(5) * R5 * mu * z
                     * (-40 * x6 - 75 * x4 * y2 + 225 * x4 * z2 - 30 * x2 * y4
                        + 210 * x2 * y2 * z2 - 156 * x2 * z4 + 5 * y6 - 15 * y4 * z2
                        - 12 * y2 * z4 + 8 * z6)
                     / r_mag15;
        G_j5(0, 1) = (63.0 / 8.0) * J(5) * R5 * mu * x * y * z
                     * (-15 * x4 - 30 * x2 * y2 + 80 * x2 * z2 - 15 * y4 + 80 * y2 * z2
                        - 48 * z4)
                     / r_mag15;
        G_j5(0, 2)
            = (21.0 / 8.0) * J(5) * R5 * mu * x
              * (5 * x6 + 15 * x4 * y2 - 120 * x4 * z2 + 15 * x2 * y4 - 240 * x2 * y2 * z2
                 + 240 * x2 * z4 + 5 * y6 - 120 * y4 * z2 + 240 * y2 * z4 - 64 * z6)
              / r_mag15;
        G_j5(1, 0) = (63.0 / 8.0) * J(5) * R5 * mu * x * y * z
                     * (-15 * x4 - 30 * x2 * y2 + 80 * x2 * z2 - 15 * y4 + 80 * y2 * z2
                        - 48 * z4)
                     / r_mag15;
        G_j5(1, 1)
            = (21.0 / 8.0) * J(5) * R5 * mu * z
              * (5 * x6 - 30 * x4 * y2 - 15 * x4 * z2 - 75 * x2 * y4 + 210 * x2 * y2 * z2
                 - 12 * x2 * z4 - 40 * y6 + 225 * y4 * z2 - 156 * y2 * z4 + 8 * z6)
              / r_mag15;
        G_j5(1, 2)
            = (21.0 / 8.0) * J(5) * R5 * mu * y
              * (5 * x6 + 15 * x4 * y2 - 120 * x4 * z2 + 15 * x2 * y4 - 240 * x2 * y2 * z2
                 + 240 * x2 * z4 + 5 * y6 - 120 * y4 * z2 + 240 * y2 * z4 - 64 * z6)
              / r_mag15;
        G_j5(2, 0)
            = (21.0 / 8.0) * J(5) * R5 * mu * x
              * (5 * x6 + 15 * x4 * y2 - 120 * x4 * z2 + 15 * x2 * y4 - 240 * x2 * y2 * z2
                 + 240 * x2 * z4 + 5 * y6 - 120 * y4 * z2 + 240 * y2 * z4 - 64 * z6)
              / r_mag15;
        G_j5(2, 1)
            = (21.0 / 8.0) * J(5) * R5 * mu * y
              * (5 * x6 + 15 * x4 * y2 - 120 * x4 * z2 + 15 * x2 * y4 - 240 * x2 * y2 * z2
                 + 240 * x2 * z4 + 5 * y6 - 120 * y4 * z2 + 240 * y2 * z4 - 64 * z6)
              / r_mag15;
        G_j5(2, 2) = (21.0 / 8.0) * J(5) * R5 * mu * z
                     * (35 * x6 + 105 * x4 * y2 - 210 * x4 * z2 + 105 * x2 * y4
                        - 420 * x2 * y2 * z2 + 168 * x2 * z4 + 35 * y6 - 210 * y4 * z2
                        + 168 * y2 * z4 - 16 * z6)
                     / r_mag15;
        G.block<3, 3>(3, 0) += R_N_B * G_j5 * R_B_N;
    }
    case 4: {
        mat3d G_j4;
        f64 r_mag13 = std::pow(r_mag, 13);
        G_j4(0, 0)
            = (15.0 / 8.0) * J(4) * R4 * mu
              * (-6 * x6 - 11 * x4 * y2 + 101 * x4 * z2 - 4 * x2 * y4 + 90 * x2 * y2 * z2
                 - 116 * x2 * z4 + y6 - 11 * y4 * z2 - 4 * y2 * z4 + 8 * z6)
              / r_mag13;
        G_j4(0, 1) = (105.0 / 8.0) * J(4) * R4 * mu * x * y
                     * (-x4 - 2 * x2 * y2 + 16 * x2 * z2 - y4 + 16 * y2 * z2 - 16 * z4)
                     / r_mag13;
        G_j4(0, 2)
            = (105.0 / 8.0) * J(4) * R4 * mu * x * z
              * (-5 * x4 - 10 * x2 * y2 + 20 * x2 * z2 - 5 * y4 + 20 * y2 * z2 - 8 * z4)
              / r_mag13;
        G_j4(1, 0) = (105.0 / 8.0) * J(4) * R4 * mu * x * y
                     * (-x4 - 2 * x2 * y2 + 16 * x2 * z2 - y4 + 16 * y2 * z2 - 16 * z4)
                     / r_mag13;
        G_j4(1, 1) = (15.0 / 8.0) * J(4) * R4 * mu
                     * (x6 - 4 * x4 * y2 - 11 * x4 * z2 - 11 * x2 * y4 + 90 * x2 * y2 * z2
                        - 4 * x2 * z4 - 6 * y6 + 101 * y4 * z2 - 116 * y2 * z4 + 8 * z6)
                     / r_mag13;
        G_j4(1, 2)
            = (105.0 / 8.0) * J(4) * R4 * mu * y * z
              * (-5 * x4 - 10 * x2 * y2 + 20 * x2 * z2 - 5 * y4 + 20 * y2 * z2 - 8 * z4)
              / r_mag13;
        G_j4(2, 0)
            = (105.0 / 8.0) * J(4) * R4 * mu * x * z
              * (-5 * x4 - 10 * x2 * y2 + 20 * x2 * z2 - 5 * y4 + 20 * y2 * z2 - 8 * z4)
              / r_mag13;
        G_j4(2, 1)
            = (105.0 / 8.0) * J(4) * R4 * mu * y * z
              * (-5 * x4 - 10 * x2 * y2 + 20 * x2 * z2 - 5 * y4 + 20 * y2 * z2 - 8 * z4)
              / r_mag13;
        G_j4(2, 2)
            = (15.0 / 8.0) * J(4) * R4 * mu
              * (5 * x6 + 15 * x4 * y2 - 90 * x4 * z2 + 15 * x2 * y4 - 180 * x2 * y2 * z2
                 + 120 * x2 * z4 + 5 * y6 - 90 * y4 * z2 + 120 * y2 * z4 - 16 * z6)
              / r_mag13;
        G.block<3, 3>(3, 0) += R_N_B * G_j4 * R_B_N;
    }
    case 3: {
        mat3d G_j3;
        f64 r_mag11 = std::pow(r_mag, 11);
        f64 r_mag4 = std::pow(r_mag, 4);
        G_j3(0, 0)
            = -5.0 / 2.0 * J(3) * R3 * mu * z
              * (r_mag2 * (3 * x2 + 3 * y2 - 4 * z2) - 21 * x2 * (x2 + y2 - 2 * z2))
              / r_mag11;
        G_j3(0, 1)
            = (105.0 / 2.0) * J(3) * R3 * mu * x * y * z * (x2 + y2 - 2 * z2) / r_mag11;
        G_j3(0, 2) = (5.0 / 2.0) * J(3) * R3 * mu * x
                     * (-3 * r_mag2 * x2 - 3 * r_mag2 * y2 + 4 * r_mag2 * z2
                        + 35 * x2 * z2 + 35 * y2 * z2 - 28 * z4)
                     / r_mag11;
        G_j3(1, 0)
            = (105.0 / 2.0) * J(3) * R3 * mu * x * y * z * (x2 + y2 - 2 * z2) / r_mag11;
        G_j3(1, 1)
            = -5.0 / 2.0 * J(3) * R3 * mu * z
              * (r_mag2 * (3 * x2 + 3 * y2 - 4 * z2) - 21 * y2 * (x2 + y2 - 2 * z2))
              / r_mag11;
        G_j3(1, 2) = (5.0 / 2.0) * J(3) * R3 * mu * y
                     * (-3 * r_mag2 * x2 - 3 * r_mag2 * y2 + 4 * r_mag2 * z2
                        + 35 * x2 * z2 + 35 * y2 * z2 - 28 * z4)
                     / r_mag11;
        G_j3(2, 0) = -3.0 / 2.0 * J(3) * R3 * mu * x
                     * (-2 * r_mag4 + 7 * r_mag2 * (x2 + y2 - 9 * z2) + 105 * z4)
                     / r_mag11;
        G_j3(2, 1) = -3.0 / 2.0 * J(3) * R3 * mu * y
                     * (-2 * r_mag4 + 7 * r_mag2 * (x2 + y2 - 9 * z2) + 105 * z4)
                     / r_mag11;
        G_j3(2, 2) = (1.0 / 2.0) * J(3) * R3 * mu * z
                     * (-54 * r_mag4 - 21 * r_mag2 * x2 - 21 * r_mag2 * y2
                        + 329 * r_mag2 * z2 - 315 * z4)
                     / r_mag11;
        G.block<3, 3>(3, 0) += R_N_B * G_j3 * R_B_N;
    }
    case 2: {
        mat3d G_j2;
        f64 r_mag9 = std::pow(r_mag, 9);
        G_j2(0, 0) = (3.0 / 2.0) * J(2) * R2 * mu
                     * (4 * x4 + 3 * x2 * y2 - 27 * x2 * z2 - y4 + 3 * y2 * z2 + 4 * z4)
                     / r_mag9;
        G_j2(0, 1) = (15.0 / 2.0) * J(2) * R2 * mu * x * y * (x2 + y2 - 6 * z2) / r_mag9;
        G_j2(0, 2)
            = (15.0 / 2.0) * J(2) * R2 * mu * x * z * (3 * x2 + 3 * y2 - 4 * z2) / r_mag9;
        G_j2(1, 0) = (15.0 / 2.0) * J(2) * R2 * mu * x * y * (x2 + y2 - 6 * z2) / r_mag9;
        G_j2(1, 1) = (3.0 / 2.0) * J(2) * R2 * mu
                     * (-x4 + 3 * x2 * y2 + 3 * x2 * z2 + 4 * y4 - 27 * y2 * z2 + 4 * z4)
                     / r_mag9;
        G_j2(1, 2)
            = (15.0 / 2.0) * J(2) * R2 * mu * y * z * (3 * x2 + 3 * y2 - 4 * z2) / r_mag9;
        G_j2(2, 0)
            = (15.0 / 2.0) * J(2) * R2 * mu * x * z * (3 * x2 + 3 * y2 - 4 * z2) / r_mag9;
        G_j2(2, 1)
            = (15.0 / 2.0) * J(2) * R2 * mu * y * z * (3 * x2 + 3 * y2 - 4 * z2) / r_mag9;
        G_j2(2, 2)
            = (3.0 / 2.0) * J(2) * R2 * mu
              * (-3 * x4 - 6 * x2 * y2 + 24 * x2 * z2 - 3 * y4 + 24 * y2 * z2 - 8 * z4)
              / r_mag9;
        G.block<3, 3>(3, 0) += R_N_B * G_j2 * R_B_N;
    }
    default: {
        mat6d G_pointmass = jacobian_tr_two_body(x_tr_rel, mu);
        G.block<3, 3>(3, 0) += G_pointmass.block<3, 3>(3, 0);
    }
    }

    return G;
}

inline mat6d jacobian_fd_od_dynamics(
    f64 t,
    const StateTr& x,
    const ODDynamicsConfig& cfg,
    f64 eps_pos = 1e-3,
    f64 eps_vel = 1e-6
) {
    mat6d G_fd = mat6d0;
    vec6d x_vec = statetr_to_vec6d(x);
    f64 eps_base = std::sqrt(std::numeric_limits<f64>::epsilon());

    for (i32 i = 0; i < 6; ++i) {
        f64 eps_i = i < 3 ? eps_pos : eps_vel;
        vec6d x_plus_vec = x_vec;
        vec6d x_minus_vec = x_vec;
        x_plus_vec(i) += eps_i;
        x_minus_vec(i) -= eps_i;

        StateTr x_plus = vec6d_to_statetr(x_plus_vec);
        StateTr x_minus = vec6d_to_statetr(x_minus_vec);
        vec6d f_plus = derivtr_to_vec6d(derivtr_od(t, x_plus, cfg));
        vec6d f_minus = derivtr_to_vec6d(derivtr_od(t, x_minus, cfg));

        G_fd.col(i) = (f_plus - f_minus) / (2.0 * eps_i);
    }

    return G_fd;
}

inline mat6d jacobian_tr_od(f64 t, const StateTr& x, const ODDynamicsConfig& cfg) {
    // jacobian dispatcher
    switch (cfg.tr_model) {
    case ODTrDynamicsModel::two_body: return jacobian_tr_two_body(x, cfg.mu);
    case ODTrDynamicsModel::zonal: {
        vec4d q_cb = cfg.q_cb0;
        if (cfg.update_body_attitude) {
            if (cfg.att_model == ODAnchorAttModel::simple_spin) {
                q_cb = step_q_simple_spin(StateAtt{.q = q_cb, .w = cfg.w_cb}, t - cfg.t0);
            }
        }
        return jacobian_tr_zonal(x, cfg.mu, cfg.R_cb_ref, cfg.zonal_degree, q_cb, cfg.J);
    }
    default: return jacobian_fd_od_dynamics(t, x, cfg);
    }
}

inline VarDerivTr deriv_var_tr_od(
    f64 t,
    const VarStateTr& y,
    const ODDynamicsConfig& cfg
) {
    mat6d J = jacobian_tr_od(t, y.x, cfg);
    VarDerivTr dy;
    dy.dx = derivtr_od(t, y.x, cfg);
    dy.dPhi = J * y.Phi;
    return dy;
}

inline VarStateTr rk4_step_var_tr_od(
    f64 t,
    const VarStateTr& y,
    f64 dt,
    const ODDynamicsConfig& cfg
) {
    VarDerivTr k1 = deriv_var_tr_od(t, y, cfg);
    VarDerivTr k2 = deriv_var_tr_od(t + dt / 2.0, y + k1 * (dt / 2.0), cfg);
    VarDerivTr k3 = deriv_var_tr_od(t + dt / 2.0, y + k2 * (dt / 2.0), cfg);
    VarDerivTr k4 = deriv_var_tr_od(t + dt, y + k3 * dt, cfg);

    VarStateTr x_next = y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    return x_next;
}

inline VarStateTr propagate_var_tr_od(
    f64 t0,
    const VarStateTr& y0,
    f64 t_interval,
    i32 n_steps,
    const ODDynamicsConfig& cfg
) {
    VarStateTr y = y0;
    if (n_steps <= 0 || t_interval == 0.0) return y;

    f64 dt_step = t_interval / static_cast<f64>(n_steps);
    f64 t = t0;
    auto f = [&cfg](f64 t, const VarStateTr& y) -> VarDerivTr {
        return deriv_var_tr_od(t, y, cfg);
    };
    for (i32 i = 0; i < n_steps; ++i) {
        auto ty
            = step_integrator<VarStateTr, VarDerivTr>(f, t, y, dt_step, cfg.integrator);
        t = ty.first;
        y = ty.second;
    }
    return y;
}

inline VarStateTr propagate_var_tr_od_adaptive(
    f64 t0,
    const VarStateTr& y0,
    f64 t_interval,
    const ODDynamicsConfig& cfg,
    const AdaptiveIntegratorConfig& adaptive_cfg
) {
    // TODO: add guards
    VarStateTr y = y0;
    if (t_interval == 0.0) return y;

    f64 t = t0;
    auto f = [&cfg](f64 t, const VarStateTr& y) -> VarDerivTr {
        return deriv_var_tr_od(t, y, cfg);
    };
    // TODO: complete this

    // propagate_dopri54<VarStateTr, VarDerivTr>(f, error_nor_fn, t, y0, t0+t_interval,
    // adaptive_cfg);

    return y;
}