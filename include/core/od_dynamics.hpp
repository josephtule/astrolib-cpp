#pragma once

#include "core/dynamics_translational.hpp"
#include "core/integrator.hpp"
#include "core/state.hpp"
#include "util/typedefs.hpp"

enum struct ODDynamicsModel : i32 {
    two_body,
    zonal,
    n_body,
};

struct ODDynamicsConfig {
    ODDynamicsModel model = ODDynamicsModel::two_body;
    f64 mu = 0.0;
    f64 body_radius = 0.0;
    i32 zonal_degree = 0;
    IntegratorType integrator = IntegratorType::rk4;
};

inline DerivTr derivtr_two_body(const StateTr& x_rel, f64 mu) {
    // x is state of target relative to the central body at the OD model origin
    DerivTr dx;
    dx.dr = x_rel.v;
    dx.dv = accel_gravity_pointmass(x_rel.r, mu);
    return dx;
}

inline DerivTr derivtr_od(f64 t, const StateTr& x, const ODDynamicsConfig& cfg) {
    switch (cfg.model) {
        // NOTE: OD only considers a two body problem, third bodies are considered
        // perturbations (no staging used)
    case ODDynamicsModel::two_body: return derivtr_two_body(x, cfg.mu);
    default: return DerivTr{};
    }
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
    f64 t_span,
    i32 n_steps,
    const ODDynamicsConfig& cfg
) {
    StateTr x = x0;
    if (n_steps <= 0 || t_span == 0) return x;

    f64 dt_step = t_span / static_cast<f64>(n_steps);
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
inline mat6d jacobian_tr_two_body(const StateTr& x_rel, f64 mu, f64 tol = tol_strict) {
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

// TODO: add jacobians for other perturbations

inline mat6d jacobian_tr_od(f64 t, const StateTr& x, const ODDynamicsConfig& cfg) {
    // jacobian dispatcher
    switch (cfg.model) {
    case ODDynamicsModel::two_body: return jacobian_tr_two_body(x, cfg.mu);
    default:
        return mat6d0;
        // TODO: populate these
        // case ODDynamicsModel::zonal:
        // case ODDynamicsModel::n_body: break;
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
    if (n_steps <= 0 || t_interval == 0) return y;

    f64 dt_step = t_interval / static_cast<f64>(n_steps);
    f64 t = t0;
    auto f = [&cfg](f64 t, const VarStateTr& y) -> VarDerivTr {
        return deriv_var_tr_od(t, y, cfg);
    };
    for (i32 i = 0; i < n_steps; ++i) {
        // y = rk4_step_var_tr_od(t, y, dt_step, cfg);
        auto ty
            = step_integrator<VarStateTr, VarDerivTr>(f, t, y, dt_step, cfg.integrator);
        t = ty.first;
        y = ty.second;
    }
    return y;
}
