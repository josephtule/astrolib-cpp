// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/integrator_tableaus.hpp"

// generic (unstaged) explicit integrator steps live here


template <typename State, typename Deriv, typename Func>
inline std::pair<f64, State> step_rk1(Func&& f, f64 t, const State& x, f64 dt) {
    Deriv k1 = f(t, x);

    State x_new = x + dt * k1;
    f64 t_new = t + dt;
    return {t_new, x_new};
}

template <typename State, typename Deriv, typename Func>
inline std::pair<f64, State> step_rk2(Func&& f, f64 t, const State& x, f64 dt) {
    Deriv k1 = f(t, x);
    Deriv k2 = f(t + dt / 2.0, x + dt * k1 / 2.0);

    State x_new = x + dt * k2;
    f64 t_new = t + dt;
    return {t_new, x_new};
}

template <typename State, typename Deriv, typename Func>
inline std::pair<f64, State> step_rk2heun(Func&& f, f64 t, const State& x, f64 dt) {
    Deriv k1 = f(t, x);
    Deriv k2 = f(t + dt, x + dt * k1);

    State x_new = x + dt * (k1 + k2) / 2.0;
    f64 t_new = t + dt;
    return {t_new, x_new};
}

template <typename State, typename Deriv, typename Func>
inline std::pair<f64, State> step_rk2ralston(Func&& f, f64 t, const State& x, f64 dt) {
    Deriv k1 = f(t, x);
    Deriv k2 = f(t + dt * 2.0 / 3.0, x + dt * k1 * 2.0 / 3.0);

    State x_new = x + dt * (k1 / 4.0 + 3.0 / 4.0 * k2);
    f64 t_new = t + dt;
    return {t_new, x_new};
}

template <typename State, typename Deriv, typename Func>
inline std::pair<f64, State> step_rk3(Func&& f, f64 t, const State& x, f64 dt) {
    Deriv k1 = f(t, x);
    Deriv k2 = f(t + dt / 2.0, x + dt * k1 / 2.0);
    Deriv k3 = f(t + dt, x - dt * k1 + 2.0 * dt * k2);

    State x_new = x + dt * (k1 + 4.0 * k2 + k3) / 6.0;
    f64 t_new = t + dt;
    return {t_new, x_new};
}

template <typename State, typename Deriv, typename Func>
inline std::pair<f64, State> step_rk4(Func&& f, f64 t, const State& x, f64 dt) {
    Deriv k1 = f(t, x);
    Deriv k2 = f(t + dt / 2.0, x + k1 * (dt / 2.0));
    Deriv k3 = f(t + dt / 2.0, x + k2 * (dt / 2.0));
    Deriv k4 = f(t + dt, x + k3 * dt);

    State x_new = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    f64 t_new = t + dt;
    return {t_new, x_new};
}

template <typename State, typename Deriv, typename Func, size_t Stages>
inline std::pair<f64, State> step_rk_tableau(
    Func&& f,
    f64 t,
    const State& x,
    f64 dt,
    const RKTableau<Stages>& tableau
) {
    GenericRKResult<State, Stages> result
        = step_generic_rk<State, Deriv>(f, t, x, dt, tableau);
    return {t + dt, result.x};
}

template <typename State, typename Deriv, typename Func>
inline std::pair<f64, State> step_integrator(
    Func&& f,
    f64 t,
    const State& x,
    f64 dt,
    IntegratorTypeFixed type = IntegratorTypeFixed::rk4
) {
    std::pair<f64, State> tx;
    switch (type) {
    case IntegratorTypeFixed::rk1: tx = step_rk1<State, Deriv>(f, t, x, dt); break;
    case IntegratorTypeFixed::rk2: tx = step_rk2<State, Deriv>(f, t, x, dt); break;
    case IntegratorTypeFixed::rk2_heun: tx = step_rk2heun<State, Deriv>(f, t, x, dt); break;
    case IntegratorTypeFixed::rk2_ralston:
        tx = step_rk2ralston<State, Deriv>(f, t, x, dt);
        break;
    case IntegratorTypeFixed::rk3: tx = step_rk3<State, Deriv>(f, t, x, dt); break;
    case IntegratorTypeFixed::rk3_ralston:
        tx = step_rk_tableau<State, Deriv>(f, t, x, dt, rk3_ralston_tableau);
        break;
    case IntegratorTypeFixed::rk4: tx = step_rk4<State, Deriv>(f, t, x, dt); break;
    case IntegratorTypeFixed::rk4_38:
        tx = step_rk_tableau<State, Deriv>(f, t, x, dt, rk4_38_tableau);
        break;
    case IntegratorTypeFixed::rk5_nystrom:
        tx = step_rk_tableau<State, Deriv>(f, t, x, dt, rk5_nystrom_tableau);
        break;
    case IntegratorTypeFixed::rk6_butcher:
        tx = step_rk_tableau<State, Deriv>(f, t, x, dt, rk6_butcher_tableau);
        break;
    default: tx = step_rk4<State, Deriv>(f, t, x, dt); break;
    }
    return tx;
}
