#pragma once

#include "util/typedefs.hpp"

// generic (unstaged) explicit integrator steps live here

enum struct IntegratorType : i32 {
    rk1,
    rk2,
    rk2_heun,
    rk2_ralston,
    rk3,
    rk4,
};
inline std::string integrator_name(IntegratorType integrator) {
    switch (integrator) {
    case IntegratorType::rk1: return "RK1";
    case IntegratorType::rk2: return "RK2";
    case IntegratorType::rk2_heun: return "RK2 Heun";
    case IntegratorType::rk2_ralston: return "RK2 Ralston";
    case IntegratorType::rk3: return "RK3";
    case IntegratorType::rk4: return "RK4";
    }
    return "Unknown";
};
inline string integrator_str(IntegratorType type){
    switch (type) {
    case IntegratorType::rk1: return "rk1";
    case IntegratorType::rk2: return "rk2";
    case IntegratorType::rk2_heun: return "rk2_heun";
    case IntegratorType::rk2_ralston: return "rk2_ralston";
    case IntegratorType::rk3: return "rk3";
    case IntegratorType::rk4: return "rk4";
    }
    return "unknown";
}

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
    Deriv k2 = f(t + dt * 3.0 / 4.0, x + dt * k1 * 3.0 / 4.0);

    State x_new = x + dt * (k1 / 3.0 + 2.0 / 3.0 * k2);
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

template <typename State, typename Deriv, typename Func>
inline std::pair<f64, State> step_integrator(
    Func&& f,
    f64 t,
    const State& x,
    f64 dt,
    IntegratorType type = IntegratorType::rk4
) {
    std::pair<f64, State> tx;
    switch (type) {
    case IntegratorType::rk1: tx = step_rk1<State, Deriv>(f, t, x, dt); break;
    case IntegratorType::rk2: tx = step_rk2<State, Deriv>(f, t, x, dt); break;
    case IntegratorType::rk2_heun: tx = step_rk2heun<State, Deriv>(f, t, x, dt); break;
    case IntegratorType::rk2_ralston:
        tx = step_rk2ralston<State, Deriv>(f, t, x, dt);
        break;
    case IntegratorType::rk3: tx = step_rk3<State, Deriv>(f, t, x, dt); break;
    case IntegratorType::rk4: tx = step_rk4<State, Deriv>(f, t, x, dt); break;
    default: tx = step_rk4<State, Deriv>(f, t, x, dt); break;
    }
    return tx;
}