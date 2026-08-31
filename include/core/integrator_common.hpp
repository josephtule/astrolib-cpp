#pragma once

#include "core/status.hpp"
#include "util/typedefs.hpp"

#include <array>
#include <cstddef>
#include <variant>

enum struct IntegratorFamily : i32 {
    fixed,
    adaptive,
};

inline string integrator_name(IntegratorFamily family) {
    switch (family) {
    case IntegratorFamily::fixed: return "Fixed-step";
    case IntegratorFamily::adaptive: return "Adaptive";
    }
    return "Unknown";
}
inline string integrator_str(IntegratorFamily family) {
    switch (family) {
    case IntegratorFamily::fixed: return "fixed";
    case IntegratorFamily::adaptive: return "adaptive";
    }
    return "unknown";
}

enum struct IntegratorTypeFixed : i32 {
    rk1,
    rk2,
    rk2_heun,
    rk2_ralston,
    rk3,
    rk3_ralston,
    rk4,
    rk4_38,
    rk5_nystrom,
    rk6_butcher,
};

inline string integrator_name(IntegratorTypeFixed type) {
    switch (type) {
    case IntegratorTypeFixed::rk1: return "RK1";
    case IntegratorTypeFixed::rk2: return "RK2";
    case IntegratorTypeFixed::rk2_heun: return "RK2 Heun";
    case IntegratorTypeFixed::rk2_ralston: return "RK2 Ralston";
    case IntegratorTypeFixed::rk3: return "RK3";
    case IntegratorTypeFixed::rk3_ralston: return "RK3 Ralston";
    case IntegratorTypeFixed::rk4: return "RK4";
    case IntegratorTypeFixed::rk4_38: return "RK4 3/8 Rule";
    case IntegratorTypeFixed::rk5_nystrom: return "RK5 Nystrom";
    case IntegratorTypeFixed::rk6_butcher: return "RK6 Butcher";
    }
    return "Unknown";
}

inline string integrator_str(IntegratorTypeFixed type) {
    switch (type) {
    case IntegratorTypeFixed::rk1: return "rk1";
    case IntegratorTypeFixed::rk2: return "rk2";
    case IntegratorTypeFixed::rk2_heun: return "rk2_heun";
    case IntegratorTypeFixed::rk2_ralston: return "rk2_ralston";
    case IntegratorTypeFixed::rk3: return "rk3";
    case IntegratorTypeFixed::rk3_ralston: return "rk3_ralston";
    case IntegratorTypeFixed::rk4: return "rk4";
    case IntegratorTypeFixed::rk4_38: return "rk4_38";
    case IntegratorTypeFixed::rk5_nystrom: return "rk5_nystrom";
    case IntegratorTypeFixed::rk6_butcher: return "rk6_butcher";
    }
    return "unknown";
}

enum struct IntegratorTypeAdaptive : i32 {
    rkf12,
    heuneuler21,
    bosha32,
    rkf54,
    cashkarp54,
    dopri54,
    rkf78,
};

inline string integrator_name(IntegratorTypeAdaptive type) {
    switch (type) {
    case IntegratorTypeAdaptive::rkf12: return "Runge-Kutta-Fehlberg 2(1)";
    case IntegratorTypeAdaptive::heuneuler21: return "Heun-Euler 2(1)";
    case IntegratorTypeAdaptive::bosha32: return "Bogacki-Shampine 3(2)";
    case IntegratorTypeAdaptive::rkf54: return "Runge-Kutta-Fehlberg 5(4)";
    case IntegratorTypeAdaptive::cashkarp54: return "Cash-Karp 5(4)";
    case IntegratorTypeAdaptive::dopri54: return "Dormand-Prince 5(4)";
    case IntegratorTypeAdaptive::rkf78: return "Runge-Kutta-Fehlberg 8(7)";
    }
    return "Unknown";
}

inline string integrator_str(IntegratorTypeAdaptive type) {
    switch (type) {
    case IntegratorTypeAdaptive::rkf12: return "rkf12";
    case IntegratorTypeAdaptive::heuneuler21: return "heun_euler21";
    case IntegratorTypeAdaptive::bosha32: return "bogacki_shampine32";
    case IntegratorTypeAdaptive::rkf54: return "rkf54";
    case IntegratorTypeAdaptive::cashkarp54: return "cash_karp54";
    case IntegratorTypeAdaptive::dopri54: return "dormand_prince54";
    case IntegratorTypeAdaptive::rkf78: return "rkf78";
    }
    return "unknown";
}

using IntegratorType = std::variant<IntegratorTypeFixed, IntegratorTypeAdaptive>;

inline IntegratorFamily integrator_family(const IntegratorType& type) {
    return std::holds_alternative<IntegratorTypeFixed>(type) ? IntegratorFamily::fixed
                                                             : IntegratorFamily::adaptive;
}

inline string integrator_name(const IntegratorType& type) {
    return std::visit([](auto method) { return integrator_name(method); }, type);
}

inline string integrator_str(const IntegratorType& type) {
    return std::visit([](auto method) { return integrator_str(method); }, type);
}

template <size_t Stages>
struct RKTableau {
    std::array<f64, Stages> c{};
    std::array<std::array<f64, Stages>, Stages> a{};
    std::array<f64, Stages> b_high{};
    std::array<f64, Stages> b_low{};

    i32 order_high = 0;
    i32 order_low = 0;
    bool fsal = false; // first-same-as-last
    bool embedded = false;
};

template <class State, class Deriv, class Func, size_t Stages>
array<Deriv, Stages> rk_generic_stages(
    Func&& f,
    f64 t,
    const State& x,
    f64 dt,
    const RKTableau<Stages>& tableau
) {
    // NOTE: unoptimized, requires a bit more multiplication then necessary, maybe
    // compiler optimizes it away anyways

    array<Deriv, Stages> k{};

    for (size_t i = 0; i < Stages; ++i) {
        State x_stage = x;

        for (size_t j = 0; j < i; ++j) {
            x_stage += dt * tableau.a[i][j] * k[j];
        }

        f64 t_stage = t + tableau.c[i] * dt;
        k[i] = f(t_stage, x_stage);
    }

    return k;
}

template <class State, size_t Stages>
struct GenericRKResult {
    StatusCode status = StatusCode::invalid_state;
    State x{};
};

template <class State, class Deriv, class Func, size_t Stages>
GenericRKResult<State, Stages> step_generic_rk(
    Func&& f,
    f64 t,
    const State& x,
    f64 dt,
    const RKTableau<Stages>& tableau
) {
    GenericRKResult<State, Stages> result{};

    array<Deriv, Stages> k = rk_generic_stages<State, Deriv>(f, t, x, dt, tableau);

    result.x = x;
    for (std::size_t i = 0; i < Stages; ++i) {
        result.x += dt * tableau.b_high[i] * k[i];
    }

    result.status = StatusCode::ok;
    return result;
}
