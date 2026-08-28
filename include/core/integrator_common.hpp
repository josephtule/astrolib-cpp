#pragma once

#include "core/status.hpp"
#include "util/typedefs.hpp"

#include <array>
#include <cstddef>

enum struct IntegratorType : i32 {
    rk1,
    rk2,
    rk2_heun,
    rk2_ralston,
    rk3,
    rk4,
};

inline string integrator_name(IntegratorType integrator) {
    switch (integrator) {
    case IntegratorType::rk1: return "RK1";
    case IntegratorType::rk2: return "RK2";
    case IntegratorType::rk2_heun: return "RK2 Heun";
    case IntegratorType::rk2_ralston: return "RK2 Ralston";
    case IntegratorType::rk3: return "RK3";
    case IntegratorType::rk4: return "RK4";
    }
    return "Unknown";
}

inline string integrator_str(IntegratorType type) {
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

template <size_t Stages>
struct RKTableau {
    std::array<f64, Stages> c{};
    std::array<std::array<f64, Stages>, Stages> a{};
    std::array<f64, Stages> b_high{};
    std::array<f64, Stages> b_low{};

    i32 order_high = 0;
    i32 order_low = 0;
    bool fsal = false;
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

    array<Deriv, Stages> k
        = rk_generic_stages<State, Deriv>(f, t, x, dt, tableau);

    result.x = x;
    for (std::size_t i = 0; i < Stages; ++i) {
        result.x += dt * tableau.b_high[i] * k[i];
    }

    result.status = StatusCode::ok;
    return result;
}

inline constexpr RKTableau<1> rk1_tableau{
    // c
    {0.0},

    // a
    {{{{0.0}}}},

    // b_high
    {1.0},

    // b_low
    {0.0},

    1,
    0,
    false
};

inline constexpr RKTableau<2> rk2_tableau{
    // c
    {0.0, 1.0 / 2.0},

    // a
    {{{{0.0, 0.0}}, {{1.0 / 2.0, 0.0}}}},

    // b_high
    {0.0, 1.0},

    // b_low
    {0.0, 0.0},

    2,
    0,
    false
};

inline constexpr RKTableau<2> rk2_heun_tableau{
    // c
    {0.0, 1.0},

    // a
    {{{{0.0, 0.0}}, {{1.0, 0.0}}}},

    // b_high
    {1.0 / 2.0, 1.0 / 2.0},

    // b_low
    {0.0, 0.0},

    2,
    0,
    false
};

inline constexpr RKTableau<2> rk2_ralston_tableau{
    // c
    {0.0, 3.0 / 4.0},

    // a
    {{{{0.0, 0.0}}, {{3.0 / 4.0, 0.0}}}},

    // b_high
    {1.0 / 3.0, 2.0 / 3.0},

    // b_low
    {0.0, 0.0},

    2,
    0,
    false
};

inline constexpr RKTableau<3> rk3_tableau{
    // c
    {0.0, 1.0 / 2.0, 1.0},

    // a
    {{{{0.0, 0.0, 0.0}}, {{1.0 / 2.0, 0.0, 0.0}}, {{-1.0, 2.0, 0.0}}}},

    // b_high
    {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},

    // b_low
    {0.0, 0.0, 0.0},

    3,
    0,
    false
};

inline constexpr RKTableau<4> rk4_tableau{
    // c
    {0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0},

    // a
    {{{{0.0, 0.0, 0.0, 0.0}},
      {{1.0 / 2.0, 0.0, 0.0, 0.0}},
      {{0.0, 1.0 / 2.0, 0.0, 0.0}},
      {{0.0, 0.0, 1.0, 0.0}}}},

    // b_high
    {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0},

    // b_low
    {0.0, 0.0, 0.0, 0.0},

    4,
    0,
    false
};

inline constexpr RKTableau<7> dopri54_tableau{
    // c
    {0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0},

    // a
    {{{{0.0}},
      {{1.0 / 5.0}},
      {{3.0 / 40.0, 9.0 / 40.0}},
      {{44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0}},
      {{19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0}},
      {{9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0}},
      {{35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0}}}},

    // b_high
    {35.0 / 384.0,
     0.0,
     500.0 / 1113.0,
     125.0 / 192.0,
     -2187.0 / 6784.0,
     11.0 / 84.0,
     0.0},

    // b_low
    {5179.0 / 57600.0,
     0.0,
     7571.0 / 16695.0,
     393.0 / 640.0,
     -92097.0 / 339200.0,
     187.0 / 2100.0,
     1.0 / 40.0},

    5,
    4,
    true
};
