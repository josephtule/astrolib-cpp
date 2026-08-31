#pragma once

#include "core/status.hpp"
#include "util/typedefs.hpp"

#include <array>
#include <cassert>
#include <cstddef>

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
    {{
        {{0.0, 0.0}},      // row 0
        {{1.0 / 2.0, 0.0}} // row 1
    }},

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
    {{
        {{0.0, 0.0}}, // row 0
        {{1.0, 0.0}}  // row 1
    }},

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
    {0.0, 2.0 / 3.0},

    // a
    {{{{0.0, 0.0}}, {{2.0 / 3.0, 0.0}}}},

    // b_high
    {1.0 / 4.0, 3.0 / 4.0},

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

inline constexpr RKTableau<2> heuneuler21_tableau{
    // c
    {0.0, 1.0},

    // a
    {{
        {0.0, 0.0}, // row 0
        {1.0, 0.0}  // row 1
    }},

    // b_high
    {1.0 / 2.0, 1.0 / 2.0},

    // b_low
    {1.0, 0.0},

    2,
    1,
    false,
    true
};

// aka rkf12
inline constexpr RKTableau<3> rkf21_tableau{
    // c
    {0.0, 1.0 / 2.0, 1.0},

    // a
    {{
        {0.0},                       // row 0
        {1.0 / 2.0},                 // row 1
        {1.0 / 256.0, 255.0 / 256.0} // row 2
    }},

    // b_high
    {1.0 / 512.0, 255.0 / 256.0, 1.0 / 512.0},

    // b_low
    {1.0 / 256.0, 255.0 / 256.0},

    2,
    1,
    false,
    true
};

// aka bosha23
inline constexpr RKTableau<4> bosha32{
    // c
    {0.0, 1.0 / 2.0, 3.0 / 4.0, 1.0},

    // a
    {{
        {0.0},                            // row 0
        {1.0 / 2.0},                      // row 1
        {0.0, 3.0 / 4.0},                 // row 2
        {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0} // row 3
    }},

    // b_high
    {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0},

    // b_low
    {7.0 / 24.0, 1.0 / 4.0, 1.0 / 3.0, 1.0 / 8.0},

    3,
    2,
    true,
    true
};

// aka rkf45
inline constexpr RKTableau<6> rkf54_tableau{
    // c
    {0.0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0},

    // a
    {{{{0.0}},
      {{1.0 / 4.0}},
      {{3.0 / 32.0, 9.0 / 32.0}},
      {{1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0}},
      {{439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0}},
      {{-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0}}}},

    // b_high
    {16.0 / 135.0, 0.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0},

    // b_low
    {25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0, 0.0},

    5,
    4,
    false,
    true
};

inline constexpr RKTableau<6> cashkarp54_tableau{
    // c
    {0.0, 1.0 / 5.0, 3.0 / 10.0, 3.0 / 5.0, 1.0, 7.0 / 8.0},

    // a
    {{
        {0.0},                                                // row 0
        {1.0 / 5.0},                                          // row 1
        {3.0 / 40.0, 9.0 / 40.0},                             // row 2
        {3.0 / 10.0, -9.0 / 10.0, 6.0 / 5.0},                 // row 3
        {-11.0 / 54.0, 5.0 / 2.0, -70.0 / 27.0, 35.0 / 27.0}, // row 4
        {1631.0 / 55296.0,
         175.0 / 512.0,
         575.0 / 13824.0,
         44275.0 / 110592.0,
         253.0 / 4096.0} // row 5
    }},

    // b_high
    {37.0 / 378.0, 0.0, 250.0 / 621.0, 125.0 / 594.0, 0.0, 512.0 / 1771.0},

    // b_low
    {2825.0 / 27648.0,
     0.0,
     18575.0 / 48384.0,
     13525.0 / 55296.0,
     277.0 / 14336.0,
     1.0 / 4.0},

    5,
    4,
    false,
    true
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
    true,
    true
};

inline RKTableau<2> genexp2_tableau(f64 alpha) {
    assert(alpha > 0);

    RKTableau<2> tableau;

    tableau.c = {0.0, alpha};
    tableau.a = {{
        {0.0},       // row 0
        {alpha, 0.0} // row 1
    }};
    tableau.b_high = {1.0 - 1.0 / (2.0 * alpha), 1.0 / (2.0 * alpha)};
    tableau.b_low = {0.0};
    tableau.order_high = 2;
    tableau.order_low = 0;
    tableau.fsal = false;
tableau.embedded = false;
    return tableau;
}

inline RKTableau<3> genexp3_tableau(f64 alpha, f64 beta) {
    assert(alpha != 0.0);
    assert(alpha != 2.0 / 3.0);
    assert(beta != 0.0);
    assert(alpha != beta);

    RKTableau<3> tableau;

    tableau.c = {0.0, alpha, beta};

    tableau.a = {{
        {0.0},        // row 0
        {alpha, 0.0}, // row 1
        {(beta / alpha) * (beta - 3.0 * alpha * (1.0 - alpha)) / (3.0 * alpha - 2.0),
         -(beta / alpha) * (beta - alpha) / (3.0 * alpha - 2.0),
         0.0} // row 2
    }};

    tableau.b_high
        = {1.0 - (3.0 * alpha + 3.0 * beta - 2.0) / (6.0 * alpha * beta),
           (3.0 * beta - 2.0) / (6.0 * alpha * (beta - alpha)),
           (2.0 - 3.0 * alpha) / (6.0 * beta * (beta - alpha))};

    tableau.b_low = {0.0, 0.0, 0.0};

    tableau.order_high = 3;
    tableau.order_low = 0;
    tableau.fsal = false;
    tableau.embedded = false;

    return tableau;
}
