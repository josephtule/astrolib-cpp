// Copyright 2025-2026 Joseph Tu Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "util/constants.hpp"
#include "util/math.hpp"
#include "util/units.hpp"

#include <algorithm>
#include <cmath>

inline bool mean_anom_to_eccen_anom(
    f64 mean_anom,
    f64 ecc,
    f64& eccen_anom,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian,
    f64 tol = tol12
) {
    if (!std::isfinite(mean_anom) || !std::isfinite(ecc) || ecc < 0.0 || ecc >= 1.0)
        return false;

    if (angle_in != UAngle::radian) {
        mean_anom = convert_angle(mean_anom, angle_in, UAngle::radian);
    }
    mean_anom = wrap_angle(mean_anom, 0.0, twopi, UAngle::radian, UAngle::radian);

    // f(E) = E - e sin(E) - M
    auto func = [&](f64 eccen_anom) -> f64 {
        return eccen_anom - ecc * std::sin(eccen_anom) - mean_anom;
    };
    // df/dE(E) = 1 - e cos(E)
    auto dfunc = [&](f64 eccen_anom) -> f64 { return 1.0 - ecc * std::cos(eccen_anom); };

    constexpr i32 max_iter = 100;
    i32 iter = 0;

    f64 x_iter = mean_anom;
    if (ecc > 0.8) {
        x_iter = pi;
    } else {
        x_iter = mean_anom + ecc * std::sin(mean_anom);
    }

    bool converged = false;
    while (iter < max_iter) {
        // scale residual so large/small angles for tolerance check
        // scale is abs largest term of func and dfunc
        f64 fx = func(x_iter);
        f64 f_scale = std::max(
            {1.0, std::abs(x_iter), std::abs(ecc * std::sin(x_iter)), std::abs(mean_anom)}
        );
        if (!std::isfinite(fx)) {
            converged = false;
            break;
        }
        if (std::abs(fx) <= tol * f_scale) {
            converged = true;
            break;
        }

        f64 dfx = dfunc(x_iter);
        f64 df_scale = std::max({1.0, std::abs(ecc * std::cos(x_iter))});
        if (!std::isfinite(dfx) || std::abs(dfx) <= tol * df_scale) {
            converged = false;
            break;
        }

        f64 x_next = x_iter - fx / dfx;
        if (!std::isfinite(x_next)) {
            converged = false;
            break;
        }

        if (std::abs(x_next - x_iter) <= tol * std::max(1.0, std::abs(x_next))) {
            f64 f_next = func(x_next);
            f64 f_next_scale = std::max(
                {1.0,
                 std::abs(x_next),
                 std::abs(ecc * std::sin(x_next)),
                 std::abs(mean_anom)}
            );
            converged = std::isfinite(f_next) && std::abs(f_next) <= tol * f_next_scale;
            x_iter = x_next;
            break;
        }

        x_iter = x_next;
        ++iter;
    }

    if (iter >= max_iter) {
        converged = false;
    }

    eccen_anom = x_iter;
    if (angle_out != UAngle::radian) {
        eccen_anom = convert_angle(eccen_anom, UAngle::radian, angle_out);
    }

    return converged;
}
