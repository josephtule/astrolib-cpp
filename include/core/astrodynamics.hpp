#pragma once

#include "util/constants.hpp"
#include "util/math.hpp"
#include "util/units.hpp"

#include <cmath>

inline bool mean_anom_to_eccen_anom(
    f64 mean_anom,
    f64 ecc,
    f64& eccen_anom,
    UAngle angle_in = UAngle::radian,
    UAngle angle_out = UAngle::radian
) {
    if (angle_in != UAngle::radian) {
        mean_anom = convert_angle(mean_anom, angle_in, UAngle::radian);
    }
    mean_anom = wrap_angle(mean_anom, 0.0, twopi, UAngle::radian, UAngle::radian);

    auto func = [&](f64 eccen_anom) -> f64 {
        return eccen_anom - ecc * std::sin(eccen_anom) - mean_anom;
    };
    auto dfunc = [&](f64 eccen_anom) -> f64 { return 1.0 - ecc * std::cos(eccen_anom); };

    // use simple newton's method to solve
    constexpr i32 max_iter = 100;
    i32 iter = 0;
    f64 x_iter = mean_anom + ecc * std::sin(mean_anom);
    bool converged = false;
    while (iter < max_iter) {
        f64 fx = func(x_iter);
        if (std::abs(fx) <= std::sqrt(eps(fx))) {
            converged = true;
            break;
        }

        f64 dfx = dfunc(x_iter);
        if (!std::isfinite(dfx) || std::abs(dfx) <= std::sqrt(eps(dfx))) {
            converged = false;
            break;
        }

        f64 x_next = x_iter - fx / dfx;

        if (std::abs(x_next - x_iter) <= std::sqrt(eps(x_next))) {
            f64 f_next = func(x_next);
            converged = std::abs(f_next) <= std::sqrt(eps(f_next));
            x_iter = x_next;
            break;
        }

        x_iter = x_next;
        ++iter;
    }

    if (iter >= max_iter) converged = false;

    eccen_anom = x_iter;
    if (angle_out != UAngle::radian) {
        eccen_anom = convert_angle(eccen_anom, UAngle::radian, angle_out);
    }

    return converged;
}
