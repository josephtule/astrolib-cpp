#pragma once

#include "core/integrator_common.hpp"
#include "core/state.hpp"
#include "core/status.hpp"
#include "util/math.hpp"
#include "util/vecdefs.hpp"

struct AdaptiveIntegratorConfig {
    f64 rel_tol = 1e-9;
    f64 abs_tol_r = 1e-9;
    f64 abs_tol_v = 1e-12;

    f64 dt_initial = 10.0;
    f64 dt_min = 1e-9;
    f64 dt_max = 600.0;

    f64 safety = 0.9;
    f64 scale_min = 0.2;
    f64 scale_max = 5.0;

    i64 max_attempts = 100000;
    i64 max_rejections = 1000;
};

inline StatusCode validate_adaptive_integrator_config(
    const AdaptiveIntegratorConfig& cfg
) {
    if (!finite_pos(cfg.rel_tol)) return StatusCode::invalid_input;
    if (!finite_pos(cfg.abs_tol_r)) return StatusCode::invalid_input;
    if (!finite_pos(cfg.abs_tol_v)) return StatusCode::invalid_input;

    if (!finite_pos(cfg.dt_initial)) return StatusCode::invalid_input;
    if (!finite_pos(cfg.dt_min)) return StatusCode::invalid_input;
    if (!finite_pos(cfg.dt_max)) return StatusCode::invalid_input;
    if (cfg.dt_min > cfg.dt_max) return StatusCode::invalid_input;

    if (!finite_pos(cfg.safety)) return StatusCode::invalid_input;
    if (!finite_pos(cfg.scale_min)) return StatusCode::invalid_input;
    if (!finite_pos(cfg.scale_max)) return StatusCode::invalid_input;
    if (cfg.scale_min > cfg.scale_max) return StatusCode::invalid_input;
    if (cfg.scale_min > 1 || cfg.scale_max < 1) return StatusCode::invalid_input;

    if (cfg.max_attempts <= 0 || cfg.max_rejections <= 0)
        return StatusCode::invalid_input;

    return StatusCode::ok;
}

struct AdaptiveIntegratorStats {
    i64 attempted_steps = 0;
    i64 accepted_steps = 0;
    i64 rejected_steps = 0;
    i64 derivative_evaluations = 0;

    f64 min_accepted_dt = 0.0;
    f64 max_accepted_dt = 0.0;
    f64 final_accepted_dt = 0.0;
};

template <typename State, typename Deriv>
struct AdaptiveTrialResult {
    StatusCode status = StatusCode::invalid_input;

    State x_high{}; // higher-order estimate

    Deriv error_delta{};
    i32 derivative_evaluations = 0;
};

template <typename State>
struct AdaptivePropagationResult {
    StatusCode status = StatusCode::invalid_input;

    f64 t = 0.0;
    State x{};

    f64 final_error_norm = 0.0;
    AdaptiveIntegratorStats stats{};
};

template <typename State, typename Deriv, typename Func>
inline AdaptiveTrialResult<State, Deriv> step_dopri54_trial(
    Func&& f,
    f64 t,
    const State& x,
    f64 dt
) {
    AdaptiveTrialResult<State, Deriv> trial{};

    constexpr size_t stages = 7;

    // TODO: optimize this to be DP54 specific
    array<Deriv, stages> k_trial
        = rk_generic_stages<State, Deriv>(f, t, x, dt, dopri54_tableau);

    State x_high = x;
    for (std::size_t i = 0; i < stages; ++i) {
        x_high += dt * dopri54_tableau.b_high[i] * k_trial[i];
    }

    Deriv error_delta{};

    for (std::size_t i = 0; i < stages; ++i) {
        f64 error_weight = dopri54_tableau.b_high[i] - dopri54_tableau.b_low[i];

        error_delta = error_delta + dt * error_weight * k_trial[i];
    }

    trial.status = StatusCode::ok;
    trial.x_high = x_high;
    trial.error_delta = error_delta;
    trial.derivative_evaluations = stages;

    return trial;
}

inline f64 adaptive_error_norm_state_tr(
    const StateTr& x,
    const StateTr& x_high,
    const DerivTr& error_delta,
    const AdaptiveIntegratorConfig& cfg
) {
    if (!finite_state_tr(x) || !finite_state_tr(x_high) || !finite_deriv_tr(error_delta))
        return inf<f64>;

    vec6d x_vec = statetr_to_vec6d(x);
    vec6d x_high_vec = statetr_to_vec6d(x_high);
    vec6d error_vec = derivtr_to_vec6d(error_delta);

    vec6d abs_tol;
    abs_tol << cfg.abs_tol_r, cfg.abs_tol_r, cfg.abs_tol_r, cfg.abs_tol_v, cfg.abs_tol_v,
        cfg.abs_tol_v;
    if (!finite_vec(abs_tol)) return inf<f64>;

    vec6d mag = x_vec.cwiseAbs().cwiseMax(x_high_vec.cwiseAbs());
    if (!finite_vec(mag)) return inf<f64>;

    vec6d scale = abs_tol + cfg.rel_tol * mag;
    if (!finite_vec(scale)) return inf<f64>;

    vec6d weighted_error = error_vec.cwiseQuotient(scale);
    return std::sqrt(weighted_error.squaredNorm() / 6.0);
}

inline f64 adaptive_step_scale(f64 error_norm, const AdaptiveIntegratorConfig& cfg) {
    if (!finite_nonneg(error_norm)) return inf<f64>;
    if (error_norm == 0.0) return cfg.scale_max;

    f64 exponent = -1.0 / static_cast<f64>(dopri54_tableau.order_low + 1);

    f64 scale = cfg.safety * std::pow(error_norm, exponent);
    scale = std::clamp(scale, cfg.scale_min, cfg.scale_max);
    return scale;
}

template <typename State, typename Deriv, typename Func, typename ErrorNormFunc>
AdaptivePropagationResult<State> propagate_dopri54(
    Func&& f,
    ErrorNormFunc&& error_norm_fn,
    f64 t0,
    const State& x0,
    f64 tf,
    const AdaptiveIntegratorConfig& cfg
) {
    AdaptivePropagationResult<State> result{};
    result.status = validate_adaptive_integrator_config(cfg);
    if (result.status != StatusCode::ok) return result;

    if (!isfinite(t0) || !isfinite(tf)) result.status = StatusCode::invalid_input;
    if (result.status != StatusCode::ok) return result;

    if (tf == t0) {
        result.status = StatusCode::ok;
        result.t = tf;
        result.x = x0;
        return result;
    }

    result.t = t0;
    result.x = x0;

    f64 dir = tf > t0 ? 1.0 : -1.0;
    f64 dt = dir * std::clamp(cfg.dt_initial, cfg.dt_min, cfg.dt_max);

    while (true) {
        if (result.stats.attempted_steps >= cfg.max_attempts) {
            result.status = StatusCode::max_steps_reached;
            return result;
        }

        f64 t_rem = tf - result.t;
        bool final_attempt = std::abs(dt) >= std::abs(t_rem);

        if (final_attempt) {
            dt = t_rem;
        }

        AdaptiveTrialResult<State, Deriv> trial
            = step_dopri54_trial<State, Deriv>(f, result.t, result.x, dt);
        ++result.stats.attempted_steps;
        result.stats.derivative_evaluations += trial.derivative_evaluations;
        if (trial.status != StatusCode::ok) {
            result.status = trial.status;
            return result;
        }

        f64 error_norm = error_norm_fn(result.x, trial.x_high, trial.error_delta);
        if (!isfinite(error_norm)) {
            result.status = StatusCode::non_finite_result;
            return result;
        }
        bool accepted = error_norm <= 1.0;

        f64 scale = adaptive_step_scale(error_norm, cfg);
        if (!accepted) {
            ++result.stats.rejected_steps;
            if (result.stats.rejected_steps >= cfg.max_rejections) {
                result.status = StatusCode::max_rejections_reached;
                return result;
            }

            if (scale > 1.0) scale = 1.0;
            dt *= scale;
            if (std::abs(dt) < cfg.dt_min) {
                result.status = StatusCode::step_size_underflow;
                return result;
            }
            continue;
        }

        result.x = trial.x_high;
        result.t += dt;
        result.final_error_norm = error_norm;

        f64 dt_abs = std::abs(dt);
        if (result.stats.accepted_steps == 0) {
            result.stats.min_accepted_dt = dt_abs;
            result.stats.max_accepted_dt = dt_abs;
        } else {
            if (dt_abs < result.stats.min_accepted_dt)
                result.stats.min_accepted_dt = dt_abs;
            if (dt_abs > result.stats.max_accepted_dt)
                result.stats.max_accepted_dt = dt_abs;
        }
        ++result.stats.accepted_steps;
        result.stats.final_accepted_dt = dt_abs;

        if (final_attempt) {
            result.t = tf;
            result.status = StatusCode::ok;
            return result;
        }

        f64 dt_magnitude = std::clamp(abs(dt) * scale, cfg.dt_min, cfg.dt_max);
        dt = dir * dt_magnitude;
    }
}
