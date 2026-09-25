// Copyright 2025-2026 Joseph Tu Le
// SPDX-License-Identifier: Apache-2.0

#include "core/propagation_comparison.hpp"
#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/integrator_common.hpp"
#include "core/status.hpp"
#include "core/world_stepper.hpp"
#include "util/math.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <utility>
#include <variant>

static WorldStepResult dispatch_selected_backend(
    World& world,
    f64 dt_step,
    const WorldStepperConfig& stepper_cfg,
    WorldStepperWorkspace& workspace,
    PropagationBackend backend
) {
    switch (backend) {
    case PropagationBackend::world_stepper_legacy:
        return step_world_legacy(world, dt_step, stepper_cfg, workspace);
    case PropagationBackend::world_stepper:
        return step_world(world, dt_step, stepper_cfg, workspace);
    case PropagationBackend::compiled_problem:
        return WorldStepResult{
            .status = StatusCode::unsupported_method,
            .t = world.t_sim()
        };
    }

    return WorldStepResult{.status = StatusCode::unsupported_method, .t = world.t_sim()};
}

static StatusCode make_propagation_sample(
    const World& world,
    const umap<string, EntityId>& body_ids,
    const PropagationRunConfig& run_cfg,
    const Celestial* central_body,
    f64 requested_t,
    const string& body_config_id,
    PropagationStateSample& sample
) {
    const auto id_it = body_ids.find(body_config_id);
    if (id_it == body_ids.end()) return StatusCode::body_not_found;

    const Body* body = world.body(id_it->second);
    if (body == nullptr) return StatusCode::body_not_found;
    if (!finite_state_tr(body->x_tr)) return StatusCode::invalid_state;
    if (run_cfg.include_attitude_output && !finite_state_att(body->x_att)) {
        return StatusCode::invalid_att_state;
    }

    PropagationStateSample temp;
    temp.requested_t = requested_t;
    temp.t = world.t_sim();
    temp.body_config_id = body_config_id;
    temp.x_tr = body->x_tr;
    temp.has_attitude = run_cfg.include_attitude_output;

    if (temp.has_attitude) {
        temp.x_att = body->x_att;
        temp.quaternion_norm = body->x_att.q.norm();
    }

    if (central_body != nullptr && body->id != central_body->id) {
        const vec3d r_rel = body->x_tr.r - central_body->x_tr.r;
        const vec3d v_rel = body->x_tr.v - central_body->x_tr.v;
        const f64 r_mag = r_rel.norm();
        if (!finite_pos(r_mag) || !finite_pos(central_body->mu)) {
            return StatusCode::invalid_state;
        }

        temp.invariants.has_orbital_invariants = true;
        temp.invariants.specific_energy =
            0.5 * v_rel.squaredNorm() - central_body->mu / r_mag;
        temp.invariants.specific_angular_momentum = r_rel.cross(v_rel);
        temp.invariants.specific_angular_momentum_norm =
            temp.invariants.specific_angular_momentum.norm();

        if (!std::isfinite(temp.invariants.specific_energy)
            || !finite_vec(temp.invariants.specific_angular_momentum)) {
            return StatusCode::non_finite_result;
        }
    }

    sample = std::move(temp);
    return StatusCode::ok;
}

StatusCode run_propagation_case(
    World& world,
    const umap<string, EntityId>& body_ids,
    const WorldStepperConfig& stepper_cfg,
    const PropagationRunConfig& run_cfg,
    PropagationRunResult& out
) {
    out = PropagationRunResult{};
    out.label = run_cfg.label;
    out.frame_label = run_cfg.frame_label;
    out.backend = run_cfg.backend;
    out.integrator = run_cfg.integrator;
    out.requested_t0 = run_cfg.t0;
    out.requested_tf = run_cfg.tf;
    out.actual_t0 = world.t_sim();
    out.actual_tf = world.t_sim();

    auto prop_ok = [&](StatusCode status) {
        out.status = status;
        out.actual_tf = world.t_sim();
        if (status != StatusCode::ok) out.failure_t = world.t_sim();
        return status;
    };

    if (run_cfg.backend != PropagationBackend::world_stepper_legacy
        && run_cfg.backend != PropagationBackend::world_stepper) {
        return prop_ok(StatusCode::unsupported_method);
    }
    if (!std::isfinite(run_cfg.t0)
        || !std::isfinite(run_cfg.tf)
        || run_cfg.sample_epochs.empty()
        || run_cfg.body_config_ids.empty()
        || run_cfg.frame_label.empty()
        || !finite_pos(run_cfg.tolerance.time)) {
        return prop_ok(StatusCode::invalid_input);
    }
    if (integrator_family(run_cfg.integrator) == IntegratorFamily::fixed
        && !finite_pos(run_cfg.nominal_dt)) {
        return prop_ok(StatusCode::invalid_input);
    }
    if (!run_cfg.step_tr && !run_cfg.step_att) {
        return prop_ok(StatusCode::invalid_input);
    }
    if (std::abs(world.t_sim() - run_cfg.t0) > run_cfg.tolerance.time) {
        return prop_ok(StatusCode::time_mismatch);
    }

    const f64 interval = run_cfg.tf - run_cfg.t0;
    const f64 dir = interval > 0.0 ? 1.0 : interval < 0.0 ? -1.0 : 0.0;
    f64 previous_t = run_cfg.t0;
    for (size_t i = 0; i < run_cfg.sample_epochs.size(); ++i) {
        const f64 epoch = run_cfg.sample_epochs[i];
        if (!std::isfinite(epoch)) return prop_ok(StatusCode::invalid_input);
        if (i > 0 && dir * (epoch - previous_t) <= run_cfg.tolerance.time) {
            return prop_ok(StatusCode::non_monotonic_time);
        }
        if (dir == 0.0 && std::abs(epoch - run_cfg.t0) > run_cfg.tolerance.time) {
            return prop_ok(StatusCode::time_mismatch);
        }
        if (dir != 0.0
            && (dir * (epoch - run_cfg.t0) < -run_cfg.tolerance.time
                || dir * (run_cfg.tf - epoch) < -run_cfg.tolerance.time)) {
            return prop_ok(StatusCode::time_mismatch);
        }
        previous_t = epoch;
    }
    if (std::abs(run_cfg.sample_epochs.back() - run_cfg.tf) > run_cfg.tolerance.time) {
        return prop_ok(StatusCode::time_mismatch);
    }

    uset<string> unique_body_ids;
    for (const string& config_id : run_cfg.body_config_ids) {
        if (config_id.empty() || !unique_body_ids.insert(config_id).second) {
            return prop_ok(StatusCode::invalid_input);
        }

        const auto id_it = body_ids.find(config_id);
        if (id_it == body_ids.end() || world.body(id_it->second) == nullptr) {
            return prop_ok(StatusCode::body_not_found);
        }
    }

    const Celestial* central_body = nullptr;
    if (run_cfg.invariants.include_orbital_invariants) {
        const auto central_it = body_ids.find(run_cfg.invariants.central_body_config_id);
        if (central_it == body_ids.end()) return prop_ok(StatusCode::body_not_found);
        central_body = world.celestial(central_it->second);

        if (central_body == nullptr) return prop_ok(StatusCode::invalid_input);
        if (!finite_state_tr(central_body->x_tr) || !finite_pos(central_body->mu)) {
            return prop_ok(StatusCode::invalid_state);
        }
    }

    WorldStepperConfig execution_cfg = stepper_cfg;
    execution_cfg.integrator_tr = run_cfg.integrator;
    execution_cfg.integrator_att = run_cfg.integrator;
    execution_cfg.substeps = 1;
    execution_cfg.ticks = 1;
    execution_cfg.dt_scale = 1.0;
    execution_cfg.step_tr = run_cfg.step_tr;
    execution_cfg.step_att = run_cfg.step_att;
    execution_cfg.paused = false;
    execution_cfg.adaptive.opts = run_cfg.adaptive;
    execution_cfg.adaptive.use_substeps = false;

    StatusCode status = validate_world_stepper_config(execution_cfg);
    if (status != StatusCode::ok) return prop_ok(status);

    WorldStepperWorkspace workspace;
    rebuild_world_stepper_workspace(world, workspace);
    WorldStepperInstrumentation instrumentation;
    workspace.instrumentation = &instrumentation;
    out.metrics.workspace_rebuilds = 1;

    const auto start = std::chrono::steady_clock::now();
    auto finish_timed = [&](StatusCode result_status) {
        const auto stop = std::chrono::steady_clock::now();
        out.metrics.runtime_ms =
            std::chrono::duration<f64, std::milli>(stop - start).count();
        out.metrics.derivative_evaluations = instrumentation.derivative_evaluations;
        out.metrics.provider_translation_queries =
            instrumentation.provider_translation_queries;
        out.metrics.provider_orientation_queries =
            instrumentation.provider_orientation_queries;
        out.metrics.stage_builds = instrumentation.stage_builds;
        return prop_ok(result_status);
    };

    out.samples.reserve(run_cfg.sample_epochs.size() * run_cfg.body_config_ids.size());

    for (f64 requested_t : run_cfg.sample_epochs) {
        switch (integrator_family(run_cfg.integrator)) {
        case IntegratorFamily::fixed: {
            while (dir * (requested_t - world.t_sim()) > run_cfg.tolerance.time) {
                const f64 remaining = requested_t - world.t_sim();
                const f64 landing_slack = std::max(
                    run_cfg.tolerance.time,
                    std::max(
                        run_cfg.nominal_dt * 1e-9,
                        1024.0
                            * std::numeric_limits<f64>::epsilon()
                            * std::max(1.0, std::abs(requested_t))
                    )
                );
                const f64 dt_step =
                    std::abs(remaining) <= run_cfg.nominal_dt + landing_slack
                        ? remaining
                        : dir * run_cfg.nominal_dt;
                WorldStepResult step = dispatch_selected_backend(
                    world,
                    dt_step,
                    execution_cfg,
                    workspace,
                    run_cfg.backend
                );
                out.stats += step.stats;
                if (step.status != StatusCode::ok) {
                    return finish_timed(step.status);
                }
            }
        } break;
        case IntegratorFamily::adaptive: {
            const f64 remaining = requested_t - world.t_sim();
            if (dir * remaining > run_cfg.tolerance.time) {
                WorldStepResult step = dispatch_selected_backend(
                    world,
                    remaining,
                    execution_cfg,
                    workspace,
                    run_cfg.backend
                );
                out.stats += step.stats;
                if (step.status != StatusCode::ok) {
                    return finish_timed(step.status);
                }
            }
        } break;
        }

        if (std::abs(world.t_sim() - requested_t) > run_cfg.tolerance.time) {
            return finish_timed(StatusCode::time_mismatch);
        }

        for (const string& config_id : run_cfg.body_config_ids) {
            PropagationStateSample sample;
            status = make_propagation_sample(
                world,
                body_ids,
                run_cfg,
                central_body,
                requested_t,
                config_id,
                sample
            );
            if (status != StatusCode::ok) return finish_timed(status);
            out.samples.push_back(std::move(sample));
        }
    }

    return finish_timed(StatusCode::ok);
}

StatusCode compare_propagation_runs(
    const PropagationRunResult& reference,
    const PropagationRunResult& candidate,
    const PropagationComparisonTolerance& tolerance,
    PropagationComparisonResult& out
) {
    out = PropagationComparisonResult{};
    out.reference_label = reference.label;
    out.candidate_label = candidate.label;

    if (!finite_nonneg(tolerance.time)
        || !finite_nonneg(tolerance.position)
        || !finite_nonneg(tolerance.velocity)
        || !finite_nonneg(tolerance.quaternion)
        || !finite_nonneg(tolerance.angular_velocity)
        || !finite_nonneg(tolerance.invariant)) {
        out.status = StatusCode::invalid_input;
        return out.status;
    }
    if (reference.status != StatusCode::ok) {
        out.status = reference.status;
        return out.status;
    }
    if (candidate.status != StatusCode::ok) {
        out.status = candidate.status;
        return out.status;
    }
    if (reference.frame_label != candidate.frame_label) {
        out.status = StatusCode::invalid_input;
        return out.status;
    }
    if (reference.samples.empty()) {
        out.status = StatusCode::invalid_state;
        return out.status;
    }
    if (reference.samples.size() != candidate.samples.size()) {
        out.status = StatusCode::size_mismatch;
        return out.status;
    }

    out.passed = true;
    out.samples.reserve(reference.samples.size());
    bool first_failure = true;

    for (size_t i = 0; i < reference.samples.size(); ++i) {
        const PropagationStateSample& ref = reference.samples[i];
        const PropagationStateSample& test = candidate.samples[i];

        if (ref.body_config_id != test.body_config_id) {
            out.status = StatusCode::invalid_input;
            return out.status;
        }
        if (std::abs(ref.requested_t - test.requested_t) > tolerance.time) {
            out.status = StatusCode::time_mismatch;
            return out.status;
        }
        if (ref.has_attitude != test.has_attitude
            || ref.invariants.has_orbital_invariants
                   != test.invariants.has_orbital_invariants) {
            out.status = StatusCode::invalid_state;
            return out.status;
        }

        PropagationSampleDifference difference;
        difference.requested_t = ref.requested_t;
        difference.body_config_id = ref.body_config_id;
        difference.time = std::abs(ref.t - test.t);
        difference.position = (ref.x_tr.r - test.x_tr.r).norm();
        difference.velocity = (ref.x_tr.v - test.x_tr.v).norm();

        if (ref.has_attitude) {
            difference.quaternion = std::min(
                (ref.x_att.q - test.x_att.q).norm(),
                (ref.x_att.q + test.x_att.q).norm()
            );
            difference.quaternion_norm =
                std::abs(ref.quaternion_norm - test.quaternion_norm);
            difference.angular_velocity = (ref.x_att.w - test.x_att.w).norm();
        }

        if (ref.invariants.has_orbital_invariants) {
            difference.specific_energy = std::abs(
                ref.invariants.specific_energy - test.invariants.specific_energy
            );
            difference.specific_angular_momentum =
                (ref.invariants.specific_angular_momentum
                 - test.invariants.specific_angular_momentum)
                    .norm();
        }

        if (!std::isfinite(difference.time)
            || !std::isfinite(difference.position)
            || !std::isfinite(difference.velocity)
            || !std::isfinite(difference.quaternion)
            || !std::isfinite(difference.quaternion_norm)
            || !std::isfinite(difference.angular_velocity)
            || !std::isfinite(difference.specific_energy)
            || !std::isfinite(difference.specific_angular_momentum)) {
            out.status = StatusCode::non_finite_result;
            return out.status;
        }

        difference.passed =
            difference.time <= tolerance.time
            && difference.position <= tolerance.position
            && difference.velocity <= tolerance.velocity
            && difference.quaternion <= tolerance.quaternion
            && difference.quaternion_norm <= tolerance.quaternion
            && difference.angular_velocity <= tolerance.angular_velocity
            && difference.specific_energy <= tolerance.invariant
            && difference.specific_angular_momentum <= tolerance.invariant;

        out.max_time = std::max(out.max_time, difference.time);
        out.max_position = std::max(out.max_position, difference.position);
        out.max_velocity = std::max(out.max_velocity, difference.velocity);
        out.max_quaternion = std::max(out.max_quaternion, difference.quaternion);
        out.max_quaternion_norm =
            std::max(out.max_quaternion_norm, difference.quaternion_norm);
        out.max_angular_velocity =
            std::max(out.max_angular_velocity, difference.angular_velocity);
        out.max_specific_energy =
            std::max(out.max_specific_energy, difference.specific_energy);
        out.max_specific_angular_momentum = std::max(
            out.max_specific_angular_momentum,
            difference.specific_angular_momentum
        );

        if (!difference.passed && first_failure) {
            first_failure = false;
            out.first_failed_t = difference.requested_t;
            out.first_failed_body_config_id = difference.body_config_id;
        }
        out.passed = out.passed && difference.passed;
        out.samples.push_back(std::move(difference));
    }

    out.status = StatusCode::ok;
    return out.status;
}

StatusCode save_propagation_run_csv(
    const string& filepath,
    const PropagationRunResult& result
) {
    if (filepath.empty() || result.status != StatusCode::ok) {
        return StatusCode::invalid_input;
    }

    const std::filesystem::path path{filepath};
    std::error_code error;
    if (std::filesystem::exists(path, error)) {
        return error ? StatusCode::file_open_failed : StatusCode::file_already_exists;
    }
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path(), error);
        if (error) return StatusCode::file_write_failed;
    }

    std::ofstream file(path);
    if (!file.is_open()) return StatusCode::file_open_failed;

    file << std::setprecision(17);
    file << "label,backend,integrator,frame,requested_t,actual_t,body_config_id,"
            "rx,ry,rz,vx,vy,vz,has_attitude,qx,qy,qz,qs,wx,wy,wz,quaternion_norm,"
            "has_invariants,specific_energy,hx,hy,hz,specific_angular_momentum_norm\n";

    for (const PropagationStateSample& sample : result.samples) {
        file
            << result.label
            << ','
            << propagation_backend_string(result.backend)
            << ','
            << integrator_str(result.integrator)
            << ','
            << result.frame_label
            << ','
            << sample.requested_t
            << ','
            << sample.t
            << ','
            << sample.body_config_id
            << ','
            << sample.x_tr.r(0)
            << ','
            << sample.x_tr.r(1)
            << ','
            << sample.x_tr.r(2)
            << ','
            << sample.x_tr.v(0)
            << ','
            << sample.x_tr.v(1)
            << ','
            << sample.x_tr.v(2)
            << ','
            << sample.has_attitude
            << ','
            << sample.x_att.q(0)
            << ','
            << sample.x_att.q(1)
            << ','
            << sample.x_att.q(2)
            << ','
            << sample.x_att.q(3)
            << ','
            << sample.x_att.w(0)
            << ','
            << sample.x_att.w(1)
            << ','
            << sample.x_att.w(2)
            << ','
            << sample.quaternion_norm
            << ','
            << sample.invariants.has_orbital_invariants
            << ','
            << sample.invariants.specific_energy
            << ','
            << sample.invariants.specific_angular_momentum(0)
            << ','
            << sample.invariants.specific_angular_momentum(1)
            << ','
            << sample.invariants.specific_angular_momentum(2)
            << ','
            << sample.invariants.specific_angular_momentum_norm
            << '\n';
    }

    if (!file.good()) return StatusCode::file_write_failed;
    file.close();
    if (file.fail()) return StatusCode::file_close_failed;
    return StatusCode::ok;
}
