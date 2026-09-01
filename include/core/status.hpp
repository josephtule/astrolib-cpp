// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

enum struct StatusCode {
    ok,
    invalid_input,
    validation_failed,
    unsupported_type,
    missing_reference,
    duplicate_id,
    inactive_entity,
    invalid_state,
    invalid_att_state,
    invalid_mass_properties,
    invalid_shape,
    empty_measurements,
    empty_history,
    empty_events,
    size_mismatch,
    time_mismatch,
    propagation_failed,
    singular_normal_matrix,
    prediction_only,
    max_iters_reached,
    correction_rejected,
    invalid_covariance,
    singular_innovation,
    observer_not_found,
    target_not_found,
    instrument_not_found,
    sample_not_found,
    interp_failed,

    file_not_found,
    file_write_failed,
    file_close_failed,
    file_open_failed,
    file_already_exists,
    file_overwritten,

    unsupported_method,

    anchor_not_found,
    body_not_found,
    gravity_model_not_found,
    attitude_type_not_found,
    parse_failed,
    celestial_model_not_found,
    matrix_invert_failed,

    step_size_underflow,
    max_steps_reached,
    max_rejections_reached,
    non_finite_derivative,
    non_finite_result,

    empty_ephemeris,
    invalid_ephemeris_metadata,
    non_monotonic_time,
};

inline std::string status_string(StatusCode status) {
    switch (status) {
    case StatusCode::ok: return "Ok";
    case StatusCode::invalid_input: return "Invalid Input";
    case StatusCode::validation_failed: return "Validation Failed";
    case StatusCode::unsupported_type: return "Unsupported Type";
    case StatusCode::missing_reference: return "Missing Reference";
    case StatusCode::duplicate_id: return "Duplicate ID";
    case StatusCode::inactive_entity: return "Inactive Entity";
    case StatusCode::invalid_state: return "Invalid State";
    case StatusCode::invalid_att_state: return "Invalid Attitude State";
    case StatusCode::invalid_mass_properties: return "Invalid Mass Properties";
    case StatusCode::invalid_shape: return "Invalid Shape";
    case StatusCode::anchor_not_found: return "Anchor not found";
    case StatusCode::empty_measurements: return "Empty Measurements";
    case StatusCode::size_mismatch: return "Size Mismatch";
    case StatusCode::propagation_failed: return "Propagation Failed";
    case StatusCode::singular_normal_matrix: return "Singular Normal Matrix";
    case StatusCode::prediction_only: return "Prediction Only";
    case StatusCode::max_iters_reached: return "Max Iterations reached";
    case StatusCode::correction_rejected: return "Correction Rejected";
    case StatusCode::invalid_covariance: return "Invalid Covariance";
    case StatusCode::singular_innovation: return "Singular Innovation";
    case StatusCode::observer_not_found: return "Observer not found";
    case StatusCode::target_not_found: return "Target not found";
    case StatusCode::time_mismatch: return "Time mismatch";
    case StatusCode::empty_events: return "Empty Events";
    case StatusCode::instrument_not_found: return "Instrument not found";
    case StatusCode::empty_history: return "Empty history";
    case StatusCode::sample_not_found: return "Sample not found";
    case StatusCode::interp_failed: return "Interpolation failed";
    case StatusCode::file_not_found: return "File not found";
    case StatusCode::file_write_failed: return "File write failed";
    case StatusCode::file_close_failed: return "File close failed";
    case StatusCode::unsupported_method: return "Unsupported Method";
    case StatusCode::body_not_found: return "Body not found";
    case StatusCode::gravity_model_not_found: return "Gravity Model not found";
    case StatusCode::attitude_type_not_found: return "Attitude type not found";
    case StatusCode::file_open_failed: return "File open failed";
    case StatusCode::parse_failed: return "Parse failed";
    case StatusCode::celestial_model_not_found: return "Celestial Model not found";
    case StatusCode::matrix_invert_failed: return "Matrix inversion failed";
    case StatusCode::file_already_exists: return "File already exists";
    case StatusCode::file_overwritten: return "File is overwritten";
    case StatusCode::step_size_underflow: return "Step size underflow";
    case StatusCode::max_steps_reached: return "Max steps reached";
    case StatusCode::max_rejections_reached: return "Max rejections reached";
    case StatusCode::non_finite_derivative: return "Non-finite derivative";
    case StatusCode::non_finite_result: return "Non-finite result";
    case StatusCode::empty_ephemeris: return "Empty ephemeris";
    case StatusCode::invalid_ephemeris_metadata: return "Invalid ephemeris metadata";
    case StatusCode::non_monotonic_time: return "Non-monotonic time";
    }

    return "Unknown Status";
}
