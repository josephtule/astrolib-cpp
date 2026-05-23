#pragma once

#include <string>

enum struct StatusCode {
    ok,
    invalid_input,
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
    unsupported_method,
    body_not_found,
};

inline std::string status_string(StatusCode status) {
    std::string str;
    switch (status) {
    case StatusCode::ok: {
        str = "Ok";
    } break;
    case StatusCode::invalid_input: {
        str = "Invalid Input";
    } break;
    case StatusCode::empty_measurements: {
        str = "Empty Measurements";
    } break;
    case StatusCode::size_mismatch: {
        str = "Size Mismatch";
    } break;
    case StatusCode::propagation_failed: {
        str = "Propagation Failed";
    } break;
    case StatusCode::singular_normal_matrix: {
        str = "Singular Normal Matrix";
    } break;
    case StatusCode::prediction_only: {
        str = "Prediction Only";
    } break;
    case StatusCode::max_iters_reached: {
        str = "Max Iterations reached";
    } break;
    case StatusCode::correction_rejected: {
        str = "Correction Rejected";
    } break;
    case StatusCode::invalid_covariance: {
        str = "Invalid Covariance";
    } break;
    case StatusCode::singular_innovation: {
        str = "Singular Innovation";
    } break;
    case StatusCode::observer_not_found: {
        str = "Observer not found";
    } break;
    case StatusCode::target_not_found: {
        str = "Target not found";
    } break;
    case StatusCode::time_mismatch: {
        str = "Time mismatch";
    } break;
    case StatusCode::empty_events: {
        str = "Empty Events";
    } break;
    case StatusCode::instrument_not_found: {
        str = "Instrument not found";
    } break;
    case StatusCode::empty_history: {
        str = "Empty history";
    } break;
    case StatusCode::sample_not_found: {
        str = "Sample not found";
    } break;
    case StatusCode::interp_failed: {
        str = "Interpolation failed";
    } break;
    case StatusCode::file_not_found: {
        str = "File not found";
    } break;
    case StatusCode::file_write_failed: {
        str = "File write failed";
    } break;
    case StatusCode::file_close_failed: {
        str = "File close failed";
    } break;
    case StatusCode::unsupported_method: {
        str = "Unsupported Method";
    } break;
    case StatusCode::body_not_found: {
        str = "Body not found";
    } break;
    }
    return str;
}
using StatusCode = StatusCode; // TODO: do full rename

inline bool od_status_success(StatusCode status) {
    return status == StatusCode::ok || status == StatusCode::prediction_only;
}

template <class A, class B>
inline B copy_od_ekf_result(const A& result_in) {
    B result_out;
    result_out.status = result_in.status;
    result_out.filter = result_in.filter;
    result_out.raw_residual_norm = result_in.raw_residual_norm;
    result_out.residual_norm = result_in.residual_norm;

    return result_out;
}