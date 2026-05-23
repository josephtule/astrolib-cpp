#pragma once

#include "core/observation_type.hpp"
#include <string>

enum struct ODStatus {
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

inline std::string od_status_string(ODStatus status) {
    std::string str;
    switch (status) {
    case ODStatus::ok: {
        str = "Ok";
    } break;
    case ODStatus::invalid_input: {
        str = "Invalid Input";
    } break;
    case ODStatus::empty_measurements: {
        str = "Empty Measurements";
    } break;
    case ODStatus::size_mismatch: {
        str = "Size Mismatch";
    } break;
    case ODStatus::propagation_failed: {
        str = "Propagation Failed";
    } break;
    case ODStatus::singular_normal_matrix: {
        str = "Singular Normal Matrix";
    } break;
    case ODStatus::prediction_only: {
        str = "Prediction Only";
    } break;
    case ODStatus::max_iters_reached: {
        str = "Max Iterations reached";
    } break;
    case ODStatus::correction_rejected: {
        str = "Correction Rejected";
    } break;
    case ODStatus::invalid_covariance: {
        str = "Invalid Covariance";
    } break;
    case ODStatus::singular_innovation: {
        str = "Singular Innovation";
    } break;
    case ODStatus::observer_not_found: {
        str = "Observer not found";
    } break;
    case ODStatus::target_not_found: {
        str = "Target not found";
    } break;
    case ODStatus::time_mismatch: {
        str = "Time mismatch";
    } break;
    case ODStatus::empty_events: {
        str = "Empty Events";
    } break;
    case ODStatus::instrument_not_found: {
        str = "Instrument not found";
    } break;
    case ODStatus::empty_history: {
        str = "Empty history";
    } break;
    case ODStatus::sample_not_found: {
        str = "Sample not found";
    } break;
    case ODStatus::interp_failed: {
        str = "Interpolation failed";
    } break;
    case ODStatus::file_not_found: {
        str = "File not found";
    } break;
    case ODStatus::file_write_failed: {
        str = "File write failed";
    } break;
    case ODStatus::file_close_failed: {
        str = "File close failed";
    } break;
    case ODStatus::unsupported_method: {
        str = "Unsupported Method";
    } break;
    case ODStatus::body_not_found: {
        str = "Body not found";
    } break;
    }
    return str;
}
using StatusCode = ODStatus; // TODO: do full rename

inline bool od_status_success(ODStatus status) {
    return status == ODStatus::ok || status == ODStatus::prediction_only;
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