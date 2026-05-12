#pragma once

#include <string>

enum struct ODStatus {
    ok,
    invalid_input,
    empty_measurements,
    size_mismatch,
    propagation_failed,
    singular_normal_matrix,
    max_iters_reached,
    correction_rejected,
    invalid_covariance,
    singular_innovation,
    observer_not_found,
    target_not_found,
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
    }
    return str;
}
