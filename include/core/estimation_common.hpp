// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/status.hpp"

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
