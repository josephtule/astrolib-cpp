// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/integrator_adaptive.hpp"
#include "core/integrator_fixed.hpp"

struct IntegratorType {
    IntegratorFamily family = IntegratorFamily::fixed;
    IntegratorTypeFixed fixed = IntegratorTypeFixed::rk4;
    IntegratorTypeAdaptive adaptive = IntegratorTypeAdaptive::dopri54;
};

inline string integrator_name(const IntegratorType& type) {
    switch (type.family) {
    case IntegratorFamily::fixed: return integrator_name(type.fixed);
    case IntegratorFamily::adaptive: return integrator_name(type.adaptive);
    }
    return "Unknown";
}

inline string integrator_str(const IntegratorType& type) {
    switch (type.family) {
    case IntegratorFamily::fixed: return integrator_str(type.fixed);
    case IntegratorFamily::adaptive: return integrator_str(type.adaptive);
    }
    return "unknown";
}
