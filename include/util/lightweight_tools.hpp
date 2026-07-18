// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "util/typedefs.hpp"

inline bool toggle(bool& b) {
    b = !b;
    return b;
}

inline string bool_str(const bool b) { return b ? "true" : "false"; }
inline string active_str(const bool b) { return b ? "active" : "inactive"; }