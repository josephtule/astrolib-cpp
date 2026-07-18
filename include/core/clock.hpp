// Copyright 2025-2026 Joseph Tu Le
// SPDX-License-Identifier: Apache-2.0

#pragma once 

#include "util/typedefs.hpp"

struct SimulationClock {
    f64 t_real; // realtime clock
    f64 t_sim; // simulation time
    f64 dt_tr; // translational dt
    f64 dt_att; // attitude dt
    f64 dt_cont; // control update dt
};