// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "util/typedefs.hpp"

const f64 pi = 3.1415926535897932384626433832795;
const f64 pio2 = pi / 2.;
const f64 pio3 = pi / 3.0;
const f64 pio4 = pi / 4.;
const f64 pio8 = pi / 8.;
const f64 pio16 = pi / 16.;
const f64 pio32 = pi / 32.;
const f64 pio64 = pi / 64.;
const f64 twopi = 2.0 * pi;

const f64 deg_to_rad = pi / 180.;
const f64 rad_to_deg = 180. / pi;

const f64 ft_to_m = 0.3048;
const f64 m_to_ft = 1. / ft_to_m;

const f64 mi_to_km = 1.609344;
const f64 km_to_mi = 1. / mi_to_km;

const f64 km_to_au = 6.6845871226706E-9;
const f64 au_to_km = 1. / km_to_au;

// raylib conversions
const f64 rl_to_u = 1000;
const f64 u_to_rl = 1. / rl_to_u;

// gravitational constant
const f64 G_m = 6.6743e-11;
const f64 G_km = 6.6743e-20;

// tolorances
const f64 tol3 = 1e-3;
const f64 tol6 = 1e-6;
const f64 tol9 = 1e-9;
const f64 tol12 = 1e-12;
const f64 tol16 = 1e-16;