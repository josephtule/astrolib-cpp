#pragma once

#include "util/constants.hpp"
#include "util/vecdefs.hpp"
#include <cmath>
#include <cstddef>


f64 inv_r3_safe(const ecref<vec3d> r, f64 epsilon = tol12);
vec3d accel_gravity_pointmass(
    ecref<vec3d> r_rel,
    f64 mu,
    f64 epsilon = tol12
);
vec3d accel_gravity_zonal(
    ecref<vec3d> r_rel, // relative pos of obj wrt to gravity source in body frame
    f64 mu,             // target body gravitational parameter
    f64 R_cb,           // target body radius
    i32 degree,
    vec7d J, // zonal gravity parameters [J0, J1, J2, J3, J4, J5, J6]
    f64 epsilon = tol12
);

void norm_legendre(
    f64 phi,
    i32 degree,
    i32 order,
    eref<matXd> P,
    eref<matXd> scales,
    f64 tol = tol12
);

vec3d accel_gravity_spherical_harmonics(
    ecref<vec3d> r_rel, // relative pos of obj wrt to gravity source in body frame
    f64 mu,
    f64 R_cb,
    i32 degree,
    i32 order,
    ecref<matXd> C,
    ecref<matXd> S,
    f64 tol = tol12
);

