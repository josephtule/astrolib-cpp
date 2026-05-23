#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <string>

#include "core/observations.hpp"
#include "core/state.hpp"
#include "util/constants.hpp"
#include "util/units.hpp"
#include "util/vecdefs.hpp"

enum struct IStatusCode : i32 {
    ok,
    invalid_input,
    non_coplanar,
    degenerate_geometry,
    root_solve_failed,
    unsupported_method,
};

inline std::string istatus_string(IStatusCode status) {
    switch (status) {
    case IStatusCode::ok: return "Ok";
    case IStatusCode::invalid_input: return "invalid_input";
    case IStatusCode::non_coplanar: return "non_coplanar";
    case IStatusCode::degenerate_geometry: return "degenerate_geometry";
    case IStatusCode::root_solve_failed: return "root_solve_failed";
    case IStatusCode::unsupported_method: return "unsupported_method";
    default: return "unknown";
    }
}

struct IODResult {
    bool success = false;
    IStatusCode status = IStatusCode::invalid_input;
    StateTr x;
    f64 residual_norm = 0.0;
    i32 iterations = 0;
};

struct IODAnglesObs3 {
    vec3d t; // measurement timestamps
    mat3d L; // LOS unit vectors as columns
    mat3d R; // observer inertial positions as columns
};

bool iod_vec_valid(ecref<vec3d> v, f64 tol = tol12);

bool iod_time_valid(f64 t1, f64 t2, f64 t3, f64 tol = tol12);

IODAnglesObs3 iod_angles3_from_radec(
    const std::array<f64, 3>& t,
    const std::array<vec2d, 3>& radecs,
    const std::array<vec3d, 3>& R,
    UAngle angle_in = UAngle::radian
);
IODAnglesObs3 iod_angles3_from_radec(
    const vec3d& t,
    const matd<2, 3>& radecs,
    const mat3d& R,
    UAngle angle_in = UAngle::radian
);
IODAnglesObs3 iod_angles3_from_radec(
    const vec3d& t,
    const vec3d& ra,
    const vec3d& dec,
    const mat3d& R,
    UAngle angle_in = UAngle::radian
);

IODResult iod_gauss(
    // times
    f64 t1,
    f64 t2,
    f64 t3,
    // line-of-sight unit vectors expressed in the inertial/body-centered frame
    ecref<vec3d> L1,
    ecref<vec3d> L2,
    ecref<vec3d> L3,
    // observation station vectors in body-centered inertial/body-centered frame
    ecref<vec3d> R1,
    ecref<vec3d> R2,
    ecref<vec3d> R3,
    f64 mu,
    f64 tol = tol6
);

IODResult iod_gauss(
    const std::array<f64, 3>& t,
    const std::array<vec3d, 3>& L,
    const std::array<vec3d, 3>& R,
    f64 mu,
    f64 tol = tol6
);

IODResult iod_gauss(
    const svec<f64>& t,
    const svec<vec3d>& L,
    const svec<vec3d>& R,
    f64 mu,
    f64 tol = tol6
);

IODResult iod_gauss(
    ecref<vec3d> t,
    ecref<mat3d> L,
    ecref<mat3d> R,
    f64 mu,
    f64 tol = tol6
);

IODResult iod_gauss(const IODAnglesObs3& arc, f64 mu, f64 tol = tol6);

IODResult iod_gibbs(
    ecref<vec3d> r1,
    ecref<vec3d> r2,
    ecref<vec3d> r3,
    f64 mu,
    f64 tol = tol6
);

IODResult iod_herrickgibbs(
    f64 t1,
    f64 t2,
    f64 t3,
    ecref<vec3d> r1,
    ecref<vec3d> r2,
    ecref<vec3d> r3,
    f64 mu,
    f64 tol = tol6
);

IODResult iod_laplace(
    f64 t1,
    f64 t2,
    f64 t3,
    ecref<vec3d> L1,
    ecref<vec3d> L2,
    ecref<vec3d> L3,
    ecref<vec3d> R1,
    ecref<vec3d> R2,
    ecref<vec3d> R3,
    f64 mu,
    vec3d w, // angular velocity of body
    f64 tol = tol6
);
