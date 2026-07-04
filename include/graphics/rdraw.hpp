#pragma once

#include "core/state.hpp"
#include "core/transform.hpp"
#include "graphics/raygen.hpp"
#include "raylib.h"
#include "util/math.hpp"

#define CYAN CLITERAL(Color){0, 255, 255, 255} // CYAN
#define CUSTOMGRAY CLITERAL(Color){30, 30, 30, 255}

template <class T>
inline void draw_axes(const vec3<T>& r, const vec4<T>& q, f32 scale = 1.0f) {
    vec3d x = ep_rotate_fast_active(q, axis_x) * scale;
    vec3d y = ep_rotate_fast_active(q, axis_y) * scale;
    vec3d z = ep_rotate_fast_active(q, axis_z) * scale;

    Vector3 r_rl = eig_to_rl(r);

    Vector3 x_rl = eig_to_rl(x) + r_rl;
    Vector3 y_rl = eig_to_rl(y) + r_rl;
    Vector3 z_rl = eig_to_rl(z) + r_rl;

    Color c1 = RED;
    Color c2 = GREEN;
    Color c3 = BLUE;
    if (scale < 0.0f) {
        c1 = MAGENTA;
        c2 = YELLOW;
        c3 = CYAN;
    }

    DrawLine3D(r_rl, x_rl, c1);
    DrawLine3D(r_rl, y_rl, c2);
    DrawLine3D(r_rl, z_rl, c3);
}
inline void draw_axes(const StateTr& x_tr, const StateAtt& x_att, f32 scale = 1.0f) {
    draw_axes(x_tr.r, x_att.q, scale);
}

enum class GridPlane { XY, XZ, ZY, YZ, YX, ZX };
inline Vector3 map_plane_uv(f32 u, f32 v, GridPlane plane) {
    switch (plane) {
    case GridPlane::XY: // u = x, v = y
        return Vector3{u, v, 0.0f};
    case GridPlane::XZ: // u = x, v = z
        return Vector3{u, 0.0f, v};
    case GridPlane::ZY: // u = z, v = y
        return Vector3{0.0f, v, u};
    case GridPlane::YZ: // u = y, v = z
        return Vector3{0.0f, u, v};
    case GridPlane::YX: // u = y, v = x
        return Vector3{v, u, 0.0f};
    case GridPlane::ZX: // u = z, v = x
        return Vector3{v, 0.0f, u};
    default: return Vector3{0.0f, 0.0f, 0.0f};
    }
}

inline void map_plane_color(
    Color& u_pos_color,
    Color& u_neg_color,
    Color& v_pos_color,
    Color& v_neg_color,
    GridPlane plane
) {
    const Color pos_x_color = RED;
    const Color pos_y_color = GREEN;
    const Color pos_z_color = BLUE;

    const Color neg_x_color = MAGENTA;
    const Color neg_y_color = YELLOW;
    const Color neg_z_color = CYAN;

    // default colors
    u_pos_color = RED;
    u_neg_color = MAGENTA;
    v_pos_color = GREEN;
    v_neg_color = YELLOW;

    switch (plane) {
    case GridPlane::XY: // u = x, v = y
        u_pos_color = pos_x_color;
        u_neg_color = neg_x_color;
        v_pos_color = pos_y_color;
        v_neg_color = neg_y_color;
        break;
    case GridPlane::XZ: // u = x, v = z
        u_pos_color = pos_x_color;
        u_neg_color = neg_x_color;
        v_pos_color = pos_z_color;
        v_neg_color = neg_z_color;
        break;
    case GridPlane::ZY: // u = z, v = y
        u_pos_color = pos_z_color;
        u_neg_color = neg_z_color;
        v_pos_color = pos_y_color;
        v_neg_color = neg_y_color;
        break;
    case GridPlane::YZ: // u = y, v = z
        u_pos_color = pos_y_color;
        u_neg_color = neg_y_color;
        v_pos_color = pos_z_color;
        v_neg_color = neg_z_color;
        break;
    case GridPlane::YX: // u = y, v = x
        u_pos_color = pos_y_color;
        u_neg_color = neg_y_color;
        v_pos_color = pos_x_color;
        v_neg_color = neg_x_color;
        break;
    case GridPlane::ZX: // u = z, v = x
        u_pos_color = pos_z_color;
        u_neg_color = neg_z_color;
        v_pos_color = pos_x_color;
        v_neg_color = neg_x_color;
        break;
    }
}

inline void draw_grid(
    const vec2f& min,
    const vec2f& max,
    const vec2f& inc,
    GridPlane plane = GridPlane::XY,
    bool color_axes = true,
    const Color& grid_color = {180, 180, 180, 180}
) {
    f32 umin = min(0);
    f32 vmin = min(1);
    f32 umax = max(0);
    f32 vmax = max(1);
    f32 uinc = inc(0);
    f32 vinc = inc(1);

    if (uinc <= 0.0f || vinc <= 0.0f) return;
    if (umax < umin || vmax < vmin) return;

    auto make_point
        = [plane](f32 u, f32 v) -> Vector3 { return map_plane_uv(u, v, plane); };

    Color u_pos_color = grid_color;
    Color u_neg_color = grid_color;
    Color v_pos_color = grid_color;
    Color v_neg_color = grid_color;

    if (color_axes) {
        u_pos_color = grid_color;
        u_neg_color = grid_color;
        v_pos_color = grid_color;
        v_neg_color = grid_color;
        const Color pos_x_color = RED;
        const Color pos_y_color = GREEN;
        const Color pos_z_color = BLUE;

        const Color neg_x_color = MAGENTA;
        const Color neg_y_color = YELLOW;
        const Color neg_z_color = CYAN;
        map_plane_color(u_pos_color, u_neg_color, v_pos_color, v_neg_color, plane);
    }

    // draw central axes
    // draw +-u axis
    DrawLine3D(make_point(0.0f, 0.0f), make_point(umax, 0.0f), u_pos_color);
    DrawLine3D(make_point(0.0f, 0.0f), make_point(-umax, 0.0f), u_neg_color);
    // draw +-v axis
    DrawLine3D(make_point(0.0f, 0.0f), make_point(0.0f, vmax), v_pos_color);
    DrawLine3D(make_point(0.0f, 0.0f), make_point(0.0f, -vmax), v_neg_color);

    // lines parallel to v axis
    for (f32 u = uinc; u < umax + eps(umax); u += uinc) { // positive
        DrawLine3D(make_point(u, vmin), make_point(u, vmax), grid_color);
    }
    for (f32 u = -uinc; u > umin - eps(umin); u -= uinc) { // negative
        DrawLine3D(make_point(u, vmin), make_point(u, vmax), grid_color);
    }

    // lines parallel to u axis
    for (f32 v = vinc; v < vmax + eps(vmax); v += vinc) { // positive
        DrawLine3D(make_point(umin, v), make_point(umax, v), grid_color);
    }
    for (f32 v = -vinc; v > vmin - eps(vmin); v -= vinc) { // negative
        DrawLine3D(make_point(umin, v), make_point(umax, v), grid_color);
    }
}