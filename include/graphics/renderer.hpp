#pragma once

#include "raylib.h"
#include "util/vecdefs.hpp"

using AssetId = u32;
inline constexpr AssetId kInvalidAssetId = 0;

struct RenderInstance {
    AssetId asset = kInvalidAssetId;
    Matrix transform;
    Color tint;
};

enum struct BodyRenderKind {
    Cube,       // default for satellites
    Ellipsoid,  // default for celestials
    Cylinder,   // default for stations
    CustomModel // user supplied 3d model
};

struct RenderComponent {
    BodyRenderKind kind = BodyRenderKind::Cube;

    // sphere {a,a,a}, oblate {a,a,b}, ellipsoid {a,b,c}, TODO: check this
    vec3f size = vec3f1; // axes scaling
    f32 scale = 1.0f;    // uniform scaling

    AssetId model_override = kInvalidAssetId;
    Color tint = WHITE;
};

struct BuiltInAssets {
    AssetId unit_cube = kInvalidAssetId;
    AssetId unit_sphere_low = kInvalidAssetId;
    AssetId unit_sphere_med = kInvalidAssetId;
    AssetId unit_sphere_high = kInvalidAssetId;
    AssetId unit_cylinder_low = kInvalidAssetId;
    AssetId unit_cylinder_med = kInvalidAssetId;
    AssetId unit_cylinder_high = kInvalidAssetId;
    // TODO: add custom primatives later, e.g. cube sat with antennae/solar panels
};

inline AssetId choose_asset(const RenderComponent& rc, const BuiltInAssets& builtin) {
    if (rc.model_override != kInvalidAssetId) return rc.model_override;
    switch (rc.kind) {
    case BodyRenderKind::Cube: return builtin.unit_cube;
    case BodyRenderKind::Ellipsoid: return builtin.unit_sphere_med;
    case BodyRenderKind::Cylinder: return builtin.unit_cylinder_med;
    case BodyRenderKind::CustomModel: return rc.model_override;
    }

    return builtin.unit_cube;
}
