#pragma once

#include "core/body.hpp"
#include "core/entity.hpp"
#include "raylib.h"

using AssetId = u32;
inline constexpr AssetId kInvalidAssetId = 0;

struct RenderInstance {
    Matrix transform;
    Color tint;
};

struct ModelAsset {
    // model asset = a single model multiple instances can call
    Model model;
};

enum struct BodyRenderKind {
    Cube,
    Ellipsoid,
    Cylinder,
    CustomModel
}; // Ellipsoid = sized sphere