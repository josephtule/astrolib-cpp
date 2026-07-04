#pragma once

#include "raylib.h"
#include "util/typedefs.hpp"

using AssetId = u32;
inline constexpr AssetId kInvalidAssetId = 0;

enum struct RenderQuality : i32 { low = 1, medium = 2, high = 4 };

enum struct RenderPatternDensity : i32 { low = 1, medium = 2, high = 4 };

struct RenderAssetConfig {
    i32 sphere_segments_base = 16;
    i32 cylinder_segments_base = 16;
    i32 cone_segments_base = 16;
    i32 checker_cells_base = 2;     // number of cells
    i32 checker_texture_size = 128; // resolution
    Color checker_color_a = RAYWHITE;
    Color checker_color_b = GRAY;
};

struct BuiltinRenderAssets {
    AssetId unit_cube = kInvalidAssetId;
    AssetId unit_sphere_low = kInvalidAssetId;
    AssetId unit_sphere_med = kInvalidAssetId;
    AssetId unit_sphere_high = kInvalidAssetId;
    AssetId unit_cylinder_low = kInvalidAssetId;
    AssetId unit_cylinder_med = kInvalidAssetId;
    AssetId unit_cylinder_high = kInvalidAssetId;
    AssetId unit_cone_low = kInvalidAssetId;
    AssetId unit_cone_med = kInvalidAssetId;
    AssetId unit_cone_high = kInvalidAssetId;
    // TODO: may need to reorient cone similiarly to cylinder
    // TODO: add custom models later, e.g. cube sat with antennae/solar panels, add cones
};

struct RenderAssets {
    // TODO: add next_id like in world? not needed until custom loaded assets implemented

    Model unit_cube{};

    Model unit_sphere_low{};
    Model unit_sphere_med{};
    Model unit_sphere_high{};

    Model unit_cylinder_low{};
    Model unit_cylinder_med{};
    Model unit_cylinder_high{};

    Model unit_cone_low{};
    Model unit_cone_med{};
    Model unit_cone_high{};

    Texture checker_low{};
    Texture checker_med{};
    Texture checker_high{};

    BuiltinRenderAssets builtin;
    bool loaded = false;
};

RenderAssets load_render_assets(const RenderAssetConfig& cfg);
void unload_render_assets(RenderAssets& assets);
