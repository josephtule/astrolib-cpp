#include "graphics/render_assets.hpp"

#include "raylib.h"
#include "raymath.h"

static Texture load_checker_texture(
    const RenderAssetConfig& cfg,
    const RenderPatternDensity& density
) {
    i32 cell_px = std::max(
        1,
        cfg.checker_texture_size / (cfg.checker_cells_base * static_cast<i32>(density))
    );

    Image image = GenImageChecked(
        cfg.checker_texture_size,
        cfg.checker_texture_size,
        cell_px,
        cell_px,
        cfg.checker_color_a,
        cfg.checker_color_b
    );
    Texture texture = LoadTextureFromImage(image);
    UnloadImage(image);

    return texture;
}

static Model load_unit_sphere_model(
    const RenderAssetConfig& cfg,
    const RenderQuality& quality,
    const Texture& texture
) {
    i32 segments = cfg.sphere_segments_base * static_cast<i32>(quality);
    Mesh sphere_mesh = GenMeshSphere(1.0f, segments, segments);
    Model sphere_model = LoadModelFromMesh(sphere_mesh);
    sphere_model.transform = MatrixIdentity();
    sphere_model.materials[0] = LoadMaterialDefault();
    sphere_model.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = texture;
    return sphere_model;
}

static Model load_unit_cylinder_model(
    const RenderAssetConfig& cfg,
    const RenderQuality& quality,
    const Texture& texture
) {
    i32 segments = cfg.cylinder_segments_base * static_cast<i32>(quality);
    Mesh cylinder_mesh = GenMeshCylinder(1.0f, 1.0f, segments);
    Model cylinder_model = LoadModelFromMesh(cylinder_mesh);
    cylinder_model.transform = MatrixIdentity();
    cylinder_model.materials[0] = LoadMaterialDefault();
    cylinder_model.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = texture;
    return cylinder_model;
}

static Model load_unit_cone_model(
    const RenderAssetConfig& cfg,
    const RenderQuality& quality,
    const Texture& texture
) {
    i32 segments = cfg.cone_segments_base * static_cast<i32>(quality);
    Mesh cone_mesh = GenMeshCone(1.0f, 1.0f, segments);
    Model cone_model = LoadModelFromMesh(cone_mesh);
    cone_model.transform = MatrixIdentity();
    cone_model.materials[0] = LoadMaterialDefault();
    cone_model.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = texture;
    return cone_model;
}

RenderAssets load_render_assets(const RenderAssetConfig& cfg) {
    RenderAssets assets{};

    assets.checker_low = load_checker_texture(cfg, RenderPatternDensity::low);
    assets.checker_med = load_checker_texture(cfg, RenderPatternDensity::medium);
    assets.checker_high = load_checker_texture(cfg, RenderPatternDensity::high);

    Mesh cube_mesh = GenMeshCube(1.0f, 1.0f, 1.0f);
    assets.unit_cube = LoadModelFromMesh(cube_mesh);
    assets.unit_cube.transform = MatrixIdentity();
    assets.unit_cube.materials[0] = LoadMaterialDefault();
    assets.unit_cube.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = assets.checker_low;

    assets.unit_sphere_low
        = load_unit_sphere_model(cfg, RenderQuality::low, assets.checker_low);
    assets.unit_sphere_med
        = load_unit_sphere_model(cfg, RenderQuality::medium, assets.checker_med);
    assets.unit_sphere_high
        = load_unit_sphere_model(cfg, RenderQuality::high, assets.checker_high);

    assets.unit_cylinder_low
        = load_unit_cylinder_model(cfg, RenderQuality::low, assets.checker_low);
    assets.unit_cylinder_med
        = load_unit_cylinder_model(cfg, RenderQuality::medium, assets.checker_low);
    assets.unit_cylinder_high
        = load_unit_cylinder_model(cfg, RenderQuality::high, assets.checker_low);

    assets.unit_cone_low
        = load_unit_cone_model(cfg, RenderQuality::low, assets.checker_med);
    assets.unit_cone_med
        = load_unit_cone_model(cfg, RenderQuality::medium, assets.checker_med);
    assets.unit_cone_high
        = load_unit_cone_model(cfg, RenderQuality::high, assets.checker_med);

    assets.builtin.unit_cube = 1;
    assets.builtin.unit_sphere_low = 2;
    assets.builtin.unit_sphere_med = 3;
    assets.builtin.unit_sphere_high = 4;
    assets.builtin.unit_cylinder_low = 5;
    assets.builtin.unit_cylinder_med = 6;
    assets.builtin.unit_cylinder_high = 7;
    assets.builtin.unit_cone_low = 8;
    assets.builtin.unit_cone_med = 9;
    assets.builtin.unit_cone_high = 10;

    assets.loaded = true;
    return assets;
}

void unload_render_assets(RenderAssets& assets) {
    if (!assets.loaded) return;

    UnloadModel(assets.unit_cube);
    UnloadModel(assets.unit_sphere_low);
    UnloadModel(assets.unit_sphere_med);
    UnloadModel(assets.unit_sphere_high);
    UnloadModel(assets.unit_cylinder_low);
    UnloadModel(assets.unit_cylinder_med);
    UnloadModel(assets.unit_cylinder_high);
    UnloadModel(assets.unit_cone_low);
    UnloadModel(assets.unit_cone_med);
    UnloadModel(assets.unit_cone_high);

    UnloadTexture(assets.checker_low);
    UnloadTexture(assets.checker_med);
    UnloadTexture(assets.checker_high);

    assets = RenderAssets{};
}
