#pragma once

#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/state.hpp"
#include "core/world.hpp"

#include "graphics/raygen.hpp"
#include "graphics/rdraw.hpp"
#include "graphics/render_assets.hpp"

#include "raylib.h"
#include "raymath.h"

#include "util/typedefs.hpp"
#include "util/vecdefs.hpp"

struct RenderDrawOptions {
    bool draw_grid_xy = true;
    bool draw_grid_xz = false;
    bool draw_grid_zy = false;
    f32 grid_x_min = -1.0e5f;
    f32 grid_x_max = 1.0e5f;
    f32 grid_y_min = -1.0e5f;
    f32 grid_y_max = 1.0e5f;
    f32 grid_z_min = -1.0e5f;
    f32 grid_z_max = 1.0e5f;
    f32 u_inc = 5000.0f;
    f32 v_inc = 5000.0f;
    bool color_axes = true;
    Color grid_color = {180, 180, 180, 180};

    bool draw_body_axes = false; // add option for individual bodies?
    bool draw_inertial_axes = true;
    f32 inertial_axes_scale = 10000.0f;
    bool draw_labels = false;

    Color background = CUSTOMGRAY;

    bool draw_fps = true;
};

enum struct RenderPrimitiveKind {
    cube,      // default for satellites
    ellipsoid, // default for celestials
    cylinder,  // default for stations
    cone,
    custom_model // user supplied 3d model
};

struct RenderComponent {
    RenderPrimitiveKind kind = RenderPrimitiveKind::cube;

    RenderQuality quality = RenderQuality::medium;
    RenderPatternDensity pattern_density = RenderPatternDensity::medium;

    // sphere {a,a,a}, oblate spheroid {a,a,b}, general ellipsoid {a,b,c}
    vec3f size = vec3f1; // axes scaling
    f32 scale = 1.0f;    // uniform scaling

    AssetId model_override = kInvalidAssetId;
    Color tint = WHITE;
};

struct RenderBodyInstance {
    EntityId entity_id = kInvalidEntityId;
    AssetId asset = kInvalidAssetId;
    RenderPrimitiveKind kind = RenderPrimitiveKind::cube;
    Matrix transform = MatrixIdentity();
    Color tint = WHITE;
};

struct RenderSceneSnapshot {
    svec<RenderBodyInstance> cubes;
    svec<RenderBodyInstance> ellipsoids;
    svec<RenderBodyInstance> cylinders;
    svec<RenderBodyInstance> cones;
    svec<RenderBodyInstance> custom_models;
};

inline AssetId choose_render_asset(
    const RenderComponent& rc,
    const BuiltinRenderAssets& builtin
) {
    if (rc.model_override != kInvalidAssetId) return rc.model_override;

    switch (rc.kind) {
    case RenderPrimitiveKind::cube: return builtin.unit_cube;
    case RenderPrimitiveKind::ellipsoid:
        switch (rc.quality) {
        case RenderQuality::low: return builtin.unit_sphere_low;
        case RenderQuality::medium: return builtin.unit_sphere_med;
        case RenderQuality::high: return builtin.unit_sphere_high;
        };
    case RenderPrimitiveKind::cylinder:
        switch (rc.quality) {
        case RenderQuality::low: return builtin.unit_cylinder_low;
        case RenderQuality::medium: return builtin.unit_cylinder_med;
        case RenderQuality::high: return builtin.unit_cylinder_high;
        }
    case RenderPrimitiveKind::cone:
        switch (rc.quality) {
        case RenderQuality::low: return builtin.unit_cone_low;
        case RenderQuality::medium: return builtin.unit_cone_med;
        case RenderQuality::high: return builtin.unit_cone_high;
        }
    case RenderPrimitiveKind::custom_model: return rc.model_override;
    }

    return builtin.unit_cube;
}

inline RenderComponent default_render_component_celestial(
    const Celestial& cel,
    Color tint = BLUE
) {
    RenderComponent rc;
    rc.quality = RenderQuality::medium;
    rc.pattern_density = RenderPatternDensity::high;
    rc.kind = RenderPrimitiveKind::ellipsoid;
    rc.tint = tint;
    rc.size << cel.semimajor_axis, cel.semimajor_axis, cel.semiminor_axis;
    rc.scale = 1.0f;

    return rc;
}

inline RenderComponent default_render_component_satellite(
    const Satellite& sat,
    Color tint = RED
) {
    RenderComponent rc;
    rc.pattern_density = RenderPatternDensity::low;
    rc.kind = RenderPrimitiveKind::cube;
    rc.tint = tint;
    // TODO: placeholder for visualization
    // TODO: add user defined from json scenario or ui input
    rc.size << 100.0f, 200.0f, 300.0f;
    rc.scale = 2.0f;

    return rc;
}

inline RenderComponent default_render_component_station(
    const Station& stat,
    Color tint = YELLOW
) {
    RenderComponent rc;
    rc.quality = RenderQuality::medium;
    rc.pattern_density = RenderPatternDensity::low;
    rc.kind = RenderPrimitiveKind::cylinder;
    rc.tint = tint;
    // TODO: placeholder for visualization
    // TODO: add user defined from json scenario or ui input
    rc.size << 10.0f, 10.0f, 10.0f;
    rc.scale = 25.0f;

    return rc;
}
inline mat3f rot_model_to_body(RenderPrimitiveKind kind) {
    // raylib cylinder height is local +y
    // render convention is local +z height
    if (kind == RenderPrimitiveKind::cylinder) {
        mat3f R_body_model;
        R_body_model << 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f;
        return R_body_model;
    }
    // TODO: add cone transform too, make +z = height
    return mat3f1;
}

inline RenderBodyInstance make_render_instance(
    EntityId entity_id,
    const StateTr& x_tr,
    const StateAtt& x_att,
    const RenderComponent& rc,
    const BuiltinRenderAssets& assets
) {
    RenderBodyInstance inst;

    inst.transform = make_transform(x_tr, x_att, rc.size, rc.scale);
    mat3f R_body_model = rot_model_to_body(rc.kind);
    mat3f R = get_scaledrot<f32>(inst.transform) * R_body_model;
    set_rotation<f32>(inst.transform, R);

    inst.entity_id = entity_id;
    inst.kind = rc.kind;
    inst.asset = choose_render_asset(rc, assets);
    inst.tint = rc.tint;

    return inst;
}

inline RenderSceneSnapshot build_render_scene_snapshot(
    const World& world,
    const BuiltinRenderAssets& assets
) {
    RenderSceneSnapshot scene;

    for (EntityId id : world.active_entity_ids()) {
        const Body* body = world.body(id);
        if (body == nullptr) continue;
        RenderComponent rc;
        RenderBodyInstance inst;
        switch (body->body_type) {
        case BodyType::celestial: {
            const Celestial* cel = world.celestial(id);
            if (cel == nullptr) continue;
            rc = default_render_component_celestial(*cel);
            inst = make_render_instance(id, cel->x_tr, cel->x_att, rc, assets);
        } break;
        case BodyType::satellite: {
            const Satellite* sat = world.satellite(id);
            if (sat == nullptr) continue;
            rc = default_render_component_satellite(*sat);
            inst = make_render_instance(id, sat->x_tr, sat->x_att, rc, assets);
        } break;
        case BodyType::station: {
            const Station* stat = world.station(id);
            if (stat == nullptr) continue;
            rc = default_render_component_station(*stat);
            StateTr x_tr_stat = world.stat_x_tr_inertial(id);
            StateAtt x_att_stat = world.stat_x_att_inertial(id);
            inst = make_render_instance(id, x_tr_stat, x_att_stat, rc, assets);
        } break;
        case BodyType::unknown: continue;
        }
        switch (inst.kind) {
        case RenderPrimitiveKind::cube: scene.cubes.push_back(inst); break;
        case RenderPrimitiveKind::ellipsoid: scene.ellipsoids.push_back(inst); break;
        case RenderPrimitiveKind::cylinder: scene.cylinders.push_back(inst); break;
        case RenderPrimitiveKind::cone: scene.cones.push_back(inst); break;
        case RenderPrimitiveKind::custom_model:
            scene.custom_models.push_back(inst);
            break;
        }
    }

    return scene;
}

inline void render_scene_snapshot(
    const RenderSceneSnapshot& scene,
    Model& sphere_model,
    Model& cube_model,
    Model& cylinder_model
) {
    // TODO: remove or update later this is for a first pass

    // ellipsoids
    for (const RenderBodyInstance& inst : scene.ellipsoids) {
        sphere_model.transform = inst.transform;
        DrawModel(sphere_model, rlvec30, 1.0f, inst.tint);
    }

    // cubes
    for (const RenderBodyInstance& inst : scene.cubes) {
        cube_model.transform = inst.transform;
        DrawModel(cube_model, rlvec30, 1.0f, inst.tint);
    }

    // cylinders
    for (const RenderBodyInstance& inst : scene.cylinders) {
        cylinder_model.transform = inst.transform;
        DrawModel(cylinder_model, rlvec30, 1.0f, inst.tint);
    }

    // TODO: add custom
}

inline Model* render_model_from_asset(RenderAssets& assets, AssetId id) {
    if (id == assets.builtin.unit_cube) return &assets.unit_cube;
    if (id == assets.builtin.unit_sphere_low) return &assets.unit_sphere_low;
    if (id == assets.builtin.unit_sphere_med) return &assets.unit_sphere_med;
    if (id == assets.builtin.unit_sphere_high) return &assets.unit_sphere_high;
    if (id == assets.builtin.unit_cylinder_low) return &assets.unit_cylinder_low;
    if (id == assets.builtin.unit_cylinder_med) return &assets.unit_cylinder_med;
    if (id == assets.builtin.unit_cylinder_high) return &assets.unit_cylinder_high;
    if (id == assets.builtin.unit_cone_low) return &assets.unit_cone_low;
    if (id == assets.builtin.unit_cone_med) return &assets.unit_cone_med;
    if (id == assets.builtin.unit_cone_high) return &assets.unit_cone_high;
    return nullptr;
}

inline void draw_render_instance(const RenderBodyInstance& inst, RenderAssets& assets) {
    Model* model = render_model_from_asset(assets, inst.asset);
    if (model == nullptr) return;

    model->transform = inst.transform;
    DrawModel(*model, rlvec30, 1.0f, inst.tint);
}

inline void draw_render_scene(const RenderSceneSnapshot& scene, RenderAssets& assets) {

    for (const RenderBodyInstance& inst : scene.ellipsoids) {
        draw_render_instance(inst, assets);
    }
    for (const RenderBodyInstance& inst : scene.cubes) {
        draw_render_instance(inst, assets);
    }
    for (const RenderBodyInstance& inst : scene.cylinders) {
        draw_render_instance(inst, assets);
    }
    for (const RenderBodyInstance& inst : scene.cones) {
        draw_render_instance(inst, assets);
    }
    for (const RenderBodyInstance& inst : scene.custom_models) {
        draw_render_instance(inst, assets);
    }
}