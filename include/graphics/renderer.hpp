#pragma once

#include "core/body.hpp"
#include "core/entity.hpp"
#include "core/state.hpp"
#include "core/world.hpp"
#include "graphics/raygen.hpp"
#include "graphics/rdraw.hpp"
#include "raylib.h"
#include "raymath.h"
#include "util/typedefs.hpp"
#include "util/vecdefs.hpp"

using AssetId = u32;
inline constexpr AssetId kInvalidAssetId = 0;

enum struct RenderPrimitiveKind {
    cube,        // default for satellites
    ellipsoid,   // default for celestials
    cylinder,    // default for stations
    custom_model // user supplied 3d model
};

struct RenderComponent {
    RenderPrimitiveKind kind = RenderPrimitiveKind::cube;

    // sphere {a,a,a}, oblate {a,a,b}, ellipsoid {a,b,c}
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
    svec<RenderBodyInstance> custom_models;
};

struct BuiltinRenderAssets {
    AssetId unit_cube = kInvalidAssetId;
    AssetId unit_sphere_low = kInvalidAssetId;
    AssetId unit_sphere_med = kInvalidAssetId;
    AssetId unit_sphere_high = kInvalidAssetId;
    AssetId unit_cylinder_low = kInvalidAssetId;
    AssetId unit_cylinder_med = kInvalidAssetId;
    AssetId unit_cylinder_high = kInvalidAssetId;
    // TODO: add custom models later, e.g. cube sat with antennae/solar panels
};

inline AssetId choose_render_asset(
    const RenderComponent& rc,
    const BuiltinRenderAssets& builtin
) {
    if (rc.model_override != kInvalidAssetId) return rc.model_override;

    switch (rc.kind) {
    case RenderPrimitiveKind::cube: return builtin.unit_cube;
    case RenderPrimitiveKind::ellipsoid: return builtin.unit_sphere_med;
    case RenderPrimitiveKind::cylinder: return builtin.unit_cylinder_med;
    case RenderPrimitiveKind::custom_model: return rc.model_override;
    }

    return builtin.unit_cube;
}

inline RenderComponent default_render_component_celestial(
    const Celestial& cel,
    Color tint = BLUE
) {
    RenderComponent rc;
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
    rc.kind = RenderPrimitiveKind::cube;
    rc.tint = tint;
    // TODO: placeholder for visualization
    rc.size << 100.0f, 200.0f, 300.0f;
    rc.scale = 2.0f;

    return rc;
}

inline RenderComponent default_render_component_station(
    const Station& stat,
    Color tint = YELLOW
) {
    RenderComponent rc;
    rc.kind = RenderPrimitiveKind::cylinder;
    rc.tint = tint;
    // TODO: placeholder for visualization
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