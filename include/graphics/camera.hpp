#pragma once

#include "raylib.h"

enum struct RCameraMode {
    locked,
    target,
    origin,
};

struct RCamera {
    Camera camera;
    
};