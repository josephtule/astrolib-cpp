#include "graphics/ui.hpp"

#include "imgui.h"
#include "implot.h"
#include "rlImGui.h"

bool init_render_ui() {
    rlImGuiSetup(true);
    ImPlot::CreateContext();
    return true;
}

void shutdown_render_ui() {
    ImPlot::DestroyContext();
    rlImGuiShutdown();
}

void begin_render_ui_frame() { rlImGuiBegin(); }

void end_render_ui_frame() { rlImGuiEnd(); }