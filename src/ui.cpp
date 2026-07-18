// Copyright 2025-2026 Joseph Tu Le
// SPDX-License-Identifier: Apache-2.0

#include "graphics/ui.hpp"

#include "imgui.h"
#include "implot.h"
#include "rlImGui.h"

bool init_render_ui() {
    rlImGuiSetup(true);
    ImGui::GetIO().ConfigWindowsMoveFromTitleBarOnly = true;
    ImGui::GetIO().ConfigFlags |= ImGuiConfigFlags_DockingEnable;

    ImPlot::CreateContext();
    return true;
}

void shutdown_render_ui() {
    ImPlot::DestroyContext();
    rlImGuiShutdown();
}

void begin_render_ui_frame() { rlImGuiBegin(); }

void end_render_ui_frame() { rlImGuiEnd(); }