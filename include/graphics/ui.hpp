#pragma once

#include <cstdarg>

#include "imgui.h"

#include "util/typedefs.hpp"
#include "util/vecdefs.hpp"

bool init_render_ui();
void shutdown_render_ui();
void begin_render_ui_frame();
void end_render_ui_frame();

namespace ImGui {

inline void TextSL(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    TextV(fmt, args);
    va_end(args);
    SameLine();
}

inline bool InputDouble3(
    const char* label,
    f64 v[3],
    const char* format = "%.3f",
    ImGuiInputTextFlags flags = 0
) {
    return ImGui::InputScalarN(
        label,
        ImGuiDataType_Double,
        v,
        3,
        nullptr,
        nullptr,
        format,
        flags
    );
}
inline bool InputDouble3(
    const char* label,
    vec3d& v,
    const char* format = "%.3f",
    ImGuiInputTextFlags flags = 0
) {
    return InputDouble3(label, v.data(), format, flags);
}

inline bool InputDouble4(
    const char* label,
    f64 v[4],
    const char* format = "%.3f",
    ImGuiInputTextFlags flags = 0
) {
    return ImGui::InputScalarN(
        label,
        ImGuiDataType_Double,
        v,
        4,
        nullptr,
        nullptr,
        format,
        flags
    );
}
inline bool InputDouble4(
    const char* label,
    vec4d& v,
    const char* format = "%.3f",
    ImGuiInputTextFlags flags = 0
) {
    return InputDouble4(label, v.data(), format, flags);
}

inline bool InputDouble7(
    const char* label,
    f64 v[7],
    const char* format = "%.3f",
    ImGuiInputTextFlags flags = 0
) {
    return ImGui::InputScalarN(
        label,
        ImGuiDataType_Double,
        v,
        7,
        nullptr,
        nullptr,
        format,
        flags
    );
}
inline bool InputDouble7(
    const char* label,
    vec7d& v,
    const char* format = "%.3f",
    ImGuiInputTextFlags flags = 0
) {
    return InputDouble7(label, v.data(), format, flags);
}

inline bool InputDouble3x3(
    const char* label,
    mat3d& m,
    const char* format = "%.3f",
    ImGuiInputTextFlags flags = 0
) {
    bool changed = false;

    ImGui::Text("%s", label);
    ImGui::PushID(label);
    for (i32 r = 0; r < 3; ++r) {
        ImGui::PushID(r);

        double row[3] = {m(r, 0), m(r, 1), m(r, 2)};

        if (ImGui::InputScalarN(
                "##row",
                ImGuiDataType_Double,
                row,
                3,
                nullptr,
                nullptr,
                format,
                flags
            )) {
            m(r, 0) = row[0];
            m(r, 1) = row[1];
            m(r, 2) = row[2];
            changed = true;
        }

        ImGui::PopID();
    }

    ImGui::PopID();

    return changed;
}

template <class Derived>
inline bool InputMatXd(
    const char* label,
    eig::MatrixBase<Derived>& m,
    const char* format = "%.3f",
    ImGuiInputTextFlags flags = 0
) {
    bool changed = false;

    ImGui::Text("%s", label);
    ImGui::PushID(label);

    const i32 rows = static_cast<i32>(m.rows());
    const i32 cols = static_cast<i32>(m.cols());

    for (i32 r = 0; r < rows; ++r) {
        ImGui::PushID(r);

        svec<f64> row(cols);

        for (i32 c = 0; c < cols; ++c) {
            row[c] = m(r, c);
        }

        if (ImGui::InputScalarN(
                "##row",
                ImGuiDataType_Double,
                row.data(),
                static_cast<int>(cols),
                nullptr,
                nullptr,
                format,
                flags
            )) {
            for (i32 c = 0; c < cols; ++c) {
                m(r, c) = row[c];
            }

            changed = true;
        }

        ImGui::PopID();
    }

    ImGui::PopID();

    return changed;
}

} // namespace ImGui