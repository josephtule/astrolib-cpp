# Copyright 2025-2026 Joseph Le
# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.20)

set(PROJECT_NAME astrolib)
set(TARGET_NAME astrolib)

project(${PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)
list(FILTER SOURCES EXCLUDE REGEX ".*/src/legacy/.*")

add_executable(${TARGET_NAME} ${SOURCES})


target_compile_definitions(astrolib PRIVATE
    PROJECT_ROOT="${CMAKE_SOURCE_DIR}"
)

target_include_directories(${TARGET_NAME} PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Add CSPICE (optional)
option(ASTROLIB_ENABLE_CSPICE "Enable the CSPICE source adapter" OFF)
set(
    CSPICE_ROOT
    "${CMAKE_CURRENT_SOURCE_DIR}/external/cspice"
    CACHE PATH
    "CSPICE toolkit root"
)
if(ASTROLIB_ENABLE_CSPICE)
    find_file(
        CSPICE_LIBRARY
        NAMES cspice.a cspice.lib
        PATHS "${CSPICE_ROOT}/lib"
        NO_DEFAULT_PATH
        REQUIRED
    )
    find_file(
        CSPICE_CSUPPORT_LIBRARY
        NAMES csupport.a csupport.lib
        PATHS "${CSPICE_ROOT}/lib"
        NO_DEFAULT_PATH
        REQUIRED
    )
    if(NOT EXISTS "${CSPICE_ROOT}/include/SpiceUsr.h")
        message(FATAL_ERROR "CSPICE headers not found under ${CSPICE_ROOT}")
    endif()

    target_include_directories(${TARGET_NAME} PRIVATE "${CSPICE_ROOT}/include")
    target_link_libraries(
        ${TARGET_NAME}
        PRIVATE
        "${CSPICE_CSUPPORT_LIBRARY}"
        "${CSPICE_LIBRARY}"
    )
    target_compile_definitions(${TARGET_NAME} PRIVATE ASTROLIB_HAS_CSPICE=1)
else()
    target_compile_definitions(${TARGET_NAME} PRIVATE ASTROLIB_HAS_CSPICE=0)
endif()

# Add Raylib
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/raylib)
target_link_libraries(${TARGET_NAME} PRIVATE raylib)

# Add Eigen
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/eigen)
# target_link_libraries(${TARGET_NAME} PRIVATE Eigen3::Eigen)
target_include_directories(${TARGET_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen
)

# Add ImGui/ImPlot/rlImGui
add_library(render_ui_deps STATIC
    external/imgui/imgui.cpp
    external/imgui/imgui_draw.cpp
    external/imgui/imgui_tables.cpp
    external/imgui/imgui_widgets.cpp
    external/imgui/misc/cpp/imgui_stdlib.cpp

    external/implot/implot.cpp
    external/implot/implot_items.cpp

    external/rlimgui/rlImGui.cpp
)
target_include_directories(render_ui_deps PUBLIC
    external/imgui
    external/implot
    external/rlimgui
)
target_include_directories(render_ui_deps PRIVATE
    external/raylib/src
)
target_link_libraries(${TARGET_NAME} PRIVATE render_ui_deps)

# Add ImGuiFD
add_library(imguifd STATIC
    external/imguifd/ImGuiFD.cpp
    external/imguifd/ImGuiFDNative.cpp
)
target_include_directories(imguifd PUBLIC
    external/imguifd
)
target_include_directories(imguifd PRIVATE
    external/imgui
)
target_compile_definitions(imguifd PUBLIC
    IMGUIFD_ENABLE_STL=1
)
target_link_libraries(${TARGET_NAME} PRIVATE imguifd)


include(FetchContent)
FetchContent_Declare(
    nlohmann_json
    URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
)
FetchContent_MakeAvailable(nlohmann_json)
target_link_libraries(${TARGET_NAME} PRIVATE nlohmann_json::nlohmann_json)

# Other options
target_compile_features(${TARGET_NAME} PRIVATE cxx_std_23)
target_compile_options(${TARGET_NAME} PRIVATE
    $<$<AND:$<CONFIG:Release>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>>:-O3>
)
target_compile_definitions(${TARGET_NAME} PRIVATE
    $<$<CONFIG:Release>:NDEBUG>
)
target_compile_options(${TARGET_NAME} PRIVATE
    $<$<AND:$<CONFIG:Release>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>>:-march=native>
)
