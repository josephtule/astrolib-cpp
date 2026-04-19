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

# Add Raylib
include(FetchContent)

if(EXISTS "${CMAKE_SOURCE_DIR}/external/raylib")
  set(FETCHCONTENT_SOURCE_DIR_RAYLIB "${CMAKE_SOURCE_DIR}/external/raylib" CACHE PATH "" FORCE)
endif()

FetchContent_Declare(
  raylib
  GIT_REPOSITORY https://github.com/raysan5/raylib.git
  GIT_TAG master  
#   GIT_TAG 5.5       # or pin a commit hash
)
FetchContent_MakeAvailable(raylib)

target_link_libraries(${TARGET_NAME} PRIVATE raylib)

# Add Eigen
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/eigen)
target_link_libraries(${TARGET_NAME} PRIVATE Eigen3::Eigen)

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
