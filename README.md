# astrolib-cpp

C++ astrodynamics and space-simulation sandbox for staged N-body orbit propagation, gravity modeling, station observations, orbit determination, scenario loading, and 3D visualization.

This project is currently diagnostics-driven and under active development. APIs, scenario schema, and renderer structure are expected to change.

## Current Capabilities

- Celestial bodies, satellites, and ground/free stations in a shared `World`
- Translational and attitude state propagation with RK integrators
- Point-mass, zonal, and spherical-harmonic gravity models
- EGM/GFC/NASA SHA gravity coefficient readers
- Time, Julian date, sidereal time, and Earth-orientation utilities
- Station geometry, RA/Dec, Az/El, range, and range-rate measurements
- IOD, batch least-squares, and EKF orbit-determination foundations
- TLE parsing and conversion to satellite initial states
- JSON scenario loading/building
- Primitive history/snapshotting
- Early Raylib-based 3D rendering diagnostics

### Current Task
Scenario exporting/importing and body profiles

## Repository Layout

```text
include/
  core/        astrodynamics, estimation, world, scenario, time, EOP
  graphics/    Raylib rendering helpers
  util/        math, units, typedefs, vector definitions

src/
  diagnostics.cpp    diagnostic/demo entry points
  main.cpp           currently selected diagnostic
  scenario_io.cpp    JSON scenario parser/builder
  world*.cpp         world model, stepping, history, measurements
  ingest.cpp         TLE parsing

scenarios/           JSON scenario files
assets/              local data files, not guaranteed to be committed
external/            git submodules
```

## Requirements

- CMake `>= 3.20`
- C++23-capable compiler
- Git
- Raylib submodule
- Eigen submodule
- Internet access during first configure for `nlohmann_json` FetchContent, unless already cached

On macOS:

```sh
brew install cmake ninja
```

## Clone

```sh
git clone <repo-url>
cd astrolib-cpp
git submodule update --init --recursive
```

The current submodules are:

```text
external/eigen
external/raylib
```

## Data Setup

Create the local asset folders:

```sh
mkdir -p assets/earth assets/moon assets/mars assets/mercury assets/output
```

Several diagnostics and scenarios expect local data files under `assets/`. Data files are not committed to the repository.

See [`docs/data.md`](docs/data.md) for expected filenames and public data sources.

## Build

Recommended:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --target astrolib
```

Without Ninja:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target astrolib
```

## Run

```sh
./build/astrolib
```

The executable currently runs whichever diagnostic is enabled in `src/main.cpp`.


Run a different diagnostics by commenting/uncommenting diagnostic calls in `src/main.cpp`.

## Scenario Files

Example scenarios live in:

```text
scenarios/parser_stress_demo.json
scenarios/earth_station_leo_demo.json
scenarios/earth_moon_sat_demo.json
scenarios/earth_tle_demo.json
```

The scenario loader currently supports:

- Celestial body definitions
- Satellite definitions
- Station definitions
- Instruments
- Gravity providers
- Translational state from position/velocity or orbital elements
- Attitude state from quaternion or axis-angle
- Simple spin attitude model
- Basic propagation flags

The JSON schema is still evolving.

## Rendering

Rendering is currently diagnostic-driven through Raylib. It is not yet a stable application UI.

The current render path can draw bodies, satellites, stations, and scenario-loaded objects. More structured renderer cleanup and UI integration are planned.


## License

TBD.

External data files may have their own licenses or usage restrictions. Check the data provider before redistributing downloaded assets.
