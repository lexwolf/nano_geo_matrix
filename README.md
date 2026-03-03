# nano_geo_matrix

**Header-only C++17 geometry and solver framework for quasi-static and
Mie nanophotonic models.**

`nano_geo_matrix` extracts and unifies the reusable computational core
developed across:

- single_particle
- nanoshell
- mie_single_particle

It provides a modular structure for:

- Geometry definitions (single sphere, nanoshell)
- Quasi-static models
- Full Mie formulations
- Field evaluation utilities
- Backend-selectable complex Bessel implementations

------------------------------------------------------------------------

## Scope

This library focuses on:

- Mode-resolved analysis of metallic nanoparticles
- Gain-assisted plasmonic systems
- Spaser / emission-threshold modeling
- Clean separation between geometry, physics, solver logic, and
  numerical backends

It is designed to be embedded as a subsystem inside larger simulation
projects.

------------------------------------------------------------------------

## Structure

```text
include/nano_geo_matrix/
  core/              # Linear algebra helpers, extraction utilities
  bessel/            # Backend-selectable Bessel wrappers
  mie/
    geometry/
    fields/
  quasi_static/
    geometry/
    spaser/
  cup/               # High-level ready-to-use interfaces
```

The `cup/` namespace preserves the higher-level API used in legacy
simulation codes.

------------------------------------------------------------------------

## Dependencies

Required:

- C++17
- Eigen ≥ 3.3
- Armadillo

Optional (backend dependent):

- Joey Dumont complex_bessel library
- Internal header-only Bessel library
- Fortran special-functions backend (experimental)

------------------------------------------------------------------------

## Bessel Backend Selection

The backend is selected at CMake configure time:

```bash
cmake -S . -B build -DNGM_BESSEL_BACKEND=BL
```

Available options:

- BL (default) --- Header-only backend
- COMPLEX_BESSEL --- libcomplex_bessel
- SPECIAL --- Fortran backend (gfortran required)

Example:

```bash
cmake -S . -B build -DNGM_BESSEL_BACKEND=COMPLEX_BESSEL
```

------------------------------------------------------------------------

## Build and Test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
```

Smoke tests validate:

- quasi-static single sphere
- quasi-static nanoshell
- Mie single sphere
- cup high-level interface

------------------------------------------------------------------------

## Usage as Subsystem

Add as subdirectory:

```cmake
add_subdirectory(nano_geo_matrix)
target_link_libraries(your_target PRIVATE nano_geo_matrix)
```

FetchContent example:

```cmake
FetchContent_Declare(
  nano_geo_matrix
  GIT_REPOSITORY https://github.com/lexwolf/nano_geo_matrix.git
  GIT_TAG v0.1.0
)
FetchContent_MakeAvailable(nano_geo_matrix)
```

------------------------------------------------------------------------

## Scientific Context

This code supports research on:

- Mode-dependent emission thresholds
- Geometry-matrix formulations
- Active nanoshells
- Gain-enhanced metallic nanoparticles

See related publications in associated research repositories for applied
examples.

------------------------------------------------------------------------

## License

MIT License (see LICENSE file)

------------------------------------------------------------------------

## Design Philosophy

- Header-only core for composability
- Explicit backend selection
- Strict geometry/physics separation
- No hidden macros
- No build-system side effects

This repository represents the stabilized computational core extracted
from research codebases.
