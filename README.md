# nano_geo_matrix

**Header-only C++17 geometry and solver framework for quasi-static and Mie nanophotonic models.**

`nano_geo_matrix` extracts and unifies the reusable computational core developed across:

- single_particle
- nanoshell
- mie_single_particle

It provides:

- Geometry definitions (single sphere, nanoshell)
- Quasi-static models
- Full Mie formulations
- Field evaluation utilities
- Backend-selectable complex Bessel implementations
- An optional CUP module with a ready-to-use high-level class and subroutines

---

## Scope

This library focuses on:

- Mode-resolved analysis of metallic nanoparticles
- Gain-assisted plasmonic systems
- Spaser / emission-threshold modeling
- Clean separation between geometry, physics, solver logic, and numerical backends

It is designed to be embedded as a subsystem inside larger simulation projects.

---

## Structure

```text
include/nano_geo_matrix/      # library proper (core reusable headers)
  core/
  bessel/
  mie/
    geometry/
    fields/
  quasi_static/
    geometry/
    spaser/

modules/cup/                  # optional high-level module
  cup.hpp                      # ready-to-use `nanosphere` class facade
  materials.hpp                # material/permittivity routines
  pump.hpp                     # pump profile helpers
  fourier.hpp                  # Fourier utilities
  quasi_static/                # CUP quasi-static routines
  mie/                         # CUP Mie routines
  data/materials/metals/       # Johnson-Christy data tables
```

The core library API is under `include/nano_geo_matrix`.
`modules/cup` is a dedicated module that exposes a ready-to-use class (`nanosphere`) and associated subroutines.

---

## CUP Module

The `cup` module provides a higher-level workflow layer built on top of the
core `nano_geo_matrix` solvers. It includes convenience routines for
material handling, spectral scans, and solver orchestration.

The module is located in:

```text
modules/cup/
```

and can be used independently of the rest of the project.

---

### Directory structure

```text
modules/cup/
│
├── cup.hpp
├── materials.hpp
├── config.hpp
│
├── data/
│   └── materials/
│       └── metals/
│           ├── goldJCeV.dat
│           └── silverJCeV.dat
│
├── mie/
│   └── solvers.hpp
│
└── quasi_static/
    └── solvers.hpp
```

The `data/materials/metals` directory contains Johnson–Christy optical
constants used by the material routines.

---

### Including the CUP module

To use CUP in a standalone project, include the core library and the module:

```cpp
#include <cup.hpp>
```

Your compiler must be able to find:

```text
include/
modules/cup/
```

For example:

```bash
g++ -std=c++17 \
  -I/path/to/nano_geo_matrix/include \
  -I/path/to/nano_geo_matrix/modules/cup \
  ...
```

---

### CUP data files and include paths

The CUP module loads Johnson–Christy material tables from:

```text
modules/cup/data/materials/metals/
```

The location of these files is determined at runtime by
`materials.hpp`, which derives the data path **relative to the
location of the header itself** using the `__FILE__` macro.

This allows the module to remain self-contained without requiring
external configuration.

However, this mechanism relies on the compiler expanding `__FILE__`
to a meaningful path. If the CUP include directory is passed as a
**relative include path**, some compilers may expand `__FILE__`
to a relative path, which can prevent the library from locating
its data files at runtime.

For this reason, it is recommended to pass the CUP include directory
using **absolute include paths**.

Example:

```bash
NGM_ROOT=$(realpath path/to/nano_geo_matrix)

g++ -std=c++17 \
  -I"$NGM_ROOT/include" \
  -I"$NGM_ROOT/modules/cup" \
  -I/usr/include/eigen3 \
  program.cpp -lgsl
```

rather than:

```bash
-I../nano_geo_matrix/include
-I../nano_geo_matrix/modules/cup
```

CMake users typically do not need to worry about this, since CMake
generates absolute include paths automatically.

---

### Typical workflow

A minimal workflow using CUP might look like:

```cpp
#include <cup.hpp>

int main()
{
    cup::nanosphere sphere;

    auto spectrum = sphere.steady_state(
        1,
        (char*)"drude",
        (char*)"silver",
        (char*)"air",
        1.0, 4.0, 1000
    );

    return 0;
}
```

The module then:

- loads the requested material model
- constructs the solver
- performs the spectral scan
- returns the computed spectrum

---

### Design philosophy

The CUP module intentionally trades some of the minimalism of the
core library for convenience.

In particular it provides:

- preconfigured solver workflows
- built-in material tables
- higher-level orchestration routines

while the core `nano_geo_matrix` library remains focused on the
low-level electromagnetic solvers.

Users who need maximal control can use the core solvers directly,
while CUP offers a more compact interface for typical research
workflows.

---

## Dependencies

Required:

- C++17
- Eigen ≥ 3.3
- Armadillo
- GSL

Optional (backend dependent):

- Joey Dumont complex_bessel library
- Internal header-only Bessel library
- Fortran special-functions backend (experimental)

---

## Bessel Backend Selection

The backend is selected at CMake configure time:

```bash
cmake -S . -B build -DNGM_BESSEL_BACKEND=BL
```

Available options:

- `BL` (default): header-only backend
- `COMPLEX_BESSEL`: `libcomplex_bessel`
- `SPECIAL`: Fortran backend (gfortran required)

Example:

```bash
cmake -S . -B build -DNGM_BESSEL_BACKEND=COMPLEX_BESSEL
```

---

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

---

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

---

## Scientific Context

This code supports research on:

- Mode-dependent emission thresholds
- Geometry-matrix formulations
- Active nanoshells
- Gain-enhanced metallic nanoparticles

See related publications in associated research repositories for applied examples.

---

## License

MIT License (see LICENSE file)

---

## Design Philosophy

- Header-only core for composability
- Optional high-level module (`modules/cup`) for rapid integration
- Explicit backend selection
- Strict geometry/physics separation
- No hidden macros
- No build-system side effects

This repository represents the stabilized computational core extracted from research codebases.
