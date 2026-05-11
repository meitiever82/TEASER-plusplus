# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

TEASER++ is a C++ library for fast and certifiably-robust 3D point cloud registration (estimating rotation/translation/scale between two point clouds) with Python (pybind11) and MATLAB (mex) bindings. Core algorithms are TLS (Truncated Least Squares), GNC (Graduated Non-Convexity), FGR, and Quatro, combined with PMC max-clique-based outlier pruning.

## Build & test

C++ (out-of-source build at `build/`):
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build           # run all gtest cases (timeout 120s each)
sudo cmake --install build       # installs libs + CMake config under teaserpp::
```

Run a single gtest: tests are bundled into `build/test/teaser/all_tests`. Filter with the standard gtest flag:
```bash
./build/test/teaser/all_tests --gtest_filter=RegistrationTest.*
./build/test/teaser/all_tests --gtest_list_tests
```

Python bindings (uses scikit-build-core; builds a static C++ core + `_teaserpp` extension):
```bash
pip install .
python -c "import teaserpp_python"
```

Examples build automatically with `BUILD_EXAMPLES=ON` (default). The standalone `examples/teaser_cpp_ply` and `examples/teaser_cpp_fpfh` also have their own `CMakeLists.txt` that consume the installed `teaserpp::teaser_registration` package — useful for verifying `make install` worked.

Code style: `.clang-format` (LLVM-based) for C++; `ruff` configured in `pyproject.toml` (line-length 88, double quotes) for Python.

## Key CMake options

Set with `-D<OPTION>=ON/OFF`:

- `BUILD_TEASER_FPFH` (OFF) — builds the `teaser_features` lib (FPFH descriptor + matcher); requires PCL ≥ 1.8 and Boost. Adds `feature-test.cc` / `matcher-test.cc` to tests.
- `BUILD_PYTHON_BINDINGS` (OFF; forced ON via `pyproject.toml`)
- `BUILD_MATLAB_BINDINGS` (OFF) — builds `teaser_mex`
- `BUILD_SHARED_LIBS` (ON; `pyproject.toml` overrides to OFF for wheels)
- `ENABLE_MKL` (OFF) — links Eigen against Intel MKL
- `ENABLE_DIAGNOSTIC_PRINT` (ON) — defines `TEASER_DIAG_PRINT`; disable for quiet runs
- `BUILD_WITH_MARCH_NATIVE` (OFF) — auto-enabled if PCL was built with `-march=native` to avoid Eigen alignment ABI mismatch; in that case `EIGEN_MAX_STATIC_ALIGN_BYTES=0` is also defined.

## Architecture

Two installable libraries live under `teaser/` (one source tree, multiple targets defined in `teaser/CMakeLists.txt`):

- **`teaser_io`** — `ply_io.cc`, uses fetched `tinyply`. Independent of the registration pipeline.
- **`teaser_registration`** — the core solver: `registration.cc`, `certification.cc`, `graph.cc`. Public deps: `Eigen3::Eigen`. Private deps: fetched `pmc` (parallel max clique), fetched `spectra` (eigen solver, header-only, included privately), and `OpenMP` when found.
- **`teaser_features`** (optional, `BUILD_TEASER_FPFH`) — `fpfh.cc`, `matcher.cc`, depends on PCL.

External deps are pulled via CMake `FetchContent` at configure time (PMC, spectra, tinyply, and for tests googletest from `main`). Network access is required on the first configure; subsequent builds reuse `build/_deps/`.

### The registration pipeline (`teaser::RobustRegistrationSolver`)

The public entry point is `RobustRegistrationSolver::solve(src, dst)` in `teaser/include/teaser/registration.h`. It composes several swappable solvers via abstract interfaces — when extending or debugging, work at the right layer:

1. **Scale** (`AbstractScaleSolver`): `TLSScaleSolver` if `params.estimate_scaling==true`, else `ScaleInliersSelector` (scale fixed at 1.0, used only for inlier voting).
2. **Translation-Invariant Measurements (TIMs)** + **Max-clique inlier graph** (`graph.cc`): outliers pruned via `INLIER_SELECTION_MODE` — `PMC_EXACT`, `PMC_HEU`, `KCORE_HEU`, or `NONE`. TIM graph shape controlled by `INLIER_GRAPH_FORMULATION` (`CHAIN` vs. `COMPLETE`).
3. **Rotation** (`GNCRotationSolver`): three concrete implementations selected by `ROTATION_ESTIMATION_ALGORITHM` — `GNC_TLS` (default), `FGR`, `QUATRO`. GNC parameters: `rotation_gnc_factor`, `rotation_max_iterations`, `rotation_cost_threshold`.
4. **Translation** (`TLSTranslationSolver`): component-wise adaptive voting via `ScalarTLSEstimator`.
5. **Certification** (`certification.cc`): `DRSCertifier` for optimality certificates on the rotation estimate (uses `spectra` for eigendecomposition).

All defaults live on `RobustRegistrationSolver::Params`. Reading that struct (registration.h ~L420) is the fastest way to understand what knobs exist.

### Test layout

- `test/teaser/` — gtest unit tests, one per major module (`registration-test.cc`, `tls-test.cc`, `rotation-solver-test.cc`, `certification-test.cc`, …). All compiled into a single `all_tests` binary registered with `gtest_add_tests`.
- `test/test-tools/` — shared test utilities (`test_tools` static lib).
- `test/benchmark/` — `registration-benchmark` binary (not a ctest target).
- Test data files (`*.ply`, `*.pcd`, `*.label`, `*.csv`, `*.txt`) are auto-copied into the build tree by the test `CMakeLists.txt` — paths in tests are relative to the binary's CWD.

### Bindings

- **Python** (`python/teaserpp_python/teaserpp_python.cc`): pybind11 module named `_teaserpp`, re-exported via `teaserpp_python/__init__.py`. Built target output is placed inside the Python package dir so the wheel is self-contained. Type stubs in `_teaserpp.pyi`.
- **MATLAB** (`matlab/teaser_mex.cc`): a mex wrapper plus the `teaser_solve.m` MATLAB-facing function.
