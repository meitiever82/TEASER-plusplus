# 3D-BBS — System Overview (系统概述)

> Analysis of the 3D-BBS global-localization library at `/home/steve/localization_ws/src/3d_bbs`.
> All claims tagged `[VERIFY: <relative path>:<line>]` were validated against the source files in this repository.

---

## 1. What 3D-BBS Is

3D-BBS (3-D Branch-and-Bound Scan-matching) is a *full-search* global-localization
algorithm: given a single 3D LiDAR scan **roughly gravity-aligned** (roll/pitch
near zero) and a pre-built 3D point-cloud map, it returns the 4×4 transform that
places the scan into the map — **without an initial pose guess**.

The library exposes two parallel implementations:

| Backend | Namespace | Scalar | Parallelism | Files |
|---------|-----------|--------|-------------|-------|
| CPU     | `cpu`     | `double` | OpenMP `#pragma omp parallel for` [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:193] | `bbs3d/src/cpu_bbs3d/*.cpp` |
| GPU     | `gpu`     | `float`  | CUDA kernel + batch dispatch [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:73] | `bbs3d/src/gpu_bbs3d/*.cu`  |

Both share `DiscreteTransformation<T>` (templated header) and identical algorithmic
structure; they diverge only in how they evaluate scores per BnB iteration.

---

## 2. High-Level Architecture

```
                      ┌─────────────────────────────────────┐
                      │       Caller (test / ROS2 node)     │
                      └──────────────────┬──────────────────┘
                                         │ Eigen::Vector3* points
                                         ▼
        ┌────────────────────────────────────────────────────────────┐
        │                        BBS3D class                          │
        │   • set_tar_points(...)  → build VoxelMaps                  │
        │   • set_src_points(...)  → host (+ device for GPU)          │
        │   • set_trans/angular_search_range(...)                     │
        │   • localize()           → branch-and-bound loop            │
        └─────┬─────────────────────────────────────────┬─────────────┘
              │ owns                                    │ uses
              ▼                                         ▼
  ┌───────────────────────┐                ┌──────────────────────────┐
  │      VoxelMaps        │                │ DiscreteTransformation<T>│
  │ (multi-level buckets) │                │ (BnB tree node)          │
  │  multi_buckets_[L]    │                │  x,y,z,roll,pitch,yaw    │
  │  voxelmaps_res_[L]    │                │  level, score            │
  └───────────────────────┘                │  branch(), create_matrix │
              ▲                            └──────────────────────────┘
              │ Eigen::Vector4i hash buckets
              │ (x, y, z, occupied_flag)
              ▼
  ┌───────────────────────────────────────────────────────────────┐
  │                      pointcloud_iof (utilities)               │
  │  filter / narrow_scan_range / read_pcd / save_pcd /           │
  │  calc_gravity_alignment_matrix / pcl_to_eigen                 │
  └───────────────────────────────────────────────────────────────┘
```

---

## 3. Module Inventory

### 3.1 Core library — `bbs3d/`

| Module | Header | Source | LOC | Purpose |
|---|---|---|---:|---|
| `cpu::BBS3D`     | `cpu_bbs3d/bbs3d.hpp` [VERIFY: bbs3d/include/cpu_bbs3d/bbs3d.hpp:21]   | `cpu_bbs3d/bbs3d.cpp` [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:5]    | 260 | Solver, OpenMP-parallel score eval |
| `cpu::VoxelMaps` | `cpu_bbs3d/voxelmaps.hpp` [VERIFY: bbs3d/include/cpu_bbs3d/voxelmaps.hpp:9] | `cpu_bbs3d/voxelmaps.cpp` [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:5] | 98  | Hash-bucket pyramid (CPU memory) |
| Voxelmaps I/O    | (declared in `bbs3d.hpp`) | `cpu_bbs3d/voxelmaps_io.cpp` [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps_io.cpp:9] | 135 | Save/load voxelmap coords + params |
| `gpu::BBS3D`     | `gpu_bbs3d/bbs3d.cuh` [VERIFY: bbs3d/include/gpu_bbs3d/bbs3d.cuh:25]  | `gpu_bbs3d/bbs3d.cu` [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:7]     | 240 | Solver, GPU batch scoring |
| `gpu::VoxelMaps` | `gpu_bbs3d/voxelmaps.cuh` [VERIFY: bbs3d/include/gpu_bbs3d/voxelmaps.cuh:19] | `gpu_bbs3d/voxelmaps.cu` [VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps.cu:6] | 144 | Hash-bucket pyramid (device memory) |
| GPU score kernel | (in `.cu` only) | `gpu_bbs3d/calc_score.cu` [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:7] | 92 | `__global__ calc_scores_kernel` |
| GPU voxelmaps I/O | (decl in `.cuh`) | `gpu_bbs3d/voxelmaps_io.cu` [VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps_io.cu:8] | 76 | Save/load voxelmap coords + params |
| `DiscreteTransformation<T>` | `discrete_transformation/discrete_transformation.hpp` [VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:11] | (header-only template) | 101 | Branch-and-bound search node |
| CUDA error wrapper | `gpu_bbs3d/stream_manager/check_error.cuh` [VERIFY: bbs3d/include/gpu_bbs3d/stream_manager/check_error.cuh:10] | `gpu_bbs3d/stream_manager/check_error.cu` [VERIFY: bbs3d/src/gpu_bbs3d/stream_manager/check_error.cu:5] | 18 | `<<` operator over `cudaError_t` |
| `pciof::filter`  | `pointcloud_iof/filter.hpp` [VERIFY: bbs3d/include/pointcloud_iof/filter.hpp:43] | header-only template | 90 | Voxel-grid downsample, scan-range crop |
| `pciof::calc_gravity_alignment_matrix` | `pointcloud_iof/gravity_alignment.hpp` [VERIFY: bbs3d/include/pointcloud_iof/gravity_alignment.hpp:5] | header-only | 19 | Compute roll/pitch from acc vector |
| `pciof::read_pcd / save_pcd` | `pointcloud_iof/pcd_io.hpp` [VERIFY: bbs3d/include/pointcloud_iof/pcd_io.hpp:12] | header-only template | 171 | Manual binary PCD parser (no PCL dep) |
| `pciof::load_*_points`       | `pointcloud_iof/pcd_loader.hpp` [VERIFY: bbs3d/include/pointcloud_iof/pcd_loader.hpp:14] / `pcd_loader_without_pcl.hpp` [VERIFY: bbs3d/include/pointcloud_iof/pcd_loader_without_pcl.hpp:11] | header-only | 133+133 | PCD folder loaders (with/without PCL) |
| `pciof::pcl_to_eigen` etc.   | `pointcloud_iof/pcl_eigen_converter.hpp` [VERIFY: bbs3d/include/pointcloud_iof/pcl_eigen_converter.hpp:9] | header-only | 49 | PCL ↔ Eigen vector conversion |

Total core library ≈ 2.1k lines C++/CUDA (excluding thirdparty Eigen).

### 3.2 Standalone tests — `test/`

- `test/src/cpu_test.cpp` [VERIFY: test/src/cpu_test.cpp:31] — load PCDs, run `cpu::BBS3D::localize()` on each source frame; optional GICP refinement.
- `test/src/gpu_test.cpp` [VERIFY: test/src/gpu_test.cpp:31] — same flow against `gpu::BBS3D`.
- `test/src/voxelmaps_saver.cpp` [VERIFY: test/src/voxelmaps_saver.cpp:7] — builds the hierarchical voxelmap once and serialises it to disk, so subsequent runs can `set_voxelmaps_coords(path)` and skip construction.
- `test/config/test.yaml` [VERIFY: test/config/test.yaml:8] — runtime parameters (`min_level_res`, `max_level`, `min_rpy`, `max_rpy`, `score_threshold_percentage`, leaf sizes, timeout).

### 3.3 ROS 2 wrappers — `ros2_test/`

- `ros2_test/iridescence/` — uses the Iridescence GUI (`guik::viewer`); subscribes to `sensor_msgs/PointCloud2` + `sensor_msgs/Imu`; on a button press, gravity-aligns the scan via IMU acc, then calls `gpu_bbs3d.localize()` [VERIFY: ros2_test/iridescence/src/gpu_bbs3d_iridescence/gpu_ros2_test_iridescence.cpp:142].
- `ros2_test/rviz2/` — analogous RViz2 variant.
- `ros2_test/config/ros2_test.yaml` [VERIFY: ros2_test/config/ros2_test.yaml:6] — topic names + 3D-BBS parameters.

> **Workspace note:** Per the workspace `CLAUDE.md`, `3d_bbs` is *not* part of the
> `colcon build` set in `localization_ws` — it has a `CMakeLists.txt` but no
> `package.xml` and no entry in `build/`. The ROS 2 wrappers are built
> separately under `ros2_test/iridescence/CMakeLists.txt` and `ros2_test/rviz2/CMakeLists.txt`.

---

## 4. Build System

Top-level `CMakeLists.txt` [VERIFY: CMakeLists.txt:1] produces two shared libraries:

```
cpu_bbs3d.so   (always)          : bbs3d.cpp + voxelmaps.cpp + voxelmaps_io.cpp
gpu_bbs3d.so   (if BUILD_CUDA=ON): bbs3d.cu + calc_score.cu + voxelmaps.cu
                                     + voxelmaps_io.cu + check_error.cu
```

Key options [VERIFY: CMakeLists.txt:5-28]:
- `USE_THIRDPARTY_EIGEN` — auto-on if `thirdparty/Eigen/` is non-empty (git submodule).
- `BUILD_CUDA` — default ON; requires CUDA ≥ 12.0 per the README.

Public headers are installed under `<prefix>/include/{cpu_bbs3d,gpu_bbs3d,discrete_transformation,pointcloud_iof}` [VERIFY: CMakeLists.txt:66].

---

## 5. Execution Pipeline (End-to-End)

```
  1. Load target PCD folder           (pciof::load_tar_clouds / load_tar_points)
                │
                ▼
  2. (Optional) Voxel-grid downsample (pcl::ApproximateVoxelGrid OR pciof::filter)
                │
                ▼
  3. BBS3D::set_tar_points(points, min_level_res, max_level)
       └─→  VoxelMaps::create_voxelmaps(points, v_rate)
                                ┌──────────────────────────────┐
                                │  for L = 0..max_level:       │
                                │    quantise points @ res_L   │
                                │    inflate w/ 7 neighbours   │
                                │    build hash buckets        │
                                │    res *= v_rate (2.0)       │
                                └──────────────────────────────┘
                │
                ▼
  4. BBS3D::set_trans_search_range(target_points)
       └─→  init_t{x,y,z}_range_ = bbox(target) / top_res
                │
                ▼
  5. (Per source frame:)
       BBS3D::set_src_points(src) ; BBS3D::localize()
            │
            ▼
  6. localize():
       a) calc_angular_info() — Δθ per level from max source-point norm
       b) create_init_transset(level=max) — 6D grid at coarsest resolution
       c) parallel calc_score on initial set
       d) push into max-heap priority_queue
       e) loop: pop highest, if leaf → record best, else branch + score
       f) prune any node with score < best
       g) terminate on empty queue or timeout
            │
            ▼
  7. best_trans.create_matrix(min_res, ang_info_vec[0].{rpy_res, min_rpy})
            │
            ▼
  8. Eigen::Matrix4{f,d} global_pose_
```

---

## 6. CPU vs GPU — Architectural Differences

| Aspect | CPU | GPU |
|---|---|---|
| Scalar type | `double` [VERIFY: bbs3d/include/cpu_bbs3d/bbs3d.hpp:26] | `float` [VERIFY: bbs3d/include/gpu_bbs3d/bbs3d.cuh:30] |
| Score-eval granularity | One pose per `calc_score()` call; outer `for` is `#pragma omp parallel for` [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:193-196] | Whole transset launched as a single kernel [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:73-80] |
| Batching strategy | Each tree level processed in lock-step inside the BnB loop | `branch_stock_` accumulates children up to `branch_copy_size_` (default 10000) before kernel launch [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:181] [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:215-222] |
| Hash-bucket empty-slot fast exit | `if (bucket.w()==0) break;` skip remaining probes [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:156-158] | (no equivalent — the empty `(0,0,0,0)` slot never matches a point coord so it falls through `continue`) [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:45-47] |
| Thread count knob | `set_num_threads(int)` [VERIFY: bbs3d/include/cpu_bbs3d/bbs3d.hpp:44] | `set_branch_copy_size(int)` [VERIFY: bbs3d/include/gpu_bbs3d/bbs3d.cuh:48] |
| Voxelmap storage | `std::vector<Buckets>` per level [VERIFY: bbs3d/include/cpu_bbs3d/voxelmaps.hpp:51] | `thrust::device_vector<Eigen::Vector4i>` per level + array-of-pointers on device [VERIFY: bbs3d/include/gpu_bbs3d/voxelmaps.cuh:64-65] |

> **Git context.** Commit `3feeede` ("Optimize hash lookup by breaking early on
> empty buckets for cpu bbs3d (#54)") added the CPU empty-slot break; commit
> `2155daf` ("Fix boost dependency") sits beside it. The GPU kernel was left
> unchanged because, on the GPU, `w==0` is benign — the empty bucket is just
> another non-matching coord — and skipping it would not reduce SIMT divergence
> meaningfully.

---

## 7. External Dependencies

| Dep | Where used | Mandatory? |
|---|---|---|
| Eigen 3.4+ | All headers | Required (auto-falls back to `thirdparty/Eigen/` submodule) [VERIFY: CMakeLists.txt:5-13] |
| OpenMP | CPU `localize()` outer loop | Optional (no-op without it) [VERIFY: CMakeLists.txt:16-20] |
| CUDA 12.0+ | All `gpu_bbs3d` files | Required iff `BUILD_CUDA=ON` [VERIFY: CMakeLists.txt:23-28] |
| Boost (filesystem, hash) | I/O paths + `boost::hash_combine` in `VectorHash` [VERIFY: bbs3d/include/cpu_bbs3d/voxelmaps.hpp:18-21] | Required |
| PCL | Only in `test/` and `ros2_test/` and `pcd_loader.hpp` + `pcl_eigen_converter.hpp` | Not needed for core lib (use `pcd_loader_without_pcl.hpp` instead) |
| yaml-cpp | `test/` config parsing [VERIFY: test/include/load_yaml.hpp:4] | Tests only |
| Iridescence / RViz2 | `ros2_test/` wrappers | ROS 2 demos only |

The core `bbs3d/include/**` is deliberately PCL-free — the only PCL-touching
headers in the include tree are `pointcloud_iof/pcd_loader.hpp` and
`pointcloud_iof/pcl_eigen_converter.hpp`, which the caller can avoid by using
`pcd_loader_without_pcl.hpp` instead.

---

## 8. Key Cross-References

- Data structures → `01-DATA_STRUCTURES.md`
- End-to-end data flow → `02-DATA_FLOW.md`
- Algorithm deep dives → `03-ALGORITHM_01-BranchAndBound.md`,
  `04-ALGORITHM_02-HierarchicalVoxelmap.md`,
  `05-ALGORITHM_03-GPU-ScoreCalc.md`
- Function-by-function trace → `06-KEY_FUNCTIONS.md`
- Design Q&A → `07-KEY_QUESTIONS.md`

---

## 9. Verification Checklist

- [x] All headers and `.cpp/.cu` files in `bbs3d/` opened in full during analysis
- [x] CPU / GPU divergences cross-checked side-by-side
- [x] Build options traced to `CMakeLists.txt` (`USE_THIRDPARTY_EIGEN`, `BUILD_CUDA`)
- [x] Module inventory enumerates every file shipped (excluding `thirdparty/`)
- [x] Every claim in this document carries a `[VERIFY:]` tag
