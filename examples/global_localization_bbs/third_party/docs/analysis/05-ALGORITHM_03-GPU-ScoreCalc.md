# 3D-BBS — Algorithm 3: GPU Score Computation

> The `__global__ calc_scores_kernel` is the only piece of CUDA in 3D-BBS.
> Everything else runs on the host. This document trace through the kernel,
> the dispatcher, and the host-side batching that feeds it.
>
> Companions: `03-ALGORITHM_01-BranchAndBound.md`, `04-ALGORITHM_02-HierarchicalVoxelmap.md`.

---

## 1. Kernel Source (annotated)

[VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:7-56]

```cpp
__global__ void calc_scores_kernel(
    const thrust::device_ptr<Eigen::Vector4i*>             multi_buckets_ptrs,  // (1)
    const thrust::device_ptr<VoxelMapInfo>                 voxelmap_info_ptr,   // (2)
    const thrust::device_ptr<AngularInfo>                  d_ang_info_vec_ptr,  // (3)
    thrust::device_ptr<DiscreteTransformation<float>>      trans_ptr,           // (4) in/out
    size_t                                                 index_size,          // (5)
    const thrust::device_ptr<Eigen::Vector3f>              points_ptr,          // (6)
    size_t                                                 num_points) {        // (7)

  // one thread = one candidate pose
  const size_t pose_index = threadIdx.x + blockIdx.x * blockDim.x;
  if (pose_index > index_size) return;                              // (A)

  // Pull this pose + the per-level constants for its level.
  DiscreteTransformation<float>& trans = *thrust::raw_pointer_cast(trans_ptr + pose_index);
  const VoxelMapInfo& vinfo  = *thrust::raw_pointer_cast(voxelmap_info_ptr + trans.level);
  const AngularInfo&  ainfo  = *thrust::raw_pointer_cast(d_ang_info_vec_ptr + trans.level);
  const Eigen::Vector4i* buckets = thrust::raw_pointer_cast(multi_buckets_ptrs)[trans.level];

  // Build (R, t) for this pose.
  int score = 0;
  for (size_t i = 0; i < num_points; i++) {
    const Eigen::Vector3f& p = thrust::raw_pointer_cast(points_ptr)[i];

    const Eigen::Vector3f t(trans.x * vinfo.res, trans.y * vinfo.res, trans.z * vinfo.res);   // (B)
    Eigen::Matrix3f R;
    R = Eigen::AngleAxisf(trans.yaw   * ainfo.rpy_res.z() + ainfo.min_rpy.z(), Vector3f::UnitZ())   // (C)
      * Eigen::AngleAxisf(trans.pitch * ainfo.rpy_res.y() + ainfo.min_rpy.y(), Vector3f::UnitY())
      * Eigen::AngleAxisf(trans.roll  * ainfo.rpy_res.x() + ainfo.min_rpy.x(), Vector3f::UnitX());
    const Eigen::Vector3f p_w = R * p + t;                                                    // (D)

    // Quantise and hash.
    const Eigen::Vector3i coord = (p_w.array() * vinfo.inv_res).floor().cast<int>();
    const uint32_t hash = (coord[0] * 73856093) ^ (coord[1] * 19349669) ^ (coord[2] * 83492791);

    // Open addressing.
    for (int j = 0; j < vinfo.max_bucket_scan_count; j++) {
      const uint32_t bi = (hash + j) % vinfo.num_buckets;
      const Eigen::Vector4i b = buckets[bi];
      if (b.x() != coord.x() || b.y() != coord.y() || b.z() != coord.z()) continue;          // (E)
      if (b.w() == 1) { score++; break; }                                                    // (F)
    }
  }
  trans.score = score;
}
```

### Pointer arguments

1. `multi_buckets_ptrs[L]` is a raw `Eigen::Vector4i*` into device memory
   pointing at the bucket array for level `L`. Indexed by `trans.level` — one
   indirect load per pose.
2. `voxelmap_info_ptr[L]` carries `{res, inv_res, num_buckets, max_bucket_scan_count}`
   for level `L`. All four are read by every iteration of the point loop, so
   the compiler keeps them in registers per thread.
3. `d_ang_info_vec_ptr[L]` carries `{num_division, rpy_res, min_rpy}` for level `L`.
   Only `rpy_res` and `min_rpy` are used inside the kernel; `num_division` is
   used host-side by `branch(...)`.
4. `trans_ptr` is the in-out array of `DiscreteTransformation<float>`. The kernel
   reads `level, x, y, z, roll, pitch, yaw` and writes `score`.
5. `index_size` = `transset_size - 1`, used as the upper bound for `pose_index`.
6. `points_ptr` = source cloud on device (shared by every thread).
7. `num_points` = `src_points_.size()`.

### Notable points

**(A) Bound check.** `pose_index > index_size` — `index_size` is
`transset_size - 1`, so threads with index in `[0, transset_size-1]` proceed.
A trailing block (e.g. 13 leftover threads out of 32) returns immediately
[VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:15-18] [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:78].

**(B–D) Transformation.** Note the operation order:

```
p_world = (R · p_body) + t
```

i.e. translation in the world frame, applied **after** rotation. This matches
the equivalent CPU expression `transform * point` where `transform =
Translation · Rz · Ry · Rx` (Eigen multiplies translation last):

```cpp
const Eigen::Vector3d transed_point = transform * points[i];
```
[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:148]

**(C) Rotation composition.** `R = Rz · Ry · Rx`, applied as `R · p`. This is
extrinsic XYZ Euler (or equivalently intrinsic ZYX Euler) — rotate around X
first, then Y, then Z. Same convention as CPU `create_matrix`
[VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:38].

**(E, F) Probe loop ordering.** The GPU checks `coord` first and `w` second.
On a miss it falls through (`continue`) without short-circuiting `w==0`.
See `04-ALGORITHM_02-HierarchicalVoxelmap.md` §5 for the rationale.

---

## 2. Dispatcher — `BBS3D::calc_scores`

[VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:58-92]

```cpp
std::vector<DT<float>> BBS3D::calc_scores(
    const std::vector<DT<float>>& h_transset,
    thrust::device_vector<AngularInfo>& d_ang_info_vec) {

  size_t N = h_transset.size();
  thrust::device_vector<DT<float>> d_transset(N);
  cudaMemcpyAsync(d_transset.data(), h_transset.data(),
                  sizeof(DT<float>) * N, cudaMemcpyHostToDevice, stream);

  constexpr size_t block_size = 32;                                  // (a)
  size_t num_blocks = (N + block_size - 1) / block_size;             // (b)

  calc_scores_kernel<<<num_blocks, block_size, 0, stream>>>(         // (c)
      voxelmaps_ptr_->d_multi_buckets_ptrs_.data(),
      voxelmaps_ptr_->d_voxelmaps_info_.data(),
      d_ang_info_vec.data(),
      d_transset.data(),
      N - 1,                                                         // (d) index_size
      d_src_points_.data(),
      src_points_.size());

  std::vector<DT<float>> h_output(N);
  cudaMemcpyAsync(h_output.data(), d_transset.data(),
                  sizeof(DT<float>) * N, cudaMemcpyDeviceToHost, stream);
  cudaStreamSynchronize(stream);
  return h_output;
}
```

**(a) `block_size = 32`** — one warp. This guarantees no within-block warp
divergence on the bound check. A warp is the smallest unit of CUDA scheduling
on every NVIDIA arch since Tesla, so this is a portable choice.

**(b) `num_blocks = ceil(N/32)`.** Minimal grid; the kernel internally bound-checks.

**(c) Stream-scoped launch.** The class owns `stream` from its constructor
[VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:16] and tears it down in the destructor
[VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:24]. All H↔D copies and the kernel run
on this single stream, so synchronisation reduces to one
`cudaStreamSynchronize`.

**(d) `index_size = N - 1`.** Compare to the kernel's `if (pose_index > index_size) return;` — pose indices ≤ `N-1` are processed. (Threads `pose_index == N-1` *do* run.)

---

## 3. Host-Side Batching

The CPU version scores `init_transset` once with OpenMP, then scores `children`
once per BnB iteration with OpenMP. The GPU version cannot afford that — a
kernel launch per child group (often only ~8 children for the default config)
would be dominated by launch overhead.

The GPU's `localize()` therefore **batches branching** into a host-side
`branch_stock` vector:

```cpp
std::vector<DT<float>> branch_stock;
branch_stock.reserve(branch_copy_size_);                    // default 10000

while (!trans_queue.empty()) {
  if (timeout) break;
  auto trans = trans_queue.top(); trans_queue.pop();

  // Flush stock just before queue becomes empty.
  if (trans_queue.empty() && !branch_stock.empty()) {
    auto out = calc_scores(branch_stock, d_ang_info_vec);   // (i)
    for (const auto& o : out) if (o.score >= best_score) trans_queue.push(o);
    branch_stock.clear();
  }

  if (trans.score < best_score) continue;                   // prune

  if (trans.is_leaf()) { best_trans = trans; best_score = trans.score; }
  else trans.branch(branch_stock,
                    trans.level - 1,
                    static_cast<int>(v_rate_),
                    ang_info_vec[trans.level - 1].num_division);     // (ii)

  if (branch_stock.size() >= branch_copy_size_) {           // (iii)
    auto out = calc_scores(branch_stock, d_ang_info_vec);
    for (const auto& o : out) if (o.score >= best_score) trans_queue.push(o);
    branch_stock.clear();
  }
}
```

[VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:181-223]

Key observations:

- **(i) "Flush near empty"**: when the queue is about to be exhausted but the
  stock still has un-scored children, score them now and push back to keep BnB
  alive. Without this, a small queue near the end would never trigger the
  size-threshold flush at (iii) and the algorithm would terminate
  prematurely.
- **(ii) Stock fill via in-place overload of `branch(...)`**: appends children
  to `branch_stock` rather than constructing a temporary vector
  [VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:41-63].
- **(iii) Threshold-triggered flush**: when `branch_stock.size() >=
  branch_copy_size_` (default 10000), score the batch on the GPU. The
  threshold trades launch-amortisation against queue freshness — too small,
  and the kernel cost dominates; too large, and the BnB explores in a more
  breadth-first manner that prunes less aggressively.

### 3.1 Effect on prune effectiveness

Because branches sit in `branch_stock` *unscored* until the threshold trips,
the priority queue temporarily holds fewer leaves than the CPU equivalent at
the same point in time. This means `best_score` may stay lower for longer,
and pruning at `trans.score < best_score` is less aggressive.

In other words: the GPU version trades **slightly worse pruning** for **much
better single-iteration throughput**. The net effect is favourable in
practice — the README cites ~189 ms / scan including this batching overhead
[VERIFY: README.md:59].

### 3.2 `branch_copy_size_` knob

Exposed via `set_branch_copy_size(int)`
[VERIFY: bbs3d/include/gpu_bbs3d/bbs3d.cuh:48]; default `10000`
[VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:9].

Tuning guidance derivable from the code:
- A single kernel launch dispatches `ceil(N/32)` blocks, each with one warp.
- `N = 10000 ⇒ 313` blocks, plenty of work to hide launch overhead on a 2060.
- Smaller values (e.g. 1024) would over-launch; larger values (e.g. 65536)
  would delay pruning. The default is well within the saturation regime for
  any modern GPU.

---

## 4. Memory Traffic Per Pose

The kernel per pose reads:
- `DiscreteTransformation<float>` slot: 32 bytes (8 ints) → from L2 typically.
- `VoxelMapInfo`: 16 bytes.
- `AngularInfo`: ~24 bytes (3 ints + 3 floats + 3 floats; struct padded).
- `Eigen::Vector4i*` per level (one pointer): 8 bytes.

And per point (`N_src` iters):
- 12 bytes `Vector3f` for the point.
- Up to `max_bucket_scan_count = 10` bucket reads × 16 bytes.

So **memory cost per pose ≈ `N_src · (12 + 10·16) = N_src · 172` bytes**.

For `N_src = 1000` and `transset_size = 10000`, that's
`1000 · 172 · 10000 ≈ 1.72 GB` of bucket-array reads per kernel launch —
heavily L1/L2-cache-hit because the same bucket array is read by every
thread in the warp.

---

## 5. Per-Frame Lifecycle of GPU State

```
constructor:
  cudaStreamCreate(&stream)

set_tar_points(...)            (one-time per session)
   |__ VoxelMaps::create_voxelmaps()         // host build
       |__ set_buckets_on_device()
              |__ cudaMemcpyAsync(d_multi_buckets_)
              |__ cudaMemcpyAsync(d_multi_buckets_ptrs_)
              |__ cudaMemcpyAsync(d_voxelmaps_info_)
              |__ cudaStreamSynchronize

set_src_points(src)            (per frame)
   |__ cudaMemcpyAsync(d_src_points_, src, H→D)

localize()                     (per frame)
   |__ ang_info_vec ← calc_angular_info()
   |__ cudaMemcpyAsync(d_ang_info_vec, ang_info_vec, H→D)
   |__ init_transset ← create_init_transset(...)
   |__ calc_scores(init_transset, d_ang_info_vec)
        |__ cudaMemcpyAsync(d_transset, transset, H→D)
        |__ calc_scores_kernel<<<..., stream>>>(...)
        |__ cudaMemcpyAsync(transset, d_transset, D→H)
        |__ cudaStreamSynchronize
   |__ BnB loop with branch_stock + threshold flushes
   |__ global_pose_ ← best_trans.create_matrix(...)

destructor:
  cudaStreamDestroy(stream)
```

All copies and launches share the **same stream**, so the only explicit
synchronisation points are the four `cudaStreamSynchronize` calls in:
- `set_buckets_on_device` once at setup [VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps.cu:143].
- `calc_scores` after each kernel [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:90].

---

## 6. CUDA Error Reporting

Wrapped at every call via the `check_error << cudaError_t` operator-overload
trick [VERIFY: bbs3d/include/gpu_bbs3d/stream_manager/check_error.cuh:10-13]:

```cpp
void CUDACheckError::operator<<(cudaError_t error) const {
  if (error == cudaSuccess) return;
  std::cerr << "warning: " << cudaGetErrorName(error)   << std::endl;
  std::cerr << "       : " << cudaGetErrorString(error) << std::endl;
}
```
[VERIFY: bbs3d/src/gpu_bbs3d/stream_manager/check_error.cu:5-15]

Important: errors **only print** — they don't throw or abort. A failed
`cudaMemcpyAsync` will silently produce garbage data. This is a deliberate
choice for an academic codebase but worth noting if you ever embed 3D-BBS in a
production stack.

---

## 7. Floating-Point Implications

The GPU runs **single precision** end-to-end:
- `set_tar_points` takes `std::vector<Eigen::Vector3f>` [VERIFY: bbs3d/include/gpu_bbs3d/bbs3d.cuh:30].
- `d_src_points_` is `thrust::device_vector<Eigen::Vector3f>` [VERIFY: bbs3d/include/gpu_bbs3d/bbs3d.cuh:110].
- `VoxelMapInfo::res, inv_res` are `float` [VERIFY: bbs3d/include/gpu_bbs3d/voxelmaps.cuh:14-16].
- The kernel does all `floor`, multiplication, and rotation in `float`.

For maps within ±10 km of origin and resolutions ≥ 0.5 m, this is well within
`float`'s precision floor (~2 cm absolute error at the 10 km horizon). For
maps with coordinates approaching ±1 km the user should ensure the
quantisation step doesn't fall below the precision of the dynamic range.

The CPU version uses `double` end-to-end [VERIFY: bbs3d/include/cpu_bbs3d/bbs3d.hpp:26], so this concern doesn't apply there.

---

## 8. Comparison: CPU vs GPU Scoring Strategy

| Aspect | CPU | GPU |
|---|---|---|
| Scoring unit | Per-pose call (function) | Per-pose CUDA thread |
| Parallelism | OpenMP `parallel for` over transset / children | Warp-aligned grid (block=32) |
| Memory layout | `std::vector<DT<double>>` on host | `thrust::device_vector<DT<float>>` |
| Per-pose state | Eigen affine transform (constructed once per pose) | Inline `R = Rz·Ry·Rx`, `t = idx·res` (Eigen on device) |
| Branching cadence | Per-iteration (one parent at a time) | Batched (`branch_copy_size_=10000`) |
| Hash probe early-exit | Yes (CPU optimisation `#54`) | No (kernel SIMT) |
| Tail-block handling | OpenMP handles uneven N gracefully | Kernel `if (pose_index > index_size) return` |
| Latency floor | OpenMP team start ~µs | Kernel launch + 1 H→D + 1 D→H ~10–50 µs |
| Throughput ceiling | RAM bandwidth, ~30 GB/s/socket | GPU bandwidth, ~300 GB/s (RTX 2060) |

---

## 9. Verification Checklist

- [x] Kernel signature, body, and exit conditions read line-by-line
- [x] `block_size = 32` confirmed in dispatcher
- [x] `index_size = transset_size - 1` confirmed alongside `> index_size` bound
- [x] Stream lifecycle: ctor → setup → per-frame → dtor traced
- [x] Batching threshold (`branch_copy_size_ = 10000`) and flush conditions enumerated
- [x] CUDA error wrapper found to print but not throw (potential pitfall)
- [x] FP32 precision tradeoff stated with code references
