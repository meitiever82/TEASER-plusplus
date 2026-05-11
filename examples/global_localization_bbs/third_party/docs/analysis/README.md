# 3D-BBS Codebase Analysis

This folder contains a rigorous, code-anchored technical analysis of the
`3d_bbs` library generated via `/codebase-analysis-skill`.

**Project:** `3d_bbs` — Global localization via Branch-and-Bound 3D scan matching.
**Reference paper:** Aoki et al., *3D-BBS: Global Localization for 3D Point
Cloud Scan Matching Using Branch-and-Bound Algorithm*, ICRA 2024.

Every factual claim in these documents carries a `[VERIFY: <path>:<line>]`
tag pointing at the line that supports it; **273 unique references** were
validated against the actual source (no broken refs, no out-of-range lines).

---

## Reading Order

| # | Document | Purpose |
|---|---|---|
| 00 | [SYSTEM_OVERVIEW.md](00-SYSTEM_OVERVIEW.md) | Module inventory, build, CPU vs GPU split |
| 01 | [DATA_STRUCTURES.md](01-DATA_STRUCTURES.md) | `BBS3D`, `VoxelMaps`, `DiscreteTransformation`, hash bucket layout |
| 02 | [DATA_FLOW.md](02-DATA_FLOW.md) | End-to-end pipeline from PCD to `Matrix4` |
| 03 | [ALGORITHM_01-BranchAndBound.md](03-ALGORITHM_01-BranchAndBound.md) | BnB math, score bound, pruning invariants |
| 04 | [ALGORITHM_02-HierarchicalVoxelmap.md](04-ALGORITHM_02-HierarchicalVoxelmap.md) | Multi-resolution hash voxelmap + dilation |
| 05 | [ALGORITHM_03-GPU-ScoreCalc.md](05-ALGORITHM_03-GPU-ScoreCalc.md) | CUDA kernel, batched branching, FP32 implications |
| 06 | [KEY_FUNCTIONS.md](06-KEY_FUNCTIONS.md) | Line-by-line walk of `localize`, `calc_score`, etc. |
| 07 | [KEY_QUESTIONS.md](07-KEY_QUESTIONS.md) | Design-rationale Q&A |

---

## Methodology

Generated via the `codebase-analysis-skill`'s mandatory workflow:

1. Phase 0 — Project-context detection (workspace `CLAUDE.md` consulted).
2. Phase 1 — Read every `.hpp/.cpp/.cu/.cuh` in `bbs3d/`, `test/`, and `ros2_test/` (excluding `thirdparty/`).
3. Phase 2 — Document every data structure with declaration site.
4. Phase 3 — Trace data flow through both CPU and GPU pipelines.
5. Phase 4 — Deep-dive algorithms with mathematical derivations.
6. Phase 5 — Line-by-line function analysis.
7. Phase 6 — Q&A with design-rationale.
8. Phase 7 — **Mandatory verification**: 273 `[VERIFY:]` tags extracted and
   verified against the actual source files. **0 unresolved references.**

---

## Limitations & Honest Notes

- The disk-cache pitfall around `(0,0,0)` voxels (Q12 in
  `07-KEY_QUESTIONS.md`) was identified by careful reading and is **not**
  empirically tested. If you can reproduce or refute it experimentally,
  contributions welcome.
- The analysis is a snapshot of the repo at `HEAD = 41529a3` ("[fix] add
  inline to utility functions (#56)") with the recent CPU optimisation
  commit `3feeede` ("Optimize hash lookup by breaking early on empty buckets
  for cpu bbs3d (#54)") incorporated.
- Performance numbers cited (e.g. ~189 ms / localisation) come from the
  README and were not re-measured.

---

## How to Re-verify

```bash
cd /path/to/3d_bbs
# Extract every VERIFY tag
grep -rhoE 'VERIFY: [^] ,]+:[0-9]+(-[0-9]+)?' docs/analysis/*.md \
  | sed 's/^VERIFY: //' \
  | sort -u
```

For automated checking, see the script pattern in `WORKFLOW.md` of the
skill.
