#!/usr/bin/env python3
"""Run SC + TEASER (or ISC + TEASER) over all submaps in a glim_ros session.

For each submap:
  1. Invoke `global_localization_sc --mode localize ...` against the chosen DB.
  2. Parse the refined 4x4 from the output file.
  3. Compare to the levelled-frame GT (strip roll/pitch from data.txt's T_world_origin).
  4. Record CSV row: submap, n_pts, n_corres, n_inliers, inlier_ratio, rmse,
     rot_err_deg, trans_err_m, elapsed_s, status.

Usage:
  run_sc_isc_eval.py --session DIR --db PATH --bin PATH --out CSV
                      [--algo TEASER|Quatro] [--noise-bound 0.15] [--sc-downsample 0.3]
                      [--top-k 30] [--only N] [--start N --end N]
"""
import argparse
import csv
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


def parse_gt(data_txt: Path):
    """Return T_world_origin as a 4x4 np.ndarray (full RPY+t)."""
    lines = data_txt.read_text().splitlines()
    for i, l in enumerate(lines):
        if l.strip().startswith("T_world_origin"):
            rows = []
            for r in range(4):
                vals = re.findall(r"-?\d+\.?\d*(?:[eE][+-]?\d+)?", lines[i + 1 + r])
                rows.append([float(v) for v in vals[:4]])
            return np.array(rows)
    raise ValueError(f"no T_world_origin in {data_txt}")


def parse_refined_pose(out_txt: Path):
    """Return (T_refined 4x4, inlier_ratio, rmse, elapsed_s, n_corres, n_inliers) or None on fail."""
    if not out_txt.exists():
        return None
    text = out_txt.read_text()
    # 4x4 matrix lines are after "# 4x4 matrix:" header
    m = re.search(r"#\s*4x4 matrix:\s*\n((?:[^#\n]+\n){4})", text)
    if not m:
        return None
    rows = []
    for r in m.group(1).strip().splitlines():
        vals = re.findall(r"-?\d+\.?\d*(?:[eE][+-]?\d+)?", r)
        rows.append([float(v) for v in vals[:4]])
    T = np.array(rows)
    # Quality line
    qm = re.search(r"#\s*quality:.*\n([\d.eE+\- ]+)", text)
    ir, rmse, elapsed = None, None, None
    if qm:
        nums = qm.group(1).split()
        if len(nums) >= 3:
            ir, rmse, elapsed = float(nums[0]), float(nums[1]), float(nums[2])
    return T, ir, rmse, elapsed


def parse_log_stats(log_path: Path):
    """Pull final correspondences/inliers from stdout log (winning candidate)."""
    if not log_path.exists():
        return None, None
    text = log_path.read_text()
    # The "=== Global Localization Result ===" block reports final counts.
    m_c = re.search(r"Correspondences:\s+(\d+)", text)
    m_i = re.search(r"Inliers:\s+(\d+)", text)
    return (int(m_c.group(1)) if m_c else None,
            int(m_i.group(1)) if m_i else None)


def levelled_gt(T_gt: np.ndarray) -> np.ndarray:
    """Strip roll/pitch from T_gt — keep yaw + translation (matches submap_levelled.pcd frame)."""
    R = T_gt[:3, :3]
    yaw = np.arctan2(R[1, 0], R[0, 0])
    R_y = np.array([
        [np.cos(yaw), -np.sin(yaw), 0],
        [np.sin(yaw),  np.cos(yaw), 0],
        [0,            0,           1],
    ])
    T_lev = np.eye(4)
    T_lev[:3, :3] = R_y
    T_lev[:3, 3] = T_gt[:3, 3]
    return T_lev


def pose_error(T_gt_lev, T_est):
    R_err = T_gt_lev[:3, :3].T @ T_est[:3, :3]
    cos_th = max(-1.0, min(1.0, (np.trace(R_err) - 1) / 2))
    rot = float(np.degrees(np.arccos(cos_th)))
    trans = float(np.linalg.norm(T_gt_lev[:3, 3] - T_est[:3, 3]))
    return rot, trans


def run_one(submap_dir: Path, map_path: Path, db_path: Path, bin_path: Path,
            args, work_dir: Path):
    out_file = work_dir / f"{submap_dir.name}.txt"
    log_file = work_dir / f"{submap_dir.name}.log"
    scan = submap_dir / "submap_levelled.pcd"
    data = submap_dir / "data.txt"
    if not scan.exists() or not data.exists():
        return {"status": "MISSING"}

    cmd = [
        str(bin_path), "--mode", "localize",
        "--map", str(map_path),
        "--db", str(db_path),
        "--scan", str(scan),
        "--top-k", str(args.top_k),
        "--algo", args.algo,
        "--sc-downsample", str(args.sc_downsample),
        "--noise-bound", str(args.noise_bound),
        "--output", str(out_file),
    ]
    if args.fallback_bbs:
        cmd += ["--fallback-bbs", str(args.fallback_bbs)]
    t0 = time.time()
    with log_file.open("w") as f:
        rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT).returncode
    elapsed = time.time() - t0

    parsed = parse_refined_pose(out_file)
    if parsed is None:
        return {"status": f"NO_SOLUTION_rc{rc}", "elapsed_s": elapsed}
    T_est, ir, rmse, refine_elapsed = parsed
    n_corres, n_inliers = parse_log_stats(log_file)

    T_gt_lev = levelled_gt(parse_gt(data))
    rot, trans = pose_error(T_gt_lev, T_est)

    out_text = out_file.read_text() if out_file.exists() else ""
    if rc == 0:
        ok_status = "OK_BBS" if "# winning anchor: -1" in out_text else "OK_SC"
    else:
        ok_status = f"FAIL_rc{rc}"
    return {
        "status": ok_status,
        "n_corres": n_corres,
        "n_inliers": n_inliers,
        "inlier_ratio": ir,
        "rmse": rmse,
        "rot_err_deg": rot,
        "trans_err_m": trans,
        "elapsed_s": elapsed,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--session", required=True, type=Path)
    ap.add_argument("--map", required=True, type=Path)
    ap.add_argument("--db", required=True, type=Path)
    ap.add_argument("--bin", required=True, type=Path,
                    help="Path to global_localization_sc executable")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output CSV path (dir for per-submap logs is derived from this)")
    ap.add_argument("--top-k", type=int, default=30)
    ap.add_argument("--algo", default="TEASER")
    ap.add_argument("--sc-downsample", type=float, default=0.3)
    ap.add_argument("--noise-bound", type=float, default=0.15)
    ap.add_argument("--only", type=int, help="Run only submap N")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=None, help="Inclusive end submap id")
    ap.add_argument("--fallback-bbs", default="", help="path to bbs_localize; passes --fallback-bbs through")
    args = ap.parse_args()

    work_dir = args.out.parent / (args.out.stem + "_logs")
    work_dir.mkdir(parents=True, exist_ok=True)

    submaps = sorted(p for p in args.session.iterdir()
                     if p.is_dir() and p.name.isdigit() and len(p.name) == 6)
    if args.only is not None:
        submaps = [p for p in submaps if int(p.name) == args.only]
    else:
        end = args.end if args.end is not None else int(submaps[-1].name)
        submaps = [p for p in submaps if args.start <= int(p.name) <= end]

    rows = []
    header = ["submap", "n_points", "n_corres", "n_inliers", "inlier_ratio",
              "rmse_m", "rot_err_deg", "trans_err_m", "elapsed_s", "status"]
    for sub in submaps:
        pts_path = sub / "points_compact.bin"
        n_pts = pts_path.stat().st_size // 12 if pts_path.exists() else 0
        print(f"=== submap {sub.name} (n_pts={n_pts}) ===", flush=True)
        r = run_one(sub, args.map, args.db, args.bin, args, work_dir)
        r["submap"] = sub.name
        r["n_points"] = n_pts
        rows.append(r)
        if r["status"].startswith("OK"):
            print(f"    rot={r['rot_err_deg']:.2f}°  tr={r['trans_err_m']:.3f}m"
                  f"  inliers={r.get('n_inliers')}  rmse={r.get('rmse'):.3f}")
        else:
            print(f"    {r['status']}")

    with args.out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow([r.get(k, "") for k in header])

    print(f"\nCSV: {args.out}")
    print(f"per-submap logs: {work_dir}")


if __name__ == "__main__":
    main()
