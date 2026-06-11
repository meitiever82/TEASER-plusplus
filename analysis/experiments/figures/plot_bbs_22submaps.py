#!/usr/bin/env python3
"""Generate PNG for BBS+TEASER 22-submap result vs sparsity/inliers.

Output: bbs_22submaps_diagnostic.png next to this script.
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# (submap, n_frames, n_points, inliers, rot_err_deg, trans_err_m)
DATA = [
    ("00", 128, 6788,  34, 0.62, 0.077),
    ("01",  73, 8006,  83, 0.90, 0.036),
    ("02", 104, 6986,  39, 2.27, 0.180),
    ("03",  73, 6973,  38, 3.45, 0.067),
    ("04",  60, 7940,  47, 0.92, 0.132),
    ("05",  45, 6941,  51, 1.71, 0.100),
    ("06",  65, 5675,  41, 1.27, 0.089),
    ("07",  57, 3808,  15, 5.21, 0.279),
    ("08",  59, 3245,  77, 1.51, 0.035),
    ("09", 120, 2802,  30, 0.87, 0.134),
    ("10", 113, 2447,  10, 7.95, 0.349),
    ("11",  30, 3573,  12, 2.02, 0.135),
    ("12",  35, 5952,  84, 1.53, 0.022),
    ("13",  26, 7007,  30, 1.80, 0.106),
    ("14",  82, 6167,  32, 1.33, 0.103),
    ("15",  62, 4672,  20, 2.28, 0.139),
    ("16",  44, 4870,  26, 3.44, 0.217),
    ("17",  44, 3357,   4, 6.53, 0.505),
    ("18",  42, 3097,  53, 0.68, 0.049),
    ("19",  45, 2944,  28, 3.81, 0.106),
    ("20",  35, 4677,  45, 2.40, 0.187),
    ("21",  41, 4345,  18, 1.61, 0.398),
]

SC_FAIL_IDS = {"07", "10", "17"}  # Phase 2 NO_SOLUTION

def status_color(rot):
    if rot > 5.0:   return "#d62728"  # red — fail
    if rot > 3.0:   return "#ff7f0e"  # orange — warning
    return "#2ca02c"                  # green — pass

fig, axes = plt.subplots(1, 3, figsize=(17, 6.5))
plt.subplots_adjust(left=0.05, right=0.99, top=0.84, bottom=0.10, wspace=0.22)

# --- Panel 1: inliers (x) vs rot_err (y) — the "money plot" ---
ax = axes[0]
for (sid, nf, npts, nin, rot, _t) in DATA:
    c = status_color(rot)
    marker = "X" if sid in SC_FAIL_IDS else "o"
    edge = "black" if sid in SC_FAIL_IDS else c
    ax.scatter(nin, rot, c=c, marker=marker, s=120, edgecolors=edge, linewidths=1.5, zorder=3)
    ax.annotate(sid, (nin, rot), fontsize=8.5, xytext=(5, 4), textcoords="offset points")
ax.axhline(5.0, ls="--", color="gray", lw=1)
ax.axhline(3.0, ls=":",  color="gray", lw=1)
ax.axvline(15,  ls="--", color="gray", lw=1)
ax.text(15.5, 9.7, "inliers ≥ 15 (proposed threshold)", color="gray", fontsize=8)
ax.text(85,    5.1, "rot = 5° (deploy threshold)",       color="gray", fontsize=8, ha="right")
ax.text(85,    3.1, "rot = 3° (strict threshold)",       color="gray", fontsize=8, ha="right")
ax.set_xlabel("# inliers (after TEASER GNC-TLS)", fontsize=11)
ax.set_ylabel("rotation error (deg)", fontsize=11)
ax.set_title("Inliers is the dominant predictor of rot error", fontsize=11.5)
ax.set_xlim(-3, 92)
ax.set_ylim(0, 10.5)
ax.grid(True, alpha=0.3)

# --- Panel 2: n_points (x) vs rot_err (y) — weaker correlation ---
ax = axes[1]
for (sid, nf, npts, nin, rot, _t) in DATA:
    c = status_color(rot)
    marker = "X" if sid in SC_FAIL_IDS else "o"
    edge = "black" if sid in SC_FAIL_IDS else c
    ax.scatter(npts, rot, c=c, marker=marker, s=120, edgecolors=edge, linewidths=1.5, zorder=3)
    ax.annotate(sid, (npts, rot), fontsize=8.5, xytext=(5, 4), textcoords="offset points")
ax.axhline(5.0, ls="--", color="gray", lw=1)
ax.axhline(3.0, ls=":",  color="gray", lw=1)
ax.axvline(4000, ls="--", color="gray", lw=1)
ax.text(4100, 9.7, "n_points = 4000\n(Phase 2 fail cutoff)", color="gray", fontsize=8)
ax.set_xlabel("# points in submap", fontsize=11)
ax.set_ylabel("rotation error (deg)", fontsize=11)
ax.set_title("n_points is only a weak proxy (see 08 vs 10)", fontsize=11.5)
ax.set_ylim(0, 10.5)
ax.grid(True, alpha=0.3)

# --- Panel 3: n_frames (x) vs n_points (y), color = rot status, size = inliers ---
ax = axes[2]
for (sid, nf, npts, nin, rot, _t) in DATA:
    c = status_color(rot)
    marker = "X" if sid in SC_FAIL_IDS else "o"
    edge = "black" if sid in SC_FAIL_IDS else c
    # size scales with inliers (15..600)
    s = max(40, min(600, nin * 7))
    ax.scatter(nf, npts, c=c, marker=marker, s=s,
               edgecolors=edge, linewidths=1.5, alpha=0.75, zorder=3)
    ax.annotate(sid, (nf, npts), fontsize=8.5, xytext=(5, 4), textcoords="offset points")
ax.set_xlabel("# frames accumulated in submap", fontsize=11)
ax.set_ylabel("# points in submap", fontsize=11)
ax.set_title("Frames vs points (marker size ∝ inliers)", fontsize=11.5)
ax.grid(True, alpha=0.3)
# Diagonals of constant pts/frame; place labels INSIDE the visible window only.
import numpy as np
ax.set_xlim(20, 135)
ax.set_ylim(2000, 8400)
x_lo, x_hi = ax.get_xlim()
y_lo, y_hi = ax.get_ylim()
xs = np.linspace(x_lo, x_hi, 50)
for pf in [50, 100, 200]:
    ys = xs * pf
    ax.plot(xs, ys, ls=":", color="lightblue", alpha=0.7, zorder=1)
    # find a point on the line that's inside the window for the label
    x_in_y = y_hi / pf       # x where line exits the top
    if x_in_y <= x_hi:
        ax.text(x_in_y - 1.5, y_hi - 150, f"{pf} pts/frame",
                color="steelblue", fontsize=7.5, ha="right", va="top",
                rotation=0, clip_on=True)
    else:
        # line stays in window for all x; put label near right edge
        ax.text(x_hi - 1, x_hi * pf + 50, f"{pf} pts/frame",
                color="steelblue", fontsize=7.5, ha="right", va="bottom",
                clip_on=True)

# Global legend — placed below the suptitle, above the panels
legend_handles = [
    mpatches.Patch(color="#2ca02c", label="rot < 3°  (pass strict)"),
    mpatches.Patch(color="#ff7f0e", label="3° ≤ rot < 5°  (warning)"),
    mpatches.Patch(color="#d62728", label="rot ≥ 5°  (fail)"),
    plt.Line2D([0], [0], marker="X", color="white", markerfacecolor="gray",
               markeredgecolor="black", markersize=12, lw=0,
               label="Phase 2 SC NO_SOLUTION case (07/10/17)"),
]
fig.legend(handles=legend_handles, loc="upper center",
           ncol=4, frameon=False, fontsize=10, bbox_to_anchor=(0.5, 0.92))

fig.suptitle("Phase 5: BBS + TEASER on 22 Airy submaps  (map_w2/20260511_110448, v2 in-tree build)",
             fontsize=13, y=0.97)

out = Path(__file__).parent / "bbs_22submaps_diagnostic.png"
fig.savefig(out, dpi=150)
print(f"wrote {out}  ({out.stat().st_size / 1024:.1f} KB)")
