#!/usr/bin/env bash
# Phase 5 batch eval: run BBS + TEASER on each of 22 submaps in map_w2/20260511_110448.
# Outputs per-submap result to OUT_DIR, plus a CSV summary.

set -e

SESSION="${SESSION:-/home/steve/map_data/map_w2/20260511_110448}"
BIN="${BIN:-/home/steve/Documents/GitHub/tools/TEASER-plusplus/build/examples/global_localization_bbs/bbs_localize}"
OUT_DIR="${OUT_DIR:-/tmp/bbs_22submaps}"

mkdir -p "$OUT_DIR"
CSV="$OUT_DIR/summary.csv"
echo "submap,n_points,bbs_ms,refine_s,n_corres,n_inliers,inlier_ratio,rmse_m,rot_err_deg,trans_err_m,status" > "$CSV"

for i in $(seq -f "%06g" 0 21); do
  SUBDIR="$SESSION/$i"
  SCAN="$SUBDIR/submap_levelled.pcd"
  GT="$SUBDIR/data.txt"
  OUT="$OUT_DIR/${i}.txt"
  LOG="$OUT_DIR/${i}.log"

  if [[ ! -f "$SCAN" || ! -f "$GT" ]]; then
    echo "$i,,,,,,,,,,MISSING" >> "$CSV"
    continue
  fi

  N_PTS=$(python3 -c "
import struct
with open('$SUBDIR/points_compact.bin','rb') as f:
    print(len(f.read())//12)
")
  echo "=== submap $i  (n_points=$N_PTS) ==="

  if "$BIN" \
      --map "$SESSION/global_map.pcd" \
      --scan "$SCAN" \
      --gt "$GT" \
      --output "$OUT" > "$LOG" 2>&1 ; then
    STATUS="OK"
  else
    STATUS="FAIL_$?"
  fi

  if [[ -f "$OUT" ]]; then
    # quality line is after `# quality:` header
    QLINE=$(awk '/^# quality:/ {getline; print}' "$OUT")
    ELINE=$(awk '/^# error:/   {getline; print}' "$OUT")
    BBS_MS=$(echo "$QLINE" | awk '{print $1}')
    REF_S=$(echo  "$QLINE" | awk '{print $2}')
    NC=$(echo     "$QLINE" | awk '{print $3}')
    NI=$(echo     "$QLINE" | awk '{print $4}')
    IR=$(echo     "$QLINE" | awk '{print $5}')
    RMSE=$(echo   "$QLINE" | awk '{print $6}')
    ROT=$(echo    "$ELINE" | awk '{print $1}')
    TRA=$(echo    "$ELINE" | awk '{print $2}')
    echo "$i,$N_PTS,$BBS_MS,$REF_S,$NC,$NI,$IR,$RMSE,$ROT,$TRA,$STATUS" >> "$CSV"
  else
    echo "$i,$N_PTS,,,,,,,,,$STATUS" >> "$CSV"
  fi
done

echo
echo "=== Summary ==="
column -s, -t < "$CSV"
echo
echo "CSV: $CSV"
