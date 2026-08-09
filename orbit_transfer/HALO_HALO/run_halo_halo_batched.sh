#!/bin/bash
# RUN_HALO_HALO_BATCHED  Unattended L1<->L2 halo catalog runs, hang-proof.
#
# Same pattern as DPO_tulip/run_dpo_batched.sh (see DRO_tulip/
# run_ladder_batched.sh for the rationale: solver calls are uninterruptible
# in-process, so batches run under an OS timeout; resume + attempt counters
# make restarts free).
#
# Usage:   ./run_halo_halo_batched.sh [stage] [cells_per_batch] [batch_seconds]
# Example: ./run_halo_halo_batched.sh all 8 240
#          ./run_halo_halo_batched.sh allB 8 240     (the L2->L1 direction)

set -u
STAGE=${1:-all}
CELLS=${2:-8}
BSEC=${3:-240}
KILL_AFTER=$(( BSEC + 140 ))

HERE="$(cd "$(dirname "$0")" && pwd)"
MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab
CATDIR="$HERE/direct/results/catalog"
LOG="$CATDIR/combined_progress.txt"
MLOG="$HERE/direct/results/hhcat_matlab_${STAGE}.log"
PUMPKYN=/Users/msc/Desktop/proj7/external/pumpkynPie

if [ "$STAGE" = "allB" ]; then MARKER="CATALOG B (L2->L1) ALL SHEETS COMPLETE"; else MARKER="CATALOG ALL SHEETS COMPLETE"; fi

mkdir -p "$CATDIR"
LOCK="$HERE/direct/results/.hhcat_driver_${STAGE}.lock"
if ! mkdir "$LOCK" 2>/dev/null; then
  echo "another halo-halo driver appears to be running (lock: $LOCK) -- exiting"
  exit 1
fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

if [ "$STAGE" = "allB" ]; then PFX="hhB_d"; else PFX="hh_d"; fi
count_lines() { echo $(( $(cat "$CATDIR"/${PFX}*_progress.txt 2>/dev/null | wc -l) )); }

echo "=== hhcat run started $(date) (stage $STAGE, $CELLS cells/batch, ${BSEC}s clean, ${KILL_AFTER}s kill)" >> "$LOG"
dry=0
for pass in $(seq 1 600); do
  before=$(count_lines); before=${before:-0}

  "$MATLAB" -batch "here=pwd; cd('$PUMPKYN'); startup(); cd(here); \
                  addpath('$HERE'); \
                  run_halo_halo_catalog('$STAGE', $CELLS, $BSEC);" \
          >> "$MLOG" 2>&1 &
  MPID=$!
  ( sleep $KILL_AFTER; kill -9 $MPID 2>/dev/null && \
    echo "  [pass $pass timed out -- a solver call hung; skipping ahead]" >> "$LOG" ) &
  WPID=$!
  wait $MPID 2>/dev/null; MSTATUS=$?
  pkill -P $WPID 2>/dev/null; kill $WPID 2>/dev/null
  if [ $MSTATUS -ne 0 ] && [ $MSTATUS -ne 137 ]; then
    echo "=== MATLAB exited with status $MSTATUS in pass $pass -- ABORTING" >> "$LOG"
    exit 1
  fi
  if tail -30 "$MLOG" 2>/dev/null | grep -qa "^Error in\|^{Error\|Unrecognized"; then
    echo "=== MATLAB ERROR in pass $pass -- ABORTING (see hhcat_matlab.log)" >> "$LOG"
    exit 1
  fi

  after=$(count_lines); after=${after:-0}
  echo "  [pass $pass: $after log lines total]" >> "$LOG"
  if grep -qa "$MARKER" "$LOG" 2>/dev/null; then
    echo "=== completion marker seen -- stopping" >> "$LOG"; break
  fi
  if [ "$after" -le "$before" ]; then dry=$((dry+1)); else dry=0; fi
  if [ "$dry" -ge 3 ]; then
    echo "=== log did not grow for 3 passes (true wedge) -- stopping" >> "$LOG"
    break
  fi
done
echo "=== hhcat run finished $(date)" >> "$LOG"
