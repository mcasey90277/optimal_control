#!/bin/bash
# RUN_GTO_BATCHED  Unattended GTO->tulip catalog runs in hang-proof batches.
#
# Same pattern as DRO_tulip/run_ladder_batched.sh and
# HALO_tulip/run_halo_batched.sh (see those for the rationale: solver calls
# are uninterruptible in-process, so batches run under an OS timeout;
# resume + attempt counters make restarts free). Deviation from a literal
# clone: run_gto_catalog's front door takes a SHEET SELECTOR as its first
# argument (run_gto_catalog(sheetSel, maxCells, batchSec)), not a
# 'pilot'/'all'/'report' stage string like the halo driver, so this script
# forwards a MATLAB index-vector expression instead of a stage word.
#
# Usage:   ./run_gto_batched.sh [sheetSel] [cells_per_batch] [batch_seconds]
# Example (pilot, orient=0 x Np=7 == sheet 2 of 16):
#          ./run_gto_batched.sh 2 8 240
# Example (whole catalog):
#          ./run_gto_batched.sh 1:16 8 240

set -u
SHEETSEL=${1:-1:16}
CELLS=${2:-8}
BSEC=${3:-240}
KILL_AFTER=$(( BSEC + 140 ))

HERE="$(cd "$(dirname "$0")" && pwd)"
MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab
CATDIR="$HERE/results/catalog"
LOG="$CATDIR/combined_progress.txt"
MLOG="$HERE/results/gtocat_matlab.log"
PUMPKYN=/Users/msc/Desktop/proj7/external/pumpkynPie

mkdir -p "$CATDIR"
LOCK="$HERE/results/.gtocat_driver.lock"
if ! mkdir "$LOCK" 2>/dev/null; then
  echo "another gto driver appears to be running (lock: $LOCK) -- exiting"
  exit 1
fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

count_lines() { echo $(( $(cat "$CATDIR"/gto_*_progress.txt 2>/dev/null | wc -l) )); }

echo "=== gtocat run started $(date) (sheetSel=$SHEETSEL, $CELLS cells/batch, ${BSEC}s clean, ${KILL_AFTER}s kill)" >> "$LOG"
dry=0
for pass in $(seq 1 600); do
  before=$(count_lines); before=${before:-0}

  "$MATLAB" -batch "here=pwd; cd('$PUMPKYN'); startup(); cd(here); \
                  addpath('$HERE'); \
                  run_gto_catalog($SHEETSEL, $CELLS, $BSEC);" \
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
    echo "=== MATLAB ERROR in pass $pass -- ABORTING (see gtocat_matlab.log)" >> "$LOG"
    exit 1
  fi

  after=$(count_lines); after=${after:-0}
  echo "  [pass $pass: $after log lines total]" >> "$LOG"
  if grep -qa "SHEETS COMPLETE" "$LOG" 2>/dev/null; then
    echo "=== completion marker seen -- stopping" >> "$LOG"; break
  fi
  if [ "$after" -le "$before" ]; then dry=$((dry+1)); else dry=0; fi
  if [ "$dry" -ge 3 ]; then
    echo "=== log did not grow for 3 passes (true wedge) -- stopping" >> "$LOG"
    break
  fi
done
echo "=== gtocat run finished $(date)" >> "$LOG"
