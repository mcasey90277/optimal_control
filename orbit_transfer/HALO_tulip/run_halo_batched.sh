#!/bin/bash
# RUN_HALO_BATCHED  Unattended halo->tulip catalog runs in hang-proof batches.
#
# Same pattern as DRO_tulip/run_ladder_batched.sh (see that file for the
# rationale: solver calls are uninterruptible in-process, so batches run
# under an OS timeout; resume + attempt counters make restarts free).
#
# Usage:   ./run_halo_batched.sh [cells_per_batch] [batch_seconds]
# Example: ./run_halo_batched.sh 8 240

set -u
CELLS=${1:-8}
BSEC=${2:-240}
KILL_AFTER=$(( BSEC + 140 ))

HERE="$(cd "$(dirname "$0")" && pwd)"
MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab
CATDIR="$HERE/direct/results/catalog"
LOG="$CATDIR/combined_progress.txt"
MLOG="$HERE/direct/results/halocat_matlab.log"
PUMPKYN=/Users/msc/Desktop/proj7/external/pumpkynPie

mkdir -p "$CATDIR"
LOCK="$HERE/direct/results/.halocat_driver.lock"
if ! mkdir "$LOCK" 2>/dev/null; then
  echo "another halo driver appears to be running (lock: $LOCK) -- exiting"
  exit 1
fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

count_lines() { cat "$CATDIR"/halo_*_progress.txt 2>/dev/null | wc -l; }

echo "=== halocat run started $(date) ($CELLS cells/batch, ${BSEC}s clean, ${KILL_AFTER}s kill)" >> "$LOG"
dry=0
for pass in $(seq 1 200); do
  before=$(count_lines); before=${before:-0}

  $MATLAB -batch "here=pwd; cd('$PUMPKYN'); startup(); cd(here); \
                  addpath('$HERE'); \
                  run_halo_catalog('all', $CELLS, $BSEC);" \
          >> "$MLOG" 2>&1 &
  MPID=$!
  ( sleep $KILL_AFTER; kill -9 $MPID 2>/dev/null && \
    echo "  [pass $pass timed out -- a solver call hung; skipping ahead]" >> "$LOG" ) &
  WPID=$!
  wait $MPID 2>/dev/null; MSTATUS=$?
  kill $WPID 2>/dev/null
  if [ $MSTATUS -ne 0 ] && [ $MSTATUS -ne 137 ]; then
    echo "=== MATLAB exited with status $MSTATUS in pass $pass -- ABORTING" >> "$LOG"
    break
  fi
  if tail -30 "$MLOG" 2>/dev/null | grep -qa "^Error in\|^{Error\|Unrecognized"; then
    echo "=== MATLAB ERROR in pass $pass -- ABORTING (see halocat_matlab.log)" >> "$LOG"
    break
  fi

  after=$(count_lines); after=${after:-0}
  echo "  [pass $pass: $after log lines total]" >> "$LOG"
  if grep -qa "ALL SHEETS COMPLETE" "$LOG" 2>/dev/null; then
    echo "=== completion marker seen -- stopping" >> "$LOG"; break
  fi
  if [ "$after" -le "$before" ]; then dry=$((dry+1)); else dry=0; fi
  if [ "$dry" -ge 3 ]; then
    echo "=== log did not grow for 3 passes (true wedge) -- stopping" >> "$LOG"
    break
  fi
done
echo "=== halocat run finished $(date)" >> "$LOG"
