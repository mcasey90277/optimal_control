#!/bin/bash
# RUN_LADDER_BATCHED  Long unattended runs of the costate-library pipeline.
#
# Why batches: some solver calls inside MATLAB cannot be interrupted from
# within MATLAB (an ODE integration can crawl near the lunar singularity
# without bound). Running the work in small batches under an OS-level timeout
# means a hang costs one batch instead of the campaign; the pipeline's own
# resume logic then skips whatever already finished.
#
# Usage:   ./run_ladder_batched.sh [ladder|extend] [cells_per_batch] [batch_seconds]
# Example: ./run_ladder_batched.sh extend 6 240
#
# Progress is appended to results/<stage>_progress.txt; every completed rung
# is already saved to disk, so the run can be stopped and restarted freely.

set -u
STAGE=${1:-ladder}
CELLS=${2:-8}
BSEC=${3:-180}
KILL_AFTER=$(( BSEC + 140 ))          # clean exit + MATLAB startup + grace

HERE="$(cd "$(dirname "$0")" && pwd)"
MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab
LOG="$HERE/direct/results/${STAGE}_progress.txt"
PUMPKYN=/Users/msc/Desktop/proj7/external/pumpkynPie

mkdir -p "$HERE/direct/results"
echo "=== $STAGE run started $(date) ($CELLS cells/batch, ${BSEC}s clean, ${KILL_AFTER}s kill)" >> "$LOG"

for pass in $(seq 1 200); do
  before=$(grep -ac "OK (" "$LOG" 2>/dev/null); before=${before:-0}

  $MATLAB -batch "here=pwd; cd('$PUMPKYN'); startup(); cd(here); \
                  addpath('$HERE'); \
                  run_costate_library('$STAGE', $CELLS, $BSEC);" \
          >> "$HERE/direct/results/${STAGE}_matlab.log" 2>&1 &
  MPID=$!
  ( sleep $KILL_AFTER; kill -9 $MPID 2>/dev/null && \
    echo "  [pass $pass timed out -- a solver call hung; skipping ahead]" >> "$LOG" ) &
  WPID=$!
  wait $MPID 2>/dev/null
  kill $WPID 2>/dev/null

  after=$(grep -ac "OK (" "$LOG" 2>/dev/null); after=${after:-0}
  echo "  [pass $pass: $after verified entries total]" >> "$LOG"
  if [ "$after" -eq "$before" ] && [ $pass -gt 2 ]; then
    echo "=== no progress in a full pass -- stopping ($after entries)" >> "$LOG"
    break
  fi
done
echo "=== $STAGE run finished $(date)" >> "$LOG"
