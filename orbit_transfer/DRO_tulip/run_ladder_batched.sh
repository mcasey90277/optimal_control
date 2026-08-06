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
if [ "$STAGE" = "catalog" ]; then
  LOG="$HERE/direct/results/catalog/combined_progress.txt"
  mkdir -p "$HERE/direct/results/catalog"
  ln -sf /dev/null "$LOG" 2>/dev/null; rm -f "$LOG"
  count_ok() { cat "$HERE"/direct/results/catalog/ladder_*_progress.txt 2>/dev/null | wc -l; }
else
  LOG="$HERE/direct/results/${STAGE}_progress.txt"
  count_ok() { wc -l < "$LOG" 2>/dev/null; }
fi
PUMPKYN=/Users/msc/Desktop/proj7/external/pumpkynPie

mkdir -p "$HERE/direct/results"
echo "=== $STAGE run started $(date) ($CELLS cells/batch, ${BSEC}s clean, ${KILL_AFTER}s kill)" >> "$LOG"

dry=0
for pass in $(seq 1 200); do
  before=$(count_ok); before=${before:-0}

  if [ "$STAGE" = "lowthrust" ]; then
    CMD="run_lowthrust_ladder('all', $CELLS, $BSEC);"
  elif [ "$STAGE" = "catalog" ]; then
    CMD="run_catalog_sweep('all', $CELLS, $BSEC);"
  else
    CMD="run_costate_library('$STAGE', $CELLS, $BSEC);"
  fi
  $MATLAB -batch "here=pwd; cd('$PUMPKYN'); startup(); cd(here); \
                  addpath('$HERE'); \
                  $CMD" \
          >> "$HERE/direct/results/${STAGE}_matlab.log" 2>&1 &
  MPID=$!
  ( sleep $KILL_AFTER; kill -9 $MPID 2>/dev/null && \
    echo "  [pass $pass timed out -- a solver call hung; skipping ahead]" >> "$LOG" ) &
  WPID=$!
  wait $MPID 2>/dev/null
  kill $WPID 2>/dev/null

  after=$(count_ok); after=${after:-0}
  echo "  [pass $pass: $after log lines total]" >> "$LOG"
  # explicit completion marker beats inference
  if grep -qa "ALL SHEETS COMPLETE\|ALL WORK COMPLETE" "$LOG" 2>/dev/null; then
    echo "=== completion marker seen -- stopping" >> "$LOG"; break
  fi
  # progress = the WORK LIST advancing (attempts included), not only fills;
  # a hang cluster burns attempts and retires -- that is still progress.
  if [ "$after" -le "$before" ]; then dry=$((dry+1)); else dry=0; fi
  if [ "$dry" -ge 3 ]; then
    echo "=== log did not grow for 3 passes (true wedge) -- stopping" >> "$LOG"
    break
  fi
done
echo "=== $STAGE run finished $(date)" >> "$LOG"
