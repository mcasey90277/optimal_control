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
# Usage:   ./run_gto_batched.sh [sheetSel] [cells_per_batch] [batch_seconds] [rungs]
# Example (pilot, orient=0 x Np=7 == sheet 2 of 16):
#          ./run_gto_batched.sh 2 8 240
# Example (whole catalog):
#          ./run_gto_batched.sh 1:16 8 240
# Example (deliverable-7 v1 fleet, 5-N-truncated ladder, task 4 amendment):
#          ./run_gto_batched.sh "[1 3:16]" 8 240 "[15 12 10 7 5]"
#
# rungs (optional, 4th arg): MATLAB vector expression forwarded as
# run_gto_catalog's rungsIn override. Omitted/empty -> run_gto_catalog's own
# default, which is the shipped deliverable-7 v1 fleet's 5-rung ladder
# [15 12 10 7 5] (RE-BASED 2026-08-29 review fix -- this now matches every
# on-disk sheet's Q.rungs; the front door's own rung-mismatch guard throws
# otherwise). The full 9-rung ladder [15 12 10 7 5 3 2 1.5 1] still exists
# but must be requested explicitly via this 4th arg, and only against a
# sheet that itself carries all 9 rungs.

set -u
SHEETSEL=${1:-1:16}
CELLS=${2:-8}
BSEC=${3:-240}
RUNGS=${4:-}
RUNGS_ARG=""
if [ -n "$RUNGS" ]; then RUNGS_ARG=", $RUNGS"; fi
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

# ROOT-CAUSE FIX (fix round 2, 2026-08-27): $LOG (combined_progress.txt) is
# append-only across EVERY invocation ever run -- it is never truncated.
# A bare `grep "SHEETS COMPLETE" "$LOG"` therefore matches a completion
# marker left by ANY prior run (e.g. an earlier pilot-only sheetSel that
# legitimately finished), even though THIS invocation's own sheets are far
# from done. That false positive is exactly what stopped the 15-sheet fleet
# after a single partial pass: a "SELECTED SHEETS COMPLETE (2)" line from
# the earlier pilot run was still sitting in the file, and the very first
# pass's blind grep found it immediately. Fix: snapshot the log's line
# count BEFORE this invocation writes anything, and scope every completion
# check to lines strictly AFTER that snapshot for the rest of this run.
START_LINE=$( { [ -f "$LOG" ] && wc -l < "$LOG"; } 2>/dev/null || echo 0 )
START_LINE=${START_LINE:-0}

echo "=== gtocat run started $(date) (sheetSel=$SHEETSEL, $CELLS cells/batch, ${BSEC}s clean, ${KILL_AFTER}s kill, rungs=${RUNGS:-default})" >> "$LOG"
dry=0
for pass in $(seq 1 600); do
  before=$(count_lines); before=${before:-0}

  "$MATLAB" -batch "here=pwd; cd('$PUMPKYN'); startup(); cd(here); \
                  addpath('$HERE'); \
                  run_gto_catalog($SHEETSEL, $CELLS, $BSEC$RUNGS_ARG);" \
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
  # Scoped to THIS invocation's own appended lines only -- see the
  # START_LINE comment above for why an unscoped grep is unsafe.
  if tail -n +"$((START_LINE+1))" "$LOG" 2>/dev/null | grep -qa "SHEETS COMPLETE"; then
    echo "=== completion marker seen -- stopping" >> "$LOG"; break
  fi
  if [ "$after" -le "$before" ]; then dry=$((dry+1)); else dry=0; fi
  if [ "$dry" -ge 3 ]; then
    echo "=== log did not grow for 3 passes (true wedge) -- stopping" >> "$LOG"
    break
  fi
done
echo "=== gtocat run finished $(date)" >> "$LOG"
