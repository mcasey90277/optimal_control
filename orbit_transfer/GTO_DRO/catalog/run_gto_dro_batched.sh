#!/bin/bash
# RUN_GTO_DRO_BATCHED  Unattended GTO->DRO catalog runs in hang-proof batches.
#
# Clone of GTO_tulip/catalog/run_gto_batched.sh (see it for the rationale:
# solver calls are uninterruptible in-process, so batches run under an OS
# timeout; resume + attempt counters make restarts free). The front door
# run_gto_dro_catalog takes a SHEET SELECTOR as its first argument, so this
# script forwards a MATLAB index-vector expression instead of a stage word.
#
# Usage:   ./run_gto_dro_batched.sh [sheetSel] [cells_per_batch] [batch_seconds] [rungs] [tf0]
# The front door's directGateKm default (1000 km, direct gate only) is used.
# Pilot (sheet 2 = orient 0, tau 1):      ./run_gto_dro_batched.sh 2 6 1800
# Deep probe (6 cells, 9 rungs): not made by this driver. planOnly never
#   solves, so run it by hand: a direct thrust_ladder_library call using the
#   plan's optsA (9 rungs, 6 cells) with an overridden, SEPARATE output sheet
#   file (see GTO_DRO/README.md once stage 2 is documented).
# Fleet:                                  ./run_gto_dro_batched.sh "[1 3:16]" 6 1800 "" 0.30
#
# Batch budget: a GTO -> DRO cell at 15 N costs 100 to 520 s per rung
# attempt; 6 cells per 1800 s batch with the OS kill at BSEC + KILL_MARGIN
# (default 3600 s, i.e. 5400 s here).
#
# rungs (4th arg, optional): MATLAB vector expression forwarded as the front
# door's rungsIn. Omitted/empty -> its default [15 12 10 7 5]; the front
# door's rung-mismatch guard refuses any sheet on disk with other rungs.
# tf0 (5th arg, optional): pass A top-rung seed (ND), forwarded as
# struct('tf0', ...); default is the front door's 0.30.

set -u
SHEETSEL=${1:-1:16}
CELLS=${2:-6}
BSEC=${3:-1800}
RUNGS=${4:-}
TF0=${5:-}
RUNGS_ARG=""
TF0_ARG=""
if [ -n "$RUNGS" ]; then RUNGS_ARG=", $RUNGS"; fi
if [ -n "$TF0" ]; then
  TF0_ARG=", struct('tf0', $TF0)"
  if [ -z "$RUNGS" ]; then RUNGS_ARG=", []"; fi
fi
# KILL_MARGIN (s): OS-kill slack beyond the clean-exit budget BSEC. The engine
# checks batchSec only BETWEEN cells and counts the attempt before solving, so
# a cell started just before the budget runs on to completion; a full healthy
# 5-rung DRO cell (each rung capped at 600 s IPOPT CPU, plus refinement) must
# be allowed to finish, or the kill lands on healthy solves, burns one of the
# cell's two attempts and can lose its cold-start pass B. A true hang now
# costs up to BSEC + KILL_MARGIN s before the kill. The tulip driver's 140 s
# (sized for 9-58 s cells) is WRONG for this campaign. Override from the
# environment: KILL_MARGIN=7200 ./run_gto_dro_batched.sh ...
KILL_MARGIN=${KILL_MARGIN:-3600}
KILL_AFTER=$(( BSEC + KILL_MARGIN ))

HERE="$(cd "$(dirname "$0")" && pwd)"
MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
CATDIR="$HERE/../results/catalog"
LOG="$CATDIR/combined_progress.txt"
MLOG="$HERE/../results/gtodrocat_matlab.log"
PUMPKYN=/Users/msc/Desktop/proj7/external/pumpkynPie

mkdir -p "$CATDIR"
LOCK="$HERE/../results/.gtodrocat_driver.lock"
if ! mkdir "$LOCK" 2>/dev/null; then
  echo "another gto-dro driver appears to be running (lock: $LOCK) -- exiting"
  exit 1
fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

count_lines() { echo $(( $(cat "$CATDIR"/gto_dro_*_progress.txt 2>/dev/null | wc -l) )); }

# $LOG (combined_progress.txt) is append-only across EVERY invocation, so a
# bare grep for a completion marker would match one left by a prior run.
# Snapshot the line count before this invocation writes anything and scope
# every completion check to lines strictly AFTER it.
START_LINE=$( { [ -f "$LOG" ] && wc -l < "$LOG"; } 2>/dev/null || echo 0 )
START_LINE=${START_LINE:-0}

echo "=== gtodrocat run started $(date) (sheetSel=$SHEETSEL, $CELLS cells/batch, ${BSEC}s clean, ${KILL_AFTER}s kill, rungs=${RUNGS:-default}, tf0=${TF0:-default})" >> "$LOG"
dry=0
for pass in $(seq 1 600); do
  before=$(count_lines); before=${before:-0}

  "$MATLAB" -batch "here=pwd; cd('$PUMPKYN'); startup(); cd(here); \
                  addpath('$HERE'); addpath('$HERE/..'); \
                  run_gto_dro_catalog($SHEETSEL, $CELLS, $BSEC$RUNGS_ARG$TF0_ARG);" \
          >> "$MLOG" 2>&1 &
  MPID=$!
  ( sleep $KILL_AFTER; kill -9 $MPID 2>/dev/null && \
    echo "  [pass $pass timed out -- a solver call hung; skipping ahead]" >> "$LOG" ) &
  WPID=$!
  wait $MPID 2>/dev/null; MSTATUS=$?
  pkill -P $WPID 2>/dev/null; kill $WPID 2>/dev/null
  # LAST LINE OF DEFENCE against a warning flood (pilot 2026-11-10: 1.96 GB
  # of ode45 tolerance warnings): above 200 MB keep the first and last 2 MB.
  MSIZE=$(wc -c < "$MLOG" 2>/dev/null | tr -d ' '); MSIZE=${MSIZE:-0}
  if [ "$MSIZE" -gt 209715200 ]; then
    { head -c 2097152 "$MLOG"; \
      printf '\n..... [trimmed %s by run_gto_dro_batched.sh: %s bytes, kept first and last 2 MB] .....\n' "$(date)" "$MSIZE"; \
      tail -c 2097152 "$MLOG"; } > "$MLOG.trim" && mv "$MLOG.trim" "$MLOG"
    echo "  [pass $pass: MATLAB log was $MSIZE bytes -- trimmed to first/last 2 MB]" >> "$LOG"
  fi
  if [ $MSTATUS -ne 0 ] && [ $MSTATUS -ne 137 ]; then
    echo "=== MATLAB exited with status $MSTATUS in pass $pass -- ABORTING" >> "$LOG"
    exit 1
  fi
  if tail -30 "$MLOG" 2>/dev/null | grep -qa "^Error in\|^{Error\|Unrecognized"; then
    echo "=== MATLAB ERROR in pass $pass -- ABORTING (see gtodrocat_matlab.log)" >> "$LOG"
    exit 1
  fi

  after=$(count_lines); after=${after:-0}
  echo "  [pass $pass: $after log lines total]" >> "$LOG"
  if tail -n +"$((START_LINE+1))" "$LOG" 2>/dev/null | grep -qa "SHEETS COMPLETE"; then
    echo "=== completion marker seen -- stopping" >> "$LOG"; break
  fi
  if [ "$after" -le "$before" ]; then dry=$((dry+1)); else dry=0; fi
  if [ "$dry" -ge 3 ]; then
    echo "=== log did not grow for 3 passes (true wedge) -- stopping" >> "$LOG"
    break
  fi
done
echo "=== gtodrocat run finished $(date)" >> "$LOG"
