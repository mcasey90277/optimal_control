#!/bin/bash
# RUN_DRO_ELFO_BATCHED  Unattended DRO -> ELFO / DRO -> tulip torus fleet in
# hang-proof batches.
#
# Copy of GTO_DRO/catalog/run_gto_dro_batched.sh (see GTO_tulip/catalog/
# run_gto_batched.sh for the rationale: solver calls are uninterruptible
# in-process, so batches run under an OS timeout; resume + attempt counters
# make restarts free). The front door run_dro_elfo_torus(sel, opts) takes a
# SHEET SELECTOR first, forwarded here as a MATLAB-quoted index string.
#
# Usage:   ./run_dro_elfo_batched.sh [sel] [cells_per_batch] [batch_seconds] [rungs] [tf0]
#   sel    index string into the front door's opts.sheets, e.g. "[1 3:5]";
#          empty or omitted = every sheet (default: the 12 ELFO sheets,
#          droTau 1, raanDeg 0:30:330)
#   rungs  MATLAB vector, e.g. "[15 12 10 7 5]"; empty = the front door's
#          default [15 12 10 7 5 3 2 1.5 1 0.75 0.5]; its rung guard refuses
#          any sheet on disk with other rungs
#   tf0    pass A top-rung seed (ND); empty = the front door's 0.3
# Environment:
#   OUTDIR       sheet folder [DRO_ELFO/results/catalog]; the driver lock,
#                the MATLAB log and combined_progress.txt (with the
#                worker's completion marker) live there too
#   KILL_MARGIN  OS-kill slack (s) beyond batch_seconds [10800]
#   EXTRA_OPTS   MATLAB expression giving a struct whose fields override the
#                front door options, e.g. a 6 x 6 tulip sheet:
#                EXTRA_OPTS="struct('nD',6,'nA',6,'sheets',struct('target','tulip','droTau',1,'raanDeg',NaN))"
#                (trusted input: it is pasted into the MATLAB command) [struct()]
#   MATLAB_BIN   MATLAB executable [/Applications/MATLAB_R2026a.app/bin/matlab]
#                (a stub here lets the driver loop be smoke-tested)
# Pilot (sheet 1 = raanDeg 0):       ./run_dro_elfo_batched.sh 1 6 1800
# Fleet:                             nohup ./run_dro_elfo_batched.sh "" 6 1800 > /dev/null 2>&1 &
#
# COMPLETION: the driver passes its own run id (opts.runId = <pid>-<epoch>);
# the front door logs "SHEETS SELECTED [run <id>]: <tags>" each pass and,
# when no selected sheet has an open cell, "... SHEETS COMPLETE [run <id>]:
# <tags>". The driver stops only on a marker of ITS run whose tag list
# equals ITS selection, so a foreign marker in a shared folder is ignored.
#
# NO-PROGRESS RULE (kept as in the GTO model): stop after 3 passes in which
# the per-sheet progress logs did not grow. NOTE: the engine writes a header
# line ("THRUST LADDER: ...") on every call, so those logs grow on every pass
# that reaches the engine; in practice it is the attempt counters (every
# cell gets at most 2 attempts, counted before solving, so even a pass killed
# on a hang retires its cell) that end the loop through the completion
# marker. The rule only catches passes that never reach the engine.
#
# SIGNALS: TERM, HUP and INT (and any exit) stop the MATLAB child (TERM, then
# KILL after a 10 s grace) and the watchdog, then remove the lock, so a
# killed driver never leaves a MATLAB writing the sheets unlocked. The lock
# directory holds a pid file ("<pid> <start time>"); a second driver prints
# that pid, whether it is alive, and how to clear a stale lock.

set -u
SEL=${1:-}
CELLS=${2:-6}
BSEC=${3:-1800}
RUNGS=${4:-}
TF0=${5:-}
EXTRA_OPTS=${EXTRA_OPTS:-"struct()"}

# Argument screens: these strings are pasted into a MATLAB command line.
if ! [[ "$SEL" =~ ^[][0-9\ :,\;]*$ ]]; then echo "sel '$SEL' is not an index list"; exit 2; fi
if ! [[ "$CELLS" =~ ^[0-9]+$ ]]; then echo "cells_per_batch '$CELLS' is not an integer"; exit 2; fi
if ! [[ "$BSEC" =~ ^[0-9]+$ ]]; then echo "batch_seconds '$BSEC' is not an integer"; exit 2; fi
if ! [[ "$RUNGS" =~ ^[][0-9.eE\ ,+-]*$ ]]; then echo "rungs '$RUNGS' is not a numeric vector"; exit 2; fi
if ! [[ "$TF0" =~ ^[0-9.eE+-]*$ ]]; then echo "tf0 '$TF0' is not a number"; exit 2; fi

OPTS_SET="o.maxCells = $CELLS; o.batchSec = $BSEC;"
if [ -n "$RUNGS" ]; then OPTS_SET="$OPTS_SET o.rungs = $RUNGS;"; fi
if [ -n "$TF0" ]; then OPTS_SET="$OPTS_SET o.tf0 = $TF0;"; fi

# KILL_MARGIN (s): OS-kill slack beyond the clean-exit budget BSEC. The engine
# checks batchSec only BETWEEN cells and counts the attempt before solving, so
# a cell started just before the budget runs on to completion; a full healthy
# cell must be allowed to finish, or the kill lands on a healthy solve and
# burns one of the cell's two attempts (and can lose its pass B). This fleet
# walks 11 rungs per cell (GTO_DRO: 5, margin 3600 s); each rung is capped
# at 600 s IPOPT CPU plus a 300 s fenced refinement, so one healthy cell
# started just before the budget can need 11 x 900 = 9,900 s: 10800 s. A
# true hang costs up to BSEC + KILL_MARGIN before the kill. The PILOT should
# measure the slowest healthy cell and tighten this. Override from the
# environment: KILL_MARGIN=7200 ./run_dro_elfo_batched.sh ...
KILL_MARGIN=${KILL_MARGIN:-10800}
if ! [[ "$KILL_MARGIN" =~ ^[0-9]+$ ]]; then echo "KILL_MARGIN '$KILL_MARGIN' is not an integer"; exit 2; fi
KILL_AFTER=$(( BSEC + KILL_MARGIN ))

HERE="$(cd "$(dirname "$0")" && pwd)"
MATLAB=${MATLAB_BIN:-/Applications/MATLAB_R2026a.app/bin/matlab}
OUTDIR=${OUTDIR:-"$HERE/../results/catalog"}
case "$OUTDIR" in *\'*) echo "OUTDIR must not contain a single quote"; exit 2;; esac
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd -P)"
LOG="$OUTDIR/combined_progress.txt"
MLOG="$OUTDIR/dro_elfo_matlab.log"
PUMPKYN=/Users/msc/Desktop/proj7/external/pumpkynPie

RUN_ID="$$-$(date +%s)"
LOCK="$OUTDIR/.dro_elfo_driver.lock"
MPID=""; WPID=""; OWN_LOCK=0

stop_children() {
  if [ -n "$WPID" ]; then pkill -P "$WPID" 2>/dev/null; kill "$WPID" 2>/dev/null; WPID=""; fi
  if [ -n "$MPID" ] && kill -0 "$MPID" 2>/dev/null; then
    kill -TERM "$MPID" 2>/dev/null
    for _ in 1 2 3 4 5 6 7 8 9 10; do kill -0 "$MPID" 2>/dev/null || break; sleep 1; done
    kill -KILL "$MPID" 2>/dev/null
  fi
  MPID=""
}
cleanup() {
  stop_children
  if [ "$OWN_LOCK" = 1 ]; then rm -rf "$LOCK"; OWN_LOCK=0; fi
}
on_signal() {
  echo "=== driver received SIG$1 -- stopping the MATLAB child, releasing the lock" >> "$LOG"
  cleanup
  trap - EXIT
  exit "$2"
}

if ! mkdir "$LOCK" 2>/dev/null; then
  OTHER=$(cat "$LOCK/pid" 2>/dev/null)
  OPID=${OTHER%% *}
  if [ -n "$OPID" ] && kill -0 "$OPID" 2>/dev/null; then STATE="ALIVE"; else STATE="NOT running (stale lock)"; fi
  echo "another dro-elfo driver holds $LOCK: pid ${OPID:-unknown} (started ${OTHER#* }) is $STATE -- exiting"
  echo "if that driver is gone, clear the stale lock with: rm -rf '$LOCK'"
  exit 1
fi
OWN_LOCK=1
echo "$$ $(date '+%Y-%m-%d %H:%M:%S')" > "$LOCK/pid"
trap cleanup EXIT
trap 'on_signal TERM 143' TERM
trap 'on_signal HUP 129' HUP
trap 'on_signal INT 130' INT

count_lines() { echo $(( $(cat "$OUTDIR"/dro_*_progress.txt 2>/dev/null | wc -l) )); }

# $LOG is append-only across EVERY invocation, so a bare grep for the
# completion marker would match one left by a prior run. Snapshot the line
# count before this invocation writes anything and scope every completion
# check to lines strictly AFTER it.
START_LINE=$( { [ -f "$LOG" ] && wc -l < "$LOG"; } 2>/dev/null || echo 0 )
START_LINE=${START_LINE:-0}

echo "=== droelfo run $RUN_ID started $(date) (sel=${SEL:-all}, $CELLS cells/batch, ${BSEC}s clean, ${KILL_AFTER}s kill, rungs=${RUNGS:-default}, tf0=${TF0:-default}, extra=$EXTRA_OPTS)" >> "$LOG"
dry=0
for pass in $(seq 1 600); do
  before=$(count_lines); before=${before:-0}

  "$MATLAB" -batch "here=pwd; cd('$PUMPKYN'); startup(); cd(here); \
                  addpath('$HERE/..'); addpath('$HERE'); \
                  o = struct('outDir', '$OUTDIR', 'runId', '$RUN_ID'); $OPTS_SET \
                  x = $EXTRA_OPTS; fx = fieldnames(x); \
                  for q = 1:numel(fx), o.(fx{q}) = x.(fx{q}); end; \
                  run_dro_elfo_torus('$SEL', o);" \
          >> "$MLOG" 2>&1 &
  MPID=$!
  ( sleep $KILL_AFTER; kill -9 $MPID 2>/dev/null && \
    echo "  [pass $pass timed out -- a solver call hung; skipping ahead]" >> "$LOG" ) &
  WPID=$!
  wait $MPID 2>/dev/null; MSTATUS=$?
  pkill -P $WPID 2>/dev/null; kill $WPID 2>/dev/null
  MPID=""; WPID=""
  # LAST LINE OF DEFENCE against a warning flood (GTO -> DRO pilot
  # 2026-11-10: 1.96 GB of ode45 tolerance warnings): above 200 MB keep the
  # first and last 2 MB.
  MSIZE=$(wc -c < "$MLOG" 2>/dev/null | tr -d ' '); MSIZE=${MSIZE:-0}
  if [ "$MSIZE" -gt 209715200 ]; then
    { head -c 2097152 "$MLOG"; \
      printf '\n..... [trimmed %s by run_dro_elfo_batched.sh: %s bytes, kept first and last 2 MB] .....\n' "$(date)" "$MSIZE"; \
      tail -c 2097152 "$MLOG"; } > "$MLOG.trim" && mv "$MLOG.trim" "$MLOG"
    echo "  [pass $pass: MATLAB log was $MSIZE bytes -- trimmed to first/last 2 MB]" >> "$LOG"
  fi
  if [ $MSTATUS -ne 0 ] && [ $MSTATUS -ne 137 ]; then
    echo "=== MATLAB exited with status $MSTATUS in pass $pass -- ABORTING (see dro_elfo_matlab.log)" >> "$LOG"
    exit 1
  fi
  if tail -30 "$MLOG" 2>/dev/null | grep -qa "^Error in\|^{Error\|Unrecognized"; then
    echo "=== MATLAB ERROR in pass $pass -- ABORTING (see dro_elfo_matlab.log)" >> "$LOG"
    exit 1
  fi

  after=$(count_lines); after=${after:-0}
  echo "  [pass $pass: $after log lines total]" >> "$LOG"
  RECENT=$(tail -n +"$((START_LINE+1))" "$LOG" 2>/dev/null)
  SELTAGS=$(printf '%s\n' "$RECENT" | grep -aF "SHEETS SELECTED [run $RUN_ID]: " | tail -1 | sed 's/^.*\]: //')
  DONETAGS=$(printf '%s\n' "$RECENT" | grep -aF "SHEETS COMPLETE [run $RUN_ID]: " | tail -1 | sed 's/^.*\]: //')
  if [ -n "$SELTAGS" ] && [ "$SELTAGS" = "$DONETAGS" ]; then
    echo "=== completion marker of run $RUN_ID seen ($DONETAGS) -- stopping" >> "$LOG"; break
  fi
  if [ "$after" -le "$before" ]; then dry=$((dry+1)); else dry=0; fi
  if [ "$dry" -ge 3 ]; then
    echo "=== log did not grow for 3 passes (true wedge) -- stopping" >> "$LOG"
    break
  fi
done
echo "=== droelfo run finished $(date)" >> "$LOG"
