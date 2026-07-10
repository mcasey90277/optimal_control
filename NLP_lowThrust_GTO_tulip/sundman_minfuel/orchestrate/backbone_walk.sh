#!/bin/zsh
# backbone_walk.sh — energy-backbone continuation walk across t_f factors.
#
# Usage:   ./backbone_walk.sh <seed_factor> <factor1> [factor2 ...]
# Example: ./backbone_walk.sh 1.15 1.20 1.25 1.30      (walk up from 1.15x)
#          ./backbone_walk.sh 1.15 1.13 1.11 1.09      (walk down)
#
# Chains energy_step: each step seeds from the last SUCCESSFULLY produced
# energy file (a failed/timed-out step is skipped and the walk continues from
# the previous anchor — validated 2026-07-09 when 1.80x timed out and
# 1.75->1.85 succeeded directly). Each step runs in its own MATLAB process
# with a watchdog: the sporadic CasADi/IPOPT MEX FATAL crash (~1 in 10
# solves) kills only that step, never the walk.
#
# The convex eps=1 energy problem is the ONLY primitive that continues
# crash-free across t_f (bang-bang continuation basin-drifts and MEX-crashes)
# — see LOW_THRUST_MINFUEL_CAMPAIGN.md "Down-sweep CRACKED".
#
# ZSH GOTCHA (cost 13 result files on 2026-07-09): never write
#   local f=$1 out=$DIR/x_$f.mat
# on ONE line — every word expands BEFORE the assignments execute, so $f is
# empty in $out. Keep each local on its own line.

set -u
DIR=$(cd "$(dirname "$0")/.." && pwd)          # sundman_minfuel library root
MAT=${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}
WATCHDOG=${WATCHDOG_S:-600}                    # per-step kill timeout [s]
LOGDIR=$DIR/results/logs
mkdir -p "$LOGDIR"
LOG=$LOGDIR/backbone_walk_$(date +%Y%m%d_%H%M%S).log

[ $# -lt 2 ] && { echo "usage: $0 <seed_factor> <factor1> [factor2 ...]"; exit 1; }
seedf=$1; shift
# canonical energy layout: results/energy/energy_f<milli>.mat (2026-07-09 migration)
efile() { printf '%s/results/energy/energy_f%04.0f.mat' "$DIR" $(($1 * 1000)); }
mkdir -p "$DIR/results/energy"
seed=$(efile $seedf)
[ ! -f "$seed" ] && { echo "seed $seed not found" | tee -a "$LOG"; exit 1; }

for f in "$@"; do
  out=$(efile $f)
  before=0; [ -f "$out" ] && before=$(stat -f %m "$out")
  echo "=== step -> ${f}x  (seed $(basename "$seed")) ===" >> "$LOG"
  "$MAT" -batch "cd('$DIR'); energy_step('$seed',$f,'$out')" >> "$LOG" 2>&1 &
  pid=$!
  waited=0
  while kill -0 $pid 2>/dev/null; do
    sleep 10; waited=$((waited+10))
    if [ $waited -ge $WATCHDOG ]; then
      echo "  TIMEOUT ${f}x -> kill" >> "$LOG"; kill -9 $pid 2>/dev/null; break
    fi
  done
  wait $pid 2>/dev/null
  after=0; [ -f "$out" ] && after=$(stat -f %m "$out")
  if [ "$after" -gt "$before" ]; then
    seed=$out; echo "  OK -> seed now ${f}x" >> "$LOG"
  else
    echo "  FAIL ${f}x (no fresh output) -> keep seed $(basename "$seed")" >> "$LOG"
  fi
done
echo "BACKBONE_WALK_DONE" >> "$LOG"
echo "done; log: $LOG"
