#!/bin/zsh
# psr_batch.sh — CRASH-ROBUST PSR pipeline sweep across t_f factors.
#
# Runs the full PSR pipeline (direct energy->fuel solve, PSR refinement, data
# export, IPOPT local-min certificate, optional movie) for many t_f factors, one
# per t_f. Each factor runs in its OWN `matlab -batch` process, so the sporadic
# UNCATCHABLE CasADi/IPOPT MEX FATAL crash (~1 in 10 solves, kills the whole
# MATLAB process -- a try/catch cannot catch it) takes down only THAT factor;
# the sweep keeps going. This is the crash-survivable sibling of run_psr_batch.m
# (which runs every factor in one process and so dies entirely on a MEX crash).
#
# ============================ HOW TO RUN ============================
# Open a terminal and:
#
#   cd /Users/msc/Desktop/optimal_control/orbit_transfer/GTO_tulip/direct/PSR
#   chmod +x psr_batch.sh            # once, to make it executable
#   ./psr_batch.sh <epsMin> <factor1> [factor2 ...]
#   ./psr_batch.sh <epsMin> energy   # run EVERY factor that has an energy seed
#
# Examples:
#   ./psr_batch.sh 0   1.12 1.13 1.15 1.20      # bang-bang for these factors
#   ./psr_batch.sh 0   energy                   # bang-bang, all seeded factors (1.12..1.95)
#   ./psr_batch.sh 0.1 energy                   # smooth eps=0.1, all seeded factors
#
# Optional environment variables (prefix them on the command line):
#   MOVIE=movie|preview|none   per-factor control movie      [default none]
#   RUNVERIFY=1                also append first-order PMP verify (slow) [0]
#   WATCHDOG_S=1800            per-factor kill timeout, seconds [default 1800]
#   MATLAB_BIN=/path/to/matlab                                 [R2025b]
# e.g.:   MOVIE=movie WATCHDOG_S=3600 ./psr_batch.sh 0 energy
#
# Run it in the background so a closed terminal doesn't kill it, and watch:
#   nohup ./psr_batch.sh 0 energy >/dev/null 2>&1 &
#   tail -f results/logs/psr_batch_*.log        # newest log
#
# RESUMABLE: psr_run_one skips any stage whose output already exists, so if the
# sweep is interrupted just re-run the SAME command -- finished factors are
# skipped instantly. Per-factor result rows are written to
# PSR_data/psr_result_f####_minEps#.mat; a summary table is printed and saved to
# PSR_data/psr_batch_summary_minEps#_<insertionLabel>.mat at the end (also
# rebuildable any time with:  matlab -batch "cd('PSR'); psr_collect_summary(<epsMin>)").
# ===================================================================

set -u
DIR=$(cd "$(dirname "$0")" && pwd)                       # the PSR folder
MAT=${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}
WATCHDOG=${WATCHDOG_S:-1800}                             # per-factor kill timeout [s]
MOVIE=${MOVIE:-none}
RUNVERIFY=${RUNVERIFY:-0}
ENERGYDIR=$DIR/../sundman_minfuel/results/energy
LOGDIR=$DIR/results/logs
mkdir -p "$LOGDIR"
LOG=$LOGDIR/psr_batch_$(date +%Y%m%d_%H%M%S).log

[ $# -lt 2 ] && { echo "usage: $0 <epsMin> <factor1> [factor2 ...] | <epsMin> energy"; exit 1; }
epsMin=$1
shift

# --- resolve the factor list ------------------------------------------------
if [ "$1" = "energy" ]; then
  factors=()
  for f in "$ENERGYDIR"/energy_f*.mat; do
    [ -e "$f" ] || continue
    b=$(basename "$f")
    n=${b#energy_f}
    n=${n%.mat}                                          # milli-factor, e.g. 1150
    factors+=( $(awk "BEGIN{printf \"%.3f\", $n/1000}") )
  done
  factors=( ${(on)factors} )                             # numeric sort (zsh)
else
  factors=( "$@" )
fi
[ ${#factors[@]} -eq 0 ] && { echo "no factors to run"; exit 1; }

echo "=== PSR BATCH  epsMin=$epsMin  movie=$MOVIE  runVerify=$RUNVERIFY  watchdog=${WATCHDOG}s ===" | tee -a "$LOG"
echo "factors (${#factors[@]}): ${factors[*]}" | tee -a "$LOG"

# --- sweep: one matlab process per factor, watchdog, continue on crash ------
opts="struct('epsMin',$epsMin,'movieMode','$MOVIE','runVerify',logical($RUNVERIFY))"
for f in "${factors[@]}"; do
  echo "" | tee -a "$LOG"
  echo "=== factor $f  ($(date +%H:%M:%S)) ===" | tee -a "$LOG"
  "$MAT" -batch "cd('$DIR'); psr_run_one($f, $opts);" >> "$LOG" 2>&1 &
  pid=$!
  waited=0
  while kill -0 $pid 2>/dev/null; do
    sleep 10
    waited=$((waited+10))
    if [ $waited -ge $WATCHDOG ]; then
      echo "  TIMEOUT factor $f after ${WATCHDOG}s -> kill; continuing" | tee -a "$LOG"
      kill -9 $pid 2>/dev/null
      break
    fi
  done
  wait $pid 2>/dev/null
  rc=$?
  # rc != 0 (incl. a MEX FATAL / segfault) just means this factor died -- the
  # loop moves to the next one. psr_run_one wrote a result row iff it finished.
  if [ $rc -eq 0 ]; then
    echo "  factor $f: OK" | tee -a "$LOG"
  else
    echo "  factor $f: process exited $rc (crash/kill) -- continuing" | tee -a "$LOG"
  fi
done

# --- summary from the per-factor result rows --------------------------------
echo "" | tee -a "$LOG"
echo "=== building summary ===" | tee -a "$LOG"
"$MAT" -batch "cd('$DIR'); psr_collect_summary($epsMin);" 2>&1 | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "log: $LOG"
