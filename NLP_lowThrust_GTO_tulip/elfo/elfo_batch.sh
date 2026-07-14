#!/bin/zsh
# elfo_batch.sh — CRASH-ROBUST GTO->ELFO min-fuel sweep across t_f factors.
#
# Runs the energy->fuel homotopy (gen_elfo_minfuel, eps 1 -> epsMin) for many t_f
# factors, ONE per factor, each in its OWN `matlab -batch` process so the sporadic
# UNCATCHABLE CasADi/IPOPT MEX FATAL crash (~1 in 10 solves; a try/catch cannot
# catch it) takes down only THAT factor -- the sweep continues. The ELFO sibling
# of PSR/psr_batch.sh.
#
# ============================ HOW TO RUN ============================
#   cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo
#   chmod +x elfo_batch.sh                        # once
#   ./elfo_batch.sh <epsMin> <factor1> [factor2 ...]
#   ./elfo_batch.sh <epsMin> energy               # every factor with an energy seed
#
# Examples:
#   ./elfo_batch.sh 0   1.15 1.20 1.25            # bang-bang fuel for these factors
#   ./elfo_batch.sh 0   energy                    # bang-bang, all seeded factors
#   ./elfo_batch.sh 0.5 energy                    # smooth eps=0.5, all seeded factors
#
# Optional env vars (prefix on the command line):
#   WATCHDOG_S=1800   per-factor kill timeout, seconds        [default 1800]
#   MATLAB_BIN=/path/to/matlab                                 [R2025b]
#
# Background + watch:
#   nohup ./elfo_batch.sh 0 energy >/dev/null 2>&1 &
#   tail -f results/logs/elfo_batch_*.log
#
# RESUMABLE: elfo_run_one skips any factor whose result row already exists, so
# re-running the SAME command finishes only the missing factors. A summary is
# printed and saved to results/elfo_batch_summary_minEps#.mat (rebuildable with
#   matlab -batch "cd('elfo'); elfo_collect_summary(<epsMin>)").
# ===================================================================

set -u
DIR=$(cd "$(dirname "$0")" && pwd)                        # the elfo folder
MAT=${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}
WATCHDOG=${WATCHDOG_S:-1800}
SEEDDIR=$DIR/results
LOGDIR=$DIR/results/logs
mkdir -p "$LOGDIR"
LOG=$LOGDIR/elfo_batch_$(date +%Y%m%d_%H%M%S).log

[ $# -lt 2 ] && { echo "usage: $0 <epsMin> <factor1> [factor2 ...] | <epsMin> energy"; exit 1; }
epsMin=$1
shift

# --- resolve the factor list ------------------------------------------------
if [ "$1" = "energy" ]; then
  factors=()
  for f in "$SEEDDIR"/energy_elfo_f[0-9]*.mat(N); do
    [ -e "$f" ] || continue
    b=$(basename "$f")
    n=${b#energy_elfo_f}
    n=${n%.mat}                                           # milli-factor, e.g. 1200
    factors+=( $(awk "BEGIN{printf \"%.3f\", $n/1000}") )
  done
  factors=( ${(on)factors} )                              # numeric sort (zsh)
else
  factors=( "$@" )
fi
[ ${#factors[@]} -eq 0 ] && { echo "no factors to run"; exit 1; }

echo "=== ELFO BATCH  epsMin=$epsMin  watchdog=${WATCHDOG}s ===" | tee -a "$LOG"
echo "factors (${#factors[@]}): ${factors[*]}" | tee -a "$LOG"

# --- sweep: one matlab process per factor, watchdog, continue on crash ------
opts="struct('epsMin',$epsMin)"
for f in "${factors[@]}"; do
  echo "" | tee -a "$LOG"
  echo "=== factor $f  ($(date +%H:%M:%S)) ===" | tee -a "$LOG"
  "$MAT" -batch "cd('$DIR'); elfo_run_one($f, $opts);" >> "$LOG" 2>&1 &
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
  if [ $rc -eq 0 ]; then
    echo "  factor $f: OK" | tee -a "$LOG"
  else
    echo "  factor $f: process exited $rc (crash/kill) -- continuing" | tee -a "$LOG"
  fi
done

# --- summary from the per-factor result rows --------------------------------
echo "" | tee -a "$LOG"
echo "=== building summary ===" | tee -a "$LOG"
"$MAT" -batch "cd('$DIR'); elfo_collect_summary($epsMin);" 2>&1 | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "log: $LOG"
