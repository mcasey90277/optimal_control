#!/bin/zsh
# elfo_energy_sweep.sh — run the GTO->ELFO min-ENERGY tf-band sweep, crash-retried.
#
# Phase 1 of the min-fuel tf-grid campaign: sweeps the PINNED transfer time UP and
# DOWN from the 1.20x base seed, banking a converged energy seed
# results/energy_elfo_f<NNNN>.mat at each grid factor (NNNN=round(1000*factor)),
# plus the band summary results/energy_elfo_tfgrid_<insertionLabel>.mat.
#
# Unlike the fuel batch (one process per factor), the energy sweep is a SINGLE
# warm-started continuation, so it runs in ONE matlab process. To survive the
# sporadic UNCATCHABLE CasADi/IPOPT MEX FATAL crash (~1 in 10 solves), this wrapper
# RE-RUNS the sweep on a non-zero exit; gen_elfo_energy_tfsweep is resumable, so
# each retry picks up from the furthest already-banked seed instead of the base.
# A natural energy-band-edge stop exits 0 (the sweep step-halves and breaks
# cleanly), so retries fire ONLY on crashes, not on a legitimate early stop.
#
# ============================ HOW TO RUN ============================
#   cd /Users/msc/Desktop/optimal_control/orbit_transfer/GTO_ELFO/direct/elfo
#   chmod +x elfo_energy_sweep.sh            # once
#   ./elfo_energy_sweep.sh <factorLo> <factorHi> <factorStep>
#
# Examples:
#   ./elfo_energy_sweep.sh 1.11 2.00 0.08    # full band (the campaign default)
#   ./elfo_energy_sweep.sh 1.15 1.30 0.05    # a narrow band
#
# Optional env vars (prefix on the command line):
#   MAXTRIES=6        max sweep restarts on crash                 [default 6]
#   MATLAB_BIN=/path/to/matlab                                    [R2025b]
#
# Background + watch:
#   nohup ./elfo_energy_sweep.sh 1.11 2.00 0.08 >/dev/null 2>&1 &
#   tail -f results/logs/elfo_energy_sweep_*.log
#
# When done, feed the banked seeds to the fuel batch:  ./elfo_batch.sh 0 energy
# ===================================================================

set -u
DIR=$(cd "$(dirname "$0")" && pwd)                        # the elfo folder
MAT=${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}
MAXTRIES=${MAXTRIES:-6}
LOGDIR=$DIR/results/logs
mkdir -p "$LOGDIR"
LOG=$LOGDIR/elfo_energy_sweep_$(date +%Y%m%d_%H%M%S).log

[ $# -lt 3 ] && { echo "usage: $0 <factorLo> <factorHi> <factorStep>"; exit 1; }
fLo=$1;  fHi=$2;  fStep=$3

echo "=== ELFO ENERGY SWEEP  band [$fLo, $fHi] step $fStep  (<=$MAXTRIES tries) ===" | tee -a "$LOG"
opts="struct('factorLo',$fLo,'factorHi',$fHi,'factorStep',$fStep)"

ok=0
try=1
while [ $try -le $MAXTRIES ]; do
  echo "" | tee -a "$LOG"
  echo "--- attempt $try/$MAXTRIES  ($(date +%H:%M:%S)) ---" | tee -a "$LOG"
  "$MAT" -batch "cd('$DIR'); setup_paths; gen_elfo_energy_tfsweep($opts);" >> "$LOG" 2>&1
  rc=$?
  if [ $rc -eq 0 ]; then
    echo "  sweep completed (exit 0) on attempt $try" | tee -a "$LOG"
    ok=1
    break
  fi
  echo "  attempt $try exited $rc (likely MEX FATAL) -- resuming from furthest banked seed" | tee -a "$LOG"
  try=$((try+1))
done

echo "" | tee -a "$LOG"
echo "log: $LOG"
if [ $ok -eq 1 ]; then
  exit 0
else
  echo "  GAVE UP after $MAXTRIES attempts (a deterministic crash point? inspect the log)" | tee -a "$LOG"
  exit 1
fi
