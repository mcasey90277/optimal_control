#!/bin/zsh
# sharpen_batch.sh — sharpen energy-backbone factors to certified min-fuel,
# in parallel batches, via the canonical driver minfuel_at_tf.
#
# Usage:   ./sharpen_batch.sh <batch_size> <factor1> [factor2 ...]
# Example: ./sharpen_batch.sh 2 1.30 1.35 1.40 1.45
#
# Each factor is INDEPENDENT (tight re-clean + fine energy->fuel homotopy from
# its own energy_<f>.mat), so factors run concurrently up to <batch_size>
# (keep it small if another MATLAB campaign shares the machine). Per-solve:
# own process + watchdog + ONE retry — the sporadic uncatchable CasADi/IPOPT
# MEX FATAL crash (~1 in 10 solves) almost always clears on a fresh retry.
#
# Output: results/minfuel/minfuel_f####_en.mat (provenance-stamped by
# minfuel_at_tf). Aggregate + certify + plot with aggregate_front.
#
# ZSH GOTCHA (cost 13 result files on 2026-07-09): never write
#   local f=$1 out=$DIR/x_$f.mat
# on ONE line — every word expands BEFORE the assignments execute, so $f is
# empty in $out. Keep each local on its own line.

set -u
DIR=$(cd "$(dirname "$0")/.." && pwd)          # sundman_minfuel library root
MAT=${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}
WATCHDOG=${WATCHDOG_S:-1800}                   # per-solve kill timeout [s]
LOGDIR=$DIR/results/logs
mkdir -p "$LOGDIR" "$DIR/results/minfuel"
LOG=$LOGDIR/sharpen_batch_$(date +%Y%m%d_%H%M%S).log

[ $# -lt 2 ] && { echo "usage: $0 <batch_size> <factor1> [factor2 ...]"; exit 1; }
BATCH=$1; shift

expected_out() {
  # mirror minfuel_config's fname: minfuel_f<milli>_en.mat
  local f=$1
  local milli
  milli=$(printf '%04.0f' $((f * 1000)))     # zsh does float arithmetic natively
  echo "$DIR/results/minfuel/minfuel_f${milli}_en.mat"
}

sharpen_one() {
  local f=$1
  local out
  out=$(expected_out $f)
  local milli
  milli=$(printf '%04.0f' $((f * 1000)))
  if [ ! -f "$DIR/energy_$f.mat" ] && [ ! -f "$DIR/results/energy/energy_f$milli.mat" ]; then
    echo "  SKIP ${f}x (no energy backbone file)" >> "$LOG"; return
  fi
  local try
  for try in 1 2; do
    rm -f "$out"
    "$MAT" -batch "cd('$DIR'); minfuel_at_tf($f,'seed','energy')" \
        > "$LOGDIR/sharpen_${f}_try$try.log" 2>&1 &
    local pid=$!
    local waited=0
    while kill -0 $pid 2>/dev/null; do
      sleep 15; waited=$((waited+15))
      if [ $waited -ge $WATCHDOG ]; then
        echo "  TIMEOUT ${f}x try$try -> kill" >> "$LOG"; kill -9 $pid 2>/dev/null; break
      fi
    done
    wait $pid 2>/dev/null
    if [ -f "$out" ] && grep -q "MINFUEL_AT_TF done" "$LOGDIR/sharpen_${f}_try$try.log"; then
      echo "  ${f}x try$try: $(grep 'MINFUEL_AT_TF done' "$LOGDIR/sharpen_${f}_try$try.log" | tail -1)" >> "$LOG"
      return
    fi
    echo "  ${f}x try$try FAILED (likely MEX crash) -> retry" >> "$LOG"
  done
  echo "  ${f}x GAVE UP after 2 tries" >> "$LOG"
}

# run in batches of $BATCH
pids=()
count=0
for f in "$@"; do
  echo "=== launch ${f}x ===" >> "$LOG"
  sharpen_one $f &
  pids+=($!)
  count=$((count+1))
  if [ $count -ge $BATCH ]; then
    for p in "${pids[@]}"; do wait $p; done
    pids=(); count=0
  fi
done
for p in "${pids[@]}"; do wait $p; done
echo "SHARPEN_BATCH_DONE" >> "$LOG"
echo "done; log: $LOG"
