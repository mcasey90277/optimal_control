#!/bin/zsh
# RUN_CAMPAIGN_WORKERS  Launch N campaign workers against a work queue,
# VERIFY that each one actually started, and watch them for INACTIVITY.
#
# Every launch bug this campaign hit is guarded here:
#   * a launcher that silently does nothing. A shell loop that split a
#     bracketed column list left three workers unlaunched, and nothing said
#     so for hours. This waits for each worker's HEARTBEAT and reports
#     loudly if it does not appear or if the process exits first.
#   * a tag that carried the whole argument. Arguments are positional,
#     validated, quoted, and echoed back.
#   * a watchdog that measured LIFETIME. A timer started at launch killed
#     healthy workers inside their fourth column. The watchdog here reads
#     the worker's heartbeat AGE and kills only a worker that has reported
#     nothing for WATCHDOG seconds -- with per-solve beats that is a hang,
#     not a long column. (Astra chain review 2026-09-13.)
#   * relaunch reusing w1..wN and truncating their logs while the old
#     workers still ran. Every launch gets an id; tags and logs carry it.
#
# Usage:
#   run_campaign_workers.sh <jobScript> <nWorkers> <unitSec> <outDir>
# The job script must read WORKER_TAG from the environment and call
# campaign_worker with it; heartbeats go to <outDir>/hb -- the same place
# the generated job writes them, so launcher and worker cannot disagree.
set -u
usage() { echo "usage: $0 <jobScript> <nWorkers> <unitSec> <outDir>"; exit 2; }
[ $# -eq 4 ] || usage
JOB=$1; N=$2; UNIT_SEC=$3; OUT=$4

# ---- validate before anything is spawned -----------------------------------
[ -r "$JOB" ] || { echo "REFUSING: job script not readable: $JOB"; exit 2; }
[[ "$N" == <1-> ]] || { echo "REFUSING: nWorkers must be a positive integer (got '$N')"; exit 2; }
[[ "$UNIT_SEC" == <1-> ]] || { echo "REFUSING: unitSec must be a positive integer (got '$UNIT_SEC')"; exit 2; }
mkdir -p "$OUT" 2>/dev/null && [ -w "$OUT" ] || { echo "REFUSING: outDir not writable: $OUT"; exit 2; }
MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
[ -x "$MATLAB" ] || { echo "REFUSING: MATLAB not found at $MATLAB"; exit 2; }
HB="$OUT/hb"; mkdir -p "$HB"
WATCHDOG=$(( UNIT_SEC * 4 ))
[ $WATCHDOG -ge 1800 ] || WATCHDOG=1800          # never shorter than the stale lease
LAUNCH=$(date +%Y%m%d_%H%M%S)
LOG="$OUT/launcher.log"
echo "launch $LAUNCH: $N worker(s), unit ~${UNIT_SEC}s, inactivity watchdog ${WATCHDOG}s, heartbeats in $HB" | tee -a "$LOG"

# ---- spawn -----------------------------------------------------------------
typeset -a TAGS PIDS
for i in $(seq 1 $N); do
  TAG="w${i}_${LAUNCH}"
  TAGS+=("$TAG")
  WORKER_TAG="$TAG" nohup "$MATLAB" -batch "run('$JOB')" > "$OUT/worker_$TAG.out" 2>&1 &
  MPID=$!
  PIDS+=("$MPID")
  echo "  spawned $TAG (pid $MPID)" | tee -a "$LOG"
  # INACTIVITY watchdog: kill only when the heartbeat has not moved for
  # WATCHDOG seconds; exit quietly when the worker exits on its own
  ( while kill -0 $MPID 2>/dev/null; do
      sleep 60
      HBF="$HB/$TAG.hb"
      if [ -f "$HBF" ]; then
        AGE=$(( $(date +%s) - $(stat -f %m "$HBF") ))
        if [ $AGE -gt $WATCHDOG ] && kill -0 $MPID 2>/dev/null; then
          kill $MPID; sleep 5; kill -9 $MPID 2>/dev/null
          echo "watchdog killed $TAG (pid $MPID): no heartbeat for ${AGE}s $(date)" >> "$LOG"
          break
        fi
      fi
    done ) > /dev/null 2>&1 &
  sleep 2
done

# ---- verify: a worker that did not start must not look like one that did ----
echo "verifying every worker checked in ..."
typeset -a FAILED
for k in $(seq 1 $N); do
  TAG=${TAGS[$k]}; MPID=${PIDS[$k]}
  ok=0
  for t in $(seq 1 90); do                        # up to 180 s each, concurrent starts
    if [ -f "$HB/$TAG.hb" ]; then ok=1; break; fi
    if ! kill -0 $MPID 2>/dev/null; then break; fi   # exited before its first beat
    sleep 2
  done
  if [ $ok -eq 1 ]; then echo "  $TAG: checked in"
  else FAILED+=("$TAG")
       if kill -0 $MPID 2>/dev/null; then echo "  $TAG: *** NO HEARTBEAT YET (process alive, pid $MPID) -- NOT verified ***"
       else echo "  $TAG: *** EXITED BEFORE ITS FIRST HEARTBEAT ***"; fi
       echo "      last output:"; tail -5 "$OUT/worker_$TAG.out" 2>/dev/null | sed 's/^/      /'
  fi
done
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "LAUNCH INCOMPLETE: ${#FAILED[@]} of $N worker(s) not verified: ${FAILED[*]}" | tee -a "$LOG"
  exit 1
fi
echo "LAUNCH OK: all $N workers checked in (tags: ${TAGS[*]})" | tee -a "$LOG"
