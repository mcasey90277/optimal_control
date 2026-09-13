#!/bin/zsh
# RUN_CAMPAIGN_WORKERS  Launch N campaign workers against a work queue, and
# VERIFY that each one actually started.
#
# Every launch bug this campaign hit is guarded here:
#   * a launcher that silently does nothing. A shell loop that split a
#     bracketed column list left three workers unlaunched, and nothing said
#     so for hours. This waits for each worker's HEARTBEAT to appear and
#     reports loudly if it does not.
#   * a tag that carried the whole argument. Arguments are passed
#     positionally and echoed back, so a mangled tag is visible immediately.
#   * a watchdog shorter than a unit of work, which destroyed nine hours of
#     walking twice. WATCHDOG defaults to 4x the measured unit time and the
#     script refuses to start if it is shorter than UNIT_SEC.
#
# Usage:
#   run_campaign_workers.sh <jobScript> <nWorkers> <unitSec> [outDir]
# The job script must read WORKER_TAG from the environment and call
# campaign_worker with it.
set -u
JOB=${1:?job script}; N=${2:?n workers}; UNIT_SEC=${3:?seconds for one unit}
OUT=${4:-$(dirname $JOB)/campaign}
MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
[ -x $MATLAB ] || MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab
HB=$OUT/hb
WATCHDOG=$(( UNIT_SEC * 4 ))
mkdir -p $HB

if [ $WATCHDOG -le $UNIT_SEC ]; then
  echo "REFUSING: watchdog ${WATCHDOG}s is not longer than one unit (${UNIT_SEC}s)."
  echo "A watchdog inside a unit destroys that unit's work instead of protecting it."
  exit 1
fi
echo "launching $N worker(s), unit ~${UNIT_SEC}s, watchdog ${WATCHDOG}s, heartbeats in $HB"

TAGS=()
for i in $(seq 1 $N); do
  TAG="w$i"
  TAGS+=($TAG)
  rm -f $HB/$TAG.hb
  WORKER_TAG=$TAG $MATLAB -batch "run('$JOB')" > $OUT/worker_$TAG.out 2>&1 &
  MPID=$!
  ( sleep $WATCHDOG; if kill -0 $MPID 2>/dev/null; then kill $MPID
      echo "watchdog killed $TAG $(date)" >> $OUT/launcher.log; fi ) &
  echo "  spawned $TAG (pid $MPID)" >> $OUT/launcher.log
  sleep 2
done

# VERIFY: a worker that did not start must not look like one that did.
echo "verifying every worker checked in ..."
FAILED=()
for TAG in $TAGS; do
  ok=0
  for t in $(seq 1 60); do
    if [ -f $HB/$TAG.hb ]; then ok=1; break; fi
    sleep 2
  done
  if [ $ok -eq 1 ]; then echo "  $TAG: checked in"
  else FAILED+=($TAG); echo "  $TAG: *** NO HEARTBEAT AFTER 120 s -- DID NOT START ***"
       echo "      last output:"; tail -5 $OUT/worker_$TAG.out 2>/dev/null | sed 's/^/      /'
  fi
done
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "LAUNCH INCOMPLETE: ${#FAILED[@]} of $N worker(s) never started: ${FAILED[*]}"
  exit 1
fi
echo "LAUNCH OK: all $N workers running"
