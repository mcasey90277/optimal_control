#!/bin/zsh
# RUN_CAMPAIGN_WORKERS  Launch N campaign workers against a work queue,
# VERIFY that each one attached to the queue, and watch each for INACTIVITY.
#
# Every launch bug this campaign hit is guarded here:
#   * a launcher that silently does nothing: the loop is arithmetic (no
#     seq), the number of children is asserted, and each worker must write
#     a READY heartbeat (it does so only after opening the queue) within
#     STARTUP_SEC or it is killed and the launch reported incomplete.
#   * a tag that carried the whole argument: arguments are positional,
#     validated, quoted, echoed back; the job path reaches MATLAB through
#     the environment, never interpolated into MATLAB code.
#   * a watchdog that measured LIFETIME and killed healthy workers inside
#     their fourth column: this one reads the heartbeat AGE and kills only
#     a worker silent for HANG_SEC. Ownership is a kernel lock held by the
#     process (unit_lock), so killing is the ONLY way a unit changes hands
#     -- the watchdog is the supervisor the queue relies on.
#   * relaunch reusing w1..wN and truncating live logs: every launch gets an
#     exclusively created directory (mktemp), and tags and logs carry it.
#
# Usage:
#   run_campaign_workers.sh <jobScript> <nWorkers> <hangSec> <outDir>
#     hangSec   seconds of heartbeat silence after which a worker is killed
#               (the entry script sizes it from the solver's own wall cap)
# The job script reads WORKER_TAG from the environment and calls
# campaign_worker; heartbeats go to <outDir>/hb, the same place the
# generated job writes them.
set -u
usage() { echo "usage: $0 <jobScript> <nWorkers> <hangSec> <outDir>"; exit 2; }
[ $# -eq 4 ] || usage
JOB=$1; N=$2; HANG_SEC=$3; OUT=$4
STARTUP_SEC=240

# ---- validate before anything is spawned -----------------------------------
[ -f "$JOB" ] && [ -r "$JOB" ] || { echo "REFUSING: job script not a readable file: $JOB"; exit 2; }
[[ "$N" == <1-64> ]] || { echo "REFUSING: nWorkers must be an integer in 1..64 (got '$N')"; exit 2; }
[[ "$HANG_SEC" == <60-> ]] || { echo "REFUSING: hangSec must be an integer >= 60 (got '$HANG_SEC')"; exit 2; }
mkdir -p "$OUT/hb" 2>/dev/null && [ -w "$OUT" ] && [ -w "$OUT/hb" ] || { echo "REFUSING: outDir not writable: $OUT"; exit 2; }
MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
[ -x "$MATLAB" ] || { echo "REFUSING: MATLAB not found at $MATLAB"; exit 2; }
HB="$OUT/hb"
LDIR=$(mktemp -d "$OUT/launch_$(date +%Y%m%d_%H%M%S)_XXXX") || { echo "REFUSING: cannot create a launch directory in $OUT"; exit 2; }
LAUNCH=${LDIR:t}
LOG="$OUT/launcher.log"
echo "$LAUNCH: $N worker(s), inactivity watchdog ${HANG_SEC}s, startup deadline ${STARTUP_SEC}s, heartbeats in $HB" | tee -a "$LOG"

# ---- spawn -----------------------------------------------------------------
typeset -a TAGS PIDS
for (( i = 1; i <= N; i++ )); do
  TAG="w${i}_${LAUNCH#launch_}"
  WORKER_TAG="$TAG" CAMPAIGN_JOB="$JOB" nohup "$MATLAB" -batch "run(getenv('CAMPAIGN_JOB'))" \
      > "$LDIR/worker_$TAG.out" 2>&1 < /dev/null &
  MPID=$!
  TAGS+=("$TAG"); PIDS+=("$MPID")
  echo "  spawned $TAG (pid $MPID)" | tee -a "$LOG"
  # SUPERVISOR for this child: startup deadline, then inactivity. It acts on
  # the pid only while that pid is still a MATLAB batch process (a recycled
  # pid is left alone), and logs the child's exit when it sees it.
  ( t0=$(date +%s); HBF="$HB/$TAG.hb"
    isours() { ps -o command= -p $MPID 2>/dev/null | grep -q -- "-batch"; }   # wrapper script or binary
    while isours; do
      sleep 15
      now=$(date +%s)
      if [ ! -f "$HBF" ]; then
        if [ $(( now - t0 )) -gt $STARTUP_SEC ]; then
          isours && kill $MPID; sleep 5; isours && kill -9 $MPID
          echo "$(date +%H:%M:%S) supervisor killed $TAG (pid $MPID): no READY heartbeat in ${STARTUP_SEC}s" >> "$LOG"; break
        fi
        continue
      fi
      AGE=$(( now - $(stat -f %m "$HBF") ))
      if [ $AGE -gt $HANG_SEC ]; then
        isours && kill $MPID; sleep 5; isours && kill -9 $MPID
        echo "$(date +%H:%M:%S) supervisor killed $TAG (pid $MPID): no heartbeat for ${AGE}s" >> "$LOG"; break
      fi
    done
    echo "$(date +%H:%M:%S) $TAG (pid $MPID) is no longer running" >> "$LOG"
  ) > /dev/null 2>&1 < /dev/null &
  sleep 1
done
[ ${#PIDS[@]} -eq $N ] || { echo "LAUNCH FAILED: spawned ${#PIDS[@]} of $N" | tee -a "$LOG"; exit 1; }

# ---- verify: READY means the worker opened the queue -------------------------
echo "verifying every worker attached to the queue ..."
typeset -a FAILED
deadline=$(( $(date +%s) + STARTUP_SEC ))
for (( k = 1; k <= N; k++ )); do
  TAG=${TAGS[$k]}; MPID=${PIDS[$k]}; ok=0
  while [ $(date +%s) -lt $deadline ]; do
    # READY is a persistent file the worker writes once it has opened the
    # queue; the heartbeat itself is overwritten within milliseconds
    if [ -f "$HB/$TAG.ready" ]; then ok=1; break; fi
    if [ -f "$HB/$TAG.hb" ] && grep -q "^fail|" "$HB/$TAG.hb"; then break; fi
    if ! kill -0 $MPID 2>>"$LOG"; then
      # EVIDENCE, not a verdict: say what the process table shows for that pid
      echo "$(date +%H:%M:%S) verify: kill -0 $MPID failed for $TAG; ps says: $(ps -o pid=,ppid=,stat=,command= -p $MPID 2>&1 | cut -c1-100)" >> "$LOG"
      break
    fi
    sleep 2
  done
  if [ $ok -eq 1 ]; then echo "  $TAG: READY"
  else FAILED+=("$TAG")
       if kill -0 $MPID 2>/dev/null; then echo "  $TAG: *** NOT READY (process alive, pid $MPID; its supervisor will kill it at ${STARTUP_SEC}s) ***"
       else echo "  $TAG: *** EXITED BEFORE READY ***"; fi
       echo "      last output:"; tail -5 "$LDIR/worker_$TAG.out" 2>/dev/null | sed 's/^/      /'
  fi
done
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "LAUNCH INCOMPLETE: ${#FAILED[@]} of $N worker(s) not ready: ${FAILED[*]}" | tee -a "$LOG"
  exit 1
fi
echo "LAUNCH OK: all $N workers attached (tags: ${TAGS[*]}; logs in $LDIR)" | tee -a "$LOG"
