#!/bin/zsh
# RUN_CAMPAIGN_WORKERS  Launch N campaign workers against a work queue, each
# under a SUPERVISOR that is its PARENT: it verifies READY, watches for
# inactivity, kills a hung or unready worker, confirms the exit, reaps it,
# and records the exit code and reason.
#
# Every launch bug this campaign hit is guarded here:
#   * a launcher that silently does nothing: arithmetic loop, child count
#     asserted, and each worker must write a READY marker (only after it
#     has opened the queue) within STARTUP_SEC or its supervisor kills it.
#   * a tag that carried the whole argument: positional, validated, quoted;
#     the job path reaches MATLAB through the environment.
#   * a watchdog that measured LIFETIME: this one measures heartbeat AGE,
#     and only after READY. Ownership is a kernel lock the process holds
#     (unit_lock), so a kill is the ONLY way a unit changes hands, and the
#     supervisor confirms the exit before it records anything.
#   * a sibling watchdog acting on a pid it did not own (pass 3): the
#     supervisor is the parent, so the pid cannot be recycled under it
#     before `wait` reaps it, and the exit code is real.
#   * relaunch reusing tags and truncating live logs: every launch gets an
#     exclusively created directory; tags, logs, pid and exit files live in it.
#
# Usage:
#   run_campaign_workers.sh <jobScript> <nWorkers> <hangSec> <outDir>
#     hangSec   seconds of heartbeat silence after which a worker is killed;
#               the entry script sizes it from the largest capped STAGE
# Env: STARTUP_SEC [240] spawn-to-READY deadline.
# Per launch directory: worker_<tag>.out, pid_<tag>, exit_<tag> ("<rc> <reason>").
set -u
usage() { echo "usage: $0 <jobScript> <nWorkers> <hangSec> <outDir>"; exit 2; }
[ $# -eq 4 ] || usage
JOB=$1; N=$2; HANG_SEC=$3; OUT=$4
STARTUP_SEC=${STARTUP_SEC:-240}

# ---- validate before anything is spawned -----------------------------------
[ -f "$JOB" ] && [ -r "$JOB" ] || { echo "REFUSING: job script not a readable file: $JOB"; exit 2; }
[[ "$N" == <1-64> ]] || { echo "REFUSING: nWorkers must be an integer in 1..64 (got '$N')"; exit 2; }
[[ "$HANG_SEC" == <60-> ]] || { echo "REFUSING: hangSec must be an integer >= 60 (got '$HANG_SEC')"; exit 2; }
[[ "$STARTUP_SEC" == <30-> ]] || { echo "REFUSING: STARTUP_SEC must be an integer >= 30 (got '$STARTUP_SEC')"; exit 2; }
mkdir -p "$OUT/hb" 2>/dev/null && [ -w "$OUT" ] && [ -w "$OUT/hb" ] || { echo "REFUSING: outDir not writable: $OUT"; exit 2; }
MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
[ -x "$MATLAB" ] || { echo "REFUSING: MATLAB not found at $MATLAB"; exit 2; }
HB="$OUT/hb"
LDIR=$(mktemp -d "$OUT/launch_$(date +%Y%m%d_%H%M%S)_XXXX") || { echo "REFUSING: cannot create a launch directory in $OUT"; exit 2; }
LAUNCH=${LDIR:t}
LOG="$OUT/launcher.log"
log() { echo "$(date +%H:%M:%S) $*" >> "$LOG"; }
echo "$LAUNCH: $N worker(s), inactivity ${HANG_SEC}s after READY, startup deadline ${STARTUP_SEC}s, heartbeats in $HB" | tee -a "$LOG"

# ---- the supervisor: parent of exactly one worker -----------------------------
supervise() {
  local TAG=$1
  WORKER_TAG="$TAG" CAMPAIGN_JOB="$JOB" "$MATLAB" -batch "run(getenv('CAMPAIGN_JOB'))" \
      > "$LDIR/worker_$TAG.out" 2>&1 < /dev/null &
  local MPID=$!
  echo $MPID > "$LDIR/pid_$TAG"
  log "spawned $TAG (pid $MPID)"
  local t0=$(date +%s) now reason="" AGE
  local HBF="$HB/$TAG.hb" RDY="$HB/$TAG.ready"
  while kill -0 $MPID 2>/dev/null; do
    sleep 5
    now=$(date +%s)
    if [ ! -f "$RDY" ]; then
      # the startup deadline applies until READY, whatever else is written
      if [ $(( now - t0 )) -gt $STARTUP_SEC ]; then reason="no READY within ${STARTUP_SEC}s"; break; fi
    elif [ -f "$HBF" ]; then
      AGE=$(( now - $(stat -f %m "$HBF" 2>/dev/null || echo $now) ))
      if [ $AGE -gt $HANG_SEC ]; then reason="no heartbeat for ${AGE}s"; break; fi
    fi
  done
  if [ -n "$reason" ]; then
    local KIDS=$(pgrep -P $MPID 2>/dev/null | tr '\n' ' ')
    log "supervisor killing $TAG (pid $MPID): $reason"
    kill $MPID 2>/dev/null
    local k
    for k in 1 2 3 4 5 6; do sleep 5; kill -0 $MPID 2>/dev/null || break; done
    if kill -0 $MPID 2>/dev/null; then kill -9 $MPID 2>/dev/null; log "$TAG did not exit on TERM; sent KILL"; fi
    if [ -n "$KIDS" ]; then sleep 2; kill -9 ${=KIDS} 2>/dev/null; log "$TAG: killed descendants $KIDS (its pool)"; fi
  fi
  wait $MPID
  local rc=$?
  echo "$rc ${reason:-not killed by supervisor}" > "$LDIR/exit_$TAG"
  log "$TAG (pid $MPID) exited rc=$rc (${reason:-not killed by supervisor})"
}

# ---- spawn -----------------------------------------------------------------
typeset -a TAGS
for (( i = 1; i <= N; i++ )); do
  TAG="w${i}_${LAUNCH#launch_}"
  TAGS+=("$TAG")
  ( supervise "$TAG" ) > /dev/null 2>&1 < /dev/null &
  sleep 1
done
[ ${#TAGS[@]} -eq $N ] || { echo "LAUNCH FAILED: spawned ${#TAGS[@]} of $N" | tee -a "$LOG"; exit 1; }

# ---- verify: READY means the worker opened the queue, and it is THIS worker ----
echo "verifying every worker attached to the queue ..."
typeset -a FAILED
deadline=$(( $(date +%s) + STARTUP_SEC + 30 ))
for (( k = 1; k <= N; k++ )); do
  TAG=${TAGS[$k]}; ok=0; why=""
  while [ $(date +%s) -lt $deadline ]; do
    if [ -f "$LDIR/pid_$TAG" ] && [ -f "$HB/$TAG.ready" ]; then
      MPID=$(cat "$LDIR/pid_$TAG"); RPID=$(cut -d'|' -f3 "$HB/$TAG.ready")
      if [ "$RPID" = "$MPID" ]; then
        if [ -f "$HB/$TAG.hb" ] && grep -q "^fail|" "$HB/$TAG.hb"; then why="attached, then FAILED: $(cut -d'|' -f4- "$HB/$TAG.hb" | cut -c1-120)"; break; fi
        ok=1; break
      else
        why="READY marker carries pid $RPID, not this worker's $MPID"; break
      fi
    fi
    if [ -f "$LDIR/exit_$TAG" ]; then why="exited before READY: $(cat "$LDIR/exit_$TAG")"; break; fi
    sleep 2
  done
  if [ $ok -eq 1 ]; then echo "  $TAG: READY (pid $(cat "$LDIR/pid_$TAG"))"
  else FAILED+=("$TAG"); echo "  $TAG: *** NOT READY -- ${why:-no READY marker before the deadline; its supervisor will kill it} ***"
       echo "      last output:"; tail -5 "$LDIR/worker_$TAG.out" 2>/dev/null | sed 's/^/      /'
  fi
done
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "LAUNCH INCOMPLETE: ${#FAILED[@]} of $N worker(s) not ready: ${FAILED[*]}" | tee -a "$LOG"
  exit 1
fi
echo "LAUNCH OK: all $N workers attached (tags: ${TAGS[*]}; pids, logs and exit records in $LDIR)" | tee -a "$LOG"
