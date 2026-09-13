#!/bin/zsh
# CAMPAIGN_SUPERVISOR  Keep a campaign at its requested worker count until
# its queue is finished, then run its finalizer once. The long-lived
# controller a multi-day unattended run needs (Astra pass 2/3: "no durable
# supervisor replaces lost capacity or triggers finalization").
#
# Every minute:
#   * count live workers from the launchers' pid files (kill -0);
#   * decide FINISHED from the queue's own files: every expected output
#     exists, or the missing ones have spent their attempts (.att >= maxAtt)
#     -- the same rule work_queue applies, read without starting MATLAB;
#   * if not finished and live < N, launch the deficit through
#     run_campaign_workers.sh (a new launch directory, unique tags), within a
#     restart budget so a worker that dies on startup cannot loop forever;
#   * when finished, run the finalizer (a MATLAB job script) exactly once
#     and exit with its status. If the remaining units are all retired, the
#     campaign is BLOCKED: no finalizer, exit 4, verdict written.
# One instance per campaign directory (a lock directory).
#
# Usage:
#   campaign_supervisor.sh <jobScript> <nWorkers> <hangSec> <outDir> <queueDir> <maxAtt> [finalizeJob]
set -u
usage() { echo "usage: $0 <jobScript> <nWorkers> <hangSec> <outDir> <queueDir> <maxAtt> [finalizeJob]"; exit 2; }
[ $# -ge 6 ] || usage
JOB=$1; N=$2; HANG_SEC=$3; OUT=$4; QDIR=$5; MAXATT=$6; FINAL=${7:-}
HERE=${0:A:h}
LAUNCHER=$HERE/run_campaign_workers.sh
MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
POLL=${POLL_SEC:-60}
MAXLAUNCH=$(( 3 * N ))          # total worker launches allowed for the campaign
LOG=$OUT/supervisor.log
VERDICT=$OUT/SUPERVISOR_VERDICT.txt
log() { echo "$(date +%Y-%m-%d\ %H:%M:%S) $*" >> "$LOG"; }
verdict() { { echo "=== $(date) ==="; echo "$*"; } >> "$VERDICT"; log "VERDICT: $*"; }

[ -r "$JOB" ] || { echo "REFUSING: job not readable: $JOB"; exit 2; }
[[ "$N" == <1-64> ]] || { echo "REFUSING: nWorkers 1..64 (got '$N')"; exit 2; }
[ -s "$QDIR/expected_units.txt" ] || { echo "REFUSING: no $QDIR/expected_units.txt (run the entry script once)"; exit 2; }
[ -z "$FINAL" ] || [ -r "$FINAL" ] || { echo "REFUSING: finalize job not readable: $FINAL"; exit 2; }
mkdir -p "$OUT"
if ! mkdir "$OUT/supervisor.lock" 2>/dev/null; then echo "REFUSING: another supervisor holds $OUT/supervisor.lock"; exit 3; fi
trap 'rmdir "$OUT/supervisor.lock" 2>/dev/null' EXIT
log "supervisor start: N=$N hangSec=$HANG_SEC queue=$QDIR maxAtt=$MAXATT finalize=${FINAL:-none} launch budget=$MAXLAUNCH"

launched=0; failedLaunches=0
live_count() {
  local n=0 p
  for pf in "$OUT"/launch_*/pid_*(N); do
    p=$(cat "$pf" 2>/dev/null) || continue
    kill -0 "$p" 2>/dev/null && n=$(( n + 1 ))
  done
  echo $n
}
# queue census from files: prints "missing retired"
census() {
  local missing=0 retired=0 f id att
  while IFS= read -r f; do
    [ -f "$f" ] && continue
    missing=$(( missing + 1 ))
    id=$(basename "$f" | sed -E 's/^[^0-9]*0*([0-9]+).*$/\1/')
    att=$(cat "$QDIR/$id.att" 2>/dev/null || echo 0)
    [[ "$att" == <0-> ]] || att=999           # unreadable record = blocked
    [ "$att" -ge "$MAXATT" ] && retired=$(( retired + 1 ))
  done < "$QDIR/expected_units.txt"
  echo "$missing $retired"
}

while true; do
  read missing retired <<< "$(census)"
  live=$(live_count)
  if [ "$missing" -eq 0 ]; then
    log "every expected output exists (live workers: $live); finalizing"
    break
  fi
  if [ "$retired" -eq "$missing" ]; then
    verdict "BLOCKED: $missing unit(s) missing and all of them RETIRED after $MAXATT attempts; not finalizing"
    exit 4
  fi
  if [ "$live" -lt "$N" ]; then
    deficit=$(( N - live ))
    if [ $(( launched + deficit )) -gt $MAXLAUNCH ]; then
      deficit=$(( MAXLAUNCH - launched ))
    fi
    if [ "$deficit" -le 0 ]; then
      if [ "$live" -eq 0 ]; then
        verdict "STALLED: launch budget ($MAXLAUNCH) spent, no live worker, $missing unit(s) missing ($retired retired)"
        exit 5
      fi
    else
      log "capacity $live/$N with $missing unit(s) missing: launching $deficit worker(s) (launched so far $launched)"
      if "$LAUNCHER" "$JOB" "$deficit" "$HANG_SEC" "$OUT" >> "$LOG" 2>&1; then
        launched=$(( launched + deficit )); failedLaunches=0
      else
        launched=$(( launched + deficit )); failedLaunches=$(( failedLaunches + 1 ))
        log "launch reported INCOMPLETE ($failedLaunches consecutive)"
        if [ "$failedLaunches" -ge 3 ]; then
          verdict "STALLED: three consecutive launches failed to produce READY workers; see $LOG"
          exit 6
        fi
      fi
    fi
  fi
  sleep "$POLL"
done

if [ -n "$FINAL" ]; then
  log "running finalizer $FINAL"
  "$MATLAB" -batch "run('$FINAL')" > "$OUT/finalize_job.out" 2>&1
  rc=$?
  if [ $rc -eq 0 ]; then verdict "FINISHED: all outputs present; finalizer exited 0"
  else verdict "FINALIZER FAILED (exit $rc): see $OUT/finalize_job.out"; fi
  exit $rc
fi
verdict "FINISHED: all outputs present (no finalizer given)"
exit 0
