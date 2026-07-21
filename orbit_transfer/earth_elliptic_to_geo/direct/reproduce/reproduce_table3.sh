#!/bin/bash
# REPRODUCE_TABLE3.SH  Per-process watchdog orchestrator for the Table-3
# reproducer engine (reproduce_row.m). One MATLAB process per rung, with
# relaunch-on-death (crash OR hang), up to a per-rung attempt cap.
#
# Model: this module's own run_task9_watchdog.sh (crash-relaunch precedent,
# same "grep the log for a genuine-failure marker, else relaunch" pattern)
# and ../GTO_tulip/direct/PSR/psr_batch.sh (one-process-per-item +
# background-launch/poll-loop hang watchdog, since this Mac has neither
# `timeout` nor `gtimeout`).
#
# WHY PER-PROCESS: reproduce_row.m composes casadi_lt_mee.m / IPOPT / MUMPS
# through CasADi's MEX layer, which has a recurring UNCATCHABLE fatal crash
# at this problem scale (documented throughout this campaign -- see
# run_task9_watchdog.sh's header). A try/catch INSIDE MATLAB cannot catch
# it, so the only robust recovery is one MATLAB process per rung, watched
# from the OUTSIDE, with a relaunch on death. reproduce_row.m's own
# REPRO_-tagged driver caches (per-round/per-stage files under results/ and
# results/repro/) are resume-safe, so a relaunch picks back up close to
# where the dead process left off rather than re-paying the whole rung.
#
# HANG DETECTION: no `timeout`/`gtimeout` on this Mac, so each attempt is
# launched in the BACKGROUND and polled with `kill -0 $pid` every ~15s; if
# still alive past the per-attempt wall cap, it is `kill -9`'d and counted
# as a (hang) death, then relaunched exactly like a crash.
#
# DEATH CLASSIFICATION: a nonzero exit code alone is NOT trusted to mean
# "crashed" (a genuine solver assertion also exits nonzero, and a hard MEX
# crash can look like an ordinary nonzero exit too) -- each attempt's death
# is cross-checked against a FRESH MATLAB crash dump (~/matlab_crash_dump.*,
# MATLAB's own post-mortem file) written strictly AFTER that attempt's own
# launch time. If one exists it is logged (filename + attempt) as positive
# evidence of a MEX-level fatal crash rather than an ordinary nonzero exit.
# A GENUINE (non-crash) solver/config failure is recognized by grepping the
# recent log for reproduce_row.m's own assertion error IDs
# (anchorNotCertified/fuelAllFailed/psrNotCertified/unknownStrategy/
# prevMissing/noTolForThrust) or table3_certified.m/table3_recipes.m's
# unknownThrust -- when one of those fires, retries for that rung STOP
# (matching run_task9_watchdog.sh's "stop, don't mask, a genuine failure"
# policy) rather than burning the whole attempt budget on a deterministic
# block.
#
# WALL CAP: flat 21600 s (6 h) per ATTEMPT for every rung, not scaled by
# thrust. This build has not yet run any of the five rungs through
# reproduce_row.m start-to-finish under THIS orchestrator, so there is no
# real per-rung wall-clock data to calibrate a thrust-scaled cap against --
# guessing a tight per-rung number risks killing a legitimate deep-rung
# solve (1 N's smallN_first continuation alone has its own internal 2.5 h
# budget; 0.5 N's PSR pass adds more on top) mid-flight, which is worse
# than an occasionally-too-generous flat cap. 6 h is "clearly enough for a
# healthy solve, clearly finite so a genuine hang doesn't run forever" for
# every rung including the deepest. Override with WALLCAP_S=<seconds> once
# real timings (Task 6's controller-run) recalibrate this.
#
# USAGE:
#   cd earth_elliptic_to_geo/direct/reproduce
#   chmod +x reproduce_table3.sh          # once
#   ./reproduce_table3.sh                 # default rung list: 10 5 2.5 1 0.5
#   ./reproduce_table3.sh 10 5            # only these rungs
#   WALLCAP_S=3600 MAXATTEMPTS=4 POLL_S=5 ./reproduce_table3.sh 1
#
# Run detached so a closed terminal doesn't kill it:
#   nohup ./reproduce_table3.sh >/dev/null 2>&1 &
#   tail -f ../results/repro/reproduce_table3.log
#
# RESUMABLE: reproduce_row.m's own driver-level caches mean re-running this
# script (or a single rung) after an interruption is cheap for whatever
# already converged.
#
# NOTE: this script is validated in this build by `bash -n` (syntax) and
# inspection only -- the full multi-hour rung ladder is NOT launched here
# (that is a later, controller-run task).
set -u

DIR=$(cd "$(dirname "$0")/.." && pwd)   # module root: this script lives in reproduce/
MAT=${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}
CASADI_DIR=${CASADI_PATH:-"$HOME/casadi-3.7.0"}
WALLCAP_S=${WALLCAP_S:-21600}      # 6 h per attempt, flat (see header)
MAXATTEMPTS=${MAXATTEMPTS:-8}      # per-rung attempt cap
POLL_S=${POLL_S:-15}       # poll interval for the hang-watchdog loop, in seconds

mkdir -p "$DIR/results/repro"
LOG="$DIR/results/repro/reproduce_table3.log"

if [ $# -gt 0 ]; then
  rungs=("$@")
else
  rungs=(10 5 2.5 1 0.5)
fi

# round(10*T) -> reproduce_row.m's REPRO_row_T<...>.mat filename suffix
round10() {
  awk -v t="$1" 'BEGIN{printf "%d", int(t*10+0.5)}'
}

echo "=== REPRODUCE_TABLE3 watchdog start $(date) rungs=(${rungs[*]}) wallcap=${WALLCAP_S}s maxattempts=$MAXATTEMPTS ===" | tee -a "$LOG"

for T in "${rungs[@]}"; do
  echo "" | tee -a "$LOG"
  echo "=== RUNG T=$T N ===" | tee -a "$LOG"
  rowfile="$DIR/results/repro/REPRO_row_T$(round10 "$T").mat"
  attempt=1
  rungDone=0

  while [ "$attempt" -le "$MAXATTEMPTS" ]; do
    tsfile=$(mktemp /tmp/reproduce_table3_ts.XXXXXX)
    echo "  --- attempt $attempt/$MAXATTEMPTS at $(date) ---" | tee -a "$LOG"

    "$MAT" -batch "addpath('$CASADI_DIR'); cd('$DIR'); setup_paths; reproduce_row($T)" >> "$LOG" 2>&1 &
    pid=$!
    waited=0
    hung=0
    while kill -0 "$pid" 2>/dev/null; do
      sleep "$POLL_S"
      waited=$((waited+POLL_S))
      if [ "$waited" -ge "$WALLCAP_S" ]; then
        echo "  HANG: T=$T attempt $attempt exceeded ${WALLCAP_S}s wall cap -- kill -9 $pid" | tee -a "$LOG"
        kill -9 "$pid" 2>/dev/null
        hung=1
        break
      fi
    done
    wait "$pid" 2>/dev/null
    rc=$?

    # death classification: a fresh crash dump written after THIS attempt's
    # own launch (tsfile's mtime), independent of the exit code / hang flag
    crashdump=$(find "$HOME" -maxdepth 1 -iname 'matlab_crash_dump.*' -newer "$tsfile" 2>/dev/null | head -1)
    rm -f "$tsfile"

    if [ "$rc" -eq 0 ] && [ "$hung" -eq 0 ]; then
      echo "  attempt $attempt: exit 0" | tee -a "$LOG"
    elif [ -n "$crashdump" ]; then
      echo "  attempt $attempt: CRASH (fresh dump: $crashdump) rc=$rc hung=$hung" | tee -a "$LOG"
    elif [ "$hung" -eq 1 ]; then
      echo "  attempt $attempt: HANG (killed after ${WALLCAP_S}s) rc=$rc" | tee -a "$LOG"
    else
      echo "  attempt $attempt: nonzero exit rc=$rc (no fresh crash dump found -- may be a genuine solver failure, see log above)" | tee -a "$LOG"
    fi

    if [ -f "$rowfile" ] && [ "$rc" -eq 0 ] && [ "$hung" -eq 0 ]; then
      rungDone=1
      break
    fi

    # genuine (non-crash) solver/config failure -- stop retrying, don't mask it
    if tail -n 300 "$LOG" | grep -qE \
        "reproduce_row:(badInput|anchorNotCertified|fuelAllFailed|psrNotCertified|unknownStrategy|prevMissing|noTolForThrust)|table3_certified:unknownThrust|table3_recipes:unknownThrust"; then
      echo "  RUNG T=$T BLOCKED: a genuine solver/config assertion fired (not a crash) -- stopping retries for this rung, needs attention" | tee -a "$LOG"
      break
    fi

    attempt=$((attempt+1))
  done

  if [ "$rungDone" -eq 1 ]; then
    echo "  RUNG T=$T: DONE after $attempt attempt(s)" | tee -a "$LOG"
  else
    echo "  RUNG T=$T: NOT confirmed done after $((attempt<MAXATTEMPTS?attempt:MAXATTEMPTS)) attempt(s) -- see log" | tee -a "$LOG"
  fi
done

echo "" | tee -a "$LOG"
echo "=== building final Table 3 (reproduce_table3_collect) ===" | tee -a "$LOG"
rungsMatlab=$(printf '%s,' "${rungs[@]}")
rungsMatlab=${rungsMatlab%,}
"$MAT" -batch "cd('$DIR'); setup_paths; reproduce_table3_collect([$rungsMatlab]);" 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "log: $LOG"
