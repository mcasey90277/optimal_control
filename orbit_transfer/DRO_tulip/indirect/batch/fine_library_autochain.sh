#!/bin/zsh
# FINE_LIBRARY_AUTOCHAIN  Wait for the 24 x 24 rib campaign to drain, then run
# the finish job (package, audit, sweep, movie). Conditions, all needed:
#   * every EXPECTED rib file (ribq/expected_units.txt, written by the entry
#     script from the sheet) exists -- not a count of a glob;
#   * no old explicit-range rib worker is alive (they hold no queue lock).
# One instance at a time (a lock directory); verdicts are appended, never
# overwritten; the finish job's exit status is recorded and is non-zero
# unless the library packaged and the movie rendered.
ROOT=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect
R=$ROOT/results_fine
MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
LOG=$R/finish_autochain.log
EXP=$R/ribq/expected_units.txt
if ! mkdir "$R/finish.lock" 2>/dev/null; then
  echo "$(date) another autochain holds $R/finish.lock; not starting" >> $LOG; exit 3
fi
trap 'rmdir "$R/finish.lock" 2>/dev/null' EXIT
[ -s "$EXP" ] || { echo "$(date) BLOCKED: no expected-units list at $EXP (run run_costate_library once)" >> $LOG; exit 2; }
echo "$(date) autochain armed: waiting for $(wc -l < $EXP | tr -d ' ') expected rib files and no rib worker" >> $LOG
while true; do
  missing=0; while IFS= read -r f; do [ -f "$f" ] || missing=$(( missing + 1 )); done < "$EXP"
  live=$(pgrep -f fine_ribs_range_job | wc -l | tr -d ' ')
  if [ "$missing" -eq 0 ] && [ "$live" -eq 0 ]; then break; fi
  if [ "$live" -eq 0 ] && [ "$missing" -gt 0 ]; then
    echo "$(date) BLOCKED: no rib worker alive but $missing expected rib file(s) missing; finish NOT started" >> $LOG
    { echo "=== $(date) ==="; echo "BLOCKED: $missing expected rib file(s) missing and no worker alive"; } >> $R/FINISH_VERDICT.txt
    exit 1
  fi
  sleep 300
done
echo "$(date) every expected rib file present, no rib worker alive: launching the finish job" >> $LOG
"$MATLAB" -batch "run('$ROOT/batch/fine_library_finish_job.m')" > $R/finish_job.out 2>&1
rc=$?
echo "$(date) finish job exited $rc" >> $LOG
grep -q "FINISH JOB END" $R/FINISH_VERDICT.txt 2>/dev/null || { echo "=== $(date) ==="; echo "FINISH JOB DIED WITHOUT A VERDICT (exit $rc)"; } >> $R/FINISH_VERDICT.txt
exit $rc
