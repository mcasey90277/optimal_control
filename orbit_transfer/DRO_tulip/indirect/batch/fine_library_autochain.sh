#!/bin/zsh
# FINE_LIBRARY_AUTOCHAIN  Wait for the 24 x 24 rib campaign to drain, then run
# the finish job (package, audit, sweep, movie). Two conditions, both needed:
#   * all 19 certified columns are on disk, and
#   * no old explicit-range rib worker is alive -- they hold no queue claims,
#     so the queue cannot see them (Astra pass 2, N).
# If the workers are all gone but columns are missing, it does NOT package a
# partial library: it writes BLOCKED and stops, for a human to relaunch.
ROOT=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect
R=$ROOT/results_fine
MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
LOG=$R/finish_autochain.log
WANT=19
rm -f $R/FINISH_VERDICT.txt
echo "$(date) autochain armed: waiting for $WANT columns and no rib worker" >> $LOG
while true; do
  n=$(ls $R/fine_rib_col*.mat 2>/dev/null | wc -l | tr -d ' ')
  live=$(pgrep -f fine_ribs_range_job | wc -l | tr -d ' ')
  if [ "$n" -ge $WANT ] && [ "$live" -eq 0 ]; then break; fi
  if [ "$live" -eq 0 ] && [ "$n" -lt $WANT ]; then
    echo "$(date) BLOCKED: no rib worker alive but only $n/$WANT columns; finish NOT started" >> $LOG
    echo "BLOCKED: $n/$WANT rib columns and no worker alive $(date)" > $R/FINISH_VERDICT.txt
    exit 1
  fi
  sleep 300
done
echo "$(date) $n columns on disk, no rib worker alive: launching the finish job" >> $LOG
"$MATLAB" -batch "run('$ROOT/batch/fine_library_finish_job.m')" > $R/finish_job.out 2>&1
rc=$?
echo "$(date) finish job exited $rc" >> $LOG
[ -f $R/FINISH_VERDICT.txt ] || echo "FINISH JOB DIED WITHOUT A VERDICT (exit $rc) $(date)" > $R/FINISH_VERDICT.txt
