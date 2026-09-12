#!/bin/zsh
# RUN_RESWEEP  Launch resweep_job.m in a clean MATLAB -batch session under an
# OS watchdog. Log to files only (matlab -batch buffers stdout). Completion is
# results_resweep/VERDICT.txt, never the log. Usage: nohup batch/run_resweep.sh &
ROOT=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect
OUT=$ROOT/results_resweep
MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
[ -x $MATLAB ] || MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab
WATCHDOG_SEC=${1:-18000}                    # 5 h backstop; the run should take 2-3 h
mkdir -p $OUT
rm -f $OUT/VERDICT.txt
if pgrep -fl resweep_job > /dev/null; then echo "a resweep_job is already running; refusing"; exit 1; fi
echo "LAUNCH $(date)  watchdog ${WATCHDOG_SEC}s  matlab $MATLAB" >> $OUT/resweep_driver.log
$MATLAB -batch "run('$ROOT/batch/resweep_job.m')" > $OUT/resweep_batch.out 2>&1 &
MPID=$!
( sleep $WATCHDOG_SEC; if kill -0 $MPID 2>/dev/null; then kill $MPID; echo "WATCHDOG KILLED $(date)" >> $OUT/resweep_driver.log; echo "CHAIN KILLED BY WATCHDOG  $(date)" > $OUT/VERDICT.txt; fi ) &
WPID=$!
wait $MPID
RC=$?
kill $WPID 2>/dev/null
echo "MATLAB EXIT $RC $(date)" >> $OUT/resweep_driver.log
[ -f $OUT/VERDICT.txt ] || echo "CHAIN EXITED WITHOUT VERDICT (rc $RC)  $(date)" > $OUT/VERDICT.txt
echo "DRIVER DONE $(date)" >> $OUT/resweep_driver.log
