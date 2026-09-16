#!/bin/zsh
# RUN_FINE_RIBS  Launch fine_ribs_job.m in a clean MATLAB -batch session under an
# OS watchdog. Log to files only (matlab -batch buffers stdout). Completion is
# results_fine/RIBS_VERDICT.txt, never the log. Usage: nohup batch/run_fine_ribs.sh &
ROOT=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect
OUT=$ROOT/results_fine
MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
[ -x $MATLAB ] || MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab
WATCHDOG_SEC=${1:-34200}                    # 9.5 h: LONGER than the job's own 8 h
                                            # budget, so the clean stop wins over the kill
mkdir -p $OUT
rm -f $OUT/RIBS_VERDICT.txt
if pgrep -fl fine_ribs_job > /dev/null; then echo "a fine_ribs_job is already running; refusing"; exit 1; fi
echo "LAUNCH $(date)  watchdog ${WATCHDOG_SEC}s  matlab $MATLAB" >> $OUT/fine_ribs_driver.log
$MATLAB -batch "run('$ROOT/batch/fine_ribs_job.m')" > $OUT/fine_ribs_batch.out 2>&1 &
MPID=$!
( sleep $WATCHDOG_SEC; if kill -0 $MPID 2>/dev/null; then kill $MPID; echo "WATCHDOG KILLED $(date)" >> $OUT/fine_ribs_driver.log; echo "CHAIN KILLED BY WATCHDOG  $(date)" > $OUT/RIBS_VERDICT.txt; fi ) &
WPID=$!
wait $MPID
RC=$?
kill $WPID 2>/dev/null
echo "MATLAB EXIT $RC $(date)" >> $OUT/fine_ribs_driver.log
[ -f $OUT/RIBS_VERDICT.txt ] || echo "CHAIN EXITED WITHOUT VERDICT (rc $RC)  $(date)" > $OUT/RIBS_VERDICT.txt
echo "DRIVER DONE $(date)" >> $OUT/fine_ribs_driver.log
