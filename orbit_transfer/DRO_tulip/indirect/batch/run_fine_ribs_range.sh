#!/bin/zsh
# RUN_FINE_RIBS_RANGE  One rib worker over a column range, as its own MATLAB.
# Usage: nohup batch/run_fine_ribs_range.sh <tag> "<matlab column vector>" &
ROOT=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect
OUT=$ROOT/results_fine
MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
TAG=$1; COLS=$2; SEC=${3:-34200}
mkdir -p $OUT
rm -f $OUT/RIBS_VERDICT_$TAG.txt
if pgrep -f "RIB_TAG=$TAG" > /dev/null; then echo "worker $TAG already running"; exit 1; fi
echo "LAUNCH $TAG cols=$COLS $(date)" >> $OUT/fine_ribs_driver.log
RIB_TAG=$TAG RIB_COLS=$COLS $MATLAB -batch "run('$ROOT/batch/fine_ribs_range_job.m')" \
    > $OUT/fine_ribs_batch_$TAG.out 2>&1 &
MPID=$!
( sleep $SEC; if kill -0 $MPID 2>/dev/null; then kill $MPID; echo "WATCHDOG killed $TAG $(date)" >> $OUT/fine_ribs_driver.log; fi ) &
WPID=$!
wait $MPID; kill $WPID 2>/dev/null
echo "WORKER $TAG EXIT $? $(date)" >> $OUT/fine_ribs_driver.log
