#!/bin/zsh
# Replacement backstop for a demonstrably healthy resweep whose original
# 5 h watchdog would have fired mid-sweep. Targets the running MATLAB by PID.
MPID=$1; SEC=$2
OUT=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results_resweep
sleep $SEC
if kill -0 $MPID 2>/dev/null; then
  kill $MPID
  echo "WATCHDOG2 KILLED $MPID $(date)" >> $OUT/resweep_driver.log
  echo "CHAIN KILLED BY WATCHDOG2 (extended)  $(date)" > $OUT/VERDICT.txt
fi
