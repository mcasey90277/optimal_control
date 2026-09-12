#!/bin/zsh
# MONITOR_RESWEEP  One-shot status: done-marker, process death, progress, STALL.
OUT=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results_resweep
LOG=$OUT/resweep_diary.log
if [ -f $OUT/VERDICT.txt ]; then echo "VERDICT: $(cat $OUT/VERDICT.txt)"; fi
if pgrep -fl resweep_job > /dev/null; then echo "process: RUNNING"; else echo "process: not running"; fi
if [ -f $LOG ]; then
  age=$(( $(date +%s) - $(stat -f %m $LOG) ))
  echo "log age ${age}s, $(wc -l < $LOG) lines; stage lines:"
  grep -E "^[0-9]\. |sweep census|written back|CHAIN" $LOG | tail -8
  echo "sweep entries measured: $(grep -c "interior [0-9]*, cand" $LOG)"
  [ $age -gt 1200 ] && ! [ -f $OUT/VERDICT.txt ] && echo "*** STALLED: no log line for ${age}s ***"
else
  echo "no diary yet"
fi
tail -3 $OUT/resweep_batch.out 2>/dev/null
