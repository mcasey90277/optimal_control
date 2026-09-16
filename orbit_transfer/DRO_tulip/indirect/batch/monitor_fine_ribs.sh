#!/bin/zsh
# MONITOR_FINE_RIBS  Status of every departure-rib worker.
# LIVENESS IS LOG AGE, not a process match: the range workers carry their tag
# in the ENVIRONMENT, which never appears in the command line, so pgrep on the
# tag reported every one of them dead while all four were certifying points.
OUT=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results_fine
echo "=== $(date +%H:%M:%S) ==="
done_=$(ls $OUT/fine_rib_col*.mat 2>/dev/null | wc -l | tr -d ' ')
echo "columns complete: ${done_} / 19   (matlab procs: $(pgrep -f MATLAB_maca64 | wc -l | tr -d ' '))"
for t in main B C D; do
  if [ "$t" = "main" ]; then L=$OUT/fine_ribs_diary.log; V=$OUT/RIBS_VERDICT.txt
  else L=$OUT/fine_ribs_diary_$t.log; V=$OUT/RIBS_VERDICT_$t.txt; fi
  if [ ! -f $L ]; then echo "  $t: no log"; continue; fi
  age=$(( $(date +%s) - $(stat -f %m $L) ))
  col=$(grep -E "^=== column" $L | tail -1 | sed -E 's/=== column ([0-9]+).*/\1/')
  pts=$(grep -c ", certified" $L)
  if [ -f $V ]; then st="DONE: $(cat $V | cut -c1-60)"
  elif [ $age -gt 1800 ]; then st="*** STALLED ${age}s ***"
  else st="alive (${age}s)"; fi
  echo "  $t: col ${col:-?}, ${pts} points certified, $st"
done
