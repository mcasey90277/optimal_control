#!/bin/zsh
# MONITOR_FINE_SHEET  One-shot status of the 24-level arrival sheet build.
# Reports done-marker, process death, progress, and STALL (log age), so that
# silence can never be mistaken for success.
OUT=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results_fine
LOG=$OUT/fine_sheet_diary.log
echo "=== $(date +%H:%M:%S) ==="
[ -f $OUT/VERDICT.txt ] && echo "VERDICT: $(cat $OUT/VERDICT.txt)"
if pgrep -f fine_sheet_job > /dev/null; then echo "process: RUNNING"; else echo "process: not running"; fi
if [ -f $LOG ]; then
  age=$(( $(date +%s) - $(stat -f %m $LOG) ))
  # the diary APPENDS across runs: read only from the last JOB START, or a
  # finished run's leftovers are reported as this run's progress
  start=$(grep -n "JOB START" $LOG | tail -1 | cut -d: -f1)
  CUR=$(mktemp); sed -n "${start},\$p" $LOG > $CUR
  cert=$(grep -c "  certified$" $CUR)
  ref=$(grep -cE "crossing .* at j =" $CUR)
  echo "log age ${age}s | crossings attempted $ref | certified $cert"
  echo "failure reasons so far:"
  grep -oE "\): .*$" $CUR | sed -E 's/[0-9][0-9.e+-]*/N/g' | sort | uniq -c | sort -rn | head -5
  rm -f $CUR
  [ $age -gt 1200 ] && [ ! -f $OUT/VERDICT.txt ] && echo "*** STALLED: no log line for ${age}s ***"
else
  echo "no diary yet"
fi
tail -2 $OUT/fine_sheet_batch.out 2>/dev/null
