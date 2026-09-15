#!/bin/zsh
setopt null_glob
# Watch the round-5 entry job (sheet rebuild -> reuse -> launch): exit on LAUNCH DONE, an error, the job dying, a silent log > 40 min, or a 60-min tick.
SP=$SCRATCH
V5=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results_fine_v5
O=$SP/v5_campaign_job.out; start=$(date +%s)
while true; do
  ev=""
  grep -q "ROUND 5 LAUNCH DONE" $O 2>/dev/null && ev="LAUNCH DONE"
  grep -q "^Error\|error:\|fatal error" $O 2>/dev/null && ev="$ev; ERROR in job"
  pgrep -f "v5_campaign_job.m" >/dev/null || ev="$ev; job process gone"
  L=$(ls -t $V5/*.log 2>/dev/null | head -1); [ -n "$L" ] && [ $(( $(date +%s) - $(stat -f %m $L) )) -gt 2400 ] && ev="$ev; log silent"
  [ -n "$ev" ] && { echo "EVENT: $ev"; break; }
  [ $(( $(date +%s) - start )) -gt 3600 ] && { echo "TICK: 60 min"; break; }
  sleep 180
done
echo "=== $(date +%H:%M) ==="
tail -4 $V5/round5.log 2>/dev/null | cut -c1-150
grep -h "candidates\|certified\|ROUND 5\|CHANGED\|columns to walk\|Error" $O 2>/dev/null | tail -8 | cut -c1-150
