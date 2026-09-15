#!/bin/zsh
setopt null_glob
SP=$SCRATCH
L=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results_fine_v6/fill_holes_direct.log; O=$SP/fill_holes_job.improve4.out; start=$(date +%s)
sleep 120
while true; do
  ev=""
  grep -q "FILL HOLES DONE" $O 2>/dev/null && ev="DONE"
  grep -q "^Error\|fatal error" $O 2>/dev/null && ev="$ev; ERROR"
  pgrep -f "MATLAB.*improve_job" >/dev/null || ev="$ev; job gone"
  [ $(( $(date +%s) - $(stat -f %m $L) )) -gt 2700 ] && ev="$ev; log silent 45 min"
  [ -n "$ev" ] && { echo "EVENT: $ev"; break; }
  [ $(( $(date +%s) - start )) -gt 5400 ] && { echo "TICK: 90 min"; break; }
  sleep 180
done
echo "=== $(date +%H:%M) ==="; grep -c "CERTIFIED" $L; grep "cell (\|queued\|improve:\|DONE" $L | tail -6 | cut -c1-170; grep -B2 -A8 "^Error" $O 2>/dev/null | head -12 | cut -c1-150
