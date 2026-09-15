#!/bin/zsh
# Wait for the running (pre-chaining) improve pass to finish, relaunch the chaining version once, then watch it.
setopt null_glob
SP=$SCRATCH
L=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results_fine_v6/fill_holes_direct.log
O1=$SP/fill_holes_job.improve1.out; O2=$SP/fill_holes_job.improve2.out
until grep -q "FILL HOLES DONE" $O1 2>/dev/null || ! pgrep -f "MATLAB.*improve_job.m" >/dev/null; do sleep 60; done
echo "=== $(date +%H:%M) pass 1 ended: $(grep -h 'FILL HOLES DONE' $O1 | cut -c1-80) ==="
cd $SP; nohup /Applications/MATLAB_R2026a.app/bin/matlab -batch "run('$SP/improve_job.m')" > $O2 2>&1 &
echo "chaining improve pass pid $!"; sleep 120
start=$(date +%s)
while true; do
  ev=""
  grep -q "FILL HOLES DONE" $O2 2>/dev/null && ev="DONE"
  grep -q "^Error\|fatal error" $O2 2>/dev/null && ev="$ev; ERROR"
  pgrep -f "MATLAB.*improve_job.m" >/dev/null || ev="$ev; job gone"
  [ $(( $(date +%s) - $(stat -f %m $L) )) -gt 2700 ] && ev="$ev; log silent 45 min"
  [ -n "$ev" ] && { echo "EVENT: $ev"; break; }
  [ $(( $(date +%s) - start )) -gt 5400 ] && { echo "TICK: 90 min"; break; }
  sleep 180
done
echo "=== $(date +%H:%M) ==="; grep -c "CERTIFIED" $L; grep "cell (\|queued\|DONE" $L | tail -6 | cut -c1-170; grep -B2 -A8 "^Error" $O2 2>/dev/null | head -12 | cut -c1-150
