#!/bin/zsh
# Watch both direct18 arcs: exit when both DONE lines exist, a job dies, a log is silent > 30 min, or after 60 min.
SP=$SCRATCH
R=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results
start=$(date +%s)
while true; do
  ev=""; ndone=0
  for d in up dn; do
    grep -q "ARC direct18 $d DONE" $SP/arc_direct18_$d.out 2>/dev/null && ndone=$((ndone+1))
    pgrep -f "arc_direct18_$d.m" >/dev/null || grep -q "DONE" $SP/arc_direct18_$d.out 2>/dev/null || ev="$ev; $d job gone"
    L=$R/arrival_arc_direct18_${d}_long.log
    [ -f $L ] && [ $(( $(date +%s) - $(stat -f %m $L) )) -gt 1800 ] && ev="$ev; $d log silent"
  done
  [ $ndone -eq 2 ] && ev="both arcs DONE"
  [ -n "$ev" ] && { echo "EVENT: $ev"; break; }
  [ $(( $(date +%s) - start )) -gt 3600 ] && { echo "TICK: 60 min"; break; }
  sleep 300
done
echo "=== $(date +%H:%M) ==="
for d in up dn; do echo "--- $d ---"; tail -2 $R/arrival_arc_direct18_${d}_long.log 2>/dev/null | cut -c1-150; grep -h "DONE\|Error\|error:" $SP/arc_direct18_$d.out 2>/dev/null | tail -3 | cut -c1-150; done
