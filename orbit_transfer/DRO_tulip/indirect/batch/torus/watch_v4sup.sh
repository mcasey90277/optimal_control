#!/bin/zsh
# Watch the round-2 supervised campaign: exit on the supervisor's verdict, a
# worker log silent > 50 min (three stage caps), the supervisor dying, or a
# 60-min tick.
V=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results_fine_v4
start=$(date +%s)
age() { [ -f "$1" ] && echo $(( $(date +%s) - $(stat -f %m "$1") )) || echo -1; }
while true; do
  now=$(date +%s); ev=""
  [ -f $V/SUPERVISOR_VERDICT.txt ] && ev="supervisor verdict"
  pgrep -f campaign_supervisor.sh >/dev/null || ev="$ev; supervisor process gone"
  for h in $V/hb/*.hb(N); do
    a=$(age "$h"); if [ $a -gt 3000 ] && ! grep -q "^done|" "$h"; then ev="$ev; $(basename $h) silent ${a}s"; fi
  done
  [ -n "$ev" ] && { echo "EVENT: $ev"; break; }
  [ $(( now - start )) -gt 3600 ] && { echo "TICK: 60 min"; break; }
  sleep 180
done
echo "=== $(date +%H:%M) ==="
tail -3 $V/supervisor.log 2>/dev/null
for w in $V/worker_*.log(N); do echo "$(basename $w | cut -c8-30): $(tail -1 $w | cut -c10-120)"; done
for h in $V/hb/*.hb(N); do echo "hb $(basename $h .hb | cut -c1-22): age $(age $h)s  $(cut -d'|' -f1,4 $h | cut -c1-60)"; done
ls $V/fine_rib_col21.mat* 2>/dev/null
tail -4 $V/SUPERVISOR_VERDICT.txt 2>/dev/null
