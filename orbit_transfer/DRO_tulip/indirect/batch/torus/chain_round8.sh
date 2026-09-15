#!/bin/zsh
# Wait for improve pass 4 to finish, then re-package round 6 (round 8: with the improved direct rib), audit, sweep (merged sidecar); then watch.
setopt null_glob
SP=$SCRATCH
I=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect
O4=$SP/fill_holes_job.improve4.out
until grep -q "FILL HOLES DONE" $O4 2>/dev/null || ! pgrep -f "MATLAB.*improve_job" >/dev/null; do sleep 60; done
echo "=== $(date +%H:%M) improve pass 4 ended: $(grep -h 'FILL HOLES DONE' $O4 | cut -c1-80) ==="
mv $I/results_fine_v6/repackage_job.out $I/results_fine_v6/repackage_job.round7.out
cd $I/results_fine_v6
nohup /Applications/MATLAB_R2026a.app/bin/matlab -batch "run('$SP/v6_repackage_job.m')" > $I/results_fine_v6/repackage_job.out 2>&1 &
echo "round 8 re-package pid $!"; sleep 120; exec $SP/watch_v6repack.sh
