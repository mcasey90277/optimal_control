#!/bin/bash
# RUN_TASK9_WATCHDOG  Auto-relaunching wrapper for run_task9_deep.m.
#
# The MEX/volume-linked MATLAB crash is recurring at this problem scale
# (2 confirmed fatal crashes so far, both during the 0.5 N anchor's Stage B
# continuation). run_task9_rung.m / run_mintime_mee.m / run_transfer_mee.m /
# psr_mee_refine.m are all resume-safe (per-round/per-stage caching), so a
# crashed MATLAB process can simply be relaunched and will pick up exactly
# where it left off. This wrapper automates that relaunch instead of
# requiring a manual intervention on every single crash, while still
# stopping (not masking) on a GENUINE solver failure (an assertion like
# anchorUncertified/fuelUncertified/allStagesFailed) rather than looping
# forever on a deterministic block.
cd /Users/msc/Desktop/optimal_control/earth_elliptic_to_geo || exit 1
LOG=results/task9_deep_run.log
MAXATTEMPTS=200

for i in $(seq 1 $MAXATTEMPTS); do
    echo "=== WATCHDOG ATTEMPT $i at $(date) ===" >> "$LOG"
    /Applications/MATLAB_R2025b.app/bin/matlab -batch \
        "cd('/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo'); run_task9_deep" >> "$LOG" 2>&1
    RC=$?
    echo "=== WATCHDOG ATTEMPT $i EXIT CODE $RC at $(date) ===" >> "$LOG"

    if grep -q "TASK 9 DEEP LADDER COMPLETE" "$LOG"; then
        echo "WATCHDOG_DONE: full ladder completed after $i attempt(s)" >> "$LOG"
        break
    fi
    if tail -n 200 "$LOG" | grep -qE "anchorUncertified|fuelUncertified|allStagesFailed"; then
        echo "WATCHDOG_BLOCKED_GENUINE_FAILURE: a solver assertion fired (not a crash) -- stopping, needs attention" >> "$LOG"
        break
    fi
    sleep 5
done
echo "WATCHDOG_LOOP_EXITED after $i attempt(s)" >> "$LOG"
