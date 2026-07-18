#!/bin/bash
# RUN_TASK9_FUELFIRST_WATCHDOG  Single-retry auto-relaunching wrapper for
# run_task9_fuelfirst.m (controller directive, 2026-07-18: "relaunch once
# per crash, resume via caches"). Unlike run_task9_watchdog.sh (which loops
# up to 200 attempts for the anchor stage's recurring crash), this one is
# deliberately capped at 2 total attempts (1 initial + 1 relaunch) -- the
# fuel-first path is expected to be a single long solve, not a multi-round
# grind, so a crash-loop guard isn't needed; a second consecutive crash is
# treated as a genuine blocker to report, not retried further.
cd /Users/msc/Desktop/optimal_control/earth_elliptic_to_geo || exit 1
LOG=results/task9_fuelfirst_run.log
MAXATTEMPTS=2

for i in $(seq 1 $MAXATTEMPTS); do
    echo "=== FUELFIRST WATCHDOG ATTEMPT $i at $(date) ===" >> "$LOG"
    /Applications/MATLAB_R2025b.app/bin/matlab -batch \
        "cd('/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo'); run_task9_fuelfirst" >> "$LOG" 2>&1
    RC=$?
    echo "=== FUELFIRST WATCHDOG ATTEMPT $i EXIT CODE $RC at $(date) ===" >> "$LOG"

    if grep -q "TASK 9 FUEL-FIRST COMPLETE" "$LOG"; then
        echo "FUELFIRST_WATCHDOG_DONE: completed after $i attempt(s)" >> "$LOG"
        break
    fi
    if grep -q "BLOCKED AT FUEL STAGE" "$LOG"; then
        echo "FUELFIRST_WATCHDOG_BLOCKED_GENUINE: fuel solve did not certify -- stopping, needs attention" >> "$LOG"
        break
    fi
    if [ "$i" -ge "$MAXATTEMPTS" ]; then
        echo "FUELFIRST_WATCHDOG_EXHAUSTED: $MAXATTEMPTS attempts used, still no completion/genuine-block marker (likely repeated crash)" >> "$LOG"
        break
    fi
    sleep 5
done
echo "FUELFIRST_WATCHDOG_LOOP_EXITED after $i attempt(s)" >> "$LOG"
