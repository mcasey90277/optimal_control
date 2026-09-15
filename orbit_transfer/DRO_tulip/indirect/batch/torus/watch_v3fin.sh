#!/bin/zsh
V3=/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results_fine_v3
O=$V3/finalize_job_rerun.out; start=$(date +%s)
until grep -q "FINALIZE: state\|^Error\|error:" $O 2>/dev/null || ! pgrep -f "results_fine_v3/finalize_job.m" >/dev/null || [ $(( $(date +%s) - $(stat -f %m $O) )) -gt 1800 ] || [ $(( $(date +%s) - start )) -gt 3600 ]; do sleep 120; done
echo "=== $(date +%H:%M) v3 finalizer: age $(( $(date +%s) - $(stat -f %m $O) ))s; audited $(grep -c '^AUDIT (' $O) ==="
grep -h "3b\.\|^4\.\|^5\.\|AUDIT SUMMARY\|^7\.\|^8\.\|CHAIN\|FINALIZE\|BLOCKERS" $O | tail -6 | cut -c1-140
grep -B2 -A6 "^Error\|error:" $O | grep -v "^\s*$" | head -10 | cut -c1-140
