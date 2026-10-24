#!/bin/zsh
# Launch drone_table1 in a clean MATLAB, detached, logging to results/.
# Usage: ./run_table1.sh        (then tail -f ../results/table1_batch.log)
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$HERE/../results"
nohup /Applications/MATLAB_R2026a.app/bin/matlab -batch "run('$HERE/job_table1.m')" \
    > "$HERE/../results/table1_batch.stdout" 2>&1 &
echo "launched, pid $!; log: $HERE/../results/table1_batch.log"
