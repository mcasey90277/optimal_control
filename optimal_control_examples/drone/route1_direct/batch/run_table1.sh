#!/bin/zsh
# Launch drone_table1 in a clean MATLAB, detached, logging to results/.
# Usage: ./run_table1.sh        (then tail -f ../results/table1_batch.log)
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$HERE/../results"
MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
if [[ ! -x "$MATLAB" ]]; then
    MATLAB=$(command -v matlab)
    if [[ -z "$MATLAB" ]]; then echo "no MATLAB: R2026a path absent and matlab not on PATH" >&2; exit 1; fi
fi
nohup "$MATLAB" -batch "run('$HERE/job_table1.m')" \
    > "$HERE/../results/table1_batch.stdout" 2>&1 &
echo "launched, pid $!; log: $HERE/../results/table1_batch.log"
