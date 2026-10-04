#!/bin/zsh
# RUN_RECERTIFY  Launch N chunk processes of the status-layer backfill, one
# MATLAB -batch each, each with its own log, under an OS watchdog.
#   batch/run_recertify.sh N          -- recertify_chunk_job (CHUNK = 1..N)
#   batch/run_recertify.sh N audit    -- backfill_v2_job BACKFILL_STAGE=audit
# Optional 3rd argument: watchdog seconds per process [172800 = 48 h].
# BACKFILL_OUT (environment) is passed through; default results/library_70mN_24x48_v2.
# Completion is each chunk's ~/RECERT_<k>_VERDICT.txt (or
# ~/BACKFILL_V2_AUDIT_<k>_VERDICT.txt), never the log. Every chunk resumes.
# Usage: nohup batch/run_recertify.sh 4 &
ROOT=${0:A:h:h}
N=${1:?usage: run_recertify.sh N [recert|audit] [watchdog_sec]}
WHAT=${2:-recert}
SEC=${3:-172800}
OUT=${BACKFILL_OUT:-$ROOT/results/library_70mN_24x48_v2}
MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
[ -x $MATLAB ] || MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab
mkdir -p $OUT
if [ "$WHAT" = audit ]; then JOB=$ROOT/batch/backfill_v2_job.m; STAGE=audit; TAG=audit
else JOB=$ROOT/batch/recertify_chunk_job.m; STAGE=; TAG=recert; fi
if [ "$WHAT" != audit ] && [ ! -f $OUT/harvest.mat ]; then echo "no $OUT/harvest.mat: run backfill_v2_job (harvest) first"; exit 1; fi
if pgrep -f "${JOB:t}" > /dev/null; then echo "$TAG chunks already running; refusing"; exit 1; fi
for K in $(seq 1 $N); do
  (
    echo "LAUNCH $TAG $K/$N $(date)  matlab $MATLAB" >> $OUT/${TAG}_driver.log
    BACKFILL_TAG=$TAG BACKFILL_STAGE=$STAGE BACKFILL_OUT=$OUT CHUNK=$K NCHUNK=$N \
        $MATLAB -batch "run('$JOB')" > $OUT/${TAG}_${K}_batch.out 2>&1 &
    MPID=$!
    ( sleep $SEC; if kill -0 $MPID 2>/dev/null; then kill $MPID; echo "WATCHDOG killed $TAG $K $(date)" >> $OUT/${TAG}_driver.log; fi ) &
    WPID=$!
    wait $MPID; RC=$?
    kill $WPID 2>/dev/null
    echo "EXIT $TAG $K rc $RC $(date)" >> $OUT/${TAG}_driver.log
  ) &
  sleep 20                                   # stagger the JVM / pool start-ups
done
wait
echo "ALL $TAG CHUNKS DONE $(date)" >> $OUT/${TAG}_driver.log
