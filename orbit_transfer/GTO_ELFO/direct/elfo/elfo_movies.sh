#!/bin/zsh
# elfo_movies.sh — render control movies (or preview stills) from saved GTO->ELFO
# min-fuel solutions, POST-HOC (no re-solve).
#
# The solve batch (elfo_batch.sh) writes solution .mat files but does NOT render
# movies (to keep the long sweep fast). This renders them afterward from those
# files via elfo_render_movies -> elfo_movie. Pure MATLAB rendering (no CasADi),
# so no MEX-crash risk -- one process, looped.
#
# ============================ HOW TO RUN ============================
#   cd /Users/msc/Desktop/optimal_control/orbit_transfer/GTO_tulip/elfo
#   chmod +x elfo_movies.sh                  # once
#   ./elfo_movies.sh all                     # a movie for every banked solution
#   ./elfo_movies.sh 1.20 1.65               # just these factors
#   MODE=preview ./elfo_movies.sh all        # 3 stills/factor (seconds) instead of MP4+GIF
#
# Optional env vars:
#   MODE=movie|preview   [default movie]  (preview = quick stills to eyeball)
#   MATLAB_BIN=/path/to/matlab             [R2025b]
#
# NOTE: a full 'movie' render is ~10-15 min/factor. Best run AFTER the solve
# batch finishes (or concurrently, accepting some slowdown). 'preview' is seconds.
# ===================================================================

set -u
DIR=$(cd "$(dirname "$0")" && pwd)
MAT=${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}
MODE=${MODE:-movie}
LOGDIR=$DIR/results/logs;  mkdir -p "$LOGDIR"
LOG=$LOGDIR/elfo_movies_$(date +%Y%m%d_%H%M%S).log

[ $# -lt 1 ] && { echo "usage: $0 [all | <factor1> <factor2> ...]   (MODE=preview for stills)"; exit 1; }

if [ "$1" = "all" ]; then
  sel="'all'"
else
  sel="[$*]"                                 # numeric factor list -> MATLAB row vector
fi

echo "=== ELFO MOVIES  mode=$MODE  sel=$sel  ($(date +%H:%M:%S)) ===" | tee -a "$LOG"
"$MAT" -batch "cd('$DIR'); setup_paths; elfo_render_movies($sel, '$MODE');" 2>&1 | tee -a "$LOG"
echo "log: $LOG"
