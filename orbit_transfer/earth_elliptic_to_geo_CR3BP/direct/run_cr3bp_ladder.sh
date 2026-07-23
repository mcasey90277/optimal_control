#!/bin/zsh
# RUN_CR3BP_LADDER  Walk the thrust ladder (Moon ON) via the validated front door.
#
# Each rung runs in its OWN MATLAB process (MEX-fatal crash isolation, the
# elfo_batch lesson) with movieMode off. Ladder-level resume: a rung whose
# final products .mat already exists is skipped. A rung failure is RECORDED
# and the ladder continues (see-how-far-we-get semantics). Per-rung logs +
# a summary table land in results/ladder_logs/.
#
# Usage:  ./run_cr3bp_ladder.sh [rungs...]   (default: 5 2.5 1 0.5 0.2 0.1)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
RES="$HERE/results";  LOGD="$RES/ladder_logs";  mkdir -p "$LOGD"
MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab
RUNGS=("$@"); [ ${#RUNGS[@]} -eq 0 ] && RUNGS=(5 2.5 1 0.5 0.2 0.1)
SUMMARY="$LOGD/ladder_summary.txt"
echo "CR3BP thrust ladder $(date '+%F %T')  rungs: ${RUNGS[*]}" > "$SUMMARY"
for T in "${RUNGS[@]}"; do
  TAG="cr3bp_T$(echo $T | tr '.' 'p')N_phi0_fuel"
  if [ -f "$RES/$TAG.mat" ]; then
    echo "SKIP  T=$T N ($TAG.mat exists)" | tee -a "$SUMMARY"; continue
  fi
  echo "RUNG  T=$T N -> $TAG  ($(date '+%T'))" | tee -a "$SUMMARY"
  "$MATLAB" -batch "cd('$HERE'); LADDER_OVERRIDES=struct('thrustN',$T,'runName','$TAG','movieMode',false); run_cr3bp_geo" \
      > "$LOGD/$TAG.log" 2>&1
  RC=$?
  if [ $RC -eq 0 ] && [ -f "$RES/$TAG.mat" ]; then
    RESULT=$(grep -a "m_f=" "$LOGD/$TAG.log" | tail -1)
    echo "PASS  T=$T N  $RESULT" | tee -a "$SUMMARY"
  else
    TAILERR=$(grep -aE "Error|stalled|gate" "$LOGD/$TAG.log" | tail -2 | tr '\n' ' ')
    echo "FAIL  T=$T N  rc=$RC  $TAILERR" | tee -a "$SUMMARY"
  fi
done
echo "LADDER DONE $(date '+%F %T')" | tee -a "$SUMMARY"
