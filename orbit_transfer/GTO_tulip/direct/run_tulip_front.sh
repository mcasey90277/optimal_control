#!/bin/zsh
# RUN_TULIP_FRONT  Sweep the Delta-V / t_f front via the validated front door.
#
# The tulip analogue of run_cr3bp_ladder.sh -- but it sweeps t_f, NOT thrust,
# and that difference is deliberate, not an omission.
#
# The earth and CR3BP campaigns sweep thrust because thrust is the dimension
# that works for them. For the tulip it is not: the 20 mN pilot rung
# (pilot_rung_20mN.m) was an honest FAILURE against a fixed-tau_f topology wall
# -- a chained 25 mN winding cannot grow to the rev count a lower-thrust rung
# needs, so the eps=1 energy step will not close. That is an open campaign
# problem (see process/LADDER_PREP_PILOT_FINDINGS.md). run_gto_tulip refuses
# off-nominal thrust outright, so a thrust sweep here would be a script that
# fails on every rung.
#
# What DOES work for this campaign is the t_f factor, and the Delta-V/t_f front
# is the campaign's mapped result.
#
# THIS SWEEP PRODUCES AN UPPER BOUND ON Delta-V, NOT THE FRONT. It runs ONE
# seed route (the energy backbone) per factor. Measured 2026-07-26, distinct
# seed routes at the SAME t_f converge to genuinely different extremals -- all
# eps=0 and machine-tight -- differing by 0.02% to 9.84% in Delta-V (worst case
# factor 1.450: 43 switches vs 26). The published front is the lower envelope
# over ALL routes, which aggregate_front computes as min(dV) per factor. So:
#   * use this sweep to POPULATE the front cheaply;
#   * then run aggregate_front, which takes the best point per factor;
#   * do not quote a single-route sweep as the front. Energy backbones exist for factors
# 1.12 - 1.95; the band 1.01-1.11 approaching min-time is the other open
# problem and is deliberately NOT in the default sweep.
#
# Each factor runs in its OWN MATLAB process (crash isolation, the elfo_batch
# lesson). Sweep-level resume: a factor whose result .mat already exists is
# skipped. A failure is RECORDED and the sweep continues -- an uncertified
# factor is information about the front, not a reason to stop. Per-factor logs
# and a summary land in results/front_logs/.
#
# Usage:  ./run_tulip_front.sh [factors...]
#         ./run_tulip_front.sh                 # default sweep, min-fuel
#         ./run_tulip_front.sh 1.15 1.20       # just these two
#         EPSMIN=1 ./run_tulip_front.sh        # min-ENERGY sweep instead
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
RES="$HERE/results";  LOGD="$HERE/results/front_logs";  mkdir -p "$LOGD" "$RES"
MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab
EPSMIN="${EPSMIN:-0}"

FACTORS=("$@")
if [ ${#FACTORS[@]} -eq 0 ]; then
  FACTORS=(1.12 1.15 1.20 1.25 1.30 1.40 1.50 1.60 1.70 1.80 1.85 1.95)
fi

SUMMARY="$LOGD/front_summary.txt"
echo "tulip Delta-V/t_f front $(date '+%F %T')  epsMin=$EPSMIN  factors: ${FACTORS[*]}" > "$SUMMARY"

for F in "${FACTORS[@]}"; do
  # Canonical name, matching minfuel_config's fname(): factor in milli-units.
  # Must match run_gto_tulip's own naming: psr_direct_f####_minEps<eTag>.mat
  MILLI=$(printf '%.0f' "$(echo "$F * 1000" | bc -l)")
  ETAG=$(echo "$EPSMIN" | tr '.' 'p')
  TAG=$(printf 'psr_direct_f%04d_minEps%s' "$MILLI" "$ETAG")

  if [ -f "$RES/$TAG.mat" ]; then
    echo "SKIP  f=$F ($TAG.mat exists)" | tee -a "$SUMMARY"; continue
  fi

  echo "SOLVE f=$F -> $TAG  ($(date '+%T'))" | tee -a "$SUMMARY"
  "$MATLAB" -batch "cd('$HERE'); setup_paths; TULIP_OVERRIDES=struct('factor',$F,'epsMin',$EPSMIN,'movieMode','none'); run_gto_tulip" \
      > "$LOGD/$TAG.log" 2>&1
  RC=$?

  # A run is a PASS only if the front door certified it AND wrote the artifact.
  # run_gto_tulip refuses to write an uncertified result, so the file's
  # existence is the gate -- rc=0 alone is not enough (an uncertified run warns
  # and exits cleanly).
  if [ $RC -eq 0 ] && [ -f "$RES/$TAG.mat" ]; then
    RESULT=$(grep -aE "delta-V|switches" "$LOGD/$TAG.log" | tr -s ' ' | tr '\n' ' ')
    echo "PASS  f=$F  $RESULT" | tee -a "$SUMMARY"
  else
    WHY=$(grep -aE "Error|not certified|stalled|no energy backbone" "$LOGD/$TAG.log" | tail -2 | tr '\n' ' ')
    echo "FAIL  f=$F  rc=$RC  $WHY" | tee -a "$SUMMARY"
  fi
done

echo "FRONT SWEEP DONE $(date '+%F %T')" | tee -a "$SUMMARY"
echo
echo "Summary: $SUMMARY"
echo "Aggregate the certified rows into a front with: aggregate_front.m"
echo "NOTE: this was a SINGLE-route sweep -- an upper bound on dV, not the"
echo "      envelope. aggregate_front takes min(dV) per factor over all routes."
