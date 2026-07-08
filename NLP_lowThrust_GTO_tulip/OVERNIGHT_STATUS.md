# Overnight run — energy->fuel seeding + tf-continuation (Jul 7-8 2026)

> **SUPERSEDED — historical.** This is a point-in-time status from the
> overnight tf-continuation. The full campaign (through cone-elimination,
> Hermite-Simpson, and CasADi+IPOPT) and the current conclusions live in
> `LOW_THRUST_MINFUEL_CAMPAIGN.md`. Kept for the tf-continuation how-to below.

## The question
Can the min-ENERGY solution seed the min-FUEL direct method so it converges
on the FULL 40-rev spiral (which min-fuel could not do from any earlier warm
start)? And can tf-continuation then grow the switch structure toward the
research-grade many-switch regime?

## Result 1 (DONE): energy seeds fuel -> YES, cleanly.
Fed the min-energy full-spiral solution (tf = 1.15x min, 93.5% on / 6.5% off,
defect 3e-4) straight into the min-fuel NLP as its warm start. Min-fuel
CONVERGED:
- **flag = 2, max defect = 2.2e-15 (machine zero)** -- better than min-energy's
  own 3e-4 floor on the same mesh (the near-bang-bang min-fuel optimum has the
  throttle pinned on its bounds, so the collocation defects drive clean to 0).
- propellant 2.9503 kg, **3 switches**, s in [0.020, 1.000].
This is the clean win: the energy->fuel warm start works where the burn+coast
and time-stretch starts failed. Saved: `minfuel_from_energy_seed.mat`.
(3 switches only, because tf = 1.15x is a tight time budget -- mild bang-bang.)

## Result 2 (RUNNING overnight): tf-continuation to many switches.
`tf_continuation_minfuel.m` -- launched `tf_continuation_minfuel(2500)`.
Log: `tf_continuation.log`.  Incremental results: `tf_continuation_results.mat`.

What it does: FIXED tulip target (the true rendezvous point, NOT the
phase-shifted coast terminus), start from the min-time solution at tf = tfMin
(0 switches, always burn), step tf up through
[1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.50 1.60 1.75 1.90 2.00] x tfMin.
Each step warm-starts from the previous converged solution, so the burn/coast
switch structure grows one relaxation at a time. With a fixed endpoint and
more time, the min-fuel optimum coasts near each apogee and burns near each
perigee -> switch count should climb from 0 toward the ~80-switch
research-grade regime as tf increases.

Robustness (built for a long unattended run):
- per-step fmincon iteration cap (1500) so no single step runs away;
- defect guard: a step with defect > 5e-3 does NOT become the next warm start
  (prevents a bad step from poisoning the chain);
- incremental `save` after EVERY step -- so even a partial run leaves the
  switches-vs-tf trend in `tf_continuation_results.mat`.

## How to read the results in the morning
```matlab
load tf_continuation_results.mat        % struct array `results`
[[results.factor]' [results.switches]' [results.mProp_kg]' [results.maxDefect]']
```
Columns: tf/tfMin, switch count, propellant (kg), max defect. The story is the
switches column climbing with tf. `tail tf_continuation.log` shows live
progress (one line per factor).

## If it broke early
Check `tf_continuation.log` for the failing factor. Most likely culprit is a
step that couldn't restore feasibility after a tf jump (defect stays > 5e-3,
logged as "(loose)"). The saved `results` still holds everything up to that
point. Resuming from the last good `Zgood` (also saved) is the fix.

## Expected headline
A monotone(ish) rise in switch count with tf -- concrete evidence that (a)
min-energy is the right homotopy root to make min-fuel tractable on the full
spiral, and (b) tf-continuation walks it into the many-switch regime that
defeats a single-shot solve. That closes the loop opened when min-fuel first
failed on the full spiral.
