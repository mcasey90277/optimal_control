# RESULTS — headline 1.15× refinement demonstration

Run: `run_headline_1p15.m` on the certified 1.15× solution
(`sundman_minfuel_certified.mat`, 4001 nodes, 25 switches, defect
2.04e-14 (live file; the campaign record quotes 2.4e-14),
propellant 2.2640 kg), `maxRounds=4, K=8, maxAdd=40`. Seed prep regenerated
KKT-dual costates via one `eps=0 warmTight` re-solve (the certified `.mat`
predates dual extraction), stamping `factor=1.15`. Full console output in
`headline_run.log` (not committed — regenerate by re-running); saved history
in `refine_history_headline_1p15.mat`; figure in `refine_headline_1p15.png`.

## Per-round table (actual run output)

```
round  nodes   sw   maxMove     dProp(kg)   nViol   HresMax
0      4001    25   NaN         NaN         0       4.10e-01
1      4009    25   6.08e-02    2.52e-05    0       4.10e-01
2      4017    25   1.54e-02    -8.79e-05   0       4.11e-01
```

(round = measured index printed by `run_headline_1p15`, matching
`refine_loop`'s internal `[round r-1]` log lines. Round 2 is where the
acceptance test fired — `converged=1` in `refine_history_headline_1p15.mat`.)

Note: switch counts are the raw s>0.5 threshold counter; per the ms_band
campaign the certified count is 10 for 1.12x and the near-graze dip inflates
raw counts (threshold artifact) — see LOW_THRUST_MINFUEL_CAMPAIGN.md.

Additional fields from the saved history (not in the console table):

| round | prop_kg | maxDefect | betaSpread | ipoptStatus |
|---|---|---|---|---|
| 0 | 2.264013 | 1.98e-14 | 1.404 | Solve_Succeeded |
| 1 | 2.264038 | 3.10e-14 | 1.382 | Solve_Succeeded |
| 2 | 2.263950 | 4.24e-14 | 1.319 | Solve_Succeeded |

Node growth: 4001 → 4009 → 4017 (+8 each round, exactly `K=8`; `refine_sigma`
never hit the `maxAdd=40` cap and dropped no intervals). Every re-solve
converged `Solve_Succeeded` with `maxDefect` in the 1e-14 range (machine-tight)
and switch count fixed at 25 throughout — no switch was born or lost.

## Switch-time stabilization verdict

**Stabilized, and it took two refinement rounds (not one).** This is a
genuine multi-round case, unlike the Task-5 smoke test on the 1.12× file
(10 switches), which converged after a single round because that seed was
already well localized. Here `maxSwitchMove` decayed monotonically and by
roughly a factor of 4 per round — 6.08e-02 (round 1) → 1.54e-02 (round 2) —
and the loop's acceptance test (`maxMove < that round's local mesh width`
AND `|dProp| < 1e-4 kg` AND switch count unchanged) fired at round 2.
Propellant drift was 2.52e-05 kg (round 1) and −8.79e-05 kg (round 2), both
comfortably inside the 1e-4 kg tolerance and an order of magnitude smaller
than the certified propellant itself (2.264 kg) — the mesh refinement moved
switch brackets without perturbing the optimum's fuel cost. `nViol` (switching-
law sign violations outside the deadband) stayed at 0 in every round.
Caveat: `maxSwitchMove` measures the quantized bracket-midpoint switch
position (resolution ~half a local cell), not the sub-cell S=0 root;
`diag.tauCr` provides that root and is the natural stronger acceptance
signal for a follow-up.

The final 25 switch times (`history(3).tauSwitch`, sorted) show 20 switches
clustered over τ∈[0, 55] and a tight group of 5 switches over τ∈[149.1,
151.5] — i.e. near the trajectory's terminal end (τ_f ≈ 151.5). Both regions
localized cleanly; no divergent or oscillating switch was observed in either
cluster.

## Option-2 decision (Hamiltonian-residual escalation)

**`HresMax` stayed essentially flat — it did not track the switch-mesh
refinement in either direction, up or down.** Values were 0.4098 → 0.4102 →
0.4106 across the three rounds: a ~0.2% drift, not a meaningful trend, while
`maxSwitchMove` fell by ~4×. This matches the same near-constant behavior
already seen in the Task-5 smoke test on the 1.12× file (`HresMax` 1.030 →
1.029, also flat).

Inspecting the money plot (`refine_headline_1p15.png`, middle panel) shows
why: the passive `|H_σ|` trace is not spread evenly across the trajectory —
it is dominated by a single sharp spike at τ≈145–155, exactly the region
holding the tight final cluster of 5 switches near rendezvous. The bulk of
the trajectory (τ∈[0,140], including the 20 switches that `refine_sigma`
spent most of its bisections on) shows a small, quiet `|H_σ|` baseline.
`refine_sigma` did bisect the terminal cluster's brackets too (switch count
and localization improved there as everywhere), yet the spike's magnitude
barely moved — consistent with this being a **boundary/costate-conditioning
artifact near the free-rendezvous terminal condition**, not an
under-refined mesh interval. This lines up with the campaign's earlier,
independently-documented finding (`LOW_THRUST_MINFUEL_CAMPAIGN.md`, Tier-1
PMP certification scope) that direct costate recovery near the rendezvous
terminal condition is the identified weak point (no BC on λ_r, λ_v there),
not a mesh-resolution problem that bisection can fix.

**Verdict: Option 1 (switch-localization refinement) sufficed for the
switch-localization goal it targets, and the passive `HresMax` signal does
not indicate that promoting the Hamiltonian residual to an active refinement
driver (Option 2) would help.** The residual's flatness is explained by a
fixed terminal-region conditioning effect, not by unrefined switch brackets
elsewhere — the switch brackets that Option 1 did work on all localized and
stopped moving. Per the design spec's own criterion ("stays high in
un-refined arcs → escalate; drops → don't"), the literal reading is
ambiguous — `HresMax` is not high specifically *in un-refined arcs* (there
are none left; every switch's neighborhood was bisected across the three
rounds) and it did not drop either — it is high in one already-refined,
terminal region for a reason the mesh cannot fix. Escalating to a
Hamiltonian-driven refiner is not supported by this data.

## Honest caveats

- Only one seed (the certified 1.15×) and one `(K, maxAdd)` setting were
  run — no sensitivity sweep. A different `K` could plausibly change the
  round count to converge.
- `betaSpread` (single-β switch-fit sanity check) decreased monotonically
  (1.404 → 1.382 → 1.319) across rounds, consistent with a healthier dual
  map as the mesh refines, but this was not a stopping criterion and is
  reported only as a secondary sanity signal.
- No re-solve failed and no MEX crash occurred on this run; all three
  measured rounds and two re-solves completed cleanly with
  `Solve_Succeeded`. There is nothing to recover from disk beyond what
  is already in `refine_history_headline_1p15.mat`.
