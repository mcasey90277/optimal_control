# Thrust-Ladder Prep (P2) — Pilot Rung Findings

**Date:** 2026-07-21/22  **Package:** `docs/superpowers/plans/2026-07-21-ladder-prep.md`
**Gate:** two warm-chained 20 mN pilot rungs (one per campaign) — the package's
validation exit gate. The infrastructure they exercise (fingerprints, opt-in
boxes + boundSat, chain helpers, hardened `Solve_Succeeded` gates) all worked;
the pilots test whether a real off-nominal rung is *reachable* by it.

## ELFO 20 mN — PASS (clean)

`certified=1`, `Solve_Succeeded`, defect **1.5e-15**, ε=0 bang-bang reached.
Structure: 11 switches, 15.7% propellant, edge 99.6%, t_f 33.46 d.
- **First certified off-nominal rung** either CR3BP campaign has produced.
- `boundSatWorst=massHi` — worst-label only, no saturation warning fired (not
  binding).
- Bonus physics: certified at ~0.99× the naive 1/T-scaled min-time — first
  hint the CR3BP `T·t_f,min` analog deviates ≥1% from exact 1/T at ELFO.
- Why it worked: the freetf `cScale` slack decouples the clock, so the
  pass-through chain reuses the source nodes without a winding mismatch.

## Tulip 20 mN — HONEST FAILURE (topology wall), machinery validated

**Run 1** requested the source t_f UNCHANGED at 20 mN → 0.92× the 20 mN
min-time → genuinely infeasible (defect floor 5e-3 at every ε, then the known
libcoinmumps MEX-fatal). This was a PILOT-DESIGN bug (fixed `807c642`: t_f now
C-law rescales by T_src/T_new, holding the 1.15× factor), NOT a solver failure
— the R0-law arithmetic predicts the infeasibility exactly.

**Run 2** (t_f = 9.043 ND, factor 1.15 held): `certified=0` — no schedule step
converged tight, the ε=1 ENERGY step included (defect ~1–8e-3 throughout;
switch counts thrashing 56–90 vs the nominal 25; edge 55–78% vs 99.6%).
`boundSatWorst=vBox`.

**Diagnosis — this is the fixed-τf topology wall the spec predicted, not a bug:**
- The chained seed carries the 25 mN winding (~40 revs); the 20 mN / 1.15×
  optimum wants MORE revs (t_f,min grows ~1/T, so at held factor the absolute
  transfer is longer and lower-thrust ⇒ more revolutions). Sundman's fixed τ_f
  freezes the seed's rev topology, so the solver cannot grow the winding — the
  ε=1 energy solve, which normally certifies trivially, cannot even close.
- The `vBox` boundSat lead is a SYMPTOM measured on a non-converged iterate,
  not a diagnosis — do not over-read it. It IS a concrete, testable next step
  (below).
- Contrast with ELFO confirms the mechanism: ELFO's `cScale` gives clock/span
  slack, so no winding mismatch; the tulip engine has no such freedom.

## What the gate established

- Package MACHINERY is validated end-to-end (ELFO PASS exercises every piece;
  the tulip boundSat diagnostic correctly surfaced the binding bound; the
  hardened gates correctly REFUSED to certify a non-solution rather than
  emitting a false positive — the single most important behavior).
- The tulip near-nominal band is **narrower than ELFO's** — a single 20% rung
  down already needs a rev-count change the fixed-τf engine can't supply.
- The spec's honest scope limit (near-nominal band; deep/rev-changing rungs
  need a free-span reformulation) is now EMPIRICALLY confirmed on the tulip
  side at a smaller thrust step than expected.

## Recorded next steps (ladder-campaign scope, not prep)

1. **Cheap disambiguation probe (do first):** re-solve ONLY the tulip 20 mN
   ε=1 energy step with a widened `vBox` (opt-in arg already built, e.g. 25)
   + more iterations. If it certifies → the tulip band exists, just needs a
   wider box (a cheap prep win); if it still stalls → confirmed topology wall,
   escalate. This directly tests the `vBox` lead.
2. **If topology-wall confirmed:** the tulip ladder needs the free-span
   reformulation (θ-domain / Δθ-free spiral phase — the CR3BP analog of the
   earth MEE+ΔL rebuild; see [[direct-vs-indirect-strengths]] and the spec's
   escalation note). ELFO can ladder near-nominal now; tulip cannot without it.
3. Per-rung anchors for either campaign: direct all-burn (Route-B) or R0-law
   estimate, refined — do not rely on pumpkyn indirect min-time (stress-probe
   pending; see the tulip TODO first ladder task).
