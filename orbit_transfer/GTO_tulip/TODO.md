# GTO_tulip — TODO

Two standing goals (2026-07-21): **(a) keep perfecting the direct code,
(b) get the indirect code working.** Full history behind every item:
`process/LOW_THRUST_MINFUEL_CAMPAIGN.md`.

## (a) Direct — perfect what works

- [ ] **Close the 1.01–1.11× near-min-time band.** The certified front starts
  at 1.12×; the energy backbone itself is only generatable for
  t_f ≈ 1.12×–1.95× (the band is hard even for the SMOOTH energy problem —
  a near-min-time conditioning wall, not just bang-bang structure). Candidate
  attacks: direct continuation from the min-time anchor downward, or the
  indirect band campaign below.
- [ ] **PSR/lib de-dup.** `direct/PSR/lib/` vendors ~20 files (params, refine
  suite, sms_* set) that partially duplicate `direct/sundman_minfuel/` and
  `../cr3bp_common/`. Deliberately left intact during the 2026-07-21
  restructure (behavior risk); fold into the shared sources with a
  reproduce-the-certified-result gate when touched.
- [ ] **Front hygiene.** Keep `aggregate_front`'s honest 3-class front current
  as new t_f points land; switch counts reported as bands (mesh-sensitivity
  lesson from the earth_elliptic P0 study applies here too).
- [ ] **Thrust ladder (Table-3 analog for the tulip).** Port the
  `../earth_elliptic_to_geo/` ladder recipe — per-rung min-time anchor,
  thrust-continuation warm-chaining, certified-only caching, R0-law check —
  to the tulip problem: sweep T_max around the nominal 25 mN, map
  t_f,min(T) and m_f(T), and test whether a T·t_f,min ≈ const law holds in
  the CR3BP. Prior art: the indirect ztl campaign was a thrust-ladder
  attempt (its P0 findings — ~75 mN sweet spot, fixed-t_f ladder argument —
  constrain the design); the earth-GEO MEE ladder is the proven direct
  machinery. FIRST TASK — pumpkyn stress probe (ztl-style): run pumpkyn
  `tfMin` at 22.5 and 20 mN, cold AND warm-chained from the 25 mN solution,
  to find where the indirect single-shooting min-time substrate dies. Do not
  block on it: per-rung anchors can be DIRECT all-burn solves (ELFO Route-B
  precedent) or R0-law estimates (earth 0.5 N precedent) — min-time + high
  revs is the regime where direct shines (no switch structure) and single
  shooting is weakest (STM sensitivity grows with revs).

## 2026-07-21 review follow-ups (doc/reviews/2026-07-21_triage.md)

Fixed same day: acceptance gates now require `Solve_Succeeded` (C1) and
`certified` requires the requested homotopy endpoint (C2).

- [ ] **Mesh-band study for the flagship (C3).** Re-solve the certified 1.15×
  solution at ≥2 finer meshes; publish the switch count as a band (the
  earth-campaign P0 protocol). LIVE EVIDENCE (2026-07-21 gate-hardening
  regression): a fresh re-solve of the certified recipe landed 24 switches /
  ΔV 3.3660 vs the published 25 / 3.3696 — same mass to ~0.1%, switch integer
  basin-sensitive even at fixed mesh. Also relax `certify_minfuel_pmp`'s strict
  integer PMP-crossing match (node-grazing switches fail it spuriously).
- [ ] **Bound-saturation diagnostic + box widening (C4).** Port the earth
  solver's `boundSaturation` warning into `casadi_minfuel_sundman`; widen the
  [-12,12] velocity boxes before any ladder work.
- [ ] **Ladder-prep trio (C5, feeds the thrust-ladder goal above):** per-rung
  thrust parameterization + artifact fingerprints; adaptive `tauf0`/`cBox`/
  state boxes; phase-correct cross-rung warm starts (steering-law
  regeneration, not index-carried controls). Plus the `rF` fallback in
  `gen_tulip_energy_2p` step_solve (C6).

## (b) Indirect — get it working

Today: machinery built and validated, no certified indirect solve yet.

- [ ] **IFS to certification.** `indirect/ifs/` (direct-seeded indirect
  finishing solve) stalls at ‖R‖ ≈ 0.023 on terminal-cluster conditioning at
  the full 1.12× problem. Next move (per campaign record): retarget to a
  clean-band t_f whose switches are non-grazing, then walk back to 1.12×.
- [ ] **ms_band.** Multiple-shooting attack on the 1.01–1.11× band — blocked by
  the same near-min-time wall; revisit after IFS certifies anywhere.
- [ ] **Use the ztl P0 findings** (recorded in `indirect/ztl/`): the min-time
  substrate is dead as a ladder start; fixed-t_f ladder argument; the
  cold-landscape asymmetry ("arrives-warm easy"); ~75 mN sweet-spot signal.
  These constrain which indirect strategies are worth another attempt.
- [ ] **Success bar:** an indirect (PMP shooting) solution of the 1.15×
  min-fuel problem certified against the direct result (25 switches,
  ΔV 3.3696 km/s) — or an honest documented refutation of why single/multiple
  shooting cannot close it at this scale.

## Housekeeping

- [ ] `indirect/min_time/` serves both tulip and ELFO retargeting — if it grows,
  consider promoting it to a shared home (noted in the restructure spec).
