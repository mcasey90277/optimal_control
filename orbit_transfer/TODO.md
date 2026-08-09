# orbit_transfer — TODO (program level)

Per-campaign items live in each subfolder's `TODO.md`; this file holds only
what spans campaigns. Task numbering follows the session task list where one
exists.

## Costate-catalog program (Darin roadmap)

- [ ] **Densify the red cells** of all three tulip catalogs (halo 46 pairs,
  DPO 65 pairs, DRO s_A=0.075 row + 6 stragglers — tasks #14/#15). The
  engines revisit any no-OK cell automatically; run `densify_ladder`-style
  passes with fresh attempt budgets.
- [ ] **L1↔L2 halo-to-halo catalog** (Darin liked the movies). Blocked on a
  schema extension: the catalog format keys arrivals by petal count `Np`;
  a halo arrival needs an **arrival-period axis** — the "new parameter
  axis" step in `DRO_tulip/doc/costate_library_sdd.tex`'s extension guide.
  Only two admissible L1 members exist (τ = 1.8037, 2.7433), so the box is
  2 × 11 L2 members. Watch item: one 10 N rung showed a lone
  tfMin-acceptance NaN with a healthy 0.08 km flown miss.
- [ ] **Sub-0.5 N rungs** deferred by Darin "until customers need it."
- [ ] Min-fuel catalogs: explicitly NOT yet (Darin, Aug 2026).

## Optimal-control library (goal-oc-library)

- [x] Migrations #2–#5 (duals_to_costates, ms_bvp + conjugate test +
  golden cells, flown-control verifier, foc covector unification) — done
  2026-08-08, each with an equivalence gate.
- [ ] Promote `costate_common` to a proper MATLAB package (`+oc` or
  similar), pumpkyn house style; SDD is the architecture document.
- [ ] Version the sheet-.mat schema and compact catalog format BEFORE a
  second data consumer exists (accepted debt in the SDD; includes the
  `cat.derive` formula-registry replacement and environment pinning).
- [ ] Acceptance-gate harness: generalize the tfMin acceptance gate as a
  per-family independent-solver harness (the 1 of 11 verification checks
  that is not yet generic).
- [ ] Batched-driver template + monitor template into the library (the
  halo/DPO drivers differ only in names; monitors must watch process
  liveness, not just log errors).

## Certification (see OPTIMALITY_CERTIFICATION.md before planning any of this)

- [ ] STRICT second-order sufficiency is still nowhere; the live path is
  the STM switching-time Hessian (Part B). The conjugate-point test
  (necessary condition) now runs on catalog cells and is exposed to the
  Earth campaigns via `foc_check` opts.msInfo — wire an Earth-campaign
  multiple-shooting solution to actually exercise it.
- [ ] foc gate coverage is uneven across the four direct campaigns
  (see the Part A matrix); close the standing gaps.

## Legacy campaign items (tracked in their own TODOs)

- GTO_tulip: 1.01–1.11× near-min-time band still open; indirect certification.
- GTO_ELFO: indirect Route C not started.
- earth_elliptic_to_geo: 0.5 N min-time wall; our own MATLAB indirect.
- earth_elliptic_to_geo_CR3BP: Phase-2 indirect not started.
