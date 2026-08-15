# orbit_transfer — TODO (program level)

Per-campaign items live in each subfolder's `TODO.md`; this file holds only
what spans campaigns. Task numbering follows the session task list where one
exists.

## Costate-catalog program (Darin roadmap)

- [ ] **Densify the red cells** of all three tulip catalogs (halo 46 pairs,
  DPO 65 pairs, DRO s_A=0.075 row + 6 stragglers — tasks #14/#15). The
  engines revisit any no-OK cell automatically; run `densify_ladder`-style
  passes with fresh attempt budgets.
- [x] **L1↔L2 halo-to-halo catalog** — DONE 2026-08-09, BOTH directions,
  shipped as deliverable 6 on the schema-v2 arrival-period axis (A: 1,952
  entries / 90%; B: 2,096 / 96%; directions measurably asymmetric). Still
  open from it: the lone 10 N tfMin-acceptance-NaN watch item, and
  densification of the remaining pairs.
- [ ] **Sub-0.5 N rungs** deferred by Darin "until customers need it."
- [ ] Min-fuel catalogs: explicitly NOT yet (Darin, Aug 2026).
- [~] **Min-energy / fixed-t_f pipeline** — pilot DONE 2026-08-14 on three
  flagship DRO→tulip cells + a γ ladder (`DRO_tulip/run_minenergy_pilot`):
  fixed-tf `ms_bvp`, `ms_minenergy`, min-energy PMP field, generic
  single-shooting acceptance gate. Next: a t_f-multiplier axis in the
  catalog schema, a fixed-tf conjugate test, then energy→fuel homotopy on
  the same seeds (the min-fuel catalog route).

## Optimal-control library (goal-oc-library)

- [x] Migrations #2–#5 (duals_to_costates, ms_bvp + conjugate test +
  golden cells, flown-control verifier, foc covector unification) — done
  2026-08-08, each with an equivalence gate.
- [~] Promote to a MATLAB package — STARTED 2026-08-09 as the top-level
  cross-folder `../oclib/+oc` (moves 1–2 done: `oc.duals_to_costates`,
  `oc.fly_control`, both with bit-identical acceptance on orbit AND
  booster flagship data). Next: move 3 = `oc.ms_bvp` +
  `oc.ms_conjugate_test` promotion + cart-pole PMP-BVP integration demo.
  See `../oclib/README.md` and `../OCP_UNIFYING_MATH.md` §5.
- [x] Schema versioning — DONE 2026-08-09 (schema v2: `catalog_schema`
  validator + named derive registry + environment pinning; v1
  compatibility proven bitwise).
- [~] Acceptance-gate harness: `costate_common/ss_bvp_accept` (2026-08-14)
  is the generic single-shooting form, used by min-energy; still open:
  route the min-time tfMin gate through the same interface.
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
