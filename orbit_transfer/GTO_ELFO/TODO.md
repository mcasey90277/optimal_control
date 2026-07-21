# GTO_ELFO — TODO

Goal: **both methods working.** Direct is in good shape; the open work is
mostly on the indirect side.

## Direct — polish

- [ ] **Near-min-time end of the ΔV–t_f front.** The front is mapped and
  labeled against t_f,min = 6.0962 ND (27.02 d), but the transition region
  just above min-time deserves the same scrutiny the tulip front got (the
  tulip campaign's 1.01–1.11× band was hard even for the smooth energy
  problem — check whether ELFO shows the same wall).
- [ ] **Switch-count bands.** Report ELFO switch counts as mesh-convergence
  bands, not integers (lesson from earth_elliptic_to_geo's P0 study).
- [ ] Keep `elfo_export_data` / movies current as new front points land.
- [ ] **Thrust ladder.** Same goal as the tulip and earth-GEO campaigns: sweep
  T_max around the nominal 25 mN with per-rung min-time anchors +
  fixed-c_tf fuel solves (thrust-continuation warm-chaining, certified-only
  caching), and check the T·t_f,min ≈ const law analog for the ELFO target.
  Port the `../earth_elliptic_to_geo/` ladder recipe.

## 2026-07-21 review follow-ups (doc/reviews/2026-07-21_triage.md)

Fixed same day: factor semantics rebased on `tfMin_ELFO` with tf-nearest seed
selection (C1); resume path certification-gated + single-trajectory saves
(C2); acceptance gates require `Solve_Succeeded` (C3).

- [ ] **tf-sweep tight re-clean (C4).** `gen_elfo_energy_tfsweep` banks
  500-iter loose solves directly; add the tight re-clean before `save_point`,
  as gravhom/fuel already do.
- [ ] **Mesh-band repeat + filename hygiene (C5).** Repeat ≥1 front point at a
  refined mesh and report switch counts as bands; parameterize N with t_f;
  drop the switch integer from identity-bearing filenames at the next format
  change.
- [ ] **Ladder-prep trio (C6, feeds the thrust-ladder goal above):** per-rung
  thrust + fingerprints/namespacing for seeds/checkpoints/rows; adaptive
  `cBox`/mesh + bound-saturation diagnostics; gate the min-time anchor save
  on full certification.

## Indirect — get it working

- [ ] **Route C: ELFO min-time indirect.** The direct min-time anchor (Route B)
  is certified; the indirect counterpart was scoped but never built. Start
  from `../GTO_tulip/indirect/min_time/` (PMP always-burn shooting root,
  pumpkyn `tfMin` machinery) retargeted to `gto_elfo_endpoints`.
- [ ] **Indirect min-fuel.** After Route C: PMP bang-bang shooting seeded from
  the certified direct solution (the same direct-seeded strategy as
  `../GTO_tulip/indirect/ifs/` — reuse its lessons, including the
  terminal-cluster conditioning failure mode).
- [ ] Replace the `indirect/` README stub with real structure once the first
  solver lands.
