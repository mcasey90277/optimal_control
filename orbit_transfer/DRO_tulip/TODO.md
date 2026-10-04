# DRO_tulip — TODO

## Live (2026-09-22): the interpolation line -- see `process/INTERPOLATION_STATUS.md`

- [ ] **Junction-state interpolation** -- NOW FIRST (FINDINGS 97): along a
  departure rib the guesses are 0.04% off yet re-flying a guessed z8 lands
  ~2,000 km away (26x the spine's amplification); interpolate the
  neighbours' ms junction states instead, rescore the rib (44% covered) and
  the spine (87%).
- [ ] **Investigate the conjugate point at s_D ~0.17** (Mike 2026-10-04: conjugate
  points along a candidate optimum warrant investigation). The rib at s_A
  0.2837 stopped on a conjugate ZERO at s_D 0.1717 (FINDINGS 96); the 24 x 48
  ribs on columns 14-16 stalled at 0.164-0.172 the same way. Is it one curve
  on the torus? What root holds the other side, and do the branches fold?
- [x] **Measured the departure axis** -- DONE 2026-10-03 (FINDINGS 96-97): one rib at 96 departure phases on one
  arrival column, wrapped as a one-column catalog, scored with
  `score_interpolator` (the arrival axis is done: 1/96 gives 87% usable).
- [x] **Fix the hole filler's pool** -- DONE 2026-10-02 (FINDINGS 93):
  refilled 24 x 48 = 1,143/1,152, audit clean, 0 pool failures.
- [ ] `reproduce_library_70mN` verdict on a FINER grid than the record: say
  "compared on shared cells", not REPRODUCED = FAIL.
- [x] ADOPTED 2026-10-03 (FINDINGS 95): `results/library_70mN_24x48_merged/` is the library of record
  (FINDINGS 94: 24 x 48 + the record's faster/missing entries = 1,152/1,152,
  audit 1,152/0, no worse than the record on any shared cell). Section-83
  pattern: repoint README/status, archive the 24 x 24 record whole.
- [ ] Fold cells: arclength polish from the blended guess (8 of 60 spine
  misses are folds).
- [ ] Junction-state interpolation instead of re-flying a guessed z8.
- [x] `results/` backed up 2026-10-02: `~/Backups/DRO_tulip_results_2026-10-02.tar.gz`
  (771 MB, same disk -- copy off-machine).


- [x] **`conj_spectrum` locates, classifies and refines its candidates;
  review-2 items on the sweep-held files applied** (FINDINGS 42, 2026-09-10):
  sign changes are candidates, interior ones refined 4x (zero vs near-miss),
  lift pair tightened to [1e-12 1e-9] after seven false "uncertified"
  verdicts, H6 clearance vs Hamiltonian residual, hMax rename, lift input
  checks. Corrected re-sweep DONE 2026-09-11: 115/115, 44 interior
  candidates all near-miss (the ridge, FINDINGS 42), 0 zero, every lift
  certified (worst 28x), torus redrawn.
- [x] **`build_70mN_library.m` full-chain rerun** — DONE 2026-09-11: live sheet/package/audit into `results_rerun/`, bitwise-identical to the shipped catalog on all 115 entries, audit 115/0 (FINDINGS 43). Still unexercised: re-walking the departure ribs (`run.ribs`, hours). Original note: dry-run verified on the
  reuse path only; run with `run.sheet`/`run.package`/`run.audit` on to prove
  the live path reproduces 115 entries before the next ship.
  Prerequisites fixed 2026-09-11: packaging refuses to strip the catalog's
  second-order writeback unless the sweep stage is on, and backs up what it
  overwrites (`guard_catalog_overwrite`); the chain's sidecar pointer is the
  v2 file (the first sweep's disagrees with the catalog) and every sidecar
  record is identity-checked on resume; `chainOverrides.outDir` rebuilds
  beside the shipped files so the rerun can be compared, not trusted.

- [ ] **Red-row campaign (task #14):** the 12 s_A = 0.075 cells of the
  12×12 torus are still red; engineered seeds / neighbor continuation.
- [ ] **Six MS-refinement stragglers (task #15)** from the wave-2 refined
  library; re-solve with higher K or fresh seeds.
- [ ] **Densify the coarse catalog** where sheets have unsolved pairs
  (`extend_thrust_ladder` / `densify_ladder` already support it).
- [ ] **Period-axis growth per Darin:** "all reasonable periods/petals" —
  the measured admissible box (DRO τ 0.05–3.25, all tulips Np 3–14) is
  wider than the 4×4 catalog; grow sheets toward the box edges.
- [ ] Rewire `sweep_phasing_shoot.m` / `ms_refine_catalog.m` /
  `build_costate_lib*.m` legacy paths onto the shared library on next
  touch (migration rule); they predate `costate_common`.
- [ ] Legacy packagers (`build_costate_lib.m`, `_v2`) and the pre-catalog
  library formats are superseded by the compact catalog — archive or mark
  deprecated when Darin confirms he's fully on deliverable-3+ format.
- [ ] `doc/costate_libraries.tex` and `dro_tulip_mintime.tex` predate the
  halo/DPO campaigns — refresh or fold into the methodology doc.
- [x] **Min-energy follow-ups (pilot 2026-08-14, `run_minenergy_pilot`):**
  γ grid (2026-09-01), fixed-tf conjugate test (2026-09-02, corrected
  09-05), γ axis = catalog schema v3 (09-02), energy→fuel homotopy and the
  min-fuel catalog v3.1 (09-02..07). FINDINGS 17-29.
- [ ] Route the min-time tfMin acceptance through the `ss_bvp_accept`
  interface for one gate shape across costs (the one open piece of the item
  above).
