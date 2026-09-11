# DRO_tulip — TODO

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
