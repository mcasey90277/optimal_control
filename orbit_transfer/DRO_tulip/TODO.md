# DRO_tulip — TODO

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
- [ ] **Min-energy follow-ups (pilot 2026-08-14, `run_minenergy_pilot`):**
  γ grid per cell (J and m_f are non-monotone in γ on (2,5) — basins), a
  fixed-tf conjugate-point test (`ms_conjugate_test` is free-time only),
  a γ / t_f axis in `catalog_schema`, then the energy→fuel ε-homotopy on the
  same seeds toward min-fuel entries. Route the min-time tfMin acceptance
  through the `ss_bvp_accept` interface for one gate shape across costs.
