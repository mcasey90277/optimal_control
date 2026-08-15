# DRO_tulip — the first costate-catalog campaign (DRO → tulip, min-time)

The campaign that built the costate-library pipeline end to end, and still
its reference implementation. Product: libraries of converged min-time PMP
costates `z8 = [λ(7); tf]` over (DRO period × tulip petals × phasing torus ×
thrust 1–15 N), every entry accepted UNCHANGED by `pumpkyn.cr3bp.tfMin`.

**Docs (read these first):**

| file | what |
|---|---|
| `doc/costate_library_methodology.tex/.pdf` | Theory manual: CR3BP, min-time PMP (continuous burn is a theorem), covector mapping, ms_tfmin math, verification taxonomy incl. the conjugate-point test. Externally reviewed. |
| `doc/costate_library_sdd.tex/.pdf` | Software design document: the layers, data contracts, extension guide (new family / new parameter axis), accepted debts. |
| `doc/COSTATE_LIBRARIES.md` | Working record of the library products. |
| `doc/note_blocky_behavior.md` | Darin's "blocky tf map" artifact = solution-family walls, measured. |
| `process/COSTATE_LIBRARY_PIPELINE.md` | The three-step process record. |
| `FINDINGS.md` | Campaign findings log. |

## Layout

- **Front doors (root):** `run_costate_library.m` (12×12 thrust ladder),
  `run_catalog_sweep.m` (multi-orbit coarse sweep), `run_ladder_batched.sh`
  (hang-proof unattended batches), `run_lowthrust_ladder.m`.
- **`direct/`** — CasADi/IPOPT Hermite-Simpson min-time transcription
  (`lib/casadi_mintime_dro.m`), the phasing sweeps, `certify/`
  (`certify_dro_mintime` = G1/G1b/G2 gate stack + `dro_residual`,
  `costate_compare`), `viz/`, `results/` (incl.
  `dsweep_12x12_cells.mat` — the full-data-contract flagship torus that
  feeds `costate_common/golden_cells`).
- **`run_minenergy_pilot.m`** (root) — the first NON-min-time run of the
  pipeline (2026-08-14): fixed-t_f MIN-ENERGY (J = ∫s² dt) on flagship
  12×12 cells at t_f = γ·t_f^min — direct energy solve
  (`casadi_mintime_dro` `objective='energy'`, `tfFix`) → harvest →
  `indirect/ms_minenergy` (fixed-tf `ms_bvp`) → gates incl. the generic
  single-shooting acceptance. Records in `direct/results/minenergy_pilot.mat`.
- **`indirect/`** — `ms_tfmin.m` (multiple shooting; thin binding of
  `costate_common/ms_bvp` since migration #3), **`ms_minenergy.m`** (its
  fixed-t_f min-energy sibling, 2026-08-14; `tests/test_ms_minenergy.m`
  = synthetic known-answer BVP), `thrust_ladder_library.m`
  (THE ladder engine, family-agnostic endpoints — halo/DPO campaigns call
  it unmodified), `extend_thrust_ladder`/`densify_ladder`, packagers
  (`build_costate_catalog.m`, `build_costate_lib*.m`), pickers + examples
  (`costate_catalog_pick`, `costate_lib_describe`,
  `costate_catalog_extremes(_movies)` — demos run on whatever
  `costate_catalog_*.mat` sits in the current folder).
- **`deliverables/`** — shipped zips: `costate_lib_and_helper_deliverable2`
  (thrust-axis library + pickers), `costate_catalog_deliverable3`
  (multi-orbit catalog, 3,936 entries, 16 sheets).
- **`reviews/`** — external code/doc review artifacts.

## Key numbers

12×12 flagship torus: solved + refined + accepted (max |Δz| ~1e-9);
determinism demonstrated (bitwise-identical rerun). Coarse catalog:
16 sheets (τ ∈ {0.5, 1, 2, 3} × Np ∈ {5, 7, 9, 12}), 3,936 entries.
Seed-quality measurement that forced multiple shooting: collocation-dual
seeds miss 36,000–560,000 km when single-shot; ms residual ~1e-13.
Cheapest entry 0.984 km/s. Cross-method tf agreement median 2e-8.
