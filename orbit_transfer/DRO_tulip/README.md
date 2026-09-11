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

## The 70 mN phase sheet (2026-09-09)

A second propulsion regime beside the shipped catalog: **70 mN, Isp 900 s,
150 kg** — the cislunar-poster engine. It needs its own catalog because a
catalog carries one thruster (the shipped DRO → tulip catalog is Isp 1710 s),
but no schema work: the min-time schema already keys sheets by
`sD_frac x sA_frac`.

**One transfer, one call:**

```matlab
T = run_dro_tulip;                    % the anchor: 17.7976 d, 0.7485 km/s
T = run_dro_tulip(0, 0.1587);         % another arrival phase: 16.2256 d
T = run_dro_tulip(0.75, 0.0754, struct('movie', true));
```

It serves the pair from `dro_tulip_library` when it is already solved, and
otherwise walks there by continuation. Either route ends in the same gate
stack, so the printed time, ΔV and fuel always come with the flown miss, the
pumpkyn `tfMin` witness, the conjugate verdict and the three hypothesis gates.
A walk that runs out of budget says so rather than returning a number.

**The sheet:**

| step | unit |
|---|---|
| arrival axis (where the folds are) | `arclength_arrival` on `costate_common/arclength_ms` |
| departure axis (where steps are cheap) | `rib_from_crossing`, driven by `build_ribs` |
| certify one candidate | `certify_root` (`certify_crossing` converts the homogeneous chart first) |
| assemble | `sheet_from_arcs`, driven by `build_arrival_sheet` |
| package | `sheet_to_catalog_file` → `costate_common/build_costate_catalog_family` |
| picture | `plot_arrival_arcs` |

**What it found.** The conjugate test is the discriminator, not a formality:
the first assembled sheet held 14 candidates and certified 2, with the other
12 refuted by the conjugate test alone after passing the residual, the flown
arrival and the `tfMin` witness. Every refutation is slower than the certified
solution at its phase, and at a fold nose the test separates two roots 26
minutes apart. FINDINGS §37 (machinery), §38 (the result).

## Studying ONE transfer, and checking it

`transfer_study.m` is the front door for understanding rather than throughput:
a script with the scaffolding exposed. Eight sections -- generate the DRO from
its period, generate the tulip from its petal count, spell out the
nondimensionalisation, solve, verify independently, then check the NECESSARY
conditions and the SUFFICIENCY hypotheses **one at a time, computed in the
script**, and finally open a rotatable 3D figure.

**Standing rule (Mike, 2026-09-10): every costate library ships these two
scripts** -- `build_70mN_library.m` (the chain) and `transfer_study.m` (one
transfer, exposed). A new campaign copies both and changes the parameter
blocks; the front-door functions underneath stay the engine.

| unit | what it does |
|---|---|
| `transfer_study.m` | the script; edit the parameter blocks and run |
| `verify_with_pumpkyn.m` | hands our costates to pumpkyn's own solver and shows, component by component, that it does not move them -- with a CONTROL EXPERIMENT in its test proving the check can fail |
| `report_optimality.m` | the report in three groups (NECESSARY / SUFFICIENCY / CROSS-CHECKS) and four line states (PASS / FAIL / NOT CHECKED / UNRESOLVED); an unchecked line BLOCKS its group, because "every check that ran passed" is vacuously true when none did |
| `plot_transfer_3d.m` | the rotatable figure; drawn from the SAME flight the script measured (`opts.flight`), every annotation recomputed from it |
| `../../costate_common/validate_flight.m` | the ONE admissibility check every flight passes: reached t_f, finite, all-burn mass law, clear of both primaries |
| `../../costate_common/pmp_pointwise_checks.m` | Pontryagin on the flight: H = 0, transversality, the adjoint equations, and the EXACT minimum-principle gap of the control the propagator applied (its test injects a wrong-sign field and watches the gap open) |
| `build_70mN_library.m` | the whole library chain as a script in the same style: anchors -> arcs -> sheet -> ribs -> package -> audit -> sweep -> pictures -> deliverable, each stage a switch, each stage's file reused when off |
| `certify_root.m` / `certify_crossing.m` | the gate stack itself, fenced by hard timeouts |
| `audit_phase_catalog.m` | audits a SHIPPED catalog the way a recipient would: re-derives every entry from the catalog's own keys and flies it |
| `package_phase_catalog.m` | sheet + ribs -> a shareable catalog |
| `build_ribs.m` | departure ribs off every certified point of a sheet |

**A deliverable does not ship until its audit is clean** --
`build_dro_deliverable` enforces that rather than trusting a checklist. See
`../doc/CERTIFICATION_DISCIPLINE.md` for why, and FINDINGS 40-41 for the
reviews that produced the rules -- including the second script review, which
found the sampled minimum-principle check to be a tautology and H6 to be a
gate that was computed but never enforced.
