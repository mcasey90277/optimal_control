# costate_common — the optimal-control library seed

Family-agnostic core of the costate-catalog pipeline, and the official seed
(2026-08-06) of a pumpkyn-style **optimal-control orbit-transfer library** —
the OC companion to pumpkyn's astrodynamics. Everything here is used by at
least two campaigns; the migration rule is *code reused by a second campaign
moves here on its next touch*.

Architecture document: `../DRO_tulip/doc/costate_library_sdd.tex`.
Theory manual: `../DRO_tulip/doc/costate_library_methodology.tex`.
Both externally reviewed (GPT-5.6-terra + Gemini 3.1, 2026-08-07/08).

## Contents by layer

**Orbits / endpoints**

| file | what |
|---|---|
| `get_family_orbit.m` | THE endpoint provider: family name + params → propagated periodic orbit (dro/tulip/halo/dpo/lyapunov via pumpkyn getters + `cont_np` + `prop`). |
| `survey_family_bounds.m` | Admissibility survey per Darin's criteria (periselene ≥ 500 km, ≤ 100 Mm from the Moon), NaN-safe periodicity guard. |
| `assert_periodic_orbit.m` | The closure guard as a standalone assertion (NaN-safe `~(err<tol)` form). |

**Covector mapping / seeds**

| file | what |
|---|---|
| `duals_to_costates.m` | **DELEGATE** since 2026-08-09: the covector rules were promoted to the cross-folder library — `../../oclib/+oc/duals_to_costates` — when booster_landing's G5 gate became the second top-level consumer. This delegate keeps all costate_common callers working. |
| `harvest_ms_seed.m` | Direct solution → multiple-shooting seed (delegates the covector rules; owns only interpolation onto segment boundaries). |

**Multiple shooting / second order**

| file | what |
|---|---|
| `ms_bvp.m` | Generic multiple-shooting BVP engine: problem as three closures (prop/rhs/terminal), block Jacobian from segment STMs, trust-region-dogleg, iterate guards. `ms_tfmin` (DRO_tulip/indirect) is its CR3BP min-time binding — bitwise-equivalent to the pre-migration solver. Self-demo: linear oscillator BVP. |
| `ms_conjugate_test.m` | Free-final-time Jacobi (conjugate-point) test: det([Φ_xλ·P, f(t)]) on the costate-scaling quotient. The naive 6×6 det is IDENTICALLY singular for min-time PMP (scaling + λ_m invariances) — read the header before changing anything. Necessary condition, junction resolution. |

**Flown verification / screens**

| file | what |
|---|---|
| `flown_control_error.m` | The G1b gate: fly the reconstructed control end-to-end, report where you actually arrive. Since oclib move 2 a thin wrapper: CR3BP dynamics + scheme-matched control closure over the shared engine `oc.fly_control` (bitwise-equivalent). |
| `true_min_altitude.m` | Propagated (between-nodes) minimum lunar altitude — checks a collocation floor where it actually binds. |
| `cr3bp_thrust_rhs.m`, `ctrl_quad.m` | Shared physics + Hermite-Simpson control reconstruction for the two above. |
| `preflight_screen.m` | Node-level sanity (altitude, tf plausibility) BEFORE any integrator touches a solve. |

**Packaging / regression**

| file | what |
|---|---|
| `build_costate_catalog_family.m` | Family-agnostic compact-catalog packager (ND-only quantities, `derive` formulas, `dep_family`/`dep_params` recipes; keeps the legacy `tauDRO` field = departure period so every picker works on every catalog). |
| `golden_cells.m` + `golden_cells_data.mat` | 20-check quality regression: three engine cells (dro/halo/dpo, flown-perturbed 1 N entries + conjugate verdicts) + one harvest cell with REAL collocation duals. **Run after any change to the files above.** A quality drop (iterations, residual) is a failure even when correctness gates pass. |

## Conventions

- Pumpkyn house style: `%% Purpose / Inputs / Outputs / Revision History`
  headers, `nargin==0` self-demos, no Code Analyzer pragmas.
- Physics only through pumpkyn calls (`tfMinProp`/`tfMinEoM`/getters);
  nothing is ever written into pumpkyn.
- The five standing principles + principle 7 (defenses against silent
  quality degraders) govern all changes — see the SDD.
