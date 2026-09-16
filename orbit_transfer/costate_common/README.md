# costate_common — the optimal-control library seed

Core of the costate-catalog pipeline, and the official seed (2026-08-06) of
a pumpkyn-style **optimal-control orbit-transfer library** — the OC companion
to pumpkyn's astrodynamics. Its admission rule (changed 2026-09-16, was
*reused by a second campaign*) is *generic by construction, with a test*:
nothing in the file is specific to one campaign's orbits or bookkeeping, and
a test pins its contract. Job control is leaving for a sibling
`../campaign_common/` (TODO step 11). Problem-agnostic pieces with a second
TOP-LEVEL consumer go one level further, to `../../oclib/+oc` (admission rule
in `../../oclib/README.md`).

Architecture document: `../DRO_tulip/doc/costate_library_sdd.tex`.
Theory manual: `../DRO_tulip/doc/costate_library_methodology.tex`.
Generated function index: `../doc/library_catalog.md` (regenerated
2026-09-16; rerun `python3 orbit_transfer/doc/gen_library_catalog.py` after
adding files).

## State of the folder (measured 2026-09-16)

**58 files and 41 tests** (61 and 58 before cleanup steps 5–7, 2026-09-16). The folder did not match its
old admission rule; the table below is what the new one sorts. Caller counts below are code references from outside the
folder (tests excluded); an "entry point" has no calling code and is run by
hand.

| who uses it | files |
|---|---|
| two or more campaigns | 9 — `get_family_orbit`, `build_costate_catalog_family`, `nd_propulsion`, `ms_tfmin`, `phase_state`, `seed_from_z8`, `catalog_schema`, `ms_conjugate_test`, `duals_to_costates` |
| DRO_tulip only | 35 (`rib_targets` moved out) |
| DPO_tulip only | 1 — `survey_family_bounds` |
| no caller outside the folder | 13 — 9 internal helpers, 4 entry points (`conj_catalog_pass`, `gates_catalog_pass`, `golden_cells`, `campaign_status`) |

Known structural issues, each a step in `TODO.md` → *Cleanup plan*:

- **Inverted dependencies:** `second_order_pass` and `golden_cells` put
  `DRO_tulip/indirect` on the path, and `golden_cells` loads a gitignored
  DRO_tulip results file.
- **The most-shared engine lives in a campaign:** `DRO_tulip/indirect/thrust_ladder_library`
  (with `casadi_mintime_dro`, `certify_dro_mintime`, `dro_residual`) is called by
  HALO, DPO, HALO_HALO and GTO, each of which adds `DRO_tulip` to its path.
- **Job control lives here:** the nine orchestration files at the bottom are
  campaign infrastructure, not optimal control.
- ~~Tests of campaign code live here~~ — fixed 2026-09-16: the 15
  DRO_tulip-only tests moved to `DRO_tulip/indirect/tests`. (`test_gto_family`
  was misclassified: it tests this folder's `get_family_orbit`, so it stays.)

## Contents by layer

The **used by** column is measured, not aspirational: `DRO` etc. are campaign
folders, *internal* means called only from within this folder.


**Orbits and endpoints**

| file | what | used by |
|---|---|---|
| `get_family_orbit.m` | THE endpoint provider: family name + params → propagated periodic orbit (dro/tulip/halo/dpo/lyapunov via pumpkyn getters + `cont_np` + `prop`). | DRO, HALO, DPO, HALO_HALO, GTO |
| `survey_family_bounds.m` | Admissibility survey per Darin's criteria (periselene ≥ 500 km, ≤ 100 Mm from the Moon), NaN-safe periodicity guard. | DPO |
| `assert_periodic_orbit.m` | The closure guard as a standalone assertion (NaN-safe `~(err<tol)` form). | internal |
| `periodic_pp.m` | The C1-**periodic** cubic through one period of a closed orbit, its derivative (from the interpolant's own coefficients) and its seam: `seam.value` is the DATA's closure, which no interpolant can improve, `seam.deriv` the interpolant's derivative mismatch at the seam. An ordinary not-a-knot spline is not C1 at s = 0 — measured on the 7-petal tulip its derivative jumps by ~1e-2 there, while in VALUE the two agree below a millimetre except within a fraction of a percent of the seam, where the gap reaches metres. Extracted 2026-09-11 from two private copies that disagreed on policy (`transfer_study`'s `periodicPP` errored without the Curve Fitting Toolbox; `arclength_arrival`'s `makePP` fell back in silence) — the policy is now `opts.onMissing`, and the fallback warns. `tests/test_periodic_pp` (convergence order, both seams, the not-a-knot contrast, three refusals). | internal (`phase_state`) |
| `phase_state.m` | THE endpoint rule: one periodic orbit in, two closures out — the state at a phase FRACTION and `dx/ds` (per unit phase, so it carries the period). Built on `periodic_pp`. Consumers: every campaign engine that evaluates an endpoint at a phase — the live instruments (`second_order_pass`, `conj_catalog_pass`, `gates_catalog_pass`, `audit_phase_catalog`), the ladder engines and phase sweeps (`thrust_ladder_library` and 10 more), the GTO campaign tools, plus `arclength_arrival` (its `dstateA` is what `dR/dsA` differentiates) and `transfer_study`. Migrated 2026-09-11, 18 files; the recipient-facing deliverable helpers keep their inline `interp1` DELIBERATELY (a shipped catalog's helpers carry no dependency on this library) — see TODO. `tests/test_phase_state` — including the EQUIVALENCE GATE: four golden vectors from the pre-move implementations reproduced **bitwise**, one of them next to the seam. | DRO, GTO |

**Units**

| file | what | used by |
|---|---|---|
| `nd_propulsion.m` | THE nondimensional propulsion conversion: a thruster in N / s / kg becomes `c_nd` and the ND thrust acceleration at unit mass fraction, plus `ndT` (thrust → T_nd, for rung ladders) and `ndC` (Isp → c_nd, for a continuation over Isp). About twenty files wrote the same three expressions by hand — one of them already commented that it used "the same ND thrust/exhaust conversion thrust_ladder_library uses". The expressions here are those character for character, so every migrated site keeps its value BITWISE. Thrust and Isp may each be omitted (`[]`) when the caller has no single value or already holds `c`. `tests/test_nd_propulsion` (16 checks, incl. bitwise reproduction of the inline formulas and three refusals). | DRO, GTO, HALO_HALO |

**Dynamics and PMP fields**

| file | what | used by |
|---|---|---|
| `cr3bp_thrust_rhs.m` | CR3BP dynamics with thrust and mass flow: the shared right-hand side behind `flown_control_error` and `true_min_altitude`. | internal |
| `mintime_rhs_point.m` | The 14-state minimum-time PMP field at a point (pumpkyn `tfMinEoM`, state and costate rows). Shared by `ms_tfmin` and `ms_tfmin_hom`; extracted 2026-09-08 on its second consumer. | internal |
| `mintime_prop_seg.m` | One multiple-shooting segment of the minimum-time PMP flow via pumpkyn `tfMinProp`, with the 14 × 14 STM on request. The propagator closure of `ms_tfmin` and `ms_tfmin_hom`, and the flow `conj_spectrum` samples; extracted 2026-09-08. | internal |
| `cr3bp_minenergy_pmp.m` + `cr3bp_minenergy_prop.m` | **Min-energy PMP field** (J = ∫s² dt, Bertrand–Epenoy ε = 1 — GTO_tulip's convention; s* = clip((T/2)(‖λ_v‖/m + λ_m/c), 0, 1)) with exact AD Jacobian (CasADi SX, built once) and an STM propagator in `tfMinProp`'s shape/tolerances. Unit-tested to equal `pumpkyn.cr3bp.tfMinEoM` at saturation to 1e-13. The one pipeline dynamics that cannot be a pumpkyn call (pumpkyn has no min-energy EoM). | DRO |
| `cr3bp_minfuel_pmp.m` + `cr3bp_minfuel_prop.m` | **Smoothed energy→fuel PMP field**, THREE families selected by `smooth = {family, p[, delta]}`: `'eps'` Bertrand–Epenoy L = (1−p)s + ps² (ε = 1 is the min-energy field exactly; always arrives, slow at the floor); `'huber'` (PLQ core s = pQ below Q = 1, JUMP to 1 above — propagated **event-split with the saltation matrix** Φ⁺ = [I + (F⁺−F⁻)nᵀ/(nᵀF⁻)]Φ⁻ since 2026-09-05; the 09-02 "refuted" verdict was our missing saltation — with it Huber is 5× faster than ε and failure-free where it works, and walls at grazing-risk switch-structure changes on ~30% of cells, FINDINGS §22–24); `'huberc'` (2026-09-06: Huber core + CONTINUOUS ramp p→1 over Q ∈ [1, 1+δ], exact argmin of a convex H; passes Huber's walls, and at fixed δ ≈ 0.03 the p-walk is failure-free, FINDINGS §25–26; the only family to reach the floor at γ = 2, §28). Q̇ is throttle-independent and continuous across the switch (`cr3bp_minfuel_qdot`). Smoothing parameters are CasADi Function inputs. `ms_minfuel` (DRO_tulip/indirect) is the fixed-tf binding; `run_minfuel_race` walks (p, δ) rungs with loose/tight gates and a family per arm. | DRO |
| `cr3bp_minfuel_qdot.m` | Closed-form dQ/dt = −T λ_vᵀλ_r/(m\|λ_v\|) along the min-fuel flow (2026-09-06, Astra review #2): the throttle drops out exactly, so nᵀ(F⁺−F⁻) = 0 at every switch (measured 9e-16) — the exact transversality of a crossing, used by the Huber propagator's branch/grazing logic and by `DRO_tulip/indirect/huber_switch_diag`. | DRO |

**Shooting solvers and continuation**

| file | what | used by |
|---|---|---|
| `ms_bvp.m` | Generic multiple-shooting BVP engine: problem as three closures (prop/rhs/terminal), block Jacobian from segment STMs, trust-region-dogleg, iterate guards. `ms_tfmin` (this folder) is its CR3BP min-time binding. **`opts.fixedTf`** (2026-08-14) drops the t_f unknown for fixed-time problems (min-energy / min-fuel catalogs); `ms_minenergy` is the first binding. Newton polish after an early fsolve exit (its ‖JᵀR‖ test fires at ‖R‖~1e-10 on short arcs). Self-demos: oscillator BVP, free and fixed t_f. **Junction contract (verified 2026-09-16):** `info.Y` returns the K junction STARTS (14 × K) with `info.tGrid` the K+1 times, and a seed's `Y` is read only in columns 1..K, so a 14 × K array fed back is lossless. A seed may carry K or K+1 columns (headers say so since 2026-09-16). | DRO (+ internal) |
| `ms_tfmin.m` | Min-time wrapper around `ms_bvp` (pumpkyn tfMinProp/tfMinEoM closures, free-tf terminal set, opt-in conjugate test). MOVED here from `DRO_tulip/indirect` 2026-08-26 (used by every catalog campaign + the GTO probe); no delegate — callers self-bootstrap this folder (the `ms_bvp` precedent). Equivalence gate: `golden_cells` 20/20. | DRO, GTO |
| `ms_tfmin_hom.m` | HOMOGENEOUS-chart minimum-time multiple shooting: the objective multiplier ρ free and (ρ, λ₀) on the unit sphere, H(t_f) = ρ + λᵀf = 0. In the normal chart the fast family's \|λ₀\| runs 46 → 1449 as thrust falls toward 72 mN; on the sphere the multipliers stay bounded and ρ → 0 (loss of normality) is a finite, visible event. `tests/test_ms_tfmin_hom`. | DRO |
| `ss_bvp_accept.m` | **Generic single-shooting acceptance gate** (2026-08-14): the pipeline's third gate for costs with no pumpkyn twin — `ms_bvp` with K = 1 on the same closures, reporting the residual AT the seed, the move \|Δz\|, and `accepted = converged ∧ \|Δz\| < 1e-6`. First concrete form of the "acceptance-gate harness" TODO. | DRO |
| `arclength_ms.m` | **Generic pseudo-arclength continuation** (2026-09-09) of a root curve R(p,q)=0 in ANY parameter, for any residual that comes with its Jacobian: scaled coordinates (the caller passes per-BLOCK scales), tangent = right null vector of the FULL SVD of [R_x R_q] (a linear solve is exactly singular AT a fold, which is where the tangent matters), bordered Newton corrector with backtracking + a branch-jump cap, Newton-count step control, fold = tangent q-sign change classified by rank-one loss of R_x with [R_x R_q] regular, and EVERY crossing of a requested q level located by Newton at fixed q and KEPT (a grid point may hold several candidates). Two defects found by test and fixed the same day: a level crossed TWICE inside one step was missed entirely (the step is now split at the localized turning point — this is exactly what a fold does), and the fold's own position was only O(ds²) accurate (now an ORIENTED secant on tau_q; svd's null-vector sign is arbitrary, and an unoriented secant made it worse). `tests/test_arclength_ms` (unit circle: fold, double crossing; cubic q = x³−3x: two folds, branch-jump detection via x-monotonicity) and `tests/test_arclength_ms_thrust` — the generic engine fed the production thrust residual reproduces the archived 2026-09-08 `arclength_thrust` diagnostic arc ROOT FOR ROOT (t_f rel err 0.0, \|λ₀\| 73.8855 = 73.8855). | DRO |
| `newton_fixed_q.m` | Scaled Newton at a FIXED continuation parameter, landing exactly on a requested level. Extracted from `arclength_ms` on its second consumer (`crossings_from_arc`, 2026-09-12) so a crossing re-scanned from a saved arc is corrected by the same step the walk would have applied. A non-finite residual or step reports not converged. | DRO (+ `arclength_ms`) |

**Covector mapping and seeds**

| file | what | used by |
|---|---|---|
| `duals_to_costates.m` | **DELEGATE** since 2026-08-09: the covector rules were promoted to the cross-folder library — `../../oclib/+oc/duals_to_costates` — when booster_landing's G5 gate became the second top-level consumer. This delegate keeps all costate_common callers working. | DRO, verify_common (delegates to `oclib`) |
| `harvest_ms_seed.m` | Direct solution → multiple-shooting seed (delegates the covector rules; owns only interpolation onto segment boundaries). | DRO |
| `seed_from_z8.m` | Fly a trusted z8 with tfMinProp and cut into K+1 junction states — the seeded-at-a-root ms seed builder (1–2 Newton iters), extracted on its third appearance (`golden_cells`, `conj_catalog_pass`, GTO probe). Bit-identical to the former inline code (measured). | DRO, GTO |
| `seed_from_entry.m` | ONE library entry becomes ONE multiple-shooting seed, by whichever of the two routes the entry needs: a FILE-BACKED entry brings its own junction states (reused, with the final column appended and column 1 pinned to the ACTUAL departure state and the entry's costates), a CATALOG entry carries z8 only and is flown (`seed_from_z8`). The same construction existed three times — inline in `transfer_study`, as `run_dro_tulip`'s private `seed_of`, and inside `build_arrival_sheet`'s seeding loop. `tests/test_seed_from_entry` (6 checks, both routes bitwise against the inline versions). | DRO |
| `flight_to_junctions.m` | A flown trajectory cut into the K+1 junction states a multiple-shooting seed needs: interpolate onto a uniform grid in NORMALIZED time (duplicate propagator samples dropped) and, for a seed at a different thrust or t_f, replace the mass row by the DERIVED all-burn identity m = 1 − T t/c rather than rescaling the old one (a rescale once carried a scaling mistake through a rung change). Six engines held this inline and `seed_from_z8` a seventh copy; all seven now call it. `seed_from_z8` queried in ABSOLUTE time, which differs in the last bits — 3e-13 on a real flight, all in the costate rows, a seed perturbation gated by `golden_cells` (residual 5.6e-14 → 1.2e-13 against a 1e-10 bar, iteration counts unchanged). `tests/test_flight_to_junctions` (9 checks incl. bitwise equivalence with the engines' inline form). | DRO |

**First-order and flown verification**

| file | what | used by |
|---|---|---|
| `fly_transfer.m` | Fly converged min-time costates ONCE and hand back what every consumer recomputed: the flight, its admissibility (through `validate_flight`), the flown miss, and the mass, propellant and Delta-V that follow. It REPORTS admissibility rather than throwing, because a script asserts where a certifier returns a reason. Delta-V is the rocket equation `c log(1/mf)` — the same one `catalog_schema`'s `deltav_from_mf` applies to catalog entries, and `tests/test_fly_transfer` asserts the two agree on the real catalog (the schema's form needs a whole catalog struct, so it is not routed through it). | DRO |
| `validate_flight.m` | ONE admissibility check for a flown all-burn trajectory: reached t_f, finite, mass law m = 1 - (T/c)t, clear of both primaries. Used by the certifier, the witness flight, `verify_with_pumpkyn` and `transfer_study`, so "admissible" means one thing. | DRO |
| `pmp_pointwise_checks.m` | Pontryagin ON THE FLIGHT, complementing the shooting residual: H = 0, transversality, the adjoint equations (two coordinate-scaled FD steps, Richardson), and the EXACT minimum-principle gap of the control the propagator APPLIED, recovered as (powered - coasting) field. A sampled sphere check against the analytic minimiser was a tautology. Its test injects a wrong-sign thrust field and watches the gap open. Gated in `certify_root`. | DRO |
| `flown_control_error.m` | The G1b gate: fly the reconstructed control end-to-end, report where you actually arrive. Since oclib move 2 a thin wrapper: CR3BP dynamics + scheme-matched control closure over the shared engine `oc.fly_control` (bitwise-equivalent). | DRO |
| `true_min_altitude.m` | Propagated (between-nodes) minimum lunar altitude — checks a collocation floor where it actually binds. | DRO |
| `preflight_screen.m` | Node-level sanity (altitude, tf plausibility) BEFORE any integrator touches a solve. | DRO |
| `ctrl_quad.m` | Lagrange quadratic control reconstruction through the node / midpoint / node values of a Hermite-Simpson interval, for the flown-control check. | internal |
| `scalar_verdict.m` | A gate verdict is a real finite SCALAR or it is not a verdict. In MATLAB `[]` makes `v ~= 1` FALSE (an empty verdict acts like a pass), NaN trips no bound, `Inf > 0` reads as a satisfied gate, and `if` on a vector means ALL. `tests/test_scalar_verdict` 18 checks. | DRO |

**Second-order instruments**

| file | what | used by |
|---|---|---|
| `ms_conjugate_test.m` | Jacobi (conjugate-point) test, BOTH time conventions. Free-final-time (min-time): det([Φ_xλ·P, f(t)]) on the costate-scaling quotient — the naive 6×6 det is IDENTICALLY singular there (scaling + λ_m invariances); read the header before changing anything. **Fixed-final-time** (min-energy/min-fuel, 2026-09-02, **corrected 2026-09-05**): the SAME instrument with `freeTime=false`, `quotientDir=[]` (the running cost breaks the scaling invariance), **rows 1:7 = the full state** (a Jacobi field must vanish in r, v AND m; the earlier `[1:6 14]` block is the terminal shooting Jacobian, meaningful only at t_f), cols 8:14. Since 09-05 for BOTH conventions: samples through t_f (`ms_bvp` returns `.Yend`), equilibrated sign test, initial-coast/saturated-arc structural zeros skipped (`.firstFullRank`), last bracket counted, spec echoed. Validated on the LQ π-conjugate case + coast + final-segment fixtures (`tests/test_conj_fixedtf`, 10/10); golden 20/20. Necessary condition, junction resolution — see `doc/extremal_and_local_min_survey.md` for what sufficiency would need. | DRO, verify_common |
| `conj_spectrum.m` | **Dense singular-SPECTRUM scan** of the free-time quotiented conjugate matrix, closing both blind spots of the sampled sign test (two conjugate times inside one segment; an even-multiplicity crossing with no sign change). Built on `mintime_prop_seg`, sub-divided. It is a CANDIDATE-DETECTION scan: every sigma_6 dip and every determinant sign change is a LOCATED candidate (t/t_f, sigma_6/sigma_5) classified start / endpoint / interior, and every interior one is refined TWICE (4x, then 16x) -- a zero keeps falling at both levels or carries a sign change, a near-miss plateaus. One level was not enough: four 70 mN entries read 0.49 at the first level and 0.8-1.0 at the second (FINDINGS 42). `tests/test_conj_spectrum`: certified anchor (no interior candidate) and refuted control (interior zero, first-level ratio 0.06). | DRO |
| `conj_resolve.m` | Candidate detection and resolution for a dense conjugate-matrix scan, written as a PURE function of a matrix-valued function of time so synthetic matrices can drive it in tests; `conj_spectrum` feeds it the propagated CR3BP matrix. Column-normalised samples, sign trusted only above `signTol`, coarse local minima located and refined, a vanishing raw column counted as a candidate in its own right (Astra review #3, 2026-09-11). `tests/test_conj_resolve`. | internal (`conj_spectrum`) |
| `mintime_hypothesis_gates.m` + `gates_catalog_pass.m` | **Sufficiency-hypothesis gates for the min-time conjugate test** (2026-09-06, `doc/mintime_second_order_audit.tex`): on the dense `tfMinProp` flight of an entry, min\|λ_v\| (strong Legendre), min Q_mt = \|λ_v\|/m + λ_m/c (all-burn is the PMP control; `tfMinEoM`'s own switching function), and the abnormal-lift probe dim S = 1 (lifts of a fixed trajectory form a linear space; λ·f is a functional on it whose kernel is the abnormal lifts). Golden cells 15/15. `gates_catalog_pass` runs them catalog-wide with the `conj_catalog_pass` campaign contract (sidecar `*_gatesprog.mat`, resume, writeback). `gates_catalog_pass` is an entry point with no calling code; 2026-09-07 it passed 18,360/18,360 entries across the four catalogs. | DRO; entry point |
| `lift_space_dim.m` | The numerical-rank rule behind the abnormal-lift gate (2026-09-07): dim S = #{sv < tol}, tol = min(max(rankTol·sv₁, 10·nullResid), 1e-3·sv₁). A null space cannot be resolved finer than its known member's residual — the first catalog-wide pass with a fixed 1e-8 under-counted 512 entries as "dim S = 0" (gaps sv₇/sv₆ ~ 1e-7..2.5e-6, i.e. clean one-dimensional null spaces); the cap keeps a noisy lift from ever inflating the count. `tests/test_lift_space_dim` (5/5). | internal |
| `lift_margin.m` | The `dim S` rank statement as a **measured margin**: `dim S >= 1` is constructive (the lift is exhibited; a zero or non-finite vector is refused), `dim S <= 1` is Eckart-Young -- `sigma_6 / \|\|dC\|\|` with the error MEASURED by rebuilding C at a second integration setting. The pair matters: with a loose 1e-7 second build the "error" is that build's own error, and seven long-arc entries with sigma_6 = 0.99 read 4-9x; the sweep uses [1e-12 1e-9] (43x on the same entry). `tests/test_lift_margin`. | DRO |
| `h6_margin.m` | H6, `lambda_m(0) < c/T`, the reduced conjugate instrument's validity condition, with the clearance of h_max = -1 + (T/c) lambda_m(0) below zero judged against the arc's own Hamiltonian residual (`opts.Hresid`, wired from `mintime_hypothesis_gates`). ENFORCED in `certify_root` since 2026-09-10 -- it had been computed and ignored. | DRO |
| `second_order_pass.m` | Campaign sweep of the three second-order instruments over a catalog: sidecar after every entry (resume for free), writeback only on a complete census (`conj_interior`, `conj_interior_cand`, `conj_near_miss`, `conj_zero`, `h6_margin`, `lift_margin`). Latest run 2026-09-15: all 576 entries of the 70 mN library of record, 0 interior crossings, worst H6 4.44x, worst lift 11x. **Since 2026-09-11 every sidecar record carries the identity of the entry it measured (cell key + z8)**, and a resumed record must match the catalog (`second_order_pass:staleSidecar`): the sidecar had been positional, so a re-packaged catalog with a different entry set would have inherited other entries' measurements. A pre-key sidecar is adopted only with `adoptLegacy`, and only if every written-back value equals the catalog's. The library of record's sidecar is `second_order_progress_v3.mat`, beside its catalog; the first sweep's `second_order_progress.mat` disagrees with the catalog's lift margins by up to 1.35e4 and must not be used. `tests/test_second_order_sidecar_identity` (5 checks, mutation-tested). **Inverted dependency:** it puts `DRO_tulip/indirect` on the path and calls that campaign's `ladder_endpoints` (cleanup step 9). | DRO |
| `conj_catalog_pass.m` | Catalog-scale conjugate-point sweep: per entry, rebuild endpoints from the sheet recipes, fly the stored z8, re-solve with `ms_tfmin(conjTest)` seeded at the solution (verdict only when \|z−z8\| < 1e-6), store `conj_pass/conj_ncross/conj_atfinal` grids + `conj_test` provenance. Campaign contract: sidecar resume, attempt-before-solve, batch budget, explicit writeback with backup. First run 2026-08-23: 15,896 entries, 15,895/1/0. Entry point: no calling code, run by hand over every shipped catalog. | entry point |

**Execution fences**

| file | what | used by |
|---|---|---|
| `run_capped.m` | **Parfeval hard per-call timeout** (promoted 2026-09-01 on its second consumer): the only fence that bounds a CRAWLING integration — in-process wall checks fire between iterations, but one CR3BP segment propagation with a near-primary iterate can grind for hours inside a single function evaluation (measured twice). Worker killed + restarted on timeout. Returns `ok = false` for a timeout OR a worker error; the caller cannot tell which. Without a pool it cannot run at all -- consumers branch on `current_pool`/`capped_pool` returning `[]`. | DRO |
| `capped_pool.m` | The pool `run_capped`'s hard fence needs, with a PRIVATE JobStorageLocation so concurrent `-batch` sessions do not collide. Must be called AFTER startup has set the path -- workers inherit it at pool creation. Returns `[]` instead of throwing when no pool can be opened (no Parallel Computing Toolbox, or its licence held by another session); R2025b has no PCT licence here, R2026a does. | DRO |
| `current_pool.m` | The open parallel pool, or `[]` -- without THROWING when the Parallel Computing Toolbox is absent or its licence is held by another MATLAB session. `certify_root`, `verify_with_pumpkyn` and `second_order_pass` read the pool through it; it never creates one (that is `capped_pool`). | DRO (+ internal) |

**Catalog schema, packaging and regression**

| file | what | used by |
|---|---|---|
| `catalog_schema.m` | THE versioned schema authority (v1/v2 min-time; **v3 = objective/γ axis, 2026-09-02**: one catalog per objective, named `axis3`, stored `mf_frac` + `deltav_from_mf`, mandatory `Yj` junctions for minfuel; **v3.1, 2026-09-07: mixed continuation families** — `smoothing.family = 'mixed'` with `smoothing.codes`, per-entry `sheets.family_code` (int8) and `delta_floor` (huberc), validated; and the validator rejects any catalog whose `thruster.c_nd` disagrees with `isp_s`, the 09-02 mislabel). `DRO_tulip/build_minfuel_catalog` is the v3/v3.1 packager (two sources: `minfuel_grid.mat` + the high-γ race via `highgamma_select`). | DRO, HALO_HALO |
| `build_costate_catalog_family.m` | Family-agnostic compact-catalog packager (ND-only quantities, `derive` formulas, `dep_family`/`dep_params` recipes; keeps the legacy `tauDRO` field = departure period so every picker works on every catalog). | DRO, HALO, DPO, HALO_HALO, GTO |
| `golden_cells.m` + `golden_cells_data.mat` | 20-check quality regression: three engine cells (dro/halo/dpo, flown-perturbed 1 N entries + conjugate verdicts) + one harvest cell with REAL collocation duals. **Run after any change to the files above.** A quality drop (iterations, residual) is a failure even when correctness gates pass. **Not runnable from a fresh clone:** the harvest cell loads `DRO_tulip/direct/results/dsweep_12x12_cells.mat` (20 MB, gitignored) and the file puts `DRO_tulip/indirect` on the path (cleanup step 8). | entry point |

**Campaign orchestration (job control, not optimal control)**

| file | what | used by |
|---|---|---|
| `work_queue.m` | A disk work queue so N worker processes on one host pull the next unclaimed unit (a rib column, an entry, a rung) instead of static ranges. Ownership is a process-held kernel lock (`unit_lock`), never transferred on heartbeat age; attempts are recorded before work starts, and a unit that always fails cannot livelock the queue. `tests/test_work_queue`, `tests/test_campaign_processes`. | DRO |
| `campaign_worker.m` | One worker process: claim a unit, run it, VALIDATE its temporary output and publish it in one rename under the unit's lock, release, repeat. A unit function that returns without a valid output is a failed attempt, not a done unit. Launched N at a time by `run_campaign_workers.sh`. | DRO |
| `unit_lock.m` | A process-held file lock (java.nio `FileChannel.tryLock`, POSIX record locking) whose authority lives in a process-wide REGISTRY keyed by the file's device:inode, so a released or stale handle can never act on a newer holder's lock. Measured: refused while held, acquired on release and after `kill -9`. | DRO (+ internal) |
| `publish_atomic.m` | Move a finished file onto its final name in ONE rename(2) (java.nio `ATOMIC_MOVE`), so a reader sees the old file or the new one, never a partial or missing one. Not `movefile`: onto an existing directory it nests the source inside it. The only way a campaign artifact, queue record or heartbeat reaches its final name. | DRO (+ internal) |
| `walk_checkpoint.m` | A resumable checkpoint for a sequential walk (a rib column): the state after the last accepted point, saved atomically at `<output>.ckpt`, carrying the unit's identity and refused by name on a mismatch, so a column killed hours in resumes from its last point. | DRO |
| `campaign_heartbeat.m` | Heartbeats as the only liveness signal, and the five states a monitor must tell apart. No heartbeat is NOT evidence of health -- a missing log once read as "alive, age 0 s", and a `pgrep` on an environment-only tag reported live workers dead. | internal |
| `campaign_status.m` | What a campaign is actually doing, read from ARTIFACTS ONLY (queue outputs and heartbeats, no process matching): units not being worked on, and workers that are not working. Entry point. | entry point |
| `safe_report.m` | Run a REPORTING block so it can never fail the work it reports on: catch, say so loudly, return false. Twice a correctly saved stage was reported FAILED because its summary print threw. | DRO |
| `fmt_num.m` | Format a number at a fixed width, rendering NaN or empty as dashes of the same width instead of throwing -- the two shapes that broke a campaign print, and the reason a movie title or table column does not jump frame to frame. | DRO |

## Tests

`tests/` holds 41 tests. Classified 2026-09-16 by what they call (after
cleanup steps 5–7: `test_conjugate_pole_predict` deleted with its function,
`test_phase_lists` moved with `rib_targets`, the 15 DRO_tulip-only tests moved
to `DRO_tulip/indirect/tests`):

| kind | count | tests |
|---|---|---|
| library only | 25 | `test_arclength_ms`, `test_arclength_ms_thrust`, `test_campaign_processes`, `test_catalog_schema_v3`, `test_conj_fixedtf`, `test_conj_resolve`, `test_cr3bp_minenergy_pmp`, `test_h6_margin`, `test_huber_saltation`, `test_lift_margin`, `test_lift_space_dim`, `test_minfuel_pmp`, `test_mintime_gates`, `test_ms_bvp_extra`, `test_ms_bvp_fixedtf`, `test_ms_tfmin_hom`, `test_nd_propulsion`, `test_periodic_pp`, `test_phase_state`, `test_scalar_verdict`, `test_second_order_parallel`, `test_second_order_pass`, `test_second_order_sidecar_identity`, `test_ss_bvp_accept`, `test_validate_flight` |
| library, through campaign fixtures | 16 | `test_arclength_arrival`, `test_certify_caps`, `test_certify_enforcement`, `test_conj_coverage`, `test_conj_spectrum`, `test_dro_tulip_seed`, `test_entry_notes`, `test_flight_to_junctions`, `test_fly_transfer`, `test_gates_h6_wiring`, `test_gto_family` (GTO_tulip fixture), `test_pmp_pointwise_checks`, `test_seed_from_entry`, `test_sheet_to_catalog_file`, `test_stm_variational`, `test_work_queue` |

No direct test: `harvest_ms_seed` (covered only through `golden_cells`),
`run_capped`, `current_pool`, `flown_control_error`, `true_min_altitude`,
`preflight_screen`, `survey_family_bounds`, `newton_fixed_q`,
`cr3bp_thrust_rhs`, `ctrl_quad`, `assert_periodic_orbit`. `duals_to_costates`
is a delegate; the implementation is tested by `oclib/tests/test_duals_to_costates`. `test_ladder_endpoints` and the campaign-code tests live in `DRO_tulip/indirect/tests`.

Run the relevant tests plus `golden_cells` after touching an engine.

## Conventions

- Pumpkyn house style: `%% Purpose / Inputs / Outputs / Revision History`
  headers, no Code Analyzer pragmas, never `i`/`j` as loop variables.
  **Current state (2026-09-16):** all 58 files use the `%% Purpose` header;
  no `%#ok` pragmas and no `i`/`j` loop variables remain (the 22 `AGROW` /
  `INUSD` warnings the pragmas hid are now visible — preallocate on next
  touch). 13 of 58 files have a `nargin == 0` self-demo.
- Physics only through pumpkyn calls (`tfMinProp`/`tfMinEoM`/getters);
  nothing is ever written into pumpkyn.
- The five standing principles + principle 7 (defenses against silent
  quality degraders) govern all changes — see the SDD.
- **Principle 8 (2026-09-10): two scripts per costate library** — a chain
  script (`DRO_tulip/indirect/build_70mN_library.m` is the model: stages
  named, switchable, file-to-file) and a study script
  (`DRO_tulip/indirect/transfer_study.m`: one transfer, scaffolding exposed,
  necessary and sufficient conditions one at a time). New campaigns copy
  those two and change the parameter blocks. The pedagogical sequence
  before them is `root_origins_study` → `anchor_study` → `transfer_study`
  → `run_phase_torus`.
