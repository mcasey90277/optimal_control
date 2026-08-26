# Stage B design: GTO → tulip min-time costate catalog

**Date:** 2026-08-26 · **Status:** approved (design session, this date) ·
**Owner:** Mike + Claude session · **Deliverable target:** costate catalog
deliverable 7 (Darin / pumpkyn `tfMin`).

## 1. Goal

A schema-v2 compact costate catalog for the **GTO → south-pole tulip**
transfer — the entry problem no shipped catalog covers (all six existing
catalogs depart from periodic cislunar orbits; GTO is where a real mission
starts). A consumer supplies *(GTO orientation, petal count, departure
anomaly, tulip arrival phase, thrust)* and receives a converged min-time PMP
solution (`z8`) that `pumpkyn.cr3bp.tfMin` accepts unchanged.

Secondary accomplishments: fills the GTO row of the goal matrix
(`doc/transfer_problem_space.md`); produces the phasing × thrust min-time
map the min-fuel paper's results section wants; exercises the Stage-A-proven
shared engine (`costate_common/ms_bvp` via `ms_tfmin`) at production scale
on this pair.

## 2. Adjudicated decisions (with rationale)

| decision | choice | why |
|---|---|---|
| propulsion regime | **Darin standard: 150 kg, Isp 1710 s, rungs [15 12 10 7 5 3 2 1.5 1] N** (Mike, 2026-08-25) | few-rev regime where the pipeline is proven; the 25 mN flagship regime is a later walk-down |
| departure-phase axis | **orientation sheets + anomaly grid** (Mike, 2026-08-26) | GTO has two angles (anomaly on the ellipse; ellipse orientation vs the Earth–Moon line). Making orientation a SHEET KEY maps exactly onto the shipped catalogs' `sheets(τ_dep, Np) × grid` structure — same schema, pickers, packager |
| grid extents | **12 anomaly × 6 arrival per sheet** (Mike, 2026-08-26), 4 orientations {0°, 90°, 180°, 270°} (the −25° flagship geometry is reachable by the picker's nearest-sheet warning, and step-2's smoke test uses it directly) × Np {5, 7, 9, 12} = 16 sheets | double departure resolution where the physics is steepest (perigee vs apogee departure differ enormously); otherwise the deliverable-4/5 template |
| entry format | compact **z8** (schema v2) | legitimate at few-rev standard thrust per the 2026-08-26 identifiability adjudication. **Binding rule for the future mN extension: many-rev entries must ship full ms junction states, not bare z8** (bare z8 determines t_f only to ~1e-4 ND / arrival to ~100 km at ~40 revs) — record in the README |
| pm branch | pm = −1 only | pm = +1 is the exact z-mirror (standard catalog convention) |

Scale: 16 × 12 × 6 × 9 = 15,552 attempts; expected ~85–92% pair
solvability → **~8,000–13,000 entries**.

## 3. Architecture

### 3.1 The one new library piece: `'gto'` in `costate_common/get_family_orbit`

`get_family_orbit` is THE single point where every catalog engine, densifier,
packager, picker and `conj_catalog_pass` turns (family, params) into
endpoint states. Adding a `'gto'` pseudo-family makes the entire existing
machinery work on this pair with near-zero modification.

Contract:

- **params:** `.sma_km`, `.ecc` (flagship values: 350 × 35,786 km ellipse),
  `.orientDeg` — THE single planar orientation angle: from the Earth→Moon
  line to the ellipse's perigee direction, measured in the frame's rotation
  sense. (There is no separate argp: in the planar rotating frame the
  ellipse has one orientation angle; the flagship's "argp −25°" is exactly
  its orientation at the chosen epoch, i.e. orientDeg = −25° in this
  convention. The unit test pins this by reproducing the flagship
  `mintime_params`/`gto_tulip_endpoints` departure state from
  orientDeg = −25° at the matching anomaly.)
- **returns:** `[tau, rv]` — the **locus of departure states at FIXED
  orientation**, parametrized by anomaly fraction over one Kepler period.
  **Sampling (REVIEW FIX): uniform in ECCENTRIC anomaly, ≥2000 points**,
  with `tau` from Kepler's equation — a sparse or time-uniform sample
  under-resolves perigee on an e≈0.72 ellipse and corrupts interpolated
  departure states (GPT+Gemini, critical). `tau` is a LOCUS PARAMETER only,
  never an epoch offset (documented invariant).
  Each sample is computed **algebraically** (ellipse state at anomaly ν in
  Earth-centered inertial → rotate by `orientDeg` → convert to the rotating
  frame, subtracting frame rotation), NOT by propagation. The
  non-periodicity of GTO in the rotating frame never enters: orientation is
  a key, not a flow consequence. `info.periodND` = the Kepler period (so
  phase-fraction → days derivations work unchanged).
- The existing `interp1(tau, rv, frac*tau(end))` construction in every
  engine then does exactly the right thing: `sD_frac` = anomaly fraction
  from perigee.

### 3.2 Campaign module: `GTO_tulip/catalog/`

Sibling of the existing `direct/` and `indirect/` (the campaign folder
already exists; the catalog is a third product line inside it):

| file | role |
|---|---|
| `run_gto_catalog.m` | driver, cloned from `HALO_tulip/run_halo_catalog.m`: builds sheet specs (orient × Np), calls `thrust_ladder_library` per sheet, resumable |
| `run_gto_batched.sh` | batched shell runner (nohup + per-pass OS kill), cloned from `run_halo_batched.sh` |
| `setup_paths.m` | this dir + `costate_common` + pumpkyn bootstrap |
| `results/catalog/` | per-sheet .mat + progress files (gitignored .mat, tracked progress txt per convention) |

Engine: `thrust_ladder_library` **unchanged** (it self-bootstraps
`costate_common` and takes departure orbits from `get_family_orbit`).
Seeding: the proven ladder strategy — cold `tfMin` at the 15 N first rung
(few-rev, closes cold), warm-chain down the rungs, neighbor continuation
across the phasing grid, attempt budgets + resume per the matlab-campaign
discipline. Harvest seeding is NOT needed at these thrusts —
**`gen_tulip_mintime` dual capture is out of Stage B scope** (stays queued
for the mN walk-down).

### 3.3 Downstream (existing machinery, small touches)

- **Conjugate sweep:** `conj_catalog_pass` before packaging (standing rule,
  2026-08-23). One small edit: its v1 recipe fallback learns nothing new —
  the packager writes full v2 recipes with `dep_family = 'gto'`, which the
  sweep already reads. Verify `'gto'` flows through; extend the fallback
  only if a gap appears.
- **Packaging:** `build_costate_catalog_family` with `dep_family = 'gto'`,
  `dep_params` carrying {sma_km, ecc, orientDeg}. **`tau_dep`/`tauDRO`
  sheet KEY carries `orientDeg`** (REVIEW FIX 2026-08-26, all three
  reviewers: the picker selects sheets by nearest (tau_dep, Np); the Kepler
  period is identical across orientations and would make 3 of 4 sheets
  unreachable). The Kepler period lives in `dep_params`/`info.periodND` and
  the derive registry. Documented prominently in the README. Sheet ordering
  by (orient, Np).
- **Validation:** `catalog_schema('validate')`; `costate_lib_describe`;
  `golden_cells` untouched (regression only).
- **Deliverable:** via `/library-deliverable` when Mike says ship — zip,
  localized paths, recipient self-test. NOT part of the build definition of
  done.

## 4. Sequencing (implementation-plan order)

1. **`'gto'` family + unit test** (TDD): state correctness against an
   independent `orb2eci`/`fromPCI` construction at several (ν, orient)
   points; anomaly-fraction convention pinned by test; `check_matlab_code`
   clean; `golden_cells` still green (no engine edits expected).
2. **Endpoint smoke:** one cold `tfMin` solve at 15 N from (orient 0°,
   ν = apogee) to the Np = 7 tulip — proves the endpoint construction feeds
   the solver sensibly before any campaign scale.
3. **Pilot sheet:** orient 0° × Np 7, full 12 × 6 × 9. Measures per-cell
   cost and solvability. **Adjudication checkpoint:** if pair solvability
   < ~70% or per-cell cost is a surprise, stop and bring numbers before the
   fleet.
4. **Full campaign:** 16 sheets, batched + monitored, days-scale unattended
   (resume for free). Findings recorded as they land.
5. **Conjugate sweep** over the finished catalog; verdicts stored.
6. **Package + validate + describe**; catalog README (including the
   orientation-axis convention documented for Darin, the identifiability
   rule for the future mN extension); records wave (STATUS_AND_ROADMAP,
   TODOs, register §6 entries per experiment, memory).

## 5. Acceptance criteria (definition of done for the build)

- `'gto'` family unit test green; `golden_cells` 20/20 unchanged.
- Pilot sheet ≥ ~70% pairs solved (else adjudicated before continuing).
- Full catalog: every stored entry passes the three standard gates
  (ms residual ≤ 1e-10, flown arrival < 100 km, tfMin acceptance
  |Δz| < 1e-6) — enforced by the engine, spot-audited by
  `costate_catalog_example`-style replay on ≥ 3 entries across sheets.
- Conjugate verdicts stored for 100% of entries (pass/fail/not-run
  accounting complete, per the 2026-08-23 sweep pattern).
- `catalog_schema('validate')` clean; `costate_lib_describe` output sane
  (grids, ranges, coverage per sheet).
- Coverage, ΔV/flight-time ranges, and the hard-corner map recorded in
  `GTO_tulip/catalog/README.md` + the program records.

## 6. Risks (named in the session)

- **Perigee × long-tulip hard corner:** expect the standard
  shortest-departure × longest-arrival red corner; handled by attempt
  budgets + densify passes, not a blocker.
- **No symmetry discount across orientations:** 90°/270° sheets are
  genuinely distinct solves (unlike pm-mirroring); the 16 sheets are 16
  real campaigns.
- **Stalls:** standard bisection-continuation playbook from the DRO
  campaign; the engines already revisit no-OK cells on resume.
- **No Earth-altitude screen exists in the pipeline** (host review: both
  `preflight_screen` and `true_min_altitude` are lunar-only). GTO transfers
  depart a 350 km Earth perigee; a converged arc could pass through the
  Earth with all gates green. Mitigation (report-only burn-in, house
  pattern): minimum Earth altitude computed and recorded in the pilot
  census and the replay audit; a hard floor is added only if violations
  appear.
- **Anomaly grid spacing** (Gemini): 12 uniform TIME fractions cluster at
  apogee on an eccentric ellipse. The driver passes an explicit
  non-uniform `sD` vector — 12 values uniform in TRUE anomaly mapped to
  time fractions — if the engine accepts explicit grids; else uniform with
  the caveat recorded.
- **Pilot gate** (GPT): judged on per-rung coverage (esp. 1 N) and
  ladder-complete pairs, not any-rung; plus a deliberate kill-and-resume
  test during the pilot.
- **Declined (flagged to Mike):** adding −25° flagship-orientation sheets
  (+25% campaign) — general-geometry catalog; the flagship campaign covers
  −25°; schema v2 permits appending later.
- **Perigee departure at 1 N** may approach multi-rev behavior; if entries
  there start needing large K or showing identifiability symptoms, apply
  the mN rule early (store junction states for those cells) rather than
  bending gates.

## 7. Out of scope (explicit)

- The 25 mN flagship-regime catalog (later walk-down; needs the junction-
  state entry format and `gen_tulip_mintime` dual capture).
- Min-energy / min-fuel objectives (separate program line).
- ELFO targets, other departure ellipses (GTO shape is fixed at the
  flagship's 350 × 35,786 / argp −25°).
- Shipping the deliverable zip (separate, on Mike's word).
