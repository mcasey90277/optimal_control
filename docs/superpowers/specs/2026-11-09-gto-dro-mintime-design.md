# GTO -> DRO minimum-time transfers: design

Date: 2026-11-09. Status: approach A approved in conversation (Mike,
2026-11-09); regime decision: catalog regime. The plan argues from this
document.

## 1. Purpose

Earth -> DRO is the leg the orbital network lacks. The four certified
min-time costate catalogs (DRO/HALO/DPO -> tulip, HALO_HALO) and the
GTO -> tulip catalog (deliverable 7) never say how a spacecraft reaches
a DRO. With a GTO -> DRO leg, GTO -> tulip becomes the first pair with
both a direct edge and a relay route (GTO -> DRO -> tulip), the
comparison the orbit_network note left untestable.

This campaign builds GTO -> DRO minimum-time low-thrust transfers in the
spatial Earth-Moon CR3BP with mass loss, in the CATALOG REGIME:
m0 = 150 kg, Isp = 1710 s, thrust ladder from 15 N down. It is delivered
in three stages; THIS spec fully specifies stage 1 and fixes the shape
of stages 2 and 3.

1. **Stage 1 (this spec): one transfer, parameters exposed.** A study
   script `gto_dro_study.m` that builds the GTO departure and the DRO
   arrival from parameters, solves one cell down the ladder to the
   requested thrust through the shared engine, verifies it independently,
   computes the necessary conditions and sufficiency hypotheses with
   PASS/FAIL, draws the transfer and saves the result with metadata.
2. **Stage 2: the catalog.** A chain script over sheets
   orientDeg {0, 90, 180, 270} x DRO tau {0.5, 1, 2, 3}, the standard
   12 x 6 phasing grid, rungs 15 N to 1 N, packaged in schema v2 with the
   orientation-aware pickers.
3. **Stage 3: network integration.** The GTO -> tulip and GTO -> DRO
   catalogs enter `orbit_network/load_catalog_sheets`, and the relay
   closure gains GTO -> DRO -> tulip beside L1 -> L2 -> tulip.

## 2. Decisions

- **Approach A: drive the existing engine, write no solver.** The shared
  ladder engine `costate_common/thrust_ladder_library` already takes a GTO
  departure (`get_family_orbit('gto', struct('orientDeg', ..))`, the
  algebraic Kepler locus, 2001 points uniform in eccentric anomaly,
  spline endpoint rule) and a DRO arrival (`get_family_orbit('dro',
  struct('tau', ..))`, pumpkynPie `getDRO` seed refined by `cont_np`).
  The GTO -> tulip catalog was built with a family-name change only. The
  pipeline is the standing one (principle 6): direct solve
  (`casadi_mintime_dro`, Hermite-Simpson, N = 400) -> fly the control
  (G1b) -> `ms_tfmin` multiple shooting -> pumpkyn `tfMin` acceptance
  UNCHANGED.
- **Catalog regime, not the 25 mN flagship.** 150 kg / 1710 s / 15 N top
  rung. Thrust, mass and Isp are knobs, but the warm recipe and the gates
  are calibrated here; the mN regime is out of scope (it needs the
  Sundman transcription and ships junction states, not bare z8).
- **Warm recipe inherited from the GTO -> tulip campaign, then measured.**
  `thrLock = true` (throttle pinned to 1), `tf0 = 0.30` ND for the TOP
  rung only, `maxCpuSec = 600`. Section 6 of the study script prints the
  returned t_f against tf0 and the rung-to-rung ratio screen
  ([1/3, 3] x guess on STEPS only), so a DRO-specific recalibration is a
  measured change, not a guess.
- **Rendezvous at a chosen DRO phase.** pumpkyn `tfMin` targets a single
  state; the arrival phase fraction is a knob (default 0), as in every
  catalog. Free arrival phase is not built.
- **Departure phase is a TIME fraction since perigee** over one Kepler
  period (0 = perigee, 0.5 = apogee), the GTO catalog's convention.
  Default 0 (perigee): the GTO campaign measured an apogee cold start at
  15 N landing in a bad local minimum.
- **Two scripts per library** (principle 8): stage 1 ships the study
  script; stage 2 ships the chain script. Front-door functions remain
  the engine.
- **Style.** pumpkyn house style for every function (header quartet
  Purpose / Inputs / Outputs / Revision History with the
  `M. Casey (c) date / Copyright Coorbital Inc.` lines, `if nargin == 0`
  self-demo, no `i`/`j`, no `norm`, no `%#ok`, `=` aligned in parameter
  blocks). The study script follows the SDA driver pattern Mike pointed
  at (`SDA/run/run_ekf_tulipVelfo_cr3bp_resacc.m`): a `%% knobs` block
  first with units in trailing comments, flat structs, light-theme
  `tiledlayout` figures with explicit axes handles and `DisplayName` on
  every line, `fprintf` summaries, one results folder per run holding a
  `.mat` with a `run_meta` struct. Its engine calls go to pumpkyn-headered
  functions (`SDA/est/run_ekf_cr3bp_resacc.m` is the model for the
  function side). Study-script discipline from `DRO_tulip/indirect/
  transfer_study.m`: numbered sections, every condition computed IN the
  script with a stable ID and PASS / FAIL / NOT APPLICABLE, one
  tolerance block.

## 3. Folder

```
orbit_transfer/GTO_DRO/
  README.md                 what it is, how to run, what it does not claim
  setup_paths.m             this dir + lib/ + costate_common + oclib + casadi;
                            pumpkynPie startup if pumpkyn.cr3bp.tfMin is absent
                            (copy of GTO_tulip/catalog/setup_paths.m)
  gto_dro_study.m           STAGE 1: the single-transfer study script
  lib/gto_dro_entry.m       single-cell driver through thrust_ladder_library
  lib/earth_clearance.m     Earth-altitude floor check on a flown arc
  tests/test_gto_dro_entry.m
  tests/test_earth_clearance.m
  tests/test_gto_dro_study.m   headless run at the demo operating point
  tests/run_tests.m
  results/                  git-ignored: per-run folders, scratch sheets
```

Path rule: `costate_common`, `oclib` (parent of `+oc`), `~/casadi-3.7.0`,
pumpkyn + pumpkynPie. No DRO_tulip or GTO_tulip folder on the path. The
`.gitignore` gets `orbit_transfer/GTO_DRO/results/`.

## 4. The engine interface stage 1 relies on (verified 2026-11-09)

```
P = thrust_ladder_library(outMat, opts)
  opts.depFamily 'gto', opts.depParams struct('orientDeg', ., 'sma_km', ., 'ecc', .)
  opts.arrFamily 'dro', opts.arrParams struct('tau', .)
  opts.rungs [high .. low] N, .ispS, .m0kg, .N, .floorKm (lunar), .gateKm,
  opts.nD, .nA, .sD0, .sA0, .cells, .thrLock, .tf0, .maxCpuSec, .maxIter,
  opts.msSeg, .msWallS, .accTol, .logFile
  P.TF .FLYKM .ACCDZ .RES .WALL .OK   [nD x nA x nRung]
  P.Z8 [8 x nD x nA x nRung], P.rungs, P.sD, P.sA, P.meta
[tau, rv, info] = get_family_orbit(family, p)      % 'gto' | 'dro' | ...
[x, seam, dx]   = phase_state(t, y, opts)           % the endpoint rule
[z, info]       = ms_tfmin(rv0, rvf, seed, Tmax, c, muStar, opts)
seed            = seed_from_z8(z8, rv0, K, Tmax, c, muStar)
g               = mintime_hypothesis_gates(z8, rv0, Tmax, c, muStar, opts)
out             = oc.ms_conjugate_test(info, spec)
[zEnd, out]     = oc.fly_control(z0, tGrid, rhs, opts)
sol             = pumpkyn.cr3bp.tfMin(rv0, rvf, lambda0_tf, Tmax, c, muStar, solverType)
[tau, y, Hf]    = pumpkyn.cr3bp.tfMinProp(tf, y0, Tmax, c, muStar)
p               = cr3bp_lt_params(thrust_N, m0_kg, Isp_s)   % muStar lStar tStar g0 c Tmax
```

`get_family_orbit('gto')` defaults: sma 24446 km, ecc 0.724781150290436,
orientDeg = angle from the Earth -> Moon line to the perigee direction,
rotating-frame sense. DRO tau is the ND period (tau = 1 is 4.43 d,
Jacobi 3.058; tau = 0.5 / 2 / 3 are 2.22 / 8.87 / 13.30 d).

Known engine gap: `casadi_mintime_dro`, `preflight_screen` and
`certify_dro_mintime` check LUNAR altitude only (`floorKm`, 500 km). No
Earth-altitude constraint or check exists anywhere in the engine. Stage 1
adds the CHECK (section 7 below); stage 2 adds a constraint only if
stage 1 or the pilot sheet shows a violation.

## 5. `lib/gto_dro_entry.m`

```
function [entry, gates, P] = gto_dro_entry(orientDeg, tauDRO, depFrac, arrFrac, thrustN, opts)
```

A copy of `GTO_tulip/catalog/gto_entry.m` with the arrival family
switched to the DRO. Behaviour kept verbatim from gto_entry:

- Ladder from 15 N DOWN to `thrustN`: `rungsStd = [15 12 10 7 5 3 2 1.5 1]`,
  `rungs = [rungsStd(rungsStd > thrustN), thrustN]`; higher rungs are
  warm-start stepping stones, the returned entry is the bottom rung.
  A `thrustN` outside [1, 15] N is refused (`gto_dro_entry:thrust`):
  above 15 N there is no standard rung to anchor on; below 1 N is the
  many-rev regime that ships junction states, not bare z8 (section 9).
- Warm recipe defaults `.thrLock true`, `.tf0 0.30`, `.maxCpuSec 600`;
  `.N 400`, `.floorKm 500`, `.gateKm 100`, `.accTol 1e-6`, `.ispS 1710`,
  `.m0kg 150`, `.msSeg [12 24]`.
- 1 x 1 grid: `nD = nA = 1`, `sD0 = depFrac`, `sA0 = arrFrac`.
- `.outMat` scratch path (default `results/scratch/gto_dro_o<ddd>_tau<..>_<thrust>.mat`,
  deleted on success unless `.keepScratch`).
- `.conjTest` [false]: re-solve `ms_tfmin` seeded at the accepted z8 with
  `conjTest = true`, sets `gates.conjPass`.
- `.gtoParams` struct passthrough (`sma_km`, `ecc`) so the GTO is a
  parameter, not a constant; default = get_family_orbit's defaults.

Returns:

- `entry`: `.ok`, `.z8` [8 x 1], `.tf_nd`, `.tf_d`, `.mf_kg`, `.dV_kms`,
  `.thrustN`, `.rungsWalked`, `.tfPerRung` [1 x nRung], `.okPerRung`,
  `.rv0`, `.rvf` [1 x 6], `.depFrac`, `.arrFrac`, `.orientDeg`, `.tauDRO`,
  `.m0kg`, `.ispS`.
- `gates`: `.msNormR`, `.flownKm`, `.accDz` (tfMin acceptance |dz|),
  `.wallS`, `.conjPass` (NaN unless `.conjTest`), `.ratioScreen`
  [1 x nRung-1] t_f ratios between consecutive rungs with `.ratioPass`.
- `P`: the raw `thrust_ladder_library` output for the cell (every rung),
  so the study script can show the ladder.

Self-demo (`nargin == 0`): orientDeg 0, tau 1, depFrac 0, arrFrac 0,
15 N. Prints the entry and gates.

Not changed from gto_entry: it never cold-solves a deep rung; it never
falls back to the free-throttle defaults.

## 6. `lib/earth_clearance.m`

```
function [hMinKm, tAtMin, pass] = earth_clearance(tau, y, floorKm, p)
```

Earth altitude along a flown arc: `y` [M x >=6] rotating-frame ND states,
Earth at `[-muStar 0 0]`, radius 6378.137 km (pumpkyn's GTO recipe uses
6378 km for the GTO itself; the clearance check uses the equatorial
radius and says so). Returns the minimum altitude, its time, and
`hMinKm >= floorKm`. Default floor 200 km (below the 350 km GTO perigee
by design: the floor catches an impact-bound arc, not a perigee dip
under thrust). Pure geometry, no propagation; the caller supplies the
flown arc. Self-demo: a circular 6678 km Earth orbit converted with
`fromPCI` gives 300 km to 1e-6.

The same function serves the lunar floor when called with the Moon
centre and radius via `p.center`, `p.radiusKm` fields, so the study
script computes BOTH clearances with one instrument (the engine's own
`floorKm` lunar gate is the second opinion, N7-style).

## 7. `gto_dro_study.m`: sections and conditions

The script is a SCRIPT (not a function), runnable from MATLAB with no
arguments; `study_opts` may pre-exist in the workspace to override knobs
(the orbit_network convention: `clearvars -except study_opts; clc`).

```
%% 0  knobs
thrustN     = 15;        % returned rung [N]; the ladder walks 15 N down to it
m0kg        = 150;       % initial mass [kg]
ispS        = 1710;      % specific impulse [s]
gto.orientDeg = 0;       % Earth->Moon line to perigee [deg], rotating sense
gto.sma_km    = 24446;   % [km]
gto.ecc       = 0.724781150290436;
gto.depFrac   = 0;       % departure phase, TIME fraction since perigee (0.5 = apogee)
dro.tau       = 1;       % DRO period [ND]; 0.5 1 2 3 are the catalog periods
dro.arrFrac   = 0;       % arrival phase, time fraction along the DRO
mesh.N        = 400;     % Hermite-Simpson intervals per rung
recipe.thrLock   = true; recipe.tf0 = 0.30; recipe.maxCpuSec = 600;
tol.msNormR   = 1e-8;  tol.flownKm = 1;    tol.accDz = 1e-6;  tol.gateKm = 100;
tol.H         = 1e-6;  tol.lamM    = 1e-8; tol.earthKm = 200; tol.moonKm = 500;
tol.ratioLo   = 1/3;   tol.ratioHi = 3;    tol.gateVms = 1;   % arrival velocity gate [m/s]
out.dir       = fullfile(here, 'results', sprintf('study_o%03d_tau%g_%gN', ...));
out.visible   = 'on';  out.conjTest = true;
```

Then, each section printing its numbers and its verdicts:

1. **Paths and constants.** `setup_paths`; `p = cr3bp_lt_params(thrustN, m0kg, ispS)`;
   print muStar, lStar, tStar, Tmax (ND), c (ND), 1 TU in days.
2. **Departure: the GTO locus.** `get_family_orbit('gto', gto)`; print
   perigee and apogee altitude (km), Kepler period (d), orientation;
   V1 `gtoGeom`: perigee altitude recomputed from `rv` agrees with
   `sma(1-ecc) - R_E` to 1 km.
3. **Arrival: the DRO.** `get_family_orbit('dro', struct('tau', dro.tau))`;
   print period (d), Jacobi constant, perilune and apolune (km);
   V2 `droClosure`: |rv(end) - rv(1)| < 1e-8 ND (the closure residual the
   catalogs use).
4. **Endpoints.** `phase_state` at `depFrac` and `arrFrac`; print both
   states (ND and km, km/s via `toPCI`), the straight-line separation and
   pumpkyn `minDeltaV(rv0, rvf)` as the energy lower bound.
5. **Solve down the ladder.** `gto_dro_entry(...)` with the knobs; a
   per-rung table (thrust, t_f ND and d, m_f, dV, flown km, ms residual,
   accepted |dz|, wall s); the ratio screen.
   N0 `ladder`: every walked rung `ok`; V3 `ratioScreen`: each STEP's
   t_f within [1/3, 3] of the previous rung's.
6. **Independent flight.** Fly the accepted z8 with pumpkyn `tfMinProp`
   (the propagator every gate uses) AND, separately, fly the direct
   solution's own control with `oc.fly_control` through pumpkyn's
   field (the G1b convention, per-interval). Print both arrival misses.
   N3 `arrival`: tfMinProp flight reaches `rvf` within `tol.gateKm`
   and the velocity within the catalog's gateVms.
7. **Clearances.** `earth_clearance` on the flown arc for Earth
   (`tol.earthKm`) and for the Moon (`tol.moonKm`); print both minima
   and when they occur. N7a `earthFloor`, N7b `moonFloor`. These are the
   two gates the engine lacks for an Earth departure (Earth) or has only
   as a constraint (Moon); here they are measured on the flown arc.
8. **Necessary conditions** from the tfMinProp flight (IDs as in
   transfer_study): N1 `bvp` (ms residual), N2 `H0` (|H(t)| max over the
   arc under tol.H; autonomous free-time), N4 `lamMf` (|lam_m(t_f)| <
   tol.lamM), N5 `adjoint` (finite-difference check of lambda-dot =
   -dH/dx at sampled times against tfMinEoM), N6 `minPrinciple` (the
   applied thrust direction equals -lam_v/|lam_v| and the throttle is 1
   wherever Q_mt > 0).
9. **Sufficiency hypotheses and the conjugate test.**
   `mintime_hypothesis_gates(z8, rv0, Tmax, c, muStar)`: S1 `minLamV > 0`,
   S2 `minQmt > 0`, S3 `dimS == 1`; S4 `conjugate`: `gates.conjPass` from
   the `conjTest` re-solve (NOT APPLICABLE when `out.conjTest` is false).
10. **Acceptance by pumpkyn.** X1 `tfMinAccept`: pumpkyn `tfMin` called
    at the accepted z8 returns within `tol.accDz` of it (the gate that
    makes an entry a library entry; the engine ran it, the script runs it
    again in the open).
11. **Figures** (light theme, white background, explicit axes, saved
    PNG at 150 dpi into `out.dir`):
    - `transfer_3d.png`: rotating frame, GTO locus (blue), DRO (green),
      transfer coloured by mass fraction, Earth and Moon spheres to
      scale, L1 and L2 triangles, departure and arrival markers; a
      second panel from above (x-y) because the problem is near-planar.
    - `timeseries.png`: tiledlayout 4 x 1 over time in days: thrust
      direction components (rotating frame), |lam_v| and Q_mt, mass (kg),
      Hamiltonian residual on a log axis. Earth and Moon altitude on a
      fifth tile.
    - `ladder.png`: t_f (d) and dV (km/s) against thrust for the rungs
      walked, log x axis.
12. **Save and census.** `save(fullfile(out.dir, 'gto_dro_study.mat'),
    'entry', 'gates', 'P', 'knobs', 'verdicts', 'run_meta')` where
    `run_meta` holds date, MATLAB version, git commit of optimal_control,
    pumpkyn and pumpkynPie commits, host, elapsed time. Print the verdict
    table (ID, value, tolerance, PASS/FAIL/N.A.). The script ENDS with
    `assert(all code gates pass)` exactly as transfer_study does: hypotheses
    (S1-S4) are results, code gates (V1-V3, N0-N7, X1) throw when
    `selfCheck.strict` is true (default true when run from tests, false
    interactively, like transfer_study).

## 8. Tests

- `test_gto_dro_entry`: (a) refuses thrustN > 15 with the named error;
  (b) rung planning: thrustN = 15 gives one rung, 7 gives [15 12 10 7],
  1 gives all nine, 4 gives [15 12 10 7 5 4]; (c) the GTO params
  passthrough reaches `get_family_orbit` (spy via the returned
  `P.meta`); (d) SLOW, tagged: the self-demo cell at 15 N returns ok with
  flownKm < 100 and accDz < 1e-6. The slow test runs under `run_capped`
  with a 20 min cap and is skipped (printed, not passed) without a pool.
- `test_earth_clearance`: circular 6678 km orbit gives 300 km to 1e-6;
  an arc that crosses the Earth's radius FAILS the floor; the Moon-centre
  mode reproduces the engine's lunar altitude on a DRO (perilune of the
  tau = 1 DRO to 1 km against `get_family_orbit`'s own `info`).
- `test_gto_dro_study`: runs the script headless (`study_opts.visible =
  'off'`, `out.dir` under tempdir, `selfCheck.strict = true`) at the demo
  operating point; checks the three PNGs and the .mat exist and that
  every code gate in `verdicts` is PASS. Deliberate-failure experiment:
  rerun with `tol.earthKm = 1e6` and assert N7a reports FAIL (so the gate
  can fail).
- `run_tests`: the four suites, one exit code, in the house style.

## 9. What stage 1 does not do (refusals, recorded)

- No new solver, no change to `casadi_mintime_dro`, `thrust_ladder_library`,
  `certify_dro_mintime` or any costate_common file. If the Earth floor is
  violated at any rung, stage 1 REPORTS it; a constraint is a stage 2
  item with its own review.
- No free arrival phase, no free departure true anomaly: both are knobs,
  rendezvous at a chosen state, as in every catalog.
- No 25 mN / 15 kg regime; the script refuses thrustN below 1 N with a
  pointer to the many-rev rule (junction states, Sundman).
- No sheet writing: `gto_dro_entry` has no `writeSheet`; the catalog is
  stage 2's chain script.
- No claim of global optimality: a PASS table says "a certified local
  minimum-time extremal that pumpkyn's tfMin accepts unchanged".

## 10. Stage 2 and stage 3, fixed in shape

**Stage 2, `build_gto_dro_catalog.m`** (chain script, build_70mN_library
style, parameter blocks first, ten switchable stages, per-stage files):
pilot sheet (orient 0, tau 1, 12 x 6, rungs [15 12 10 7 5]) -> tf0
calibration from the pilot (median top-rung t_f) -> all 16 sheets through
`thrust_ladder_library` with `tauDRO` = orientDeg as the SHEET KEY (the
GTO catalog's rule) and the DRO period in `arrParams` -> deep rungs
[3 2 1.5 1] as a separate pass -> conjugate sweep (`conj_catalog_pass`)
-> hypothesis gates (`gates_catalog_pass`) -> package
(`build_costate_catalog_family`, schema v2, `costate_catalog_gto_dro.mat`)
-> audit (3-entry replay) -> pictures -> deliverable zip. Plus the
batched OS-timeout driver (`run_gto_dro_batched.sh`) per the campaign
discipline doc. Orientation-aware pickers copied from `GTO_tulip/catalog`.

**Stage 3, network integration:** `load_catalog_sheets` gains
`gto_tulip` and `gto_dro` rows (dep_family 'gto', dep key = orientDeg,
the 12 x 6 grid shape accepted beside 6 x 6), `orbit_node_table` gains
the four GTO orientations as nodes (sigma = 0 by convention, flagged
non-periodic in the rotating frame), `relay_closure` runs
GTO -> DRO -> tulip, and the note gets the first direct-versus-relay
pair, stated with the lemma's caveat.

## 11. Risks carried into the plan

| risk | where it shows | what the plan does |
|---|---|---|
| tf0 = 0.30 tuned for tulips | section 5 top rung fails or returns a far t_f | knob; the pilot in stage 2 calibrates it; the study prints t_f / tf0 |
| 180 deg orientation weak (54% for tulips) | stage 2 coverage | pilot both 0 and 180 before the fleet |
| 3 to 1 N wall (zero GTO entries for tulips) | stage 2 deep pass | separate pass, documented as open; never a floor claim |
| no Earth-altitude constraint in the engine | N7a | measured at every rung; constraint only on a measured violation |
| bad local minimum from apogee at 15 N | section 5 | default depFrac 0; apogee is a knob with the warning printed |
| the DRO is Moon-centred, the GTO Earth-centred: a long straight-line gap | minDeltaV bound in section 4 | printed as the lower bound beside the solved dV |
| fsolve licence for tfMin | X1 | `solverType` knob; `lsqnonlin` fallback printed as such |

## 12. Acceptance for stage 1

- `gto_dro_study` runs from a clean MATLAB with no arguments at the demo
  point (orient 0, tau 1, perigee, phase 0, 15 N) and every code gate is
  PASS; the hypotheses S1-S4 are printed with their values.
- The three figures read at slide size (fonts >= 12 pt, black axes).
- `tests/run_tests` passes (slow test may print SKIPPED without a pool).
- README states the operating point, the knobs, the gate table and the
  refusals of section 9.
- Committed and pushed on main, one commit per reviewed task.
