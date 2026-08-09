# Booster Landing — 3-DOF Powered-Descent Guidance + TVLQR

Falcon-9-class landing-burn campaign: solve the min-fuel powered-descent
guidance (PDG) problem two independent ways — direct Hermite-Simpson
collocation and lossless convexification — certify they agree and are
PMP-shaped, close the loop with a time-varying LQR tracker, and
Monte-Carlo the landing under dispersions.

## What this is

Phase 1 (vacuum, flat Earth, constant gravity — convexification is exactly
valid here): a 7-state (position, velocity, mass), 3-control (thrust
vector) point-mass model. Guidance solves against a de-rated thrust
ceiling (`P.etaT = 0.87`) so the tracker keeps real authority in the
[Tmin, Tmax] annulus during dispersed flight; the guidance targets a small
descent rate at touchdown (`P.vf = [0;0;-1.5]` m/s, not exactly zero)
because zero velocity at zero altitude is singular once minimum throttle
already exceeds vehicle weight. Both adjudications, with the measurements
behind them, are documented in `lib/booster_params.m`.

Phase 2 (opt-in exponential-atmosphere drag, `P.drag.on`) re-solves the
same problem with drag on, warm-started from the Phase-1 vacuum solution
(`run_booster_landing(struct('phase2', true))`) — vacuum stays the
certified, flagship configuration; Phase 2 is a comparison campaign on top
of it, not a replacement. Flagship measurement (`P.N=60`, 2026-08-09):
drag **saves** 434.7 kg of fuel versus vacuum (3535.4 -> 3100.7 kg, braking
against the atmosphere) at a cost of +0.94 s of flight time (16.595 ->
17.531 s); gates G0/G1/G2/G5 pass (G3/G4 skipped — no convex twin exists for
the drag-on problem, lossless convexification is only exact in vacuum);
wind Monte Carlo (200 runs) succeeds 99.0% of the time. Phase 2 is
certified against exactly the same thresholds as Phase 1 apart from the two
gates that need a convex twin: this campaign once loosened G5's
primer-alignment threshold from 1 deg to 10 deg under drag, but that
loosening was **withdrawn 2026-08-09** when an external code review found
its cause was a time-base bug in the gate (the segment defect dual was
compared against the *node* control instead of the segment *midpoint*
control it is the multiplier for). Corrected, the drag primer angle is
1.2e-6 deg — the `acos` machine-precision floor, identical to vacuum, at
every grid — and the gate is now 0.01 deg for both. See
`certify/certify_pdg.m`'s "G5 primer TIME BASE" note for the measurement
tables.

## How to run

Full nominal campaign (a few minutes: two solves, TVLQR, a 200-run Monte
Carlo, two plots and a movie):

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/booster_landing'); run_booster_landing"
```

No arguments needed — `run_booster_landing` self-paths via `setup_paths`
and runs from a fresh MATLAB. All knobs are an optional `cfg` struct, never
by editing files:

```matlab
R = run_booster_landing(struct( ...
    'doMovie', true, ...        % write <outdir>/landing.mp4   [def true]
    'doMC',    true, ...        % run the Monte Carlo(s)        [def true]
                                 % (gates BOTH Phase 1's R.mc and, when
                                 % phase2 is also true, Phase 2's R.mcD +
                                 % phase2_footprint.png -- same knob, both
                                 % phases)
    'Nrun',    200, ...         % MC sample count               [def 200]
    'phase2',  false, ...       % also run the drag-on re-solve stage
                                 % (R.solD/.repD/.ctrlD/.Pd, +.mcD if doMC)
                                 % and write phase2_vac_vs_drag.png +
                                 % phase2_footprint.png              [def false]
    'outdir',  '/Users/you/Desktop/optimal_control/booster_landing/results', ...
                                 % product dir; RELATIVE paths resolve
                                 % against MATLAB's CURRENT working
                                 % directory, not this campaign folder
                                 % [def <campaign>/results, an absolute
                                 % path built from mfilename('fullpath')]
    'P',       struct('N', 60, 'Nconv', 120)));  % booster_params() overrides
```

`cfg.P` is field-merged onto `booster_params()` **field by field**, so any
single parameter (grid size, thrust de-rate, boundary conditions,
dispersion seed, ...) can be overridden without touching a library file.
This is a flat overwrite, not a re-derivation: `booster_params.m` computes
a few fields FROM others (`P.Tmin = 0.40*P.Tmax`, `P.gvec = [0;0;-P.g0]`),
so overriding a parent (`cfg.P.Tmax`, `cfg.P.g0`) without also overriding
its dependent would otherwise desync the pair. `run_booster_landing`
detects exactly that case and re-derives the dependent automatically
(printing a `[cfg.P] ... re-derived ...` note when it fires) *unless* the
caller also supplies the dependent explicitly, which always wins. This
covers the two known parent/dependent pairs above; a `cfg.P` override of
any *other* field not documented as derived in `booster_params.m` is a
plain, unchecked overwrite — as always, know what you're overriding.
`R` returns everything the campaign produced:
`.P .solC .solV .rep .ctrl .out0 .when`, plus `.mc` if `cfg.doMC` is true
(Phase 1's Monte Carlo, same gate as Phase 2's `.mcD` below — `.mc` is not
unconditional), plus, when `cfg.phase2` is true, `.solD .repD .ctrlD .Pd`
(and `.mcD` too, if `cfg.doMC` is also true) — the same struct is saved to
`outdir/booster_run.mat` either way.

Fast end-to-end smoke test (coarse N=40/Nconv=90 grid, measured ~24 s
wall time including MATLAB startup — dominated by the TVLQR Riccati
integration and a 6-run Monte Carlo, not the coarse NLP solves, which
each take well under 1 s — exercises the whole pipeline including the
cfg-override contract):

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_run_front_door"
```

(The tests self-bootstrap with
`addpath(fullfile(fileparts(mfilename('fullpath')), '..')); setup_paths;`,
which is cwd-independent. Elsewhere this README shortens that to
`addpath('..')` for readability — that short form only works when MATLAB's
cwd is already `tests/`, so do not copy it into a new test.)

Full unit-test suite (all `tests/test_*.m`, one MATLAB invocation per
test): note each test must run in its OWN `-batch` call, not looped with
`run()` inside a single session -- `run()` executes a script in the
CALLER's workspace (these are scripts, not functions), so a driver script
that itself uses variables like `k`/`fs` would have them clobbered by
whatever a test happens to name its own script-level variables (found
during the task-11 close-out review while writing a suite runner):

```
for t in test_params test_dynamics_jac test_colloc_smoke \
         test_convex_lossless test_certify_nominal test_tvlqr_riccati \
         test_closed_loop_nominal test_monte_carlo_small test_viz_smoke \
         test_phase2_drag test_run_front_door; do
  /Applications/MATLAB_R2025b.app/bin/matlab -batch \
    "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; $t" \
    || echo "FAILED: $t"
done
```

## Folder map

| Folder | Contents |
|---|---|
| `run_booster_landing.m` | Front door: solve -> certify -> track -> MC -> plots/movie -> `results/booster_run.mat` |
| `setup_paths.m` | Adds campaign dirs + CasADi 3.7 to the MATLAB path |
| `lib/` | `booster_params` (single source of truth for all constants/BCs/grid/MC settings), `pdg_dynamics` (3-DOF EOM), `solve_pdg_colloc` (Hermite-Simpson NLP, CasADi+IPOPT), `solve_pdg_convex` (lossless convexification + golden-section tf search), `hs_quad_ctrl` (flyable per-segment control reconstruction, shared by certify/tvlqr/sim), `tvlqr_design` (phase-scheduled time-varying LQR), `sim_closed_loop` (truth-model dispersed sim, altitude-indexed tracking), `run_monte_carlo` (dispersed landing campaign) |
| `certify/` | `certify_pdg` (gates G0 time-base consistency, G1 discrete defect, G2 continuous residual, G2ff feedforward feasibility, G3 cross-method agreement, G4 losslessness, G5 PMP bang-bang + primer structure), `print_certify_report` (gate table) |
| `viz/` | `plot_pdg_solution` (2x2 trajectory/throttle/mass/speed comparison), `plot_footprint` (MC landing scatter + dispersion ellipse), `movie_landing` (booster + throttle-trace MP4, 1280x720 locked), `plot_vacuum_vs_drag` (Phase-2 1x3 trajectory/throttle/mass overlay, vacuum vs drag) |
| `tests/` | One `test_*.m` per unit (self-bootstrapping: `addpath(fullfile(fileparts(mfilename('fullpath')), '..')); setup_paths;` — cwd-independent, throws on failure), `test_run_front_door.m` (end-to-end contract on a fast grid), `test_phase2_drag.m` (Phase-2 drag solve + gates + `plot_vacuum_vs_drag` smoke) |
| `results/` | Generated products (git-ignored except `.gitkeep`): `booster_run.mat`, `pdg_solution.png`, `landing.mp4`, plus `footprint.png` when `cfg.doMC`; from a `cfg.phase2` run, `phase2_vac_vs_drag.png`, and `phase2_footprint.png` only when `cfg.phase2` **and** `cfg.doMC` are both true (it plots the wind Monte Carlo, which `doMC` gates) |

## Expected flagship result (baseline for a cold reproduction)

Nominal grid (`P.N=60`, `P.Nconv=120`), no-args run, MATLAB R2025b,
measured 2026-08-09 (4:52 wall time in one `-batch` call):

| Quantity | Value |
|---|---|
| `tf` (both solvers agree to ~3 ms) | 16.595 s |
| `mf` (collocation) / fuel used | 26464.6 kg / 3535.4 kg |
| Gates G0-G5 | ALL PASS |
| G3 `|dmf|` (cross-method mass agreement) | 0.43 kg (gate: < 1.0 kg) |
| G3 trajectory L-inf (cross-method position agreement) | 0.373 m (gate: < 5.0 m) |
| G5 primer angle (PMP alignment, midpoint semantics) | 1.2e-6 deg (gate: < 0.01 deg) |
| Nominal closed loop | miss 0.01 m, `vtd` 0.98 m/s, landed=true |
| Monte Carlo (200 runs) | 99.5% success (199 landed / 0 arrest / 1 depleted / 0 horizon) |

A re-run should land close to these numbers (IPOPT/BLAS threading can
shift the last digit or two, see `booster_params.m`'s `P.tf_hi` note for a
documented case of that sensitivity) — a large deviation (gates failing,
success rate materially under 95%, `tf` off by more than ~0.1 s) means
something changed and is worth investigating before trusting the run.

Phase 2 (`run_booster_landing(struct('phase2', true))`), same grid,
measured 2026-08-09 (Phase 1 numbers above unchanged — Phase 2 only adds a
drag-on re-solve on top):

| Quantity | Vacuum (Phase 1) | Drag (Phase 2) |
|---|---|---|
| `tf` | 16.595 s | 17.531 s (`dtf` = +0.94 s) |
| fuel used | 3535.4 kg | 3100.7 kg (`dfuel` = 434.7 kg SAVED by drag) |
| Gates | G0-G5 ALL PASS | G0/G1/G2/G5 PASS (G3/G4 skipped, no convex twin) |
| Wind Monte Carlo (200 runs) | 99.5% (vacuum, no wind) | 99.0% (`P.drag.on` also enables `sig.wind`) |

## Spec / plan pointers

- Design spec: `../docs/superpowers/specs/2026-08-08-booster-landing-design.md`
- Implementation plan (task-by-task): `../docs/superpowers/plans/2026-08-08-booster-landing.md`
- Per-task briefs/reports: `../booster_landing/process/task-{1..12}-{brief,report}.md`

## Adjudication summary (the two numbers a cold reader needs)

Dates below are as recorded in the source comments at adjudication time;
the design spec (`docs/superpowers/specs/2026-08-08-booster-landing-design.md`)
is the single source of truth if any date or number here and in code ever
disagree.

- **`P.etaT = 0.87`** — guidance solves against `etaT*Tmax = 735.15 kN`,
  not the engine's full 845 kN. The un-de-rated min-fuel optimum rides its
  thrust ceiling with only ~2% net-deceleration margin left over, so a
  -5% thrust dispersion was unrecoverable by any tracker (measured 67
  m/s touchdown, reproduced under two different controllers — an
  authority limit, not a tuning one). 0.87 is the smallest de-rate that
  clears a 7-case dispersion battery with margin (worst case 1.67 m/s
  against the 2.0 m/s gate), at a cost of ~82 kg (~1.9%) of the
  propellant budget. Full sweep in `lib/booster_params.m` and the
  task-7b report.
- **`P.vf = [0;0;-1.5]` m/s** (not `[0;0;0]`) — with `P.Tmin` already
  exceeding vehicle weight at every mass on this trajectory, a
  zero-velocity-at-zero-altitude terminal condition is singular: the
  closed loop cannot complete the last stretch of descent and instead
  arrests just above the pad. Targeting a small legs-absorbed descent
  rate (within the `P.vtd_max = 2.0` m/s mission gate) is what makes a
  genuine touchdown achievable at all. See `lib/booster_params.m` and the
  task-7 report.
- **G3's `|dmf| < 1.0` kg gate** is a *measurement* threshold over the
  convex solver's Taylor-linearized mass-bound model error (~0.70-0.94 kg,
  non-shrinking under refinement), not an agreement tolerance the two
  solvers should be tuned to close further. See `certify/certify_pdg.m`.
