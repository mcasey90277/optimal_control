# Task 10 Report — Front door (`run_booster_landing`) + flagship run

Commit: `698ccb8` — "booster_landing: front door + flagship run"

## Implementation

Created three files, exactly per the brief's file list:

- `run_booster_landing.m` — front door. Implements the brief's skeleton
  verbatim in structure (setup_paths -> cfg/P merge -> colloc -> convex ->
  certify+print -> TVLQR -> closed loop (+MC) -> plots/movie -> save ->
  summary), with one addition: the summary block now also prints the
  failure-mode breakdown (`landed/arrest/horizon` counts) and two
  adjudication one-liners (etaT de-rate rationale, vf!=0 rationale) per
  this task's supersession instructions ("Summary block: include the
  failure-mode breakdown and the etaT/vf adjudication one-liners").
- `tests/test_run_front_door.m` — end-to-end contract test on a fast grid.
- `README.md` — what/how-to-run/folder map/spec+plan pointers/adjudication
  summary.

### Deviation from the brief's literal test grid (documented)

The brief's test used `cfg.P = struct('N', 20, 'Nconv', 60)`. Running it
literally hit a real gate failure: `certify_pdg`'s G5 primer-alignment
check measures ~1.72 deg at N=20, over the <1 deg threshold — task-5's
report already found and documented that G5 needs N>=40 to clear that
gate (a genuine coarse-grid discretization effect, not a bug). Since the
brief pre-dates task 5 by definition (per this task's own supersession
note), I changed the test's grid to `N=40, Nconv=90` — the exact coarse
pair `test_certify_nominal.m` already uses and already knows clears all
five gates outright — updated the `size(R.solC.X,2)` assertion from 21 to
41 accordingly, and documented the change as an ADAPTATION note in the
test file itself, in the house style. This preserves the test's actual
intent (fast, and a genuine `cfg.P` override — N=40 is not any function's
internal default) without asserting a criterion the code cannot honestly
meet at N=20.

## TDD evidence

1. Wrote `tests/test_run_front_door.m` first.
2. Ran it before `run_booster_landing.m` existed:
   ```
   Unrecognized function or variable 'run_booster_landing'.
   ```
   Confirmed fail as expected (undefined front door).
3. Implemented `run_booster_landing.m`.
4. First pass attempt (before the grid fix, N=20/Nconv=60) failed
   honestly on `assert(R.rep.all_pass ...)` — G5 primer angle 1.72 deg,
   `ALL GATES: FAIL`. This is exactly the kind of failure worth showing:
   the front door itself was correct (solved, certified, printed, wrote
   products, wired cfg overrides) — the coarse test grid was the problem.
   Fixed by moving to N=40/Nconv=90 (see deviation note above).
5. Re-ran: `test_run_front_door PASS` — `size(R.solC.X,2)==41` (cfg.P.N
   reached the solver), `rep.all_pass==true` (all five gates PASS at
   N=40/Nconv=90), `booster_run.mat` and `pdg_solution.png` written.

## Flagship run — full console summary

Command: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); run_booster_landing"` — one `-batch` call, no args, fresh MATLAB. Wall time: 4:52 (well under the 600000 ms / 10 min budget; no incremental-save restructuring was needed).

```
=== [1/5] Guidance: collocation (N=60) ===
EXIT: Optimal Solution Found.
    tf=16.595 s  mf=26464.59 kg
=== [2/5] Guidance: lossless convexification (golden tf) ===
    tf=16.599 s  mf=26464.16 kg  gap=1.36e-04
=== [3/5] Certification ===

Gate                                    value              threshold   verdict
------------------------------------------------------------------------------------
G1 max HS defect                  1.39238e-07                < 1e-06   PASS
G2 pos residual [m]                0.00884687                    < 1   PASS
G2 vel residual [m/s]             0.000298943                  < 0.1   PASS
G2 mass residual [kg]             8.23742e-06                  < 0.5   PASS
  -> G2 gate                                                           PASS
G2ff below-Tmin frac              1.72212e-16                < 1e-06   PASS
G2ff above-Tmax frac              1.58356e-16                < 1e-06   PASS
  -> G2ff gate                                                         PASS
G3 |dmf| [kg]                         0.42792                    < 1   PASS
G3 |dtf| [s]                       0.00340864                  < 0.2   PASS
G3 traj Linf [m]                     0.372616                 (info)
  -> G3 gate                                                           PASS
G4 lossless gap [m/s^2]           0.000135597        < 1e-4*TmaxG/m0   PASS
G5 bound fraction                    0.991736                >= 0.95   PASS
G5 interior switches                        1                   <= 2   PASS
G5 structure (bound+switch+max-last)                                         PASS
G5 primer angle [deg]                0.574131                    < 1   PASS
  -> G5 gate                                                           PASS
------------------------------------------------------------------------------------
ALL GATES                                                              PASS

=== [4/5] Tracking: TVLQR + closed loop + Monte Carlo ===
  MC 25/200 ... MC 200/200
=== [5/5] Products -> /Users/msc/Desktop/optimal_control/booster_landing/results ===

==================== SUMMARY ====================
tf        16.595 s      fuel  3535.4 kg (mf 26464.6)
gates     ALL PASS
nom miss  0.01 m  vtd 0.98 m/s  (landed=1, stop=touchdown)
MC        200 runs, success 99.5%  (landed 199 / arrest 1 / horizon 0)
design    guidance thrust de-rate etaT=0.87 (ceiling 735.15 kN of Tmax=845.00 kN) --
          adjudicated 2026-08-09 so the tracker keeps thrust headroom against
          dispersions (see booster_params.m P.etaT note).
          terminal vf=-1.50 m/s (not 0) -- v(tf)=0 is singular under Tmin>weight;
          legs absorb touchdown speed up to the 2.0 m/s vtd_max gate
          (adjudicated 2026-08-08, see booster_params.m P.vf note).
=================================================
```

## Per-criterion verification (updated wording per this task's supersessions)

1. **Both solvers converged** — CONFIRMED. Both IPOPT logs end
   `EXIT: Optimal Solution Found.` (colloc N=60, convex golden-section
   inner solves).
2. **All gates PASS** — CONFIRMED. `rep.all_pass = true`; every row in
   the printed table above is PASS, no skips (both solvers ran, so G3/G4
   were evaluated, not skipped).
3. **Final masses agree < 1.0 kg** — CONFIRMED. `G3_dmf = 0.428 kg`
   (`|26464.59 - 26464.16|`), consistent with the adjudicated ~0.43 kg
   Taylor-bound model-error measurement cited in the supersession note,
   well under the 1.0 kg gate.
4. **Throttle plot shows bang-bang min-max with terminal max-throttle
   arc** — CONFIRMED by direct inspection. I rendered
   `results/pdg_solution.png` (2688x2250) and looked at panel (b),
   "throttle (bang-bang, both solvers)": both curves sit flat at
   `Tmin/Tmax=0.40` from t=0 to ~t=6.1s, step vertically to the
   `etaT=0.87` guidance ceiling at the single interior switch, and hold
   flat at that ceiling through touchdown at t~16.6s — a clean min-first
   / max-last two-segment bang-bang step, matching G5's own measurement
   (1 interior switch, bound_frac=0.9917, structOk=true).
5. **MC success >= 95% at 200 runs** — CONFIRMED. 99.5% (199/200):
   `n_landed=199, n_arrest=1, n_horizon=0`. Cross-checked against
   `results/footprint.png`, which independently renders
   "success rate: 99.5% (199/200), landed=199 arrest=1 horizon=0,
   vtd: mean=1.05 max=8.29 gate=2.0 m/s" and shows all 199 successes
   inside the 15 m pad-radius circle, the single failure (an arrest, not
   a miss) just inside the 3-sigma ellipse.
6. **`results/landing.mp4` plays clean** — CONFIRMED via `VideoReader`
   (not just file existence): `Width=1280 Height=720 FrameRate=30.000
   Duration=12.000 NumFrames=360` — matches the locked 1280x720
   divide-by-16 frame size (the documented fix for the H.264
   diagonal-streak artifact) and the default `opts.duration=12,
   opts.fps=30` (360 = 12*30, exact). Extracted three frames (t=0, t~6.5s
   mid-flight over the switch region, t~11.9s near touchdown at alt=0.0m)
   and visually inspected each — clean two-panel render (booster/thrust
   arrow + throttle trace with the etaT/Tmin/Tmax bound lines), no
   diagonal color shear in any frame.

## Full suite results (10 test files, all pass)

Ran every `tests/test_*.m` in one MATLAB session (fresh `clearvars` between
each, via a scratch runner script — not committed, scratchpad only):

```
test_certify_nominal           PASS
test_closed_loop_nominal       PASS
test_colloc_smoke              PASS
test_convex_lossless           PASS
test_dynamics_jac              PASS
test_monte_carlo_small         PASS
test_params                    PASS
test_run_front_door            PASS
test_tvlqr_riccati             PASS
test_viz_smoke                 PASS
ALL TESTS PASS
```

(Note: the task instructions said "all 8 tests now" — the actual count in
`tests/` is 10 files including this task's new one; all 10 are green. The
"8" figure in the instructions appears stale relative to the current repo
state; reported here as the real, current count rather than silently
matched to the stale number.)

## Files changed

- `/Users/msc/Desktop/optimal_control/booster_landing/run_booster_landing.m` (new)
- `/Users/msc/Desktop/optimal_control/booster_landing/tests/test_run_front_door.m` (new)
- `/Users/msc/Desktop/optimal_control/booster_landing/README.md` (new)
- Generated (git-ignored, not committed): `booster_landing/results/booster_run.mat`,
  `pdg_solution.png`, `footprint.png`, `landing.mp4`

## Self-review

- Front door matches the brief's skeleton structurally; the only
  substantive addition is the failure-mode breakdown + etaT/vf
  adjudication lines in the summary block, per this task's explicit
  supersession instruction — not scope creep.
- `cfg.P` override is verified two ways: the fast test's
  `size(R.solC.X,2)==41` assertion (N reaches `solve_pdg_colloc`), and
  the flagship's own `P.N=60` default banner line printed and matching
  the solved grid.
- Viz functions are called bare (no output capture) in the front door
  per the binding note — confirmed by reading the file: neither
  `plot_pdg_solution(...)`, `plot_footprint(...)`, nor
  `movie_landing(...)` calls capture a return value, so each closes its
  own figure handle (their nargout==0 contract) — no figure-handle leak
  across a long `-batch` run.
- `certify_pdg`/`print_certify_report` are called exactly as their own
  headers document (report-only, front door decides pass/fail from
  `rep.all_pass`) — no re-implementation of gate logic in the front door.
- README's adjudication summary is a condensed, cross-referenced version
  of the two `booster_params.m` notes (etaT, vf) plus the G3 gate
  reclassification from `certify_pdg.m` — not new claims, just pointers +
  the load-bearing numbers a cold reader needs before diving into the
  full comments.
- Test-grid deviation (N=20->N=40) is documented in-line in the test file
  itself (ADAPTATION-style comment matching house convention) as well as
  here, so a future reader hitting the same G5 coarse-grid gate does not
  re-discover it from scratch.
- Did not touch any file outside `booster_landing/`; `git add` was
  file-by-file (never `-A`), matching the repo's dirty working tree in
  unrelated projects.

## Concerns

- None blocking. The single MC failure (1/200, an "arrest" — a vertical
  velocity-null just above the pad rather than a true touchdown) is
  within the campaign's own documented failure taxonomy
  (`sim_closed_loop.m`'s ADAPTATION note) and does not affect the 99.5%
  success rate clearing the 95% gate.
- The task instructions' "all 8 tests" is stale (real count is 10); noted
  above so it doesn't look like an unexplained discrepancy.
- Phase 2 (drag) is scaffolded in `booster_params.m`/`pdg_dynamics.m` but
  not exercised by any campaign run yet — correctly described as such in
  the README ("not yet exercised by a campaign run — vacuum is the
  certified, flagship configuration"), not overclaimed.

---

## Fix report (2026-08-09 review) — 4 Important + minors

Review verdict was "Approved with 4 Important fixes" (the contract itself
— cfg wiring, no param re-fetch, same-trajectory pipeline — held under
direct verification). All four addressed below, plus the minors. **No
numeric result in this report or in `README.md` changed** — nothing about
the solve, the certification, the tracker, or the Monte Carlo was
touched; only checkpointing, error handling, a stale doc number, a missing
baseline block, and a parameter-consistency gap were fixed. **The full
flagship was deliberately NOT re-run** (per the coordinator's explicit
instruction — nothing numeric changed, so re-running it would only spend
another ~5 minutes to reproduce the same `tf`/`mf`/gates/MC numbers
already recorded above).

### Important 1 — checkpoint save before the viz stage; viz wrapped in try/catch

**Problem:** `booster_run.mat` was written once, at the very end, after
all three viz calls (`plot_pdg_solution`, `plot_footprint`,
`movie_landing`). A throw inside any of them — `movie_landing`'s H.264
divisible-by-16 frame-size path is this campaign's documented bite point
— discarded the entire ~5-minute solve+certify+MC run with nothing on
disk.

**Fix** (`run_booster_landing.m`): the solve+certify+MC block now saves a
checkpoint immediately (`R.P`/`R.when` set, `save(...)`, with a printed
`[checkpoint] ... -> <path>` line) BEFORE the viz stage begins. The three
viz calls are now wrapped in a single `try/catch`: on failure, a `warning`
(id `run_booster_landing:vizFailed`) reports which viz call/message failed
and reassures the caller that the checkpointed products are already on
disk, then execution falls through to a final re-save (`R.when` refreshed,
any viz products that DID complete before the failure folded in) instead
of the function erroring out and losing everything. A deliberately
triggered viz failure was not manufactured for this fix (would require
breaking `movie_landing` itself, outside this file's remit) — the change
was verified structurally (checkpoint line appears in every run's console
output below, e.g. the Tmax-override probe) and by code inspection: the
`try` block's only content is the three viz calls, and the `catch` never
rethrows, so any exception from any of the three is now caught and
reported rather than propagated.

### Important 2 — README fast-test timing corrected

**Problem:** README said "~1 s" for `test_run_front_door`; the coordinator
measured 23.1 s.

**Fix:** Re-measured directly — `time matlab -batch ... test_run_front_door`
gives `24.765s total` wall clock (MATLAB startup + N=40/Nconv=90 solves +
TVLQR Riccati integration + 6-run MC + 2 plots). README now says
"measured ~24 s wall time including MATLAB startup — dominated by the
TVLQR Riccati integration and a 6-run Monte Carlo, not the coarse NLP
solves, which each take well under 1 s" — an honest number with the
breakdown that explains why it isn't the "~1 s" the brief's language
implied.

### Important 3 — expected-baseline block added to README

**Fix:** added an "Expected flagship result (baseline for a cold
reproduction)" section with a table: `tf=16.595 s`, `mf=26464.6 kg` /
`fuel=3535.4 kg`, gates ALL PASS, `G3 |dmf|=0.43 kg` (gate <1.0), nominal
`miss=0.01 m`/`vtd=0.98 m/s`/landed=true, MC `99.5%` (199 landed / 1
arrest / 0 horizon) — exactly the numbers from this task's own flagship
run, dated, with a one-line note on what magnitude of deviation on a
re-run should prompt investigation (citing the `P.tf_hi` IPOPT/BLAS
threading-sensitivity note already in `booster_params.m` as a documented
precedent for small run-to-run digit noise).

### Important 4 — cfg.P re-derive for Tmin/gvec, README narrowed

**Problem:** `cfg.P` is merged onto `booster_params()` as a flat
field-by-field overwrite. `booster_params.m` DERIVES `P.Tmin =
0.40*P.Tmax` and `P.gvec = [0;0;-P.g0]` from their parents. Overriding
only the parent (e.g. `cfg.P.Tmax`) left the dependent stale — an
inconsistent `(Tmax, Tmin)` pair that `certify_pdg` would certify without
any complaint (nothing in the gate math checks that the two are
consistent with each other, only that the solved trajectory respects
whatever annulus it was given).

**Fix (a), `run_booster_landing.m`:** after the `cfg.P` merge, two
guarded re-derive blocks fire only when the PARENT was overridden and the
DEPENDENT was NOT also explicitly overridden (an explicit dependent
override always wins over the auto-derivation):
```
if isfield(cfg.P, 'Tmax') && ~isfield(cfg.P, 'Tmin')
    P.Tmin = 0.40 * P.Tmax;
    fprintf('    [cfg.P] Tmax overridden -> re-derived Tmin = 0.40*Tmax = %.1f N\n', P.Tmin);
end
if isfield(cfg.P, 'g0') && ~isfield(cfg.P, 'gvec')
    P.gvec = [0; 0; -P.g0];
    fprintf('    [cfg.P] g0 overridden -> re-derived gvec = [0;0;-g0] = [%.4f; %.4f; %.4f]\n', P.gvec);
end
```

**Fix (b), README:** the `cfg.P` paragraph now states explicitly that the
merge is a flat overwrite, names the two known parent/dependent pairs,
describes the auto-re-derive behavior and its precedence rule, and warns
that any OTHER field not documented as derived in `booster_params.m` is a
plain unchecked overwrite.

**Verification (deliberate probe, requested by the coordinator, full
output below).** Ran a coarse (`N=40, Nconv=90`, `doMC=false` for speed)
campaign with `cfg.P.Tmax = 900e3` (up from the default 845e3) and NO
`cfg.P.Tmin` override:
```
    [cfg.P] Tmax overridden -> re-derived Tmin = 0.40*Tmax = 360000.0 N
...
PROBE CHECK: P.Tmax=900000.0  P.Tmin=360000.0  expected Tmin=360000.0  match=1
PROBE CHECK: default P.Tmin (unmodified) would have been 338000.0 (STALE if not re-derived)
```
`R.P.Tmin` is exactly `0.40*900000 = 360000` N, matching the expected
re-derived value and NOT the stale default-Tmax-derived `338000` N a flat
overwrite alone would have left in place. The note also printed exactly
once, at the point of the merge, as designed. (This probe run's gates
predictably FAIL G3 — a 900 kN engine on a trajectory tuned around 845 kN
is a genuinely different physical problem, not a regression; the probe
was about the re-derive mechanism, not about producing a passing
alternate configuration.)

### Minors (all folded in)

- **`tests/test_run_front_door.m`:** added
  `assert(numel(R.mc.ok) == cfg.Nrun, ...)` (Nrun actually reached
  `run_monte_carlo`), `assert(~isfile(fullfile(od,'landing.mp4')), ...)`
  (cfg.doMovie=false is honored — no movie file materializes),
  `assert(isfile(fullfile(od,'footprint.png')), ...)` (footprint plot
  written). Dropped the dead `|| ischar(R.rep.G3_pass)` OR clause from the
  gate assertion — the front door always supplies `solV`, so G3 is never
  `'skipped'` there, and a live version of that clause could only ever
  mask a real `all_pass=false`, never legitimately fire.
- **README:** cfg example's `outdir` changed from the bare string
  `'results'` to an absolute path, with an explicit note that relative
  paths resolve against MATLAB's current working directory, not the
  campaign folder, and that the built-in default is already an absolute
  path built from `mfilename('fullpath')`.
- **README adjudication dates:** kept the literal dates as recorded in
  the source comments, but added a line naming the design spec as the
  single source of truth if a date or number here and in code ever
  disagree.

### Re-verification after the fixes

- `test_run_front_door` re-run clean after all four Important fixes +
  minors: `test_run_front_door PASS` (all assertions, including the three
  new ones, hold).
- Full 10-file suite re-run clean:
  ```
  test_certify_nominal           PASS
  test_closed_loop_nominal       PASS
  test_colloc_smoke              PASS
  test_convex_lossless           PASS
  test_dynamics_jac              PASS
  test_monte_carlo_small         PASS
  test_params                    PASS
  test_run_front_door            PASS
  test_tvlqr_riccati             PASS
  test_viz_smoke                 PASS
  ALL TESTS PASS
  ```
- `booster_landing/results/` (the committed flagship's generated products)
  is untouched by any of this fix work — all probes/re-runs used
  `tempdir`-scoped `outdir`s, so the flagship numbers quoted throughout
  this report and the README's new baseline table remain exactly what was
  measured in the original flagship run, not regenerated.

### Files changed (this fix pass)

- `/Users/msc/Desktop/optimal_control/booster_landing/run_booster_landing.m` (modified)
- `/Users/msc/Desktop/optimal_control/booster_landing/tests/test_run_front_door.m` (modified)
- `/Users/msc/Desktop/optimal_control/booster_landing/README.md` (modified)
