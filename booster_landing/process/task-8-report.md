# Task 8 Report: `run_monte_carlo`

## 1. Summary

Implemented `lib/run_monte_carlo.m` and `tests/test_monte_carlo_small.m` per
the task-8 brief, with the Task-7 evolution applied on top: `mc.ok` now also
requires `out.td.landed`, and the failure-mode breakdown
(`n_landed`/`n_arrest`/`n_horizon`) plus the raw drawn dispersions are
carried in the output for post-hoc attribution of any failure tail.

TDD followed: failing test written and confirmed to fail on undefined
`run_monte_carlo`, then the implementation was written and the test passes.
A 50-run probe at production settings (N=60, default 1σ dispersions) was
then run and is reported below.

## 2. Implementation

`lib/run_monte_carlo.m`:

- `opts.Nrun` (default 200), `opts.sig` (defaults: `.r0=[100;100;50]` m,
  `.v0=[10;10;10]` m/s, `.thrust=0.015`, `.isp=0.01`, `.wind=[10;10;0]` m/s,
  wind drawn only when `P.drag.on`), overridable field-by-field.
- `rng(P.seed)` called once at top; per-run draws in fixed order (`dr0`,
  `dv0`, `thrust_scale`, `isp_scale`, `[wind]`), so two calls with the same
  `Nrun`/`P.drag.on` are bit-identical, and a smaller `Nrun` call reproduces
  a prefix of a larger call's draws.
- Each run calls `sim_closed_loop(sol, ctrl, P, d)` and reads `out.td`
  directly (`.r`, `.vtd`, `.m`, `.miss`, `.landed`, `.stop`) rather than
  recomputing miss from `.land` — DRY with what `sim_closed_loop` already
  computes.
- **Task-7 evolution — `mc.ok`:** `out.td.landed && out.td.miss <
  P.pad_radius && out.td.vtd < P.vtd_max && out.td.m >= P.mdry`. An arrest
  or horizon termination is scored a failure unconditionally, because
  `sim_closed_loop`'s own header documents that an arrest's `vtd` is
  structurally optimistic (`vz=0` *is* the arrest event, so it can't also be
  evidence against the vertical-speed failure it represents).
- **Task-7 evolution — failure-mode tracking:** `mc.landed` (Nrun×1
  logical), `mc.stop` (Nrun×1 cellstr), `mc.n_landed`/`mc.n_arrest`/
  `mc.n_horizon` (counts, sum to Nrun).
- **Added for attributability** (task instruction 3): `mc.dr0`, `mc.dv0`
  (Nrun×3), `mc.thrust_scale`, `mc.isp_scale` (Nrun×1), `mc.wind` (Nrun×3,
  only when `P.drag.on`) — the raw draws themselves, so a failure
  population's draw magnitudes can be read straight off `mc` without
  re-deriving them.
- Progress print every 25 runs, as in the brief.

No deviation from the brief's core algorithm/output contract (`mc.land`,
`mc.vtd`, `mc.mprop`, `mc.ok`, `mc.success_rate` all present with the
specified shapes and semantics); the additions are strictly additive.

## 3. Wall-time check (task instruction 4)

Before running the 20-run test, measured a single dispersed `sim_closed_loop`
call at the production grid (N=P.N=60): **0.216 s**. The brief's "~2 min per
run" concern does not apply to this campaign — no `odeset` changes were made
(tolerances untouched, as instructed). A full 200-run MC at this rate is
~1-2 minutes, well inside the front door's (Task 9/10) budget.

## 4. TDD evidence

**Step 1 — failing test**, `tests/test_monte_carlo_small.m` run before
`run_monte_carlo.m` existed:

```
Error: Unrecognized function or variable 'run_monte_carlo'.
Error in test_monte_carlo_small (line 16)
```//confirmed FAIL as expected.

**Step 2 — implementation**, then re-run:

```
test_monte_carlo_small PASS  success=65%  landed=20 arrest=0 horizon=0
```

(20 runs at the test's coarse N=30 grid — all 20 genuinely landed; 7 of 20
failed the miss/vtd/mass gates, giving 65% `success_rate`. `mc1.success_rate
> 0` and the failure-mode-sum/determinism assertions all passed.)

## 5. 50-run probe (production grid, default dispersions)

Ran `run_monte_carlo(sol, ctrl, P, struct('Nrun', 50))` with `sol` from
`solve_pdg_colloc(P, struct('N', P.N))` (N=60, the production grid) and
`ctrl = tvlqr_design(sol, P)` — i.e. the actual shipped controller, not the
test's coarse N=30 shortcut.

**Result: 50/50 success (100%). 0 failures.**

| metric | value |
|---|---|
| success_rate | 100.0% (50/50) |
| n_landed | 50 |
| n_arrest | 0 |
| n_horizon | 0 |
| wall time | 77.0 s total, 1.54 s/run |
| lateral draw \|dr0_xy\| range (all 50) | 6.1 – 293.6 m (mean 135.7 m) |

Failure population: **empty** — nothing to attribute. As a substitute
sanity check (does the campaign actually get stressed by this batch), the
10 largest lateral-offset draws, all of which landed cleanly:

| run | \|dr0_xy\| (m) | dr0_z (m) | stop | miss (m) | vtd (m/s) | mprop (kg) |
|---|---|---|---|---|---|---|
| 31 | 293.6 | -11.72 | touchdown | 7.25 | 1.165 | 202.4 |
| 5 | 291.8 | 41.35 | touchdown | 8.26 | 0.893 | 209.9 |
| 8 | 283.6 | -40.13 | touchdown | 8.02 | 1.184 | 264.2 |
| 13 | 246.6 | 20.39 | touchdown | 6.88 | 1.058 | 313.5 |
| 45 | 227.5 | 95.38 | touchdown | 5.37 | 1.055 | 307.6 |
| 49 | 220.5 | -7.29 | touchdown | 4.49 | 0.941 | 313.6 |
| 7 | 217.0 | 40.02 | touchdown | 5.36 | 1.021 | 174.6 |
| 25 | 212.7 | 53.20 | touchdown | 3.68 | 1.316 | 191.9 |
| 30 | 210.4 | 8.25 | touchdown | 5.23 | 1.227 | 349.0 |
| 24 | 202.0 | -121.20 | touchdown | 6.34 | 1.000 | 523.9 |

Every one of these lands well inside the 15 m pad radius and well under the
2.0 m/s gate, including a draw at 293.6 m lateral offset (nearly 3σ) plus a
sizeable altitude offset in the same draw.

**Reconciling this with task instruction 3's "known tail" note:** the brief
told me to expect lateral outliers beyond ~100 m to fail, citing task-7's
round-5 measurement ("lateral capability envelope: clean to ~100 m of
offset; vtd starts crossing 2.0 m/s somewhere past ~100 m"). Checking the
task-7 report directly: that characterization is dated to round 5, **before**
task-7b's altitude-indexed-guidance overhaul and the `P.etaT=0.87` de-rate.
Task-7b's own report (§7, "Standing note for task 8") repeats the same
round-5 number forward as the standing expectation for task 8, but the round-5
sweep was run under the round-5 (pre-altitude-indexing, pre-de-rate)
controller, not the one this campaign actually ships. My probe uses the
current, final controller (altitude-indexed tracking + etaT de-rate,
`tvlqr_design`/`sim_closed_loop` as of the task-7b commit) and measures a
materially better lateral envelope — no failures out to ~294 m in this
50-draw sample.  **This is a genuine, measured improvement carried over from
task-7b that the standing note in the task-7 report did not re-verify**, not
a bug in this task's harness (the harness reproduces the deterministic
20-run test correctly, and the coarser N=30 test grid does show a realistic
35% failure rate, so the mechanism is exercised and working).

**Practical implication:** at `Nrun=50` and the default 1σ envelope, the
tail this task was told to expect did not appear. It may still appear at
`Nrun=200` (deeper into the tail, e.g. runs beyond 3σ, or unfavorable
combinations of thrust/Isp bias with large lateral+altitude draws), which is
exactly what the Task 9/10 front door's full run is for. I did not run
`Nrun=200` here per the instruction to leave that to the front door.

## 6. Files changed

- `lib/run_monte_carlo.m` (new)
- `tests/test_monte_carlo_small.m` (new)

## 7. Self-review

- Determinism: verified structurally (`rng` once, fixed draw order per run)
  and by the test's `isequal(mc1.land, mc2.land)` plus the added
  `n_landed`/`n_arrest`/`n_horizon` equality checks across two calls.
- `mc.ok` cannot be true when `mc.landed` is false — asserted directly in
  the test (`all(mc1.ok(~mc1.landed) == false)`), guarding against exactly
  the class of bug (arrest scored as success via optimistic vtd) that
  motivated the Task-7 evolution.
- No loop variable named `i`/`j` (`krun` per the brief, `kk`/`fn` elsewhere
  in test/probe scripts only, not in the committed library file's core
  loop — the library uses `krun` and `k` for the small field-copy loop,
  matching the brief and the rest of the codebase's style).
- Function header follows the project's documented MATLAB header format
  (purpose, INPUTS, OUTPUTS with sizes, REFERENCES).
- Tolerances in `sim_closed_loop` were not touched, per instruction.
- Only `lib/run_monte_carlo.m` and `tests/test_monte_carlo_small.m` were
  staged and committed; other unrelated modified/untracked files elsewhere
  in the repo (other projects) were left alone.

## 8. Concerns / open items for later tasks

1. **The "known tail" is stale.** The task-8 brief (and task-7's own
   standing note) expects lateral-outlier failures past ~100 m; measured
   behavior with the final task-7b controller shows none out to ~294 m in
   a 50-run sample. Task 9/10 (front door, full Nrun=200) should not assume
   the tail is where task-7 said it would be — report whatever the 200-run
   number actually is rather than anchoring to this stale expectation.
2. Because the 50-run probe had zero failures, I could not demonstrate the
   failure-population reporting machinery (`mc.dr0` etc.) against a live
   failure — it's exercised in the code path but only validated against a
   sample where nothing failed. The 20-run test at the coarser N=30 grid
   does hit failures (7/20) and passed the `landed`-gate assertions there,
   so the mechanism is real, but a first look at an actual MC failure entry
   (draw values, stop mode) will happen for the first time at Nrun=200.
3. `mc.mprop` at these draws ranges roughly 175–525 kg above `P.mdry` even
   under large dispersions — plenty of margin; not a concern, just noting
   the campaign is not propellant-limited in this dispersion regime.

## 9. Test commands (for reproduction)

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_monte_carlo_small"
```

---

## 10. Fix report (post-review, 2026-08-09)

Review verdict: Needs fixes -- one Important, three cheap minors. All four
addressed.

### Important 1: `opts.sig` field/shape validation

`run_monte_carlo.m`'s `sig` merge loop previously accepted any field name
in `opts.sig` unconditionally (a typo'd field, e.g. `thrust_sig`, was
silently absorbed as a dead field while the matching default stayed live)
and never checked shapes (a `1x3` `r0` would broadcast against a `3x1`
`randn` draw into a `3x3`, corrupting `d.dr0` and only surfacing as a
mysterious crash deep inside `sim_closed_loop`).

Fixed: `fieldnames(opts.sig)` is now checked against the known set
`{r0,v0,thrust,isp,wind}` -- any unrecognized name throws
`run_monte_carlo:unknownSigField` naming both the bad field(s) and the
valid set. Each recognized field is then shape-checked before being
merged: `r0`/`v0`/`wind` must be `3x1` (`run_monte_carlo:badSigShape`),
`thrust`/`isp` must be scalar. Both checks run before the loop, so a bad
`opts.sig` fails fast with zero wasted `sim_closed_loop` calls.

Verified with two bad-input probes in the test (see excerpt below): a
typo'd field name and a wrong-shaped `r0` both error cleanly with the new
identifiers, and the message names the valid field set.

### Minor 1: N=30 test-grid comment

`test_monte_carlo_small.m` now has an explicit comment at the
`solve_pdg_colloc(P, struct('N', 30))` call warning a future reader that
this is a coarse smoke grid, not the production grid (`P.N=60`), and that
its ~65% success rate is not a robustness claim -- pointing at the task-8
report's 50-run production-grid result (100%) for the real number.

### Minor 2: full-struct determinism check

The two `isequal(mc1.land, mc2.land)` / per-field failure-mode-count
equality assertions were replaced by one `assert(isequal(mc1, mc2), ...)`
covering every field (`land`, `vtd`, `mprop`, `landed`, `stop`, `ok`,
`dr0`, `dv0`, `thrust_scale`, `isp_scale`, counts, `success_rate`) in a
single check. Confirmed passing after the change (mc1 and mc2 are still
bit-identical structs under the same seed).

### Minor 3: unified failure-mode-count derivation

`mc.n_landed` now derives from `mc.stop` via `strcmp(mc.stop,
'touchdown')`, the same mechanism as `n_arrest`/`n_horizon`, rather than
from `mc.landed` separately (previously two different source fields fed
three counts that happened to agree by construction but were never
checked against each other). The `mc.landed`/`mc.stop` invariant
(`sim_closed_loop`'s documented contract: `.landed <=> stop=='touchdown'`)
is now explicitly guarded: `assert(mc.n_landed == nnz(mc.landed), ...,
'run_monte_carlo:landedStopMismatch', ...)` -- a future
`sim_closed_loop` change that broke that contract would now be caught
here instead of silently propagating.

### Test re-run (covering test)

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_monte_carlo_small"
```

Relevant tail of output:

```
  bad-field probe error (expected): opts.sig has unrecognized field(s) {thrust_sig} -- valid fields are {r0, v0, thrust, isp, wind} (a typo here would otherwise silently leave the default dispersion live)
  bad-shape probe error (expected): opts.sig.r0 must be a 3x1 vector (got [1 3])
test_monte_carlo_small PASS  success=65%  landed=20 arrest=0 horizon=0
```

### Files changed (fix pass)

- `lib/run_monte_carlo.m` -- `opts.sig` field/shape validation added;
  `n_landed`/`n_arrest`/`n_horizon` unified onto `mc.stop` with a guarded
  `landed`/`stop` invariant; header doc updated to describe the
  validation.
- `tests/test_monte_carlo_small.m` -- N=30 coarse-grid comment added;
  determinism check strengthened to full-struct `isequal(mc1, mc2)`;
  two bad-input probes added (unknown field, wrong shape) asserting the
  correct error identifiers.

### Self-review (fix pass)

- Both new error paths use named identifiers (`run_monte_carlo:...`) so
  callers can `catch` on identifier rather than parsing message text, and
  both were exercised by the test, not just written and assumed correct.
- Validation runs strictly before `rng(P.seed)` and the run loop, so a bad
  `opts.sig` never burns a solve or advances the RNG stream -- calling
  `run_monte_carlo` again with corrected `opts.sig` after a caught error
  still reproduces the documented deterministic sequence.
- No change to `sim_closed_loop` or its tolerances; no change to the
  `mc.ok` success criteria; no change to Nrun/sig defaults. Scope held to
  the review's four points.
