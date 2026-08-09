# Task 11: Phase 2 — atmosphere + drag — close-out

## Close-out session (2026-08-09)

**Handoff discrepancy found and resolved.** The close-out brief stated the Phase-2
flagship had already completed with fresh `results/booster_run.mat` (Phase-2
fields), `results/phase2_vac_vs_drag.png`, and `results/phase2_footprint.png` on
disk. On inspection this was **not** the disk state at session start:
`results/booster_run.mat` had only Phase-1 fields (`P ctrl mc out0 rep solC solV
when`), and neither `phase2_*.png` existed. `ps aux` showed **two duplicate**
`matlab -batch "run_booster_landing(struct('phase2', true))"` processes already
running concurrently (PIDs 39958 started 09:32, 40298 started 09:35), racing to
the same output paths — evidently a prior `rm -f booster_run.mat phase2_*.png`
+ relaunch (`phase2_flagship3.log`) stacked on top of an still-running earlier
launch (`phase2_flagship2.log`), with the "stuck monitor" from the handoff
belonging to one of these. Per the task's explicit constraint ("NO
re-running the flagship"), no third run was launched. Instead: waited
(non-blocking `ps` poll, single bounded Bash call) for both existing processes
to exit naturally — both finished within ~2 min of the check, converging on
**identical** numbers (both solved the same deterministic problem), so the race
was benign; whichever finished last (flagship3, 09:50:28) wrote the final
`booster_run.mat` and both PNGs.

## Verified Phase-2 numbers (final `results/booster_run.mat`)

- **Phase 1 unchanged:** `solC.tf = 16.595192 s`, `rep.all_pass = 1` — matches
  the prior-implementer baseline exactly.
- **Phase 2 fields present:** `solD`, `repD`, `mcD` all populated.
- `solD.tf = 17.530799 s`, `solD.mf = 26899.276 kg`.
- `P.m0 = 30000 kg` → fuelVac = 3535.411 kg, fuelDrag = 3100.724 kg.
- **dfuel = 434.687 kg** (vac − drag): drag **saves** fuel, as expected
  (atmospheric braking assists the burn).
- **dtf = +0.936 s** (drag − vac): drag-on flight takes slightly longer.
- Throttle structure (`Um`): bang-bang between `Tmin = 338000 N` and
  `etaT*Tmax = 735150 N` (etaT=0.87 de-rated ceiling) — min–max structure,
  consistent with Phase 1's adjudicated suicide-burn structure.
- `repD` gates: G1_pass=1, G2_pass=1, G2ff_pass=1, G3_pass='skipped',
  G4_pass='skipped', G5_pass=1, all_pass=1 — exactly the G1/G2/G5-only
  contract from the brief (no convex twin under drag).
- `mcD`: success_rate = 0.990 (99.0%), n_landed=200, n_arrest=0, n_horizon=0
  out of 200 wind-Monte-Carlo runs.

## Test suite

All 11 tests run individually (fresh `matlab -batch` per test, to avoid
script-workspace collisions from a shared-workspace runner): **11/11 PASS**
— `test_params`, `test_dynamics_jac`, `test_colloc_smoke`,
`test_convex_lossless`, `test_certify_nominal`, `test_tvlqr_riccati`,
`test_closed_loop_nominal`, `test_monte_carlo_small`, `test_viz_smoke`,
`test_phase2_drag` (PASS, `fuel vac=3535.4 drag=3100.9 kg` — the coarse
N=30 warm-start check, close to but not identical to the N=60 flagship
numbers above by design), `test_run_front_door` (fast N=40/Nrun=6 config,
not the flagship — confirmed before running).

## Commit

Everything for Task 11 was uncommitted at session start. Committed as
`96a1ab4` ("booster_landing: Phase 2 atmosphere -- drag re-solve, wind MC,
vacuum comparison"), staging only campaign files:
- `booster_landing/run_booster_landing.m` (modified — `cfg.phase2` branch)
- `booster_landing/tests/test_phase2_drag.m` (new)
- `booster_landing/viz/plot_vacuum_vs_drag.m` (new)
- `booster_landing/certify/certify_pdg.m` (modified — G5 primer relaxation
  under `P.drag.on`, documented in-file with measured evidence)
- `booster_landing/certify/print_certify_report.m` (modified — G5 primer
  row prints as info-only under drag)
- `booster_landing/lib/solve_pdg_colloc.m` (modified — `opts.tol` defaults
  to 1e-11 under `P.drag.on`, measured fix for a G1 mass-row defect that
  cleared the 1e-6 kg gate only at the tighter tolerance)

The three extra modified files (beyond the brief's three) are legitimate
in-flight fix-round work with detailed measured justification already
written into their own file headers — not new design work performed by
this close-out session.

## Concerns for the controller

- **Process hygiene**: two duplicate flagship processes ran concurrently
  writing to the same output paths. They happened to converge on identical
  numbers so no corruption resulted, but this was fragile (a non-deterministic
  MC seed or an interrupted write could have produced an inconsistent file).
  Worth flagging so future sessions check `ps aux` for stray MATLAB batch
  jobs before assuming "stuck" means "safe to ignore."
- No other design issues found; nothing routed back.

## Correction (review response round, 2026-08-09)

The "benign race" framing above (paragraph 1) overstated how much was left
to chance. Per the review: the two processes converging on **identical**
numbers was not luck, it was **guaranteed** — `run_monte_carlo.m` calls
`rng(P.seed)` once, deterministically, before all draws (line 83), and the
IPOPT solve itself is deterministic given identical inputs/warm-starts, so
two independent processes running the same deterministic problem from the
same params were never going to diverge numerically. The real (avoided, not
realized) risk was never "wrong numbers" — it was **interleaved writes**: two
processes calling `save(..., 'booster_run.mat')` at genuinely the same
instant could have produced a torn/corrupt `.mat` file regardless of
whether the underlying numbers agreed. That risk didn't materialize here
(the two `save` calls landed ~2 min apart per the logs), but the report
above was wrong to call the situation "benign" in a way that implied the
numeric agreement was the reassuring part.

## Fix round (review response, 2026-08-09)

Review verdict on the first close-out was **Needs fixes** — 5 Important
items plus minors. All addressed, re-verified with fresh synchronous
MATLAB runs (no flagship re-run: nothing numeric changes, see below), and
committed.

**1. `cfg.doMC` silently ignored by the Phase-2 branch (Important).** Fixed
in `run_booster_landing.m`: `R.mcD = run_monte_carlo(...)` and
`plot_footprint(R.mcD, ...)` are now both wrapped in `if cfg.doMC`; the P2
summary line and the final SUMMARY block's MC(wind) line are also guarded
(print "(MC skipped, cfg.doMC=false)" / "skipped" instead of dereferencing
a field that no longer exists). Verified with a dedicated fast probe
(`cfg.phase2=true, cfg.doMC=false`, N=30/Nconv=60, `doMovie=false`):
confirmed `R.mc` and `R.mcD` both absent, `footprint.png` and
`phase2_footprint.png` both absent, `R.solD/.repD/.ctrlD/.Pd` and
`phase2_vac_vs_drag.png` all still present. Full probe output landed in
this fix round's test run (see below).

**2. G5 primer gate was REMOVED under drag, not loosened (Important).**
This was a real bug: `rep.G5_pass = rep.G5_structOk` alone (the first
version) would have silently passed a 179-degree sign-flip regression --
exactly the failure class the sign-flip fix elsewhere in `certify_pdg.m`
exists to catch. Fixed: `rep.G5_pass = rep.G5_structOk && rep.G5_primer_deg
< 10` under `P.drag.on` (was `< 1` in vacuum, unchanged). Also fixed
`print_certify_report.m`'s G5 primer row to print as a real scored `prow`
against the 10 deg threshold (was an unscored "(info, drag-relaxed)" row
that no longer matched what `rep.G5_pass` actually checked). Header
comments in both files rewritten to say LOOSENED not removed, with the
10-deg margin rationale (>3x over the measured 2.61 deg production-grid
value). Verified: a direct `certify_pdg` + `print_certify_report` probe on
a fresh drag solve prints `G5 primer angle [deg] (drag-loosened)  5.11215
< 10   PASS`, consistent with `rep.G5_pass=1`.

**3. `plot_vacuum_vs_drag` rendering defects (Important).** Two real
defects, both caught by actually rendering and looking at the PNG (not
just re-running the tests, which only prove the file exists):
  - (a) The annotation box occluded panel (a)'s title/data. Went through
    two attempts before landing on a robust fix -- worth recording because
    both intermediate failures were only caught by re-rendering: first
    attempt hand-set `tl.Position` to carve out a footer strip, which
    instead broke tiledlayout's own automatic title-space reservation (the
    sgtitle and the top pixel row of panel (a) both got clipped by
    exportgraphics' tight crop). Second attempt tried `tl.RowHeight` for a
    short second row, which doesn't exist as a property on this MATLAB's
    `TiledChartLayout` (threw outright). Final fix: an 8-row GridSize
    tiledlayout, the 3 panels each spanning rows 1-7 via
    `nexttile(tl, col, [7 1])`, a genuine footer TILE (not a
    figure-normalized `annotation` box) at row 8 spanning all 3 columns,
    with `tl.Position` left at its automatic default throughout -- so
    tiledlayout keeps doing its own title/spacing bookkeeping and nothing
    is hand-guessed in normalized units. Re-rendered and inspected: title
    fully visible, no panel occluded, footer text in its own clear strip.
  - (b) The three `yline` calls (Tmin/Tmax, etaT, mdry) were leaking into
    their axes' legends as `data1`/`data2`. Fixed: `'HandleVisibility',
    'off'` added to all three.
  - Both confirmed by `imread`-based inspection of the actual rendered PNG
    (3068x1233 initially, now 2876x1306 with the 8-row layout), not just a
    passing test.

**4. README was false (Important).** "Phase 2 ... scaffolded ... not yet
exercised by a campaign run" replaced with the actual flagship result
(dfuel 434.7 kg, dtf +0.94 s, gates status, wind-MC 99.0%) and a pointer to
the G5 loosening note. Added `cfg.phase2` (and clarified `cfg.doMC`'s
Phase-2 scope) to the `cfg` documentation block, added `.solD/.repD/.ctrlD/
.Pd/.mcD` to the `R` output description, added a Phase-2 row to the
"Expected flagship result" baseline table, and added `plot_vacuum_vs_drag`
/ `test_phase2_drag` to the folder map. Also fixed (found independently
while building this fix round's own suite runner, not part of the
reviewer's list): the README's documented "full unit-test suite" command
looped `run(...)` over all test files inside a SINGLE `matlab -batch` call
-- `run()` executes a script in the CALLER's workspace (these are scripts,
not functions), so the driver's own loop variables get clobbered by
whatever a test happens to name its script-level variables. Replaced with
a shell loop that gives each test its own `-batch` invocation (the pattern
this close-out session actually used throughout, both times).

**5. No fast test covered the new branch/viz (Important).**
`tests/test_phase2_drag.m` now also calls
`plot_vacuum_vs_drag(solC, solD, <tempfile>)` on the solC/solD it already
solved (no extra solve), closes the figure, and asserts the PNG exists.
Re-ran: `test_phase2_drag PASS  fuel vac=3535.4 drag=3100.9 kg`, plus a
follow-up `imfinfo`/`imread` check confirmed the written PNG is a valid,
correctly-sized image.

**Minors folded in:**
- `isfield(P,'drag')` guard added in `certify_pdg.m` at both the
  `rep.drag_on` assignment and the G5 threshold branch (parity with
  `solve_pdg_colloc.m`'s existing `isfield(P,'drag') && P.drag.on` pattern)
  -- a caller with a pre-drag-field `P` no longer throws, it's just treated
  as vacuum.
- `R.ctrlD` and `R.Pd` now stored on the Phase-2 branch (parity with Phase
  1's `R.ctrl`/`R.P`) -- populates on the next full run, no re-run needed
  for this close-out.
- `plot_vacuum_vs_drag.m`'s "read off solD.P rather than assumed" claim is
  now a real `assert(isequal(...))` on Tmin/Tmax/etaT between `solC.P` and
  `solD.P`, not just header prose.
- `certify_pdg.m`'s ~48-53 header prose fixed: the original text
  attributed the G5-drag-relaxation instruction and the primer SIGN-flip
  instruction to the same brief quote ("companion instruction for the
  OTHER failure mode"); they are two separate, unrelated brief
  instructions, and conflating them was a documentation bug in the
  original fix, now corrected.
- The "benign race" correction is its own section above per the reviewer's
  request (not folded inline, to keep the original session's account
  intact and the correction clearly attributed).

**Test suite re-run (fix round):** `test_phase2_drag` PASS (with new viz
assert), `test_certify_nominal` PASS (vacuum G5 threshold/logic unchanged
by the `isfield` guard), `test_run_front_door` PASS, plus the dedicated
`doMC=false` phase2 fast probe PASS (all assertions above), plus the other
7 unit tests re-run individually: PASS. **11/11 (+1 dedicated probe, not
part of the 11) all PASS.** No flagship re-run -- nothing numeric changes
(all edits are gate-scoring/viz/docs/guard fixes on top of already-verified
solves), consistent with the task's synchronous-only, no-flagship-re-run
constraint.

**Commit:** `booster_landing: fix Phase-2 review findings -- cfg.doMC
guard, G5 primer loosened-not-removed, plot_vacuum_vs_drag rendering,
README truth, viz test coverage` (see git log for hash), touching
`run_booster_landing.m`, `certify/certify_pdg.m`,
`certify/print_certify_report.m`, `viz/plot_vacuum_vs_drag.m`,
`tests/test_phase2_drag.m`, `README.md`.
