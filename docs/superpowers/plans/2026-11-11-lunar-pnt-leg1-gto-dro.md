# Lunar PNT Mission Leg 1 (GTO to DRO) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Six 5 N transfer stages leave one inertially fixed GTO and deliver six satellites onto the six DRO slots the existing mission sim starts from, with gates, mass and power bookkeeping, a combined timeline, captured transfers and movies.

**Architecture:** New library functions and one study script beside the leg-2 code in `orbit_transfer/lunar_pnt_mission/`. Schedule on the GTO to DRO catalog table (a pure function), polish each transfer with the departure pinned and the arrival sliding along the DRO slot, gate, then build combined products from the leg-1 results and the saved leg-2 results. The shared engine and the GTO_DRO campaign are called, never edited.

**Tech Stack:** MATLAB R2026a through the MCP tools in the shared session; `costate_common` (`ms_tfmin`, `seed_from_z8`, `fly_transfer`, `phase_state`, `nd_propulsion`, `get_family_orbit`, `conj_spectrum`, `mintime_hypothesis_gates`, `thrust_ladder_library`); `GTO_DRO/lib/earth_clearance`; pumpkyn.

**Spec:** `docs/superpowers/specs/2026-11-11-lunar-pnt-leg1-gto-dro-design.md`

## Global Constraints

- All new code under `orbit_transfer/lunar_pnt_mission/` (`lib/`, `tests/`, study script at the folder root). Results under its git-ignored `results/`.
- A fleet run is writing `orbit_transfer/GTO_DRO/results/` and reads the working tree on every pass: NO edits under `orbit_transfer/costate_common`, `orbit_transfer/GTO_DRO/lib`, `orbit_transfer/GTO_DRO/catalog`; nothing written under `orbit_transfer/GTO_DRO/results`; never kill a MATLAB process. Nothing under `proj7`, `pumpkyn`, `pumpkynPie`, `oclib` or the deliverable-9 kit is modified.
- Leg-1 regime, read from the package and never retyped: thrust 5 N (rung index 5 of `[15 12 10 7 5]`), Isp 1710 s, m0 150 kg; `muStar`, `lStar = 389703.264829278 km`, `tStar = 382981.289129055 s` from the package constants. DRO: `get_family_orbit('dro', struct('tau', 1))`, period 1 ND (4.432654 d).
- Package: knob `catalogMat`, default `orbit_transfer/GTO_DRO/results/costate_catalog_gto_dro.mat` (one variable; `sheets(k)` fields `tauDRO`, `dep_params`, `arr_params`, `sD_frac [1x12]`, `sA_frac [1x6]`, `has_solution`, `tf_nd`, `conj_pass` each `[12 x 6 x 5]`, `entry_index`, `z8`). Use every sheet with `tauDRO == 1`. A cell is CLEAN iff `has_solution && conj_pass == 1` at rung 5.
- GTO: the catalog's default elements through `get_family_orbit('gto', struct('orientDeg', o, 'muStar', muStar))` and the periodic `phase_state`; `tauG` is its period in ND. Orientation of the inertial GTO at leg-1 time t: `orient0 + sgn*(180/pi)*t` deg, with `sgn` measured in Task 1, never assumed. GTO phase: `mod(gtoPhase0 + t/tauG, 1)`.
- Slots: satellite k of leg 2 is at DRO phase `(k-1)/6` at leg-2 time 0 (`mission_inputs().droPhase0`). In leg-1 time slot j is at `mod(delta + (j-1)/6 + t/tauD, 1)`; T0 is the smallest t >= last arrival with `mod(delta + t/tauD, 1/6) == 0`, and satellites are labelled so that the one at phase 0 at T0 is leg-2 spacecraft 1.
- Polish thresholds, from leg 2: `ms_tfmin` residual < 1e-8; flown miss < 100 km and < 10 m/s to accept an inner solve; the arrival gate reports achieved maxima. Clearance floors: Earth 200 km, Moon 500 km, on the located minimum.
- A satellite whose transfer fails any gate is re-scheduled on its next candidate or reported failed with a reason. An unpolished catalog entry is never substituted.
- Mass bookkeeping (closed form): `f = 1 - mf/m0` of the transfer; stage propellant `mp = f*(150 + d)/(1 - f)` with `d = dryFrac*mp`, i.e. `mp = 150*f/(1 - f - dryFrac*f)`; stack wet `150 + d + mp`; thrust `5*wet/150` N; jet power `thrust*Isp*g0/2`. `dryFrac = 0.15` knob.
- House MATLAB style: header block with Purpose / Inputs with sizes / Outputs with sizes / Revision History (`M. Casey (c) 11/11/2026`, `Copyright Coorbital Inc.`); `=` aligned knob blocks; `if nargin == 0` self-demo on library functions; never `i`/`j` indices; never `norm`; no `%#ok`; ASCII; no emojis; no signature line in the README; `datetime`, not `datestr`.
- Tests are plain functions returning a logical, registered in `tests/run_tests.m`; existing leg-2 tests stay green. Every gate gets a deliberate-failure experiment.
- MATLAB only through the MCP tools in the shared session; never `matlab -batch`; never `restoredefaultpath`; no `clear all`. Subagents never commit.
- Every saved results file carries `run_meta` (date, git HEAD with dirty flag, MATLAB version, knobs, package identity: file SHA-256, `created`, `n_entries`, sheet orientations used).
- Verified interfaces for `ms_tfmin`, `seed_from_z8`, `fly_transfer`, `phase_state`, `nd_propulsion`, `mintime_hypothesis_gates`, the conjugate re-solve and pumpkyn `tfMin`: `.superpowers/sdd/2026-11-10-lunar-pnt-mission/interfaces.md` sections A3 to A5. Leg-2 patterns to copy: `lib/schedule_rendezvous.m`, `lib/polish_rendezvous.m`, `lib/export_mission_transfers.m`, `lib/mission_movie.m`.

## Review Focus

1. The orientation walk loses the solution (inner solve fails or jumps to another branch part way): the polish must return `ok = false` with the orientation reached and a reason, never a transfer from a different departure state. Pinned in Task 2.
2. A polished transfer passes arrival but dives below the Earth floor between nodes at the first perigee pass: the clearance gate must use the located minimum and fail it. Pinned in Task 4 (deliberate failure with a pilot-demoted trajectory shape).
3. No complete schedule exists on the available sheets within the orientation bound: the scheduler throws `schedule_gto_dro:noAssignment` with the per-satellite candidate counts, never a partial assignment. Pinned in Task 3.
4. The package file is replaced by the fleet's between a leg-1 run and a combined-products build, or the leg-2 results come from a different DRO: combined products refuse by identity mismatch. Pinned in Task 6.
5. A re-flown polished transfer is compared to the slot at the wrong epoch (leg-1 vs leg-2 clock, or the 1/6 relabelling): the hand-off gate recomputes every satellite's DRO phase at T0 from its flown arrival state, independent of the schedule's bookkeeping. Pinned in Task 4.

---

## File map

| File | Responsibility |
|---|---|
| `lib/leg1_inputs.m` | package, clean table, GTO and slot state handles, phys, identity |
| `lib/polish_gto_dro.m` | orientation walk + pinned-departure arrival root-find for one transfer |
| `lib/schedule_gto_dro.m` | pure scheduler on tables |
| `lib/leg1_gates.m` | clearance, one-GTO, hand-off phase, second-order verdicts, mass and power |
| `lib/solve_gto_dro_offgrid.m` | fallback: one off-grid cell through the engine with the fleet recipe |
| `lunar_pnt_leg1_study.m` | numbered study script: schedule, polish, gates, save, capture, figures |
| `lib/export_leg1_transfers.m` | sidecar of the polished transfers |
| `lib/combined_timeline.m` | leg 1 + leg 2 timeline figure and table |
| `lib/leg1_movie.m` | leg-1 movie; `lib/mission_movie.m` gains `opts.dayOffset` |
| `lib/full_mission_movie.m` | joins the two movies |
| `tests/test_*.m` | one per file above |

---

### Task 1: `leg1_inputs`

**Files:** Create `lib/leg1_inputs.m`, `tests/test_leg1_inputs.m`; modify `tests/run_tests.m`, `setup_lunar_pnt_paths.m` (add `GTO_DRO/lib` for `earth_clearance` with a `which` assertion, as the file already does for other folders).

**Interfaces — Produces:** `L = leg1_inputs(opts)` with `opts.catalogMat`; fields: `muStar, lStar, tStar, thrustN (5), ispS, m0kg, phys` (as `mission_inputs().phys`, built with `nd_propulsion`), `tauG, tauD, sgn`, `sheetOrient [1 x nS]`, `sD [1x12]`, `sA [1x6]`, `tf [12 x 6 x nS]` (NaN where not clean), `clean [12 x 6 x nS]` logical, `z8 [8 x 12 x 6 x nS]`, `gtoState = @(orientDeg, phase) [6x1]`, `slotState = @(phase) [6x1]` (DRO), `identity` struct (`sha256, created, n_entries, orientations`), `checks` cell table as in `mission_inputs`.

- [ ] Write `test_leg1_inputs`: asserts regime values (5 N, 1710 s, 150 kg), `abs(tauD - 1) < 1e-9`, `tauG` within 1e-6 of the Kepler period of the default elements, 34 clean cells on the 0 deg sheet of the pilot package (skip this count with a printed note when the identity is not the pilot's), `gtoState(o, s)` equals the catalog's own departure endpoint for a sheet orientation, and the SIGN test: an inertial state on the GTO converted with `pumpkyn.cr3bp.fromPCI` at times 0 and `dt` equals `gtoState(o0 + sgn*(180/pi)*dt/tStar, s)` to 1e-9 ND for `sgn` as stored and fails for `-sgn` (deliberate failure). One clean entry re-flown with `fly_transfer(z8, gtoState(...), slotState(sA), phys)` misses by < 1 km.
- [ ] Run it; expect failure (`leg1_inputs` undefined).
- [ ] Implement `leg1_inputs`. The sign is computed by the same fromPCI comparison at load time and stored.
- [ ] Run the test and `run_tests`; all pass. `check_matlab_code` clean.

### Task 2: probe, then `polish_gto_dro`

**Files:** Create `lib/polish_gto_dro.m`, `tests/test_polish_gto_dro.m`.

**Interfaces — Consumes:** `L` from Task 1. **Produces:** `P = polish_gto_dro(L, cand, opts)` where `cand` has `iD, iA, iS` (sheet), `orientDeg` (true orientation at departure), `tDep`, `slotPhaseAtDep` (slot phase at `tDep`). `P` fields: `ok`, `reason`, `z8 [8x1]`, `tf`, `tDep`, `tArr`, `rv0 [6x1]`, `rvf [6x1]`, `orientDeg`, `sD`, `arrPhase`, `missKm`, `missMps`, `res`, `nEval`, `walk` (orientations visited, tf at each), `orientReached`.

Method: (a) WALK: from the sheet orientation to `cand.orientDeg` at the cell's own arrival phase, in steps of `opts.stepDeg` (default 2), each step `ms_tfmin` seeded by `seed_from_z8` of the previous solution flown from its own departure state; a failed step is halved up to 3 times, then `ok = false`, reason `'orientation walk stalled'`, `orientReached` set. A step whose `tf` changes by more than `opts.jumpFrac` (default 0.05) of the previous is a branch jump: same halving, then reason `'branch jump'`. (b) ROOT: departure fixed; `g(tf) = tf - T(tf)` with target `slotState(mod(cand.slotPhaseAtDep + tf/tauD, 1))`; secant with bracketing and the discontinuity rule of `polish_rendezvous`; stop at `|g| < opts.tolG` (default 1e-10).

- [ ] PROBE first (scratch script in the scratchpad, not committed): take the clean cell `(iD 1, iA 1)` of the 0 deg sheet, walk to +10 deg and to -25 deg, then root-find onto a slot 0.03 of a period past the cell's arrival phase. Record per step: iterations, residual, tf, flown miss, wall seconds. If the walk cannot reach +10 deg with any step size down to 0.25 deg, STOP and report BLOCKED with the table (the plan's fallback then becomes the main route and the controller rules).
- [ ] Write `test_polish_gto_dro`: (1) identity: `cand` at the sheet orientation and the cell's own arrival phase returns the catalog `tf` to 1e-8 and `z8` to 1e-6; (2) +10 deg and a slot offset of 0.03: `ok`, flown miss < 1 km and < 1 cm/s (loosen only to a MEASURED value, reported), `rv0` equals `L.gtoState(cand.orientDeg, sD)` to 1e-12, arrival phase satisfies `mod(slotPhaseAtDep + tf/tauD, 1)`; (3) Review Focus 1: `opts.stepDeg = 60, opts.maxHalvings = 0` on a 60 deg target returns `ok = false`, a reason, `orientReached` at the sheet orientation, and empty `z8`; (4) deliberate failure of the miss gate through `opts.gateOverride`.
- [ ] Run; expect failure. Implement. Run test and `run_tests`; pass. `check_matlab_code` clean.

### Task 3: `schedule_gto_dro`

**Files:** Create `lib/schedule_gto_dro.m`, `tests/test_schedule_gto_dro.m`.

**Interfaces — Produces:** `S = schedule_gto_dro(tf, sD, sA, sheetOrient, tauG, tauD, sgn, opts)`; `tf [nD x nA x nS]` with NaN for unusable cells; `opts`: `orient0Grid` (default `-180:5:175`), `gtoPhase0Grid` (default `(0:11)/12`), `deltaGrid` (default `(0:11)/72`), `nMaxRevs` (default 12), `tolPhase` (default 0.02), `orientBound` (default 25), `propKgPerNd`, `exclude` (`[m x 4]` rows `iD iA iS n` barred, for re-scheduling). Output: `orient0, gtoPhase0, delta, T0, lastArr`, and per satellite `k = 1:6` (leg-2 numbering): `iD, iA, iS, n, tDep, tArr, tf, orientDeg, slotPhaseAtDep, resid`, plus `cands{k}` (all admissible candidates of the chosen global triple, sorted by `tArr`) and `nCand [1x6]`.

Candidate rule for a global triple: `tDep = tauG*(mod(sD(iD) - gtoPhase0, 1) + n)`, `orientDeg = orient0 + sgn*(180/pi)*tDep`, admissible sheet iff `abs(wrap180(orientDeg - sheetOrient(iS))) <= orientBound`, `tArr = tDep + tf`, `resid = wrap(sA(iA) - mod(delta + (j-1)/6 + tArr/tauD, 1))`, admissible iff `abs(resid) <= tolPhase`. Per triple: each slot j takes its earliest admissible candidate (two stages may share a departure instant; the six slots are independent once the triple is fixed), and the triple is infeasible if any slot has none. Objective over triples: minimise T0 (definition in Global Constraints), ties within 1e-12 to less total propellant, then to the smaller `[orient0 gtoPhase0 delta]` lexicographically. Relabel slots to leg-2 numbering at T0.

- [ ] Write `test_schedule_gto_dro` on synthetic tables small enough to work by hand: (1) a 2-phase, 6-arrival table with constant `tf` where the answer (`T0`, each `tDep`) is computed in the test from the formulas; (2) orientation bound: a table whose only admissible candidates need `orient0` in a known window, assert the chosen `orient0` lies in it; (3) relabelling: satellite 1 is at DRO phase 0 at `T0` (recomputed in the test from `tArr`, `sA`, `tauD`); (4) Review Focus 3: all-NaN column for one slot throws `schedule_gto_dro:noAssignment` and the message contains the six candidate counts; (5) `exclude` removes the previous choice and returns the next; (6) tie-break to less propellant.
- [ ] Run; expect failure. Implement as a pure function (no file access). Run test and `run_tests`; pass. `check_matlab_code` clean.

### Task 4: `leg1_gates`

**Files:** Create `lib/leg1_gates.m`, `tests/test_leg1_gates.m`.

**Interfaces — Consumes:** `L`, a polished `P` (Task 2), schedule fields. **Produces:** `G = leg1_gates(L, P, ctx, opts)` with `ctx.orient0, ctx.gtoPhase0, ctx.T0, ctx.kSat`; `G` fields: `arrival` (`missKm, missMps, pass`), `oneGto` (`dStateNd, pass`: `P.rv0` against `L.gtoState` recomputed from `ctx` and `P.tDep` alone), `handoff` (`phaseAtT0, errPhase, pass`: DRO phase of the FLOWN arrival state, found by nearest point on the DRO, advanced by `(T0 - tArr)/tauD`, against `(kSat-1)/6`), `earth` and `moon` (`hMinKm, tAtMin, pass`; located minimum: refine each sample minimum of `earth_clearance` with the root of `(r - r_b).v` on a dense re-flight), `second` (`conjPass, conjCovered, minLamV, minQmt, dimS, pass`, instruments as in interfaces A4), `mass` (`mfFrac, propKgPer150, stagePropKg, stageDryKg, wetKg, thrustN, jetKW, dvKms`), `pass` (all of arrival, oneGto, handoff, earth, moon), `refuted` (second-order fail; reported, not a blocker).

- [ ] Write `test_leg1_gates`: on the Task 2 polished transfer all gates pass and the mass block satisfies the closed form to 1e-12 (hand-worked numbers for `f = 0.42`, `dryFrac = 0.15`); deliberate failures, one per gate: arrival (perturb `z8`), oneGto (`ctx.orient0 + 1`), handoff (`ctx.T0 + 0.01`), Earth (Review Focus 2: a synthetic arc whose samples straddle a perigee at 150 km with both neighbours above 200 km must fail on the located minimum while the sample minimum passes), Moon (floor raised above the measured minimum); Review Focus 5: handoff computed with `kSat` shifted by one fails.
- [ ] Run; expect failure. Implement. Run test and `run_tests`; pass. `check_matlab_code` clean.

### Task 5: `solve_gto_dro_offgrid` (fallback)

**Files:** Create `lib/solve_gto_dro_offgrid.m`, `tests/test_solve_gto_dro_offgrid.m`.

**Interfaces — Produces:** `E = solve_gto_dro_offgrid(L, orientDeg, sD, sA, opts)`: one cell at an off-grid orientation through `thrust_ladder_library` with the fleet recipe read from `GTO_DRO/catalog/run_gto_dro_catalog.m` (`thrLock`, `tf0 = 0.30`, `directGateKm = inf`, `gateKm = 100`, `msHardCapSec = 300` with `requireFence`, floors Moon 500 km and Earth 200 km, N = 400, rungs `[15 12 10 7 5]`, `keepDirect = false`), sheet written to a temp folder under `lunar_pnt_mission/results/offgrid/`. `E`: `ok`, `reason`, `z8`, `tf`, `rungsOK [1x5]`, `wallS`, `sheetMat`. Read the option names from the front door and the engine header; do not guess them.

- [ ] Write `test_solve_gto_dro_offgrid`: orientation 12 deg, `sD = 0.5`, `sA = 0`: `ok`, the 5 N entry re-flies to < 1 km, `tf` within 15 percent of the 0 deg sheet's cell (a plausibility bound, loosened only to a measured value); refusal when no pool can be opened and `requireFence` (inject with `opts.forceNoPool`). Expected wall time about 5 to 15 min; run it as a background MCP call.
- [ ] Run; expect failure. Implement. Run test and `run_tests`; pass. `check_matlab_code` clean.

### Task 6: study script, capture, combined timeline

**Files:** Create `lunar_pnt_leg1_study.m`, `lib/export_leg1_transfers.m`, `lib/combined_timeline.m`, `tests/test_lunar_pnt_leg1_study.m`, `tests/test_export_leg1_transfers.m`, `tests/test_combined_timeline.m`.

**Interfaces — Consumes:** Tasks 1 to 5. **Produces:** `results/lunar_pnt_leg1_study.mat` with `S` (schedule), `P1 [1x6]` (polish output plus `G` gates, flown history `t [1xM]` in leg-1 ND, `X [7xM]`, `alpha [3xM]`, `route` in `{'catalog','offgrid'}`, `attempts`), `T0`, `knobs`, `run_meta`, `failed` (indices and reasons); `results/leg1_transfers_gto_dro_5N.mat` (sidecar); figures `leg1_schedule.png`, `leg1_3d.png`, `leg1_mass.png`, `mission_timeline.png`.

Study sections, in the style of `lunar_pnt_mission_study.m` (knob block first; each gate computed in the script and printed PASS/FAIL): 0 knobs; 1 inputs and identity; 2 schedule; 3 polish each satellite, with the loop: on polish failure or a failed blocking gate, add the candidate to `exclude`, re-schedule that satellite only among `S.cands{k}` (same global triple), up to `maxAttempts` (default 6), then try `solve_gto_dro_offgrid` once, then mark failed; 4 gates table; 5 mass and power table followed by the caveat sentence from the spec; 6 save; 7 capture; 8 figures. `export_leg1_transfers` writes one struct per accepted transfer: `z8, tf_nd, rv0, rvf, thrustN, ispS, m0kg, orientDeg, sD, arrPhase, tDep, tArr, seededFrom (iD iA iS)`, the verdicts with what was and was not checked, and the package identity; it never writes outside `results/`. `combined_timeline(leg1Mat, leg2Mat)` checks compatibility first (Review Focus 4: leg-2 DRO params and `droPhase0` against leg 1, and the leg-1 package identity against the current file unless `opts.allowStalePackage`), then draws per satellite: GTO wait, leg 1, DRO coast to T0, DRO wait, leg 2, on station, with leg-2 times shifted by T0, and returns the table.

- [ ] Write the three tests. Study test: runs the saved results (not the solve) and asserts structure, six entries or a non-empty `failed` with reasons, every accepted `G.pass`, hand-off phases, `T0 >= max(tArr)`; headline numbers are pinned AFTER the first measured run (step below). Export test: round trip, each stored transfer re-flies to its stored miss, refusal to write outside `results/`. Timeline test: synthetic leg-1 and leg-2 structs with hand-worked bar ends; identity mismatch refuses.
- [ ] Run; expect failure. Implement the three files.
- [ ] Run the study (background MCP call; expected tens of minutes). Report every printed number. Pin the headline numbers (T0 in days, the six flight times, total propellant, minimum Earth and Moon altitudes, achieved arrival maxima) in the study test.
- [ ] If fewer than six satellites close: do not loosen anything; report which and why, with the candidate counts, and finish the task with DONE_WITH_CONCERNS.
- [ ] `run_tests` passes; `check_matlab_code` clean.

### Task 7: movies

**Files:** Create `lib/leg1_movie.m`, `lib/full_mission_movie.m`, `tests/test_leg1_movie.m`; modify `lib/mission_movie.m` (new `opts.dayOffset`, default 0, added to every displayed day; default output byte-for-byte unchanged in its frames, checked by the existing `test_mission_movie`).

**Interfaces — Produces:** `info = leg1_movie(leg1Mat, opts)` writing `results/lunar_pnt_leg1.mp4` (1280 x 720, frame size a multiple of 16); `info = full_mission_movie(leg1Mat, leg2Mat, opts)` writing `results/lunar_pnt_full_mission.mp4`.

Leg-1 movie, following `mission_movie` (reuse its helpers by copying the local functions it needs; do not restructure that file): dark frame; Earth-Moon rotating frame; textured Moon as in `mission_movie`, Earth drawn as a sphere; DRO outline always on; the GTO outline always on and redrawn each frame at `orient(t)`; each stage a marker riding the GTO until its departure, then its flown path with the yellow-orange plume along the exhaust direction while thrusting, then on the DRO; day counter and days-per-second readout; camera with fixed view direction, zoom from the data so every satellite stays in view (wide while any stage is near Earth, closing on the Moon as slots fill), first frame wide on Earth and GTO, last frame the END framing of `mission_movie`'s first frame; bottom strip: satellites in the DRO versus mission day. `full_mission_movie` renders leg 2 with `dayOffset = T0` in days and joins the two (one writer, or `ffmpeg` concat if present; frame size and rate must match).

- [ ] Write `test_leg1_movie`: a 2 s preview render (`opts.preview`) exists, is 1280 x 720, has the expected frame count; the time map is monotone and reaches the last arrival; at 20 sampled frames every satellite's projected position is inside the scene box (the in-view rule, computed from `info`); `dayOffset = 3` shifts the first displayed day of `mission_movie` by 3 (from `info`), default unchanged.
- [ ] Run; expect failure. Implement. Render both movies. Extract six still frames (departures, mid-transfer, first arrival, last arrival, hand-off, a leg-2 frame) into `results/movie_preview/` and LOOK at them: no clipped satellite, outlines present, plume behind the motion, counters legible.
- [ ] `run_tests` passes; `check_matlab_code` clean.

### Task 8: README and close-out

**Files:** Modify `orbit_transfer/lunar_pnt_mission/README.md`, `tests/run_tests.m` (final table), `orbit_transfer/README.md` row if the folder has one for this study.

- [ ] README section "Leg 1: GTO to DRO" with the measured numbers (schedule, flight times, gates achieved, clearances, second-order verdicts, mass and power table), the caveat in the spec's words, what the gates do and do not cover, a deliberate-failure table for the leg-1 gates, how to run (`lunar_pnt_leg1_study`, `combined_timeline`, `leg1_movie`, `full_mission_movie`), the Jacobi constants of the DRO (3.0576) and tulip (3.1866), and the statement of which package identity the numbers come from.
- [ ] `run_tests` full pass, output pasted in the report.
