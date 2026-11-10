# Lunar PNT Mission Simulation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One study script that schedules six spacecraft from a DRO onto the six optimally phased tulip slots by certified 70 mN minimum-time transfers (earliest complete constellation), polishes each transfer onto its exact slot state, and scores the assembled constellation for south-pole GDOP and the paper's R1 to R3 requirements over 30 days.

**Architecture:** A new campaign folder `orbit_transfer/lunar_pnt_mission/` with four library functions behind a numbered study script. The transfer side reads the 70 mN library of record and uses the `costate_common` instruments unchanged; the geometry side calls proj7's `gdop_grid_from_states` once per interval over which the on-station set is constant and applies `run_requirements_table`'s recipe on proj7's scoring utilities. Nothing outside the new folder is modified.

**Tech Stack:** MATLAB R2026a through the MCP tools in the shared session; pumpkyn / pumpkynPie (`pumpkyn.cr3bp.tfMin`, `tfMinProp`); `orbit_transfer/costate_common` (`get_family_orbit`, `phase_state`, `fly_transfer`, `seed_from_z8`, `ms_tfmin`, `pmp_pointwise_checks`, `mintime_hypothesis_gates`); the deliverable-9 kit (`transfers_at`, `catalog_transfer`, the catalog `.mat`); proj7 `pipelines/gdop` (`tulip_phasing`, `gdop_grid_from_states`, `max_min_access_bounds`, `max_min_access_block_bounds`); CasADi 3.7 (needed by `mintime_hypothesis_gates`).

**Spec:** `docs/superpowers/specs/2026-11-10-lunar-pnt-mission-design.md`. Verified interfaces of every existing function named here (signatures, shapes, units, traps, quoted with line numbers): `.superpowers/mission_sim_interfaces.md` — every implementer reads the sections for the functions their task calls before writing code.

## Global Constraints

- Regime, read from the catalog and never retyped: thrust 0.07 N, Isp 900 s, m0 150 kg; `muStar = 0.012150585609624`, `lStar = 389703.264829278 km`, `tStar = 382981.289129055 s`; DRO period 1 ND (4.432654 d); tulip Np 7, pm -1, period `(5/6)*2*pi` ND (23.209321 d).
- Catalog: `orbit_transfer/DRO_tulip/deliverables/costate_library_70mN_user_deliverable9/costate_catalog_dro_tulip_70mN.mat`, variable `costate_catalog_dro_tulip_70mN`, single sheet; `sD_frac = (0:23)/24`, `sA_frac = 0.0129 + (0:47)/48`; `tf_nd [24 x 48]`, `z8 [8 x 1152]` indexed by `entry_index(iD, iA)`. Primaries only (spec decision 5).
- Phasing only from `tulip_phasing()`: `dtau6 = [0.778317343114462 1.14555587028854 2.51741969625761 2.89299534327950 4.26677262887274 4.63164623641942]`. Slot j at clock time t has tulip phase `mod((dtau6(j) + t)/tauT, 1)`. Spacecraft k at clock time t has DRO phase `mod((k-1)/6 + t/tauD, 1)`.
- One orbit definition: `get_family_orbit` evaluated through the periodic `phase_state` (columns `[6 x 1]`), for DRO states, slot states and transfer endpoints alike.
- Schedule: tolerance 0.005 of a tulip period, `nMaxRevs = 10`, objective = minimise the last arrival, tie-break total propellant.
- Acceptance thresholds for a polished transfer (from `DRO_tulip/indirect/transfer_study.m`): `ms_tfmin` residual < 1e-8; flown miss < 100 km and < 10 m/s; `PW.Hmax < 1e-6`, `PW.lamMf < 1e-6`, `PW.adjErr < 1e-7`, `max(PW.fullGap, PW.fieldGap) <= 1e-12`, `PW.throttleErr < 1e-10`, `PW.minQmt >= 0`; gates `minLamV > 0`, `minQmt > 0`, `dimS == 1`; conjugate `info.conj.pass == 1 && info.conj.covered`.
- GDOP: elevation mask 0 deg. Paper grid = five chunks with `lat_range` pairs `[-89.9 -87] [-86.5 -84] [-83.5 -81] [-80.5 -78] [-77.5 -75]`, `lat_res 0.5`, `lon_range [0 358]`, `lon_res 1` (10,770 points), cadence 348.1 s (`time_samples = round(span_s/348.1)` per interval). Coarse grid = `lat_range [-89.5 -75.5]`, `lat_res 2`, plus a pole chunk `[-89.9 -89.9]`, `lon_range [0 358]`, `lon_res 2`, cadence 300 s.
- Requirement recipe (verbatim from `run_requirements_table`): `rak = struct('dt', dtSec, 'bin_duration', 86400)`; R1 = `max_min_access_bounds(visibility_count >= 2, rak)`, pass >= 0.90; R2 = `max_min_access_bounds(gdop < 6, rak)`, pass >= 0.40; R3 = `max_min_access_block_bounds(tau_sec, gdop < 6, struct('min_block_duration', 4*3600, 'total_duration', 9*3600))/3600`, pass >= 9 h; floors are the minimum over grid points. R3 on `gdop < 7` is an extra row.
- Validation reference: `proj7/figures/requirements_table_results.mat`, `S2(1).r`: R1 floor 1.000000, R2 floor 0.754032, R3 floor 17.987224 h. Never run `run_requirements_table` (it overwrites that file).
- Nothing under `proj7`, `pumpkyn`, `pumpkynPie`, `costate_common`, `oclib` or the kit folder is modified. `addpath`, never `genpath`.
- House MATLAB style: pumpkyn header quartet with sizes, `=` aligned knob blocks, `if nargin == 0` self-demo on library functions, never `i`/`j` as indices, never `norm` (use `vecnorm` or `sqrt(sum(...))`), no `%#ok`, no emojis, ASCII README, no author signature line in the README. `M. Casey (c) 11/10/2026`, `Copyright Coorbital Inc.` in headers.
- Results go to `orbit_transfer/lunar_pnt_mission/results/` (git-ignored via a folder `.gitignore`), each saved file carrying `run_meta` (date, git HEAD, MATLAB version, knobs).
- MATLAB only through the MCP tools in the shared session; never `matlab -batch`; never `restoredefaultpath`; no `clear all`. The GTO -> DRO fleet runs in a separate background MATLAB at the same time: do not kill MATLAB processes.
- Tests are plain functions returning a logical; `tests/run_tests.m` is a function that prints a table and throws on failure (pattern: `GTO_DRO/tests/run_tests.m`).
- Subagents never commit; the controller commits after each clean review, with `/Library/Developer/CommandLineTools/usr/bin/git`.

## Review Focus

1. No complete assignment exists (some spacecraft cannot reach any free slot within `nMaxRevs` at the tolerance): the scheduler must refuse with `schedule_rendezvous:noAssignment` and still return the 6 x 6 table in the error path's log, never return a partial assignment as if complete. Pinned in Task 2.
2. A polish whose root needs a departure before clock time 0, or which fails a gate: the study must take the scheduler's next-best candidate for that (spacecraft, slot) pair and say so, never substitute the unpolished cell. Pinned in Task 3 (`P.ok = false` with a reason) and Task 6 (fallback loop).
3. The on-station set changes between intervals: stitching per-interval GDOP calls must neither duplicate nor drop an epoch, and an interval with fewer than four satellites must yield NaN GDOP and a correct visibility count, not an error. Pinned in Task 5.
4. A function resolving from the wrong folder (the session path already carries `DRO_tulip/indirect` copies and a different `costate_catalog_pick`): the path setup's `which` assertions must fail loudly. Pinned in Task 1.
5. The proj7 reference file is missing or its phasing no longer matches `tulip_phasing()`: the validation gate must report UNRESOLVED and the study must stop before claiming mission results, not pass silently. Pinned in Task 5.

---

## File map

```
orbit_transfer/lunar_pnt_mission/
  setup_lunar_pnt_paths.m          Task 1
  lib/mission_inputs.m             Task 1   catalog + constants + orbit handles + phasing, verified
  lib/schedule_rendezvous.m        Task 2
  lib/polish_rendezvous.m          Task 3
  lib/assemble_histories.m         Task 4
  lib/mission_grid.m               Task 5   'paper' | 'coarse' grid definitions
  lib/mission_gdop.m               Task 5
  lunar_pnt_mission_study.m        Task 6
  tests/run_tests.m                Task 1 (created), extended by every task
  tests/test_mission_inputs.m      Task 1
  tests/test_schedule_rendezvous.m Task 2
  tests/test_polish_rendezvous.m   Task 3
  tests/test_assemble_histories.m  Task 4
  tests/test_mission_gdop.m        Task 5
  tests/test_lunar_pnt_mission_study.m  Task 6
  README.md                        Task 6 (skeleton), Task 7 (measured results)
  .gitignore                       Task 1   "results/"
orbit_transfer/README.md           Task 7   one row for the new folder
```

---

### Task 1: Paths and verified inputs

**Files:**
- Create: `orbit_transfer/lunar_pnt_mission/setup_lunar_pnt_paths.m`, `lib/mission_inputs.m`, `tests/run_tests.m`, `tests/test_mission_inputs.m`, `.gitignore`

**Interfaces:**
- Consumes: `get_family_orbit(family, p)`, `phase_state(t, y)`, `tulip_phasing()`, the catalog file (see the interfaces file, sections A1, A2, A3, B6).
- Produces:
  - `setup_lunar_pnt_paths()` (no arguments, no outputs): adds this folder, `lib`, the kit folder, `costate_common`, `oclib`, `$HOME/casadi-3.7.0`, `proj7/pipelines/gdop`; runs pumpkynPie's `startup` only if `pumpkyn.cr3bp.tfMin` is not found; then asserts with `which` that `transfers_at` and `catalog_transfer` resolve inside the kit folder, `fly_transfer`, `ms_tfmin`, `get_family_orbit`, `phase_state` inside `costate_common`, and `gdop_grid_from_states`, `tulip_phasing`, `max_min_access_bounds`, `max_min_access_block_bounds` inside `proj7/pipelines/gdop`; error id `setup_lunar_pnt_paths:shadow` naming the function and the folder it resolved to.
  - `M = mission_inputs()` returning a struct with: `.cat` (catalog struct), `.sh` (its sheet), `.muStar .lStar .tStar` (lStar in km, tStar in s), `.thrustN .ispS .m0kg .Tnd .cnd`, `.phys` (struct with `Tnd cnd muStar lStar tStar m0kg` for `fly_transfer`), `.tauD .tauT` (ND), `.sD [1 x 24] .sA [1 x 48] .tf [24 x 48]`, `.depState` and `.arrState` (function handles, phase fraction -> `[6 x 1]`), `.dtau6 [1 x 6]`, `.droPhase0 [1 x 6] = (0:5)/6`, `.slotPhase0 [1 x 6] = dtau6/tauT`, `.propKgPerNd` (scalar: propellant in kg per ND of burn, `m0kg*Tnd/cnd`), `.checks` (struct array `.name .value .pass` of the verifications below). `mission_inputs` throws `mission_inputs:check` if any check fails.

- [ ] **Step 1: Write the failing test** `tests/test_mission_inputs.m`

```matlab
function ok = test_mission_inputs()
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));
M = mission_inputs();
assert(isequal(size(M.tf), [24 48]) && all(isfinite(M.tf(:))), 'tf table');
assert(abs(M.thrustN - 0.07) < 1e-15 && M.ispS == 900 && M.m0kg == 150, 'regime');
assert(abs(M.tauD - 1) < 1e-12 && abs(M.tauT - (5/6)*2*pi) < 1e-12, 'periods');
assert(abs(M.sA(1) - 0.0129) < 1e-9 && abs(M.sD(2) - 1/24) < 1e-12, 'catalog grids');
assert(max(abs(M.slotPhase0 - [0.1486 0.2188 0.4808 0.5525 0.8149 0.8846])) < 1e-4, 'slot phases');
assert(isequal(size(M.depState(0.3)), [6 1]) && isequal(size(M.arrState(0.3)), [6 1]), 'state handles');
% the catalog's own endpoints: flying cell (3,5) from depState(sD(3)) lands on arrState(sA(5))
z8 = M.sh.z8(:, M.sh.entry_index(3, 5));
F  = fly_transfer(z8, M.depState(M.sD(3)), M.arrState(M.sA(5)), M.phys);
assert(F.flyKm < 100 && F.flyVms < 10, 'cell (3,5) flown miss %.3g km', F.flyKm);
assert(abs(M.propKgPerNd*M.tf(3,5) - F.propellantKg) < 1e-6, 'propellant rate');
assert(all([M.checks.pass]), 'input checks');
ok = true;  fprintf('test_mission_inputs: PASS\n');
end
```

- [ ] **Step 2: Run it** with `run_matlab_file`. Expected: fails, `setup_lunar_pnt_paths.m` not found.

- [ ] **Step 3: Implement `setup_lunar_pnt_paths` and `mission_inputs`.** Checks `mission_inputs` records and enforces: the catalog has one sheet and all 1,152 `has_solution` and `status == 4`; `sh.tau_dep == tauD` and `sh.tau_arr == tauT` to 1e-12 against `get_family_orbit`'s periods; `tulip_phasing`'s `tau0 == tauT`, `Np == sh.Np`, `pm == sh.pm`; `Tnd` and `cnd` computed from the catalog's constants equal `cat_.thruster.c_nd` and the README formula to 1e-12 relative.

- [ ] **Step 4: Write `tests/run_tests.m`** in the `GTO_DRO/tests/run_tests.m` pattern (names list starts with `'test_mission_inputs'`; error id `lunar_pnt_mission:run_tests:fail`) and `.gitignore` containing `results/`.

- [ ] **Step 5: Run the test and `run_tests`.** Expected: `test_mission_inputs: PASS`, 1 of 1.

- [ ] **Step 6: Shadow experiment** (the deliberate failure for Review Focus 4): in the shared session, temporarily `addpath` `orbit_transfer/DRO_tulip/indirect` at the FRONT, call `setup_lunar_pnt_paths`, record whether it throws `setup_lunar_pnt_paths:shadow` or re-fronts the kit; `rmpath` afterwards. The function must end with the kit copy winning or a throw; report which and make the test assert it.

- [ ] **Step 7: Commit** (controller)

```bash
git add orbit_transfer/lunar_pnt_mission
git commit -m "lunar_pnt_mission: path setup with which-assertions and mission_inputs (70 mN catalog, orbit handles, tulip phasing, verified)"
```

---

### Task 2: `schedule_rendezvous`

**Files:**
- Create: `lib/schedule_rendezvous.m`, `tests/test_schedule_rendezvous.m`
- Modify: `tests/run_tests.m`

**Interfaces:**
- Consumes: nothing from the catalog directly; it is a pure function of tables so it can be tested by hand.
- Produces: `S = schedule_rendezvous(tf, sD, sA, tauD, tauT, droPhase0, slotPhase0, opts)` with `tf [nD x nA]` ND, `sD [1 x nD]`, `sA [1 x nA]` phase fractions, `droPhase0 [1 x K]`, `slotPhase0 [1 x K]`, `opts.tol` [0.005], `opts.nMaxRevs` [10], `opts.propKgPerNd` [1]. Candidate for spacecraft k, slot j, cell (iD, iA), wait n in `0:nMaxRevs`: `tDep = tauD*mod(sD(iD) - droPhase0(k), 1) + n*tauD`; `tArr = tDep + tf(iD, iA)`; `resid = wrap(sA(iA) - mod(slotPhase0(j) + tArr/tauT, 1))` into (-0.5, 0.5]; admissible iff `abs(resid) <= tol`.
  Output fields: `.assign [1 x K]` (slot of each spacecraft), `.iD .iA .nRev .tDep .tArr .resid .tf .propKg` (each `[1 x K]`, the chosen candidate per spacecraft), `.tLast` (scalar), `.totalPropKg`, `.best` (struct with `[K x K]` matrices `tArr iD iA nRev resid`, NaN where no admissible candidate), `.ranked` (`{K x K}` cell, each an `[m x 5]` matrix `[tArr iD iA nRev resid]` of ALL admissible candidates sorted by tArr, for the polish fallback), `.nAdmissible [K x K]`, `.opts`. Assignment: minimise `max(tArr)` over permutations using each pair's earliest candidate; ties (within 1e-12) broken by smaller total propellant, then by lexicographic permutation. Error `schedule_rendezvous:noAssignment` (message lists the pairs with no candidate) when no permutation is complete.

- [ ] **Step 1: Write the failing test** `tests/test_schedule_rendezvous.m` with four checks:
  1. **Analytic case.** `tauD = 1, tauT = 4, sD = [0 0.5], sA = [0 0.25 0.5 0.75], tf = [2 2 2 2; 3 3 3 3], droPhase0 = [0 0.5], slotPhase0 = [0 0.5], tol = 0.01, nMaxRevs = 3`. Expected: spacecraft 1 best arrival 2 for both slots (cell iD 1, n 0; iA 3 for slot 1, iA 1 for slot 2); spacecraft 2 best arrival 3 for both slots (iD 2, n 0; iA 4 for slot 1, iA 2 for slot 2); `S.tLast == 3`; every `abs(S.resid) < 1e-12`.
  2. **Brute force.** A seeded random table (`rng(7)`, `nD = 6, nA = 12, K = 3`, `tf = 2 + rand(6, 12)`, `tauD = 1`, `tauT = 5.2`, `droPhase0 = (0:2)/3`, `slotPhase0 = [0.11 0.47 0.82]`, `tol = 0.02`, `nMaxRevs = 6`): the test enumerates every candidate and every permutation itself with plain loops and asserts the function's `tLast`, `assign` and `best.tArr` equal the brute-force values to 1e-12.
  3. **Tolerance respected.** In check 2's result every chosen `abs(resid) <= 0.02`, and `S.nAdmissible` equals the brute-force count per pair.
  4. **Refusal.** Check 1's inputs with `tol = 1e-6` and `sA = [0.1 0.35 0.6 0.85]` (no arrival phase can match): expects `schedule_rendezvous:noAssignment`.

- [ ] **Step 2: Run it.** Expected: `Unrecognized function 'schedule_rendezvous'`.

- [ ] **Step 3: Implement.** Vectorise over (iD, iA, n) per pair; permutations by `perms(1:K)` (720 rows for K = 6). Header states the clock and phase conventions exactly as in Global Constraints. Self-demo (`nargin == 0`): the real catalog through `mission_inputs`, printing the 6 x 6 best-arrival table in days.

- [ ] **Step 4: Run the test, then the self-demo.** Expected: PASS; the demo prints a complete assignment (record `tLast` in days, the assignment, the admissible-candidate counts per pair, and the wall time in the report).

- [ ] **Step 5: Add to `run_tests`, run it** (2 of 2).

- [ ] **Step 6: Commit** (controller): `lunar_pnt_mission: schedule_rendezvous (catalog-cell candidates with whole-revolution waits, phase tolerance, bottleneck assignment; tested against brute force)`.

---

### Task 3: `polish_rendezvous`

**Files:**
- Create: `lib/polish_rendezvous.m`, `tests/test_polish_rendezvous.m`
- Modify: `tests/run_tests.m`

**Interfaces:**
- Consumes: `M` from Task 1; a candidate row from Task 2's `S.ranked{k, j}`; `pumpkyn.cr3bp.tfMin(rv0, rvf, z8, Tnd, cnd, muStar)` (prints an iteration table: call it inside `evalc`; returns `[8 x 1]`, no exit flag); `fly_transfer`, `seed_from_z8(z8, rv0, 24, Tnd, cnd, muStar)`, `ms_tfmin(rv0, rvf, seed, Tnd, cnd, muStar, struct('tolR', 3e-11, 'wallSec', 600, 'conjTest', true))`, `pmp_pointwise_checks(t, Y, Tnd, cnd, muStar)`, `mintime_hypothesis_gates(z8, rv0, Tnd, cnd, muStar, struct('relTol', 1e-12))`.
- Produces: `P = polish_rendezvous(M, k, j, cand, opts)` with `cand = struct('iD', 'iA', 'nRev', 'tDep', 'tArr')`, `opts.tolG` [1e-9 ND], `opts.maxEval` [25], `opts.familyFrac` [0.02], `opts.gates` [true] (false skips the hypothesis gates and the conjugate solve, for fast tests).
  The arrival epoch is pinned: `tArr = cand.tArr`, `rvf = M.arrState(mod(M.slotPhase0(j) + tArr/M.tauT, 1))`. For a trial `tf`: `tDep = tArr - tf`, `rv0 = M.depState(mod(M.droPhase0(k) + tDep/M.tauD, 1))`, `T(tf)` = the eighth component of `tfMin` solved from `rv0` to `rvf` warm-started from the previous solution (first from the cell's z8), and `g(tf) = tf - T(tf)`. Root of g by secant from `tf0 = M.tf(iD, iA)` and `T(tf0)`, falling back to bisection once a sign change is bracketed; stop at `abs(g) < tolG`.
  Output fields: `.ok` (logical) `.reason` (char; one of `'accepted'`, `'departure before epoch 0'`, `'root not found'`, `'inner solve not converged'`, `'gate failed: <names>'`), `.z8 [8 x 1] .tf .tDep .tArr .rv0 [6 x 1] .rvf [6 x 1]`, `.nEval .gResid`, `.tfCell .familyChanged` (`abs(tf - tfCell)/tfCell > familyFrac`), `.flight` (the `fly_transfer` struct of the final z8), `.gates` (struct: `.resid .flyKm .flyVms .passFly .pw .passPW .hyp .passHyp .conjPass .conjCovered .conjReason .passConj`; NaN / `'skipped'` when `opts.gates` is false), `.propKg .dvKms .wallS`.
  Rules: an inner `tfMin` result counts as converged only if `fly_transfer` on it gives `flyKm < 100 && flyVms < 10`; otherwise the evaluation is retried once through `ms_tfmin` seeded by `seed_from_z8` of the previous z8, and if that fails `.reason = 'inner solve not converged'`. `tDep < 0` at the root gives `.ok = false`, `.reason = 'departure before epoch 0'` (no gates run). The final certification is one `ms_tfmin(..., conjTest)` from `seed_from_z8(P.z8, ...)`, counted only when it reproduces the root (`info.converged && max(abs(z - P.z8)) < 1e-6`).

- [ ] **Step 1: Write the failing test** `tests/test_polish_rendezvous.m`:
  1. Take the real schedule (`S = schedule_rendezvous(M.tf, M.sD, M.sA, M.tauD, M.tauT, M.droPhase0, M.slotPhase0, struct('propKgPerNd', M.propKgPerNd))`), spacecraft `k = 1`, its slot `j = S.assign(1)`, its chosen candidate. `P = polish_rendezvous(M, 1, j, cand, struct('gates', false))`. Assert: `P.ok`; `abs(P.gResid) < 1e-9`; `P.flight.flyKm < 100 && P.flight.flyVms < 10`; `abs((P.tDep + P.tf) - cand.tArr) < 1e-12` (pinned epoch); the departure state equals `M.depState` at spacecraft 1's DRO phase at `P.tDep` to 1e-12; `P.rvf` equals the slot state at `tArr` to 1e-12; `abs(P.tf - P.tfCell)/P.tfCell < 0.05` (a cell-sized gap is a small perturbation; if the measured value is larger, report it rather than loosening); `P.nEval <= 25`.
  2. The same call with `gates = true`: every gate field finite, `P.gates.passFly`, and the verdicts printed (assert `passPW` and `passHyp`; print `conjPass`, assert only that it is 0, 1 or NaN with a reason).
  3. Departure-before-epoch refusal: a hand-made candidate with `tArr = 0.5*M.tf(1,1)` (arrival earlier than any transfer could allow) returns `P.ok == false` and `P.reason` starting `'departure before epoch 0'`.
  4. Deliberate gate failure: `polish_rendezvous` exposes its gate evaluation as a local function reachable through `opts.gateOverride = struct('flyKmMax', 1e-9)`; with it the result is `ok == false`, reason starting `'gate failed'` (proves the gate can fail).

- [ ] **Step 2: Run it.** Expected: `Unrecognized function 'polish_rendezvous'`.

- [ ] **Step 3: Implement.** Self-demo: spacecraft 1 of the real schedule with gates on, printing the table row.

- [ ] **Step 4: Run the test.** Record in the report: `nEval`, `tf` cell vs polished (days), phase gap closed (hours), flown miss, every gate value, conjugate verdict, wall time per polish (this sizes Task 6: six polishes plus fallbacks).

- [ ] **Step 5: Add to `run_tests`, run it** (3 of 3).

- [ ] **Step 6: Commit** (controller): `lunar_pnt_mission: polish_rendezvous (pinned arrival epoch, root-find on t_f with the departure state moving, tfMin inner solves judged by the flown miss, acceptance gates and conjugate test on the polished extremal)`.

---

### Task 4: `assemble_histories`

**Files:**
- Create: `lib/assemble_histories.m`, `tests/test_assemble_histories.m`
- Modify: `tests/run_tests.m`

**Interfaces:**
- Consumes: `M`; a `[1 x 6]` cell or struct array of accepted `P` (fields `.tDep .tArr .flight.t .flight.Y`, the flight time axis starting at 0 at departure), the slot of each spacecraft `assign [1 x 6]`.
- Produces: `H = assemble_histories(M, P, assign, opts)` with `opts.dtSec` [300], `opts.serviceDays` [30], `opts.tEnd` (ND; default `max(tArr) + serviceDays*86400/tStar`). Fields: `.t [N x 1]` ND uniform from 0, `.r [N x 3 x 6]` ND barycentric rotating, `.phase [N x 6]` int8 (1 DRO coast, 2 transfer, 3 on station), `.tDep .tArr [1 x 6]`, `.tFirst .tLast`, `.assign`, `.jumpKm` (struct `.dep [1 x 6]`, `.arr [1 x 6]`: position difference between the DRO state and the flight's first sample at `tDep`, and between the flight's last sample and the slot state at `tArr`), `.dtSec`. Also `H = assemble_histories(M, [], 1:6, struct('onSlots', true, 'tEnd', tEnd))`: all six satellites on their slots from t = 0 (no transfers; `.phase` all 3, `.tDep = .tArr = zeros`), used by the validation gate.
  Sampling rule: DRO coast `M.depState(mod(droPhase0(k) + t/tauD, 1))`; transfer by `interp1(tDep + flight.t, flight.Y(:, 1:3), t, 'pchip')` for `tDep <= t < tArr`; on station `M.arrState(mod(slotPhase0(j) + t/tauT, 1))` for `t >= tArr`.

- [ ] **Step 1: Write the failing test** `tests/test_assemble_histories.m`:
  1. On-slots mode over one tulip period: `size(H.r) == [N 3 6]`, all `H.phase == 3`, satellite j at `t = 0` equals `M.arrState(M.slotPhase0(j))` positions to 1e-12, and at `t = H.t(end)` every satellite's distance from the Moon (`[1 - muStar, 0, 0]`) is between 1,900 km and 120,000 km.
  2. One real transfer (spacecraft 1 polished with `gates = false`, the other five given the SAME `P` shifted is not physical, so build the six-element input with spacecraft 1 polished and spacecraft 2 to 6 as `onSlots`-style entries via `opts.onStation = [false true true true true true]`): `H.phase(:, 1)` goes 1 -> 2 -> 3 monotonically with the switches at `tDep` and `tArr` (to one sample); `H.jumpKm.dep(1) < 1` and `H.jumpKm.arr(1) < 100`; the largest sample-to-sample position step of spacecraft 1 is below `3 km/s * dtSec` (no teleport at a hand-over).
  3. Deliberate failure: feeding a `P` whose `tArr` is shifted by +0.05 ND makes `H.jumpKm.arr(1) > 100` (the continuity number can fail).

- [ ] **Step 2: Run it.** Expected: `Unrecognized function 'assemble_histories'`.

- [ ] **Step 3: Implement** (`opts.onStation [1 x 6]` logical marks satellites that are on their slot from t = 0 and need no `P`).

- [ ] **Step 4: Run the test; add to `run_tests`; run it** (4 of 4).

- [ ] **Step 5: Commit** (controller): `lunar_pnt_mission: assemble_histories (DRO coast, flown transfer, tulip slot on one clock; hand-over continuity measured)`.

---

### Task 5: `mission_grid` and `mission_gdop`

**Files:**
- Create: `lib/mission_grid.m`, `lib/mission_gdop.m`, `tests/test_mission_gdop.m`
- Modify: `tests/run_tests.m`

**Interfaces:**
- Consumes: `H` from Task 4; `gdop_grid_from_states(tau, rTgt, mu, lStar, tStar, params)` (half-open uniform resampling over `[tau(1), tau(end))` with `time_samples` points; outputs `.gdop_time_series [Nt x Npts]`, `.visibility_count [Nt x Npts]`, `.lat_vec .lon_vec` rad, `.tau_i`); `max_min_access_bounds(index, rak)`; `max_min_access_block_bounds(tauSec, index, req)`; the reference file `proj7/figures/requirements_table_results.mat`.
- Produces:
  - `g = mission_grid(name)`, `name` in `{'paper', 'coarse'}`: `.name`, `.chunks` (cell of `[1 x 2]` `lat_range` pairs), `.lat_res .lon_range .lon_res .dtSec .mask`, with the values in Global Constraints.
  - `G = mission_gdop(H, M, window, grid, opts)`: `window = [t0 t1]` ND clock times; `grid` from `mission_grid`; `opts.score` [true] (R1 to R3), `opts.chunkSel` [all chunks] (indices, for tests), `opts.keepSeries` [false] (keep the full GDOP series; only sensible on the coarse grid), `opts.logFile` [''].
    Method: build the uniform epoch list `tq = t0 + (0:Nt-1)*dt` with `dt = grid.dtSec/tStar` and `Nt = round((t1 - t0)/dt)`; split it into maximal runs over which the on-station set (`H.phase == 3` at the epoch) is constant; for each run with M_on >= 1 satellites call `gdop_grid_from_states` on `tau = [tq(run); tq(run(end)) + dt]` (one extra epoch so the half-open resampling returns exactly the run's epochs) with `rTgt` interpolated from `H` for the on-station satellites only and `time_samples = numel(run)`; runs with no satellite on station get NaN GDOP and zero visibility without a call. Assemble per chunk, score per grid point, discard the series unless `keepSeries`.
    Output fields: `.t [Nt x 1]` ND, `.tSec`, `.nOn [Nt x 1]`, `.lat .lon [1 x Npts]` deg, `.poleGdop [Nt x 1]` (the grid point of lowest latitude at longitude 0), `.r1 .r2 .r3 .r3g7 [1 x Npts]` (per-point scores: fractions, fractions, hours, hours), `.r1Floor .r2Floor .r3Floor .r3g7Floor`, `.r1Nfail .r2Nfail .r3Nfail` (thresholds 0.90, 0.40, 9 h), `.pass` (struct of three logicals), `.worst` (struct: lat/lon of each floor), `.runs` (`[nRuns x 3]`: first index, last index, M_on), `.grid`, `.window`, `.wallS`, and `.gdop [Nt x Npts]` when `keepSeries`.
  - `V = mission_gdop('validate', M, opts)`: builds the on-slots history over one tulip period, runs the paper grid (or `opts.chunkSel`), loads `S2(1).r` from the reference file and returns `.status` (`'PASS' | 'FAIL' | 'UNRESOLVED'`), `.reason`, `.ours` and `.ref` (the three floors each), `.tol` (R1 0.005, R2 0.005, R3 0.10 h), `.perPointMaxDiff` when the reference stores per-point arrays for the chunk. `UNRESOLVED` when the file is missing, has no `S2`, or its stored phasing/sample count cannot be matched to `tulip_phasing()` and Nt 5760.

- [ ] **Step 1: Write the failing test** `tests/test_mission_gdop.m`:
  1. `mission_grid('paper')`: five chunks, the pairs in Global Constraints; total points 10,770 when expanded with `lat_res 0.5` and `lon 0:1:358`. `mission_grid('coarse')`: cadence 300 s.
  2. **Stitch.** On-slots history over 2 days, modified so satellites 5 and 6 come on station at t = 0.7 d and 1.3 d (`H.phase` edited in the test). `G = mission_gdop(H, M, [0 2 d], coarse, struct('keepSeries', true, 'score', false))`: `numel(G.t) == round(2*86400/300)`; `all(abs(diff(G.t) - 300/M.tStar) < 1e-12)` (no duplicated or dropped epoch); `G.nOn` steps 4 -> 5 -> 6 at the right epochs; `G.runs` has three rows.
  3. **Equivalence at a seam-free call.** The same history with all six on station: `mission_gdop`'s `G.gdop` equals a direct `gdop_grid_from_states` call on the same epochs to 1e-9 relative on the finite entries (the wrapper adds nothing).
  4. **Fewer than four.** `H.phase` edited so only three are on station for the first half day: GDOP NaN there, `G.nOn == 3`, no error.
  5. **Validation seam, one chunk.** `V = mission_gdop('validate', M, struct('chunkSel', 5))` (latitudes -77.5 to -75, where the reference's R2 worst point (-75, 151) lies): `V.status` is `'PASS'` or the per-point comparison is reported; assert `abs(V.ours.r2Floor - 0.754032) < 0.005` for this chunk if the reference's worst point is in it (the test reads the reference's worst-point latitude and skips the assertion with a printed reason otherwise).
  6. **UNRESOLVED.** `mission_gdop('validate', M, struct('refFile', fullfile(tempdir, 'nope.mat')))` returns status `'UNRESOLVED'`.
  7. **Deliberate failure.** `mission_gdop('validate', M, struct('chunkSel', 5, 'dtauOverride', M.dtau6 + [0 0 0 0 0 0.3]))` (one satellite moved off its slot) returns `'FAIL'`.

- [ ] **Step 2: Run it.** Expected: `Unrecognized function 'mission_grid'`.

- [ ] **Step 3: Implement both.** Before the full implementation, time one paper chunk over one day through `gdop_grid_from_states` plus the three scorers and write the measured seconds and peak memory (from `memory` is unavailable on macOS: use array sizes) into the report; `mission_gdop` logs per-chunk wall time. If one-period validation of a single chunk exceeds 5 minutes, check 5 uses a reduced longitude step and the full comparison moves to the study (state this in the report).

- [ ] **Step 4: Run the test; add to `run_tests`; run it** (5 of 5).

- [ ] **Step 5: Commit** (controller): `lunar_pnt_mission: mission_gdop (per-interval on-station GDOP through proj7's engine, R1-R3 by run_requirements_table's recipe on the proj7 scorers, validation against the saved Table 1 column) and mission_grid`.

---

### Task 6: The study script

**Files:**
- Create: `lunar_pnt_mission_study.m`, `tests/test_lunar_pnt_mission_study.m`, `README.md`
- Modify: `tests/run_tests.m`

**Interfaces:**
- Consumes: everything above.
- Produces: `lunar_pnt_mission_study.m`, a SCRIPT in the style of `orbit_transfer/GTO_DRO/gto_dro_study.m` (read it first: numbered `%%` sections, a knob block in section 0, a `study_opts` struct in the caller's workspace overriding knobs, a `gate(name, pass, detail)` table printed at the end, a VERDICT line, results saved with `run_meta`). Knobs: `tol` 0.005, `nMaxRevs` 10, `dtSec` 300, `serviceDays` 30, `scoreGrid` `'paper'`, `figureGrid` `'coarse'`, `makeFigures` true, `makeMovie` false, `headless` false, `outDir` `results/`, `maxFallback` 3 (candidates tried per spacecraft), `chunkSel` [] (all), `validate` true.
  Sections and their gates:
  0. Knobs, paths.
  1. Inputs: `M = mission_inputs()`; gate I1 = all input checks.
  2. Schedule: `S`; print the 6 x 6 best-arrival table (days) and the assignment; sensitivity rows for `tol` in `[0.0025 0.005 0.01]` and `nMaxRevs` in `[5 10 20]` (last-arrival day each); gate C1 = a complete assignment exists.
  3. Polish: for each spacecraft, walk `S.ranked{k, assign(k)}` in arrival order up to `maxFallback` until `P.ok`; print one row per spacecraft (cell, wait revs, t_dep d, t_arr d, t_f cell and polished, evaluations, flown miss, propellant kg, dV, family changed, which candidate rank was used); gates P1..P6 = each spacecraft accepted, PF = flown-miss gate on all, PP = pointwise, PH = hypothesis, PC = conjugate (PASS / FAIL / UNRESOLVED as in `transfer_study`). If a fallback moved a spacecraft's arrival, recompute and print the realised last-arrival time.
  4. Histories: `H`; gates H1 = every `jumpKm.dep < 1`, H2 = every `jumpKm.arr < 100`.
  5. Validation: `V = mission_gdop('validate', M, ...)` on the full paper grid; gate G0 = `V.status` PASS. If G0 is not PASS the script prints the VERDICT as `STOPPED AT VALIDATION`, saves, and returns before section 6.
  6. GDOP: assembly window `[H.tFirst, H.tLast]` on the coarse grid with `keepSeries` (pole GDOP series, satellites on station, first time pole GDOP < 6); service window `[H.tLast, H.tLast + 30 d]` on the paper grid scored; the one-period on-slots numbers beside the 30-day ones; gates R1, R2, R3 from `G.pass`; the GDOP < 7 R3 row printed.
  7. Figures (each to `outDir` as PNG, slide-size fonts, legends on every axes): `schedule.png` (per-spacecraft bars: DRO wait, transfer, on station, arrivals marked), `assembly_3d.png` (rotating frame, DRO, tulip, six transfers, Moon), `pole_gdop.png` (pole GDOP vs days with arrival lines and the GDOP = 6 line; satellites-on-station on a second axis), `service_map.png` (per-point R2 fraction over the paper grid in a polar south-pole view), `propellant.png` (table figure). Figures are created invisible when `headless`.
  8. Movie: when `makeMovie`, in the style of `pumpkynPie/Demos/lowThrust_constellation_animation.m`, frames resized to 1280 x 720; otherwise skipped with a printed line.
  9. Gate table, VERDICT line (`VERDICT: PASS` only if every gate is PASS; otherwise the failing gate names), save `results/lunar_pnt_mission_study.mat` with `M`-free content: `S, P, H` (positions as single), `V`, `Gassembly`, `Gservice`, `gates`, `run_meta`.

- [ ] **Step 1: Write the failing test** `tests/test_lunar_pnt_mission_study.m`: sets `study_opts = struct('headless', true, 'outDir', fullfile(tempdir, 'lunar_pnt_test'), 'chunkSel', 5, 'serviceDays', 3, 'makeMovie', false)` and runs the script; asserts the `.mat` exists with `run_meta`; a gate row exists for every name above; I1, C1, PF, H1, H2 are PASS; the VERDICT line was printed (capture with `evalc`); the five PNGs exist. With `chunkSel = 5` and 3 service days the R gates are computed on a reduced problem and are asserted present, not asserted PASS (R3 needs at least one full day bin: assert it is finite). A second run with `study_opts.tol = 1e-7` asserts gate C1 FAIL and a VERDICT naming it (the schedule gate can fail).

- [ ] **Step 2: Run it.** Expected: script not found.

- [ ] **Step 3: Write the script and the README skeleton** (purpose in three sentences, how to run, the knob table, the gate table with one line per gate, the refusals from spec section 9, a "Measured results" heading left for Task 7).

- [ ] **Step 4: Run the test.** Record the wall time of each section.

- [ ] **Step 5: Add to `run_tests`, run it** (6 of 6).

- [ ] **Step 6: Commit** (controller): `lunar_pnt_mission: lunar_pnt_mission_study (schedule, exact-slot polish with fallback, histories, validation gate, assembly and 30-day service GDOP with R1-R3, figures), headless test, README`.

---

### Task 7: The full run and the record (controller task)

**Files:**
- Modify: `orbit_transfer/lunar_pnt_mission/README.md` ("Measured results (2026-11-..)"), `orbit_transfer/README.md` (one row)
- Results (not committed): `results/lunar_pnt_mission_study.mat`, the five PNGs

- [ ] **Step 1: Size the full run** from Task 5's measured chunk time: full paper grid, one period (validation) + 30 days (service). If the estimate exceeds 2 hours, run it as a nohup'd script file with a log (matlab-campaign discipline) in a separate `matlab -batch`, otherwise in the shared session.
- [ ] **Step 2: Run the study with default knobs.** Read the gate table.
- [ ] **Step 3: If G0 (validation) is not PASS:** stop; report the three floors against the reference and the per-point differences; do not edit tolerances to pass.
- [ ] **Step 4: README "Measured results":** the assignment, each spacecraft's wait, departure and arrival day, polished t_f, propellant and dV; the last-arrival day; the sensitivity rows; the validation numbers against the reference; the 30-day R1, R2, R3 floors (and the one-period ones beside them); every gate verdict, including any spacecraft that needed a fallback candidate or failed; the figures by file name. Numbers only as measured.
- [ ] **Step 5: `orbit_transfer/README.md` row**, memory note, commit and push.

---

## Self-review (done while writing)

- **Spec coverage:** section 3 folder -> Tasks 1 to 6 file map; section 4 -> Task 2; section 5 -> Task 3 (root-find, pinned epoch, departure >= 0, family flag, four instruments, fallback order in Task 6); section 6 histories -> Task 4, geometry / two windows / validation gate -> Task 5 and Task 6 sections 5 and 6; section 7 -> Task 6; section 8 tests -> one test per task plus the deliberate-failure experiments (Task 1 step 6, Task 3 check 4, Task 4 check 3, Task 5 check 7, Task 6 second run); section 9 refusals -> README (Task 6); section 10 risks: shadowing (Task 1), compute sizing (Task 5 step 3, Task 7 step 1), basin loss (Task 3 inner-solve rule, Task 6 fallback), no admissible candidate (Task 2 refusal, Task 6 sensitivity rows); section 11 acceptance -> Task 7. The alternatives-sensitivity sentence of decision 5 is NOT planned (it needs off-grid phase binning the scheduler does not have); recorded here as a deliberate omission to raise with Mike rather than a silent gap.
- **Types:** `M` fields used later (`depState`, `arrState`, `slotPhase0`, `droPhase0`, `phys`, `propKgPerNd`, `tf`, `sD`, `sA`, `tauD`, `tauT`, `dtau6`, `tStar`, `lStar`, `muStar`) are all produced by Task 1; `S.ranked{k, j}` rows `[tArr iD iA nRev resid]` (Task 2) are what Task 6 turns into `cand` structs for Task 3; `P.tDep/.tArr/.flight` (Task 3) are what Task 4 consumes; `H.phase/.t/.r/.tFirst/.tLast` (Task 4) are what Task 5 consumes; `G.pass`, `V.status` (Task 5) feed Task 6's gates.
- **Review Focus:** each line names its pinning task (2, 3 + 6, 5, 1, 5).
