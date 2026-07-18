# Table-3 Reproducer Engine — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** A click-go, unattended, resumable engine that re-solves our version of Gergaud Table 3 **from scratch** (10/5/2.5/1/0.5 N), capturing every per-rung recipe as runnable code, and built so the same engine extends to 0.2/0.1 N later.

**Architecture:** One parameterized `reproduce_row(T)` composes the EXISTING certified drivers (`run_mintime_mee`, `run_transfer_mee`, `psr_mee_refine`) with a per-rung **recipe registry** and a **pluggable anchor strategy** (`coldB`/`chain`/`smallN_first`/`R0law`). Isolation from the campaign caches is by `REPRO_`-prefixed tags + a `results/repro/` namespace — **no solver/driver edits**. A shell watchdog runs each rung in its own process, relaunching on the uncatchable MEX crashes. Design: `docs/superpowers/specs/2026-07-18-table3-reproducer-engine-design.md`.

**Tech Stack:** MATLAB R2025b, CasADi 3.7.0 + IPOPT, the `earth_elliptic_to_geo/` module + its drivers.

## Global Constraints

- **MATLAB R2025b ONLY** — `/Applications/MATLAB_R2025b.app/bin/matlab`. Run tests via `matlab -batch "run('/abs/path/test_X.m')"`; add `addpath(fullfile(getenv('HOME'),'casadi-3.7.0'))` at the top of any test that solves.
- **Never use `i`/`j` as loop/index variables** (imaginary unit) — use `k`, `kk`, `idx`.
- **Every MATLAB function needs a full header** (purpose, INPUTS with sizes, OUTPUTS with sizes, REFERENCES).
- **NO solver-core or driver edits.** The engine COMPOSES `run_mintime_mee`/`run_transfer_mee`/`psr_mee_refine`/`casadi_lt_mee`/`interp_warmstart`/`mee_seed` through their existing options only. If a task appears to need a driver edit, STOP and report it.
- **From-scratch isolation:** all engine solves use `REPRO_`-prefixed tags (`cfg.tag`, `cfg.fuelTag`, PSR `opts.tag`) and, where supported, `outDir`/output under `results/repro/`. This forces a genuine cold re-solve (Stage-A lookup never finds a `REPRO_` fuel cache on first run) and never clobbers the campaign `MEE_M2_*`/`MEE_mintime_*` caches. A `reuseCampaignCache` flag (default false) may instead read the campaign caches for a fast "present" mode.
- **Honesty:** a stage that does not certify (defect/termErr gate) THROWS; a row is never fabricated. 0.2/0.1 N are seeded in the registry but NOT executed in this build.
- **Verify, don't assume:** every reproduced row asserts against `table3_certified` (§ Task 1); a run that lands elsewhere is a loud failure.
- **Lib-extraction-friendly:** a code-tidy + generic-library refactor is the NEXT effort. Keep every new function generic, single-responsibility, with clean dependencies (no reaching into another function's cache files by hand; pass data explicitly) so it lifts cleanly into a future `lib/`.

Certified numbers (source: `CAMPAIGN.md`): 10 N → 1377.10 kg / 19 sw / 7.326 rev / tfmin 22.2206 ND; 5 N → 1364.54 / 32 / 14.157 / 44.6796; 2.5 N → 1369.79 / 76 / 27.841 / 89.253; 1 N → 1371.44 / 171 / 69.152 / 223.808; 0.5 N → 1375.28 / 362 / 138.597 / 446.28 (R0-law estimate, `anchorSource='R0law'`). R0 constant = 223.14 ND.

---

## File Structure

- `earth_elliptic_to_geo/table3_certified.m` — CREATE: certified numbers table (pure).
- `earth_elliptic_to_geo/table3_recipes.m` — CREATE: recipe registry (pure).
- `earth_elliptic_to_geo/verify_row.m` — CREATE: assert a reproduced row vs certified (pure).
- `earth_elliptic_to_geo/reproduce_row.m` — CREATE: the engine (composes drivers; coldB/chain/R0law).
- `earth_elliptic_to_geo/anchor_smallN_first.m` — CREATE: the 1 N anchor strategy (harvested).
- `earth_elliptic_to_geo/reproduce_table3_collect.m` — CREATE: assemble + print the reproduced table.
- `earth_elliptic_to_geo/reproduce_table3.m` — CREATE: thin in-process wrapper (top rungs).
- `earth_elliptic_to_geo/reproduce_table3.sh` — CREATE: watchdog orchestrator.
- `earth_elliptic_to_geo/run_task9_rung.m` — MODIFY: thin deprecated shim delegating to `reproduce_row`.
- Tests: `test_table3_recipes.m`, `test_verify_row.m`, `test_reproduce_row_smoke.m`, `test_anchor_smallN_first.m`.

Interfaces (established here):
- `cert = table3_certified(T)` → struct `.thrustN .m_f_kg .switches .revs .tfmin .anchorSource`.
- `recipe = table3_recipes(T)` → struct `.thrustN .anchor(.strategy,.npr,.nprLo,.nprHi,.mtMaxIter,.warmFrom) .fuel(.npr,.seedThr,.maxIter,.warmFrom) .psr([]|struct .maxRounds,.nbr,.globalEvery,.globalFactor) .tfmin_or_R0`.
- `pass = verify_row(row, cert, tol)` → logical + throws on mismatch (tol has `.m_f_kg`, `.revsRel`, `.switchesAbs`).
- `row = reproduce_row(T, opts)` → row struct (from `gergaud_row`) + writes `results/repro/REPRO_row_T<10T>.mat` holding `.row .anchor .fuel(.X,.U,.dL,.sigma) .certified`.
- `anchor = anchor_smallN_first(T, par, warmAnchor, opts)` → struct matching `run_mintime_mee` output (`.tfmin .tfmin_h .dL_mt .revs .N .solverOut .certified`).

---

## Task 1: Foundation — certified table, recipe registry, verify (pure, no solve)

**Files:** Create `table3_certified.m`, `table3_recipes.m`, `verify_row.m`; Test `test_table3_recipes.m`, `test_verify_row.m`.

**Interfaces:** Produces the three pure functions above. Consumes nothing (data + logic only).

- [ ] **Step 1: Write failing tests**

```matlab
% test_table3_recipes.m
here = fileparts(mfilename('fullpath')); cd(here);
c10 = table3_certified(10);
assert(abs(c10.m_f_kg-1377.10)<1e-2 && c10.switches==19 && abs(c10.revs-7.326)<1e-3);
c05 = table3_certified(0.5);
assert(strcmp(c05.anchorSource,'R0law') && abs(c05.tfmin-446.28)<0.1);
r1 = table3_recipes(1);
assert(strcmp(r1.anchor.strategy,'smallN_first') && r1.anchor.nprLo==15 && r1.anchor.nprHi==25);
assert(~isempty(r1.psr) && r1.psr.maxRounds==2);
r05 = table3_recipes(0.5);
assert(strcmp(r05.anchor.strategy,'R0law') && r05.fuel.npr==12 && r05.psr.maxRounds==5);
r10 = table3_recipes(10);
assert(strcmp(r10.anchor.strategy,'coldB') && isempty(r10.psr));
% seeded deep rungs exist but are flagged not-run
r02 = table3_recipes(0.2);  assert(r02.seeded==true);
fprintf('test_table3_recipes PASSED\n');
```
```matlab
% test_verify_row.m
here = fileparts(mfilename('fullpath')); cd(here);
cert = table3_certified(10);
tol  = struct('m_f_kg',0.5,'revsRel',0.01,'switchesAbs',0);
good = struct('thrustN',10,'m_f_kg',1377.0,'switches',19,'revs',7.33);
assert(verify_row(good,cert,tol)==true);
bad  = struct('thrustN',10,'m_f_kg',1300,'switches',19,'revs',7.33);
threw=false; try, verify_row(bad,cert,tol); catch, threw=true; end
assert(threw,'verify_row must throw on a mass mismatch');
fprintf('test_verify_row PASSED\n');
```

- [ ] **Step 2: Run tests — expect FAIL** (functions absent).

- [ ] **Step 3: Implement**

`table3_certified.m` — a switch/lookup returning the numbers above; error on an unknown T. `table3_recipes.m` — return the per-rung recipe struct per the registry table in the spec (§6), with `.seeded=true` and `.psr.maxRounds>=5` for 0.2/0.1 and `.seeded=false` for 10..0.5. `verify_row.m` — compare `m_f_kg` (abs `tol.m_f_kg`), `revs` (relative `tol.revsRel`), `switches` (abs `tol.switchesAbs`); on any breach `error('verify_row:mismatch', ...)` naming the field, expected, got; else return true. Full headers.

- [ ] **Step 4: Run tests — expect PASS.**

- [ ] **Step 5: Commit** `feat(repro): certified table + recipe registry + verify_row (pure, tested)`.

---

## Task 2: `reproduce_row` engine — coldB / chain / R0law + live 10 N re-solve

**Files:** Create `reproduce_row.m`; Test `test_reproduce_row_smoke.m`.

**Interfaces:** Consumes `table3_recipes`, `table3_certified`, `verify_row` (Task 1); `run_mintime_mee`, `run_transfer_mee`, `psr_mee_refine` (existing). Produces `reproduce_row(T, opts)` + `results/repro/REPRO_row_T<10T>.mat`.

- [ ] **Step 1: Write the failing live smoke test**

```matlab
% test_reproduce_row_smoke.m — 10 N re-solved FROM SCRATCH (coldB), asserted vs certified
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
here = fileparts(mfilename('fullpath')); cd(here);
row = reproduce_row(10);          % coldB: no prev rung needed; writes results/repro/
cert = table3_certified(10);
assert(abs(row.m_f_kg-cert.m_f_kg)<0.5, '10 N m_f reproduced');
assert(row.switches==cert.switches, '10 N switches reproduced');
assert(abs(row.revs-cert.revs)/cert.revs<0.01, '10 N revs reproduced');
assert(row.certified==true);
assert(isfile(fullfile(here,'results','repro','REPRO_row_T100.mat')));
fprintf('test_reproduce_row_smoke PASSED\n');
```

- [ ] **Step 2: Run — expect FAIL** (`reproduce_row` absent). (This test SOLVES; ~minutes.)

- [ ] **Step 3: Implement `reproduce_row.m`**

Structure (compose drivers with `REPRO_` tags; write summaries to `results/repro/`):
```
function row = reproduce_row(T, opts)
%   opts: .reuseCampaignCache(false) .m0kg(1500) .ispS(2000)
recipe = table3_recipes(T);  cert = table3_certified(T);  par = kepler_lt_params(T,m0,isp);
reproDir = results/repro; ensure exists.
tagFuel = sprintf('REPRO_MEE_M2_%s', ttag(T));  tagMt = sprintf('REPRO_MEE_mintime_%s', ttag(T));
% ---- ANCHOR (strategy dispatch) ----
switch recipe.anchor.strategy
  case 'coldB'   : anchor = run_mintime_mee(T, recipe.anchor.npr, ...
                     struct('tag',tagMt,'fuelTag','REPRO_none','maxIter',3000));  % Stage A skipped
  case 'chain'   : prev = load_prev(T, recipe.anchor.warmFrom, reproDir);          % prev repro anchor
                   dLGuess = prev.anchor.dL_mt*(recipe.anchor.warmFrom/T);
                   anchor = run_mintime_mee(T, recipe.anchor.npr, struct('tag',tagMt, ...
                     'fuelTag','REPRO_none','maxIter',recipe.anchor.mtMaxIter, ...
                     'nRevSeed',max(1,round(dLGuess/(2*pi))), 'warmStartAnchor', ...
                     struct('X',prev.anchor.solverOut.X,'U',prev.anchor.solverOut.U,'dL',dLGuess,'N',prev.anchor.N)));
  case 'R0law'   : anchor = struct('tfmin', 223.14/T, 'anchorSource','R0law', ...);  % no solve
  case 'smallN_first': anchor = anchor_smallN_first(T, par, load_prev(...).anchor, recipe.anchor);  % Task 3
end
tfMinAnchor = anchor.tfmin;
% ---- FUEL ----
fuelCfg = struct('thrustN',T,'ctf',1.5,'tfMinAnchor',tfMinAnchor,'tag',tagFuel, ...
   'seedThr',recipe.fuel.seedThr,'betaMode','tangential','nodesPerRev',recipe.fuel.npr, ...
   'maxIter',recipe.fuel.maxIter,'m0kg',m0,'ispS',isp);
if ~isempty(recipe.fuel.warmFrom)
   prevF = load_prev(T, recipe.fuel.warmFrom, reproDir);
   dLGuessFuel = prevF.fuel.dL*(recipe.fuel.warmFrom/T);
   fuelCfg.warmStart = struct('sigma',prevF.fuel.sigma,'X',prevF.fuel.X,'U',prevF.fuel.U,'dL',dLGuessFuel);
end
res = run_transfer_mee(fuelCfg);  assert(res.report.certified, ...);   % THROW if not certified
sol = struct('X',res.fuel.X,'U',res.fuel.U,'dL',res.fuel.dL,'sigma',res.sigma);
rep = res.report;
% ---- PSR (optional) ----
if ~isempty(recipe.psr)
   psrOut = psr_mee_refine(res, struct('tag',[tagFuel '_PSR'],'outDir',reproDir, ...
       'maxRounds',recipe.psr.maxRounds,'nbr',recipe.psr.nbr, ...
       'globalEvery',getf(recipe.psr,'globalEvery',3),'globalFactor',getf(recipe.psr,'globalFactor',1.3), ...
       'maxIter',recipe.fuel.maxIter));
   sol = struct('X',psrOut.finalOut.X,'U',psrOut.finalOut.U,'dL',psrOut.finalOut.dL,'sigma',psrOut.finalSigma);
   rep = psr_report(psrOut);   % m_f_kg/switches/revs/edge/incDeg/defect/certified from finalOut
end
% ---- ROW + VERIFY + SAVE ----
row = gergaud_row(struct('thrustN',T,'tfmin_ND',tfMinAnchor,'ctf',1.5,'tf_ND',1.5*tfMinAnchor, ...
     'm_f_kg',rep.m_f_kg,'switches',rep.switches,'revs',rep.revs,'edge',rep.edge, ...
     'incl_deg',rep.incDeg,'defect',rep.defect,'certified',true, ...
     'note',ternary(strcmp(anchor.anchorSource,'R0law'),'0.5 N: R0-law tfmin estimate (anchor-free)','')));
verify_row(row, cert, defaultTol(T));   % exact sw for 10/5/2.5; +-few for 1/0.5 (Task-1/§11)
save(fullfile(reproDir, sprintf('REPRO_row_T%d.mat', round(10*T))), 'row','anchor','sol','rep');
fprintf(gergaud_row_str(row));
```
`load_prev(T, prevT, reproDir)` loads `REPRO_row_T<10*prevT>.mat` and errors with "run rung <prevT> first" if absent. `defaultTol(T)`: `switchesAbs=0` for T in {10,5,2.5}, `=3` (logged) for {1,0.5}; `m_f_kg=0.5`, `revsRel=0.02`. Write the full header; keep `load_prev`/`psr_report`/`defaultTol`/`ttag` as local helpers (generic, lib-friendly).

- [ ] **Step 4: Run the smoke test — expect PASS** (10 N re-solves from scratch, asserts vs certified, writes `results/repro/REPRO_row_T100.mat`).

- [ ] **Step 5: Commit** `feat(repro): reproduce_row engine (coldB/chain/R0law) + 10 N live re-solve smoke`.

---

## Task 3: `anchor_smallN_first` — the 1 N anchor strategy (harvested)

**Files:** Create `anchor_smallN_first.m`; wire it into `reproduce_row`'s `smallN_first` case; Test `test_anchor_smallN_first.m`.

**Interfaces:** Consumes `run_mintime_mee` (for the `nprLo` seed+continuation), `interp_warmstart`, `casadi_lt_mee`, `kepler_lt_params`. Produces `anchor_smallN_first(T, par, warmAnchor, aopts)` returning a struct matching `run_mintime_mee`'s output.

Source of the recipe (harvest verbatim logic, do not invent): `results/task7c_step1_manual.m` (manual relaxed-stall continuation: warmTight iff status∈{Solve_Succeeded,Solved_To_Acceptable_Level,Maximum_Iterations_Exceeded} AND defect<1e-4; maxIter 75, raised to 150 once defect<1e-3; NO decadeMin floor; always advance the warm start; keep best; per-round save) and `results/task7c_step1b_refine.m` (interp_warmstart from `nprLo` to `nprHi`, single `warmTight=true`, `maxIter=1500`).

- [ ] **Step 1: Write the failing test** (unit-level; no full 6-min solve required)

```matlab
% test_anchor_smallN_first.m — logic + resume-cache behavior, not a full solve
here = fileparts(mfilename('fullpath')); cd(here);
% (a) the continuation's warmTight predicate:
assert(smallN_warmtight('Maximum_Iterations_Exceeded', 1e-5)==true);
assert(smallN_warmtight('Restoration_Failed', 1e-9)==false);
assert(smallN_warmtight('Solve_Succeeded', 1e-3)==false);   % defect too high
% (b) the maxIter ramp: 75 default, 150 once defect<1e-3
assert(smallN_maxiter(1e-2)==75 && smallN_maxiter(5e-4)==150);
% (These small pure helpers are nested in / exported by anchor_smallN_first.m.)
fprintf('test_anchor_smallN_first PASSED\n');
```
(Expose `smallN_warmtight` and `smallN_maxiter` as testable helpers — either nested-function handles returned for test, or small local functions the test can reach; controller will confirm the mechanism at dispatch.)

- [ ] **Step 2: Run — expect FAIL.**

- [ ] **Step 3: Implement `anchor_smallN_first.m`**

Stage 1 (`nprLo`, default 15): seed by `interp_warmstart` of `warmAnchor` (the previous rung's converged anchor) onto the `nprLo` grid; run the manual relaxed-stall continuation loop harvested from `task7c_step1_manual.m` (per-round cache `REPRO_<tag>_smallN_round%02d.mat` under `results/repro/`, resume on restart, wall-budget guard, `isGood = success && maxDefect<1e-8 && termErr<1e-8`). Stage 2 (`nprHi`, default 25): `interp_warmstart` the certified `nprLo` anchor onto the `nprHi` grid, one `casadi_lt_mee(...,'mode','mintime','warmTight',true,'maxIter',1500)`; assert certified. Return `struct('tfmin',o.tf,'tfmin_h',o.tf*par.TU_s/3600,'dL_mt',o.dL,'revs',o.dL/(2*pi),'N',Nhi,'solverOut',o,'certified',true,'anchorSource','smallN_first')`. Full header citing the two source scripts.

- [ ] **Step 4: Run the unit test — expect PASS.** Then a controller-run live check (dispatch note): `reproduce_row(1)` reaches a certified anchor (the ~6-min `nprLo` stage + refine) — this is exercised in the Task-6 background validation, not required to fully run inside this task's review.

- [ ] **Step 5: Commit** `feat(repro): anchor_smallN_first 1 N strategy (harvested from task7c hand-scripts)`.

---

## Task 4: Orchestration — collect, in-process wrapper, watchdog

**Files:** Create `reproduce_table3_collect.m`, `reproduce_table3.m`, `reproduce_table3.sh`.

**Interfaces:** Consumes `reproduce_row`, `table3_certified`, `gergaud_row_str`, `results/repro/REPRO_row_T*.mat`.

- [ ] **Step 1: Write the failing test** (collect assembles a table from existing repro rows)

```matlab
% (append to test_table3_recipes.m or a new test_reproduce_collect.m)
here = fileparts(mfilename('fullpath')); cd(here);
% requires results/repro/REPRO_row_T100.mat from Task 2's smoke:
if isfile(fullfile(here,'results','repro','REPRO_row_T100.mat'))
  tbl = reproduce_table3_collect([10]);   % returns a struct array + prints
  assert(abs(tbl(1).m_f_kg-1377.10)<0.5 && tbl(1).thrustN==10);
  fprintf('test_reproduce_collect PASSED\n');
else
  fprintf('test_reproduce_collect SKIPPED (no repro row yet)\n');
end
```

- [ ] **Step 2: Run — expect FAIL/SKIP.**

- [ ] **Step 3: Implement**

`reproduce_table3_collect.m(thrustList)` — load each `REPRO_row_T*.mat`, print the fixed-width Table-3 block (reuse `gergaud_row_str`) + the R0-law spread, and return the row struct array; a missing rung prints `MISSING` (not an error). `reproduce_table3.m(thrustList)` — thin in-process wrapper: `for T = thrustList, reproduce_row(T); end; reproduce_table3_collect(thrustList)` (convenient for the crash-free top rungs; documents that deep rungs should use the watchdog). `reproduce_table3.sh` — model on `PSR/psr_batch.sh` + `run_task9_watchdog.sh`: for each `T` in `10 5 2.5 1 0.5`, loop `matlab -batch "cd('.../earth_elliptic_to_geo'); addpath(getenv('HOME')+'/casadi-3.7.0'); reproduce_row($T)"` in its own process; on nonzero exit, classify via a fresh `~/matlab_crash_dump.*`, log the attempt, and relaunch (caches resume) up to a per-rung cap (e.g. 8); then `reproduce_table3_collect`. Log to `results/repro/reproduce_table3.log`.

- [ ] **Step 4: Run the collect test — expect PASS.** Shellcheck `reproduce_table3.sh` (or a dry `bash -n`).

- [ ] **Step 5: Commit** `feat(repro): watchdog orchestrator + in-process wrapper + table collector`.

---

## Task 5: Subsume `run_task9_rung` + docs

**Files:** Modify `run_task9_rung.m` (thin shim), `README.md`, `TODO.md`.

- [ ] **Step 1:** Reduce `run_task9_rung.m` to a thin deprecated shim: a header noting it is superseded by `reproduce_row`/`table3_recipes`, delegating the deep-rung body to the engine's `chain`+coarse-fuel+PSR path (or, if a clean delegation is not 1:1, keep it callable but add a prominent "DEPRECATED: use reproduce_row" banner and a pointer). Do not silently break its existing signature if anything references it (grep first).

- [ ] **Step 2:** `README.md` — add a "Reproducing from scratch" subsection: `reproduce_table3.sh` (watchdog) / `reproduce_table3.m` (in-process) / `reproduce_row(T)`, the `results/repro/` namespace, the recipe registry, and the verify-against-certified guarantee. `TODO.md` — add under "Done" the reproducer engine; add a "Not a goal / future" note capturing the upcoming **code-tidy + generic `lib/` refactor** (so it is not lost) and that 0.2/0.1 N recipes are seeded-not-run.

- [ ] **Step 3: Commit** `docs(repro): document the reproducer engine; note the coming code-tidy/lib refactor`.

---

## Task 6: Full-ladder validation (controller-run, after the SDD tasks + final review)

Not a subagent code task — the controller executes the delivered engine:
- [ ] Run `reproduce_row(10)`, `(5)`, `(2.5)` live (in order; ~tens of minutes); confirm each asserts against certified.
- [ ] Launch `reproduce_table3.sh` (or the 1 N + 0.5 N rungs) under a **background watchdog**; let it run to completion (hours, resumable, auto-relaunch on crash); collect and report the reproduced 1 N and 0.5 N rows vs certified.
- [ ] Record the full reproduced table + any crash/relaunch counts in the ledger.

---

## Self-Review (completed)

- **Spec coverage:** registry/certified/verify (T1); engine coldB/chain/R0law + from-scratch isolation + 10 N live proof (T2); 1 N harvested smallN_first (T3); watchdog + collect + in-process (T4); subsume run_task9_rung + docs incl. lib-refactor note (T5); full background validation (T6). ✓
- **Placeholder scan:** code steps carry real code or precise composition; 0.2/0.1 are intentionally seeded-not-run. No TBD in executable paths.
- **Type consistency:** `anchor_smallN_first` returns the `run_mintime_mee` output shape so the fuel stage consumes anchors uniformly; `reproduce_row` writes a single `REPRO_row_T*.mat` schema the collector + `load_prev` both read.
- **No-driver-edit invariant** is a Global Constraint and is honored by composing through `tag`/`fuelTag`/`warmStart`/`outDir` only (confirmed those knobs exist).
