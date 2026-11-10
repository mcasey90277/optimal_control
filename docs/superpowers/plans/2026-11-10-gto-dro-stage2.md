# GTO -> DRO Min-Time Catalog (Stage 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the GTO -> DRO min-time costate catalog (16 sheets, orientDeg {0,90,180,270} x DRO period tau {0.5,1,2,3} ND, 12 x 6 phasing, thrust ladder from 15 N) through the shared ladder engine, with a pilot sheet that measures the three things stage 1 left open (top-rung basin, perigee column, deep rungs) before the fleet, an Earth-clearance sweep that closes the certifier gap, and the standard chain (census -> sweeps -> repair -> package -> conjugate -> gates -> audit -> pictures -> deliverable).

**Architecture:** `GTO_DRO/catalog/run_gto_dro_catalog.m` is a clone of `GTO_tulip/catalog/run_gto_catalog.m` with the DRO arrival, the Earth floor and a `planOnly` mode; `run_gto_dro_batched.sh` is its OS-timeout driver. Two new engine options (`seedTop`, `redo`, default off, byte-identical otherwise) let a campaign tool re-solve one cell's top rung from an external seed, which is how the perigee column (phase continuation from the neighbour cell) and refuted top-rung basins (reverse warm start from the cell's own 12 N entry) get repaired. `build_gto_dro_catalog.m` runs the post-fleet stages in the `build_70mN_library` style, each fenced, with a blocker list that refuses the deliverable.

**Tech Stack:** MATLAB R2026a (shared session via the MCP tools for development; `matlab -batch` only inside the nohup'd batched driver), CasADi 3.7 + IPOPT through `costate_common/casadi_mintime_dro`, pumpkyn/pumpkynPie, `costate_common` (thrust_ladder_library, ms_tfmin, seed_from_z8, fly_transfer, conj_catalog_pass, gates_catalog_pass, build_costate_catalog_family, nd_propulsion, phase_state, get_family_orbit), `GTO_DRO/lib` (gto_dro_entry, earth_clearance).

**Spec:** `docs/superpowers/specs/2026-11-09-gto-dro-mintime-design.md` (section 10 fixes the shape of this stage; sections 2, 4, 9, 11 carry the decisions and risks it argues from). Stage 1's measured findings are in `orbit_transfer/GTO_DRO/README.md` ("Findings (stage 1, 2026-11-09)") and must be read by every implementer whose task mentions a basin, the perigee or the Earth floor.

## Global Constraints

- Catalog regime: `m0kg = 150`, `ispS = 1710`, rungs from 15 N; warm recipe `thrLock = true`, `tf0 = 0.30` (top rung only; the pilot may re-base it), `maxCpuSec = 600`, `maxIter = 6000`, `N = 400`, `floorKm = 500` (Moon), `floorEarthKm = 200` (Earth), `gateKm = 100`, `accTol = 1e-6`.
- Sheet key rule (spec section 2, GTO catalog): `opts.tauDRO = orientDeg` for every sheet; the DRO period rides in `arrParams.tau` and in `meta.arrParams`, never in the key. Sheet file tag `gto_dro_o%03d_tau%g` (e.g. `gto_dro_o000_tau1`, `gto_dro_o090_tau0.5`).
- Phasing grid: `nD = 12`, `nA = 6`, `sD0 = 0`, `sA0 = 0` (departure phase = time-since-perigee fraction of the Kepler period; arrival phase = fraction of the DRO period, phase 0 = `get_family_orbit`'s first state). The grid is uniform in time fraction (engine limitation, recorded in `run_gto_catalog.m`'s header).
- Engine changes are options with a default that leaves every other campaign byte-identical, pinned by a test, and documented in `costate_common/README.md`. No campaign folder goes on `costate_common`'s path; `costate_common/tests` must run without `GTO_DRO` on the path.
- House MATLAB style: header quartet (Purpose / Inputs with sizes / Outputs with sizes / Revision History `M. Casey (c) MM/DD/YYYY`, `Copyright Coorbital Inc.`), `=` aligned in knob blocks, `if nargin == 0` self-demo on library functions, never `i` or `j` as indices, never `norm` (use `vecnorm` or `sqrt(sum(...))`), no `%#ok`, no emojis, ASCII-only READMEs, no author signature lines in READMEs.
- Every script that solves writes to a git-ignored `results/` folder with a `run_meta` (date, git HEAD, MATLAB version, knobs); nothing under `results/` is committed.
- Long runs follow the matlab-campaign discipline: nohup'd script files, log to files, OS-level timeout per batch, attempt counters on disk before a solve, completion marker in a log file, stall-aware monitors.
- Subagents never commit; the controller commits after each clean review with a message ending `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`; git via `/Library/Developer/CommandLineTools/usr/bin/git`; never stage the unrelated `missiles/`, `papers/`, `.DS_Store` or `min_fuel_papers` deletions.
- MATLAB during development runs only through the MCP tools (`evaluate_matlab_code`, `run_matlab_file`, `run_matlab_test_file`) in the shared R2026a session; never `matlab -batch` from a subagent.

## Review Focus

1. A `redo` call on a cell whose repaired top rung FAILS must leave the sheet with the cell's original entries, not an empty cell (the engine clears the cell before re-solving; the repair tool restores from its backup). Pinned in Task 4's test (`test_repair_gto_dro_cells`, the bogus-seed case).
2. `seedTop` with a wrong mesh (`X` not `[7 x N+1]`, `U` not `[4 x N+1]`) must refuse with `thrust_ladder_library:seedTop` before any solve, never reach CasADi with a shape mismatch. Pinned in Task 1's test.
3. A sheet on disk whose `rungs` differ from the call's must refuse (`run_gto_dro_catalog:rungMismatch`) before the engine's resume path can silently reset it. Pinned in Task 3's test (planOnly cannot exercise the on-disk guard, so the test writes a fake sheet `.mat` with `rungs = [15 12]` and calls the real front door on that one sheet with `maxCells = 0`, expecting the refusal).
4. The Earth sweep must demote an entry whose flown arc dips below the floor to `OK = false` in the sheet AND keep a backup of the pre-sweep sheet, so a wrong floor is recoverable. Pinned in Task 2's test (impossible floor 1e6 km).
5. `planOnly` output must carry `floorEarthKm = 200`, `tauDRO = orientDeg`, `arrParams.tau = tau` and 16 distinct tags; a consumer that trusts the plan must not get a tulip arrival by inheritance from the cloned front door. Pinned in Task 3's test.

---

## File map

```
orbit_transfer/costate_common/
  thrust_ladder_library.m              MODIFY: opts.seedTop, opts.redo (default off)
  tests/test_ladder_seedtop.m          CREATE
  README.md                            MODIFY: the two options, one paragraph
orbit_transfer/GTO_DRO/
  lib/seed_top_from_entry.m            CREATE: sheet entry -> engine seedTop struct
  lib/earth_clearance_sweep.m          CREATE: per-sheet flown Earth altitude, sidecar, writeback
  lib/repair_gto_dro_cells.m           CREATE: 'unsolved' (phase continuation) and 'basin' (reverse warm start) modes
  catalog/run_gto_dro_catalog.m        CREATE: the fleet front door (clone of run_gto_catalog with DRO arrival)
  catalog/run_gto_dro_batched.sh       CREATE: OS-timeout driver
  catalog/build_gto_dro_catalog.m      CREATE: post-fleet chain (census .. deliverable)
  catalog/audit_gto_dro_entries.m      CREATE: 3-entry replay with the Earth floor as a gate
  catalog/costate_catalog_pick.m       COPY from GTO_tulip/catalog (vendored, deliverable convention)
  catalog/costate_lib_describe.m       COPY from GTO_tulip/catalog
  tests/fixture_demo_sheet.m           CREATE: cached 1 x 1 demo sheet (15, 12 N) for the tests
  tests/test_seed_top_from_entry.m     CREATE
  tests/test_earth_clearance_sweep.m   CREATE
  tests/test_run_gto_dro_catalog.m     CREATE
  tests/test_repair_gto_dro_cells.m    CREATE
  tests/test_build_gto_dro_catalog.m   CREATE
  tests/run_tests.m                    MODIFY: add the five names
  README.md                            MODIFY: stage 2 section, pilot findings, fleet status
  reviews/                             CREATE (Task 7): Astra review + dispositions
orbit_transfer/README.md               MODIFY (Task 8): GTO_DRO row mentions the catalog
```

Results layout (git-ignored, under `GTO_DRO/results/`): `catalog/<tag>.mat` + `<tag>_progress.txt` + `combined_progress.txt` (sheets, engine-native), `catalog/<tag>_earth.mat`, `catalog/<tag>_repairs.mat` (sidecars), `catalog/gto_dro_o000_tau1_deep.mat` (the deep probe), `costate_catalog_gto_dro.mat` (the package), `costate_catalog_gto_dro_conjprog.mat` / `_gatesprog.mat` (sweep sidecars), `audit_gto_dro.mat`, `pictures/`, `deliverable/`, `test_fixture/demo_sheet.mat`.

---

### Task 1: Engine options `seedTop` and `redo` in `thrust_ladder_library`

**Files:**
- Modify: `orbit_transfer/costate_common/thrust_ladder_library.m` (header Inputs block after `.maxCpuSec`; the todo/resume block at lines 138-172; the per-cell seed line at 195)
- Modify: `orbit_transfer/costate_common/README.md` (one paragraph under the ladder engine's entry)
- Test: `orbit_transfer/costate_common/tests/test_ladder_seedtop.m`

**Interfaces:**
- Consumes: `casadi_mintime_dro(rv0, rvf, Tmax, c, muStar, N, X0, U0, tf0, opts)` with `X0 [7 x N+1] = [r; v; m]`, `U0 [4 x N+1] = [alpha; thr]`.
- Produces: `opts.seedTop` (struct with `.X [7 x N+1]`, `.U [4 x N+1]`, `.tf` scalar ND; applied to the TOP rung of every cell this call solves; `[]` = off) and `opts.redo` (logical, default false; when true, every cell in `opts.cells` is re-solved even if it already has OK rungs: its attempt counter is reset to 0 and its TF/FLYKM/ACCDZ/RES/WALL/OK/Z8 slices are cleared before the ladder starts). `redo` without `cells` refuses with `thrust_ladder_library:redoCells`. Neither option changes the saved variable list or the resume size check.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_ladder_seedtop()
% TEST_LADDER_SEEDTOP  opts.seedTop / opts.redo of thrust_ladder_library.
%
% Four checks on the stage 1 GTO -> DRO demo cell (orient 0, DRO tau 1,
% departure phase 0.10, arrival phase 0, Earth floor 200 km), the one cell
% whose two top-rung basins are MEASURED (GTO_DRO/README.md, stage 1):
%   (a) a bad seedTop mesh refuses before any solve
%   (b) redo without cells refuses
%   (c) a normal call solves rungs [15 12]; the 15 N cold start lands near
%       0.3040 ND and the 12 N rung near 0.3070 ND
%   (d) redo + seedTop built from the 12 N entry re-solves the cell and the
%       15 N rung lands on the faster basin (t_f < 0.29 ND), OK, with the
%       attempt counter reset to 1
% The option-off identity is covered by the existing bitwise regression
% DRO_tulip/direct/certify/tests/test_minenergy_objective (no seedTop,
% no redo -> same code path as before this change).
%
% INPUTS:  none
% OUTPUTS: ok - logical
%
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
if isempty(which('pumpkyn.cr3bp.tfMin'))
    od = cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(od);
end
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));

outDir = fullfile(tempdir, 'test_ladder_seedtop');
if isfolder(outDir), rmdir(outDir, 's'); end
mkdir(outDir);
outMat = fullfile(outDir, 'cell.mat');

N    = 400;
base = struct('rungs', [15 12], 'ispS', 1710, 'm0kg', 150, 'N', N, ...
    'floorKm', 500, 'floorEarthKm', 200, 'maxIter', 6000, 'maxCpuSec', 600, ...
    'gateKm', 100, 'accTol', 1e-6, 'nD', 1, 'nA', 1, 'sD0', 0.10, 'sA0', 0, ...
    'depFamily', 'gto', 'depParams', struct('orientDeg', 0), ...
    'arrFamily', 'dro', 'arrParams', struct('tau', 1), 'tauDRO', 0, ...
    'thrLock', true, 'tf0', 0.30, 'resume', false, ...
    'logFile', fullfile(outDir, 'cell.log'));

% (a) bad seedTop mesh refuses before any solve
o = base;  o.cells = [1 1];  o.redo = true;
o.seedTop = struct('X', zeros(7, N), 'U', zeros(4, N+1), 'tf', 0.3);
try
    thrust_ladder_library(outMat, o);
    error('test:noThrow', 'bad seedTop mesh did not refuse');
catch err
    assert(strcmp(err.identifier, 'thrust_ladder_library:seedTop'), ...
        'wrong error for a bad seedTop: %s', err.identifier);
end
assert(~isfile(outMat), 'a refused call must not write the sheet');

% (b) redo without cells refuses
o = base;  o.redo = true;
try
    thrust_ladder_library(outMat, o);
    error('test:noThrow', 'redo without cells did not refuse');
catch err
    assert(strcmp(err.identifier, 'thrust_ladder_library:redoCells'), ...
        'wrong error for redo without cells: %s', err.identifier);
end

% (c) the normal ladder: 15 N cold start, 12 N warm from it
P = thrust_ladder_library(outMat, base);
assert(all(P.OK(1, 1, :)), 'both rungs must be OK on the demo cell');
tf15 = P.TF(1, 1, 1);  tf12 = P.TF(1, 1, 2);
fprintf('  cold ladder: tf(15 N) = %.5f, tf(12 N) = %.5f ND\n', tf15, tf12);
assert(abs(tf15 - 0.3040) < 0.005, '15 N cold start expected near 0.3040 ND, got %.5f', tf15);
assert(abs(tf12 - 0.3070) < 0.005, '12 N rung expected near 0.3070 ND, got %.5f', tf12);

% (d) redo + seedTop from the 12 N entry: the reverse warm start
ndp  = nd_propulsion(12, 1710, 150, P.meta.lStar, P.meta.tStar);
[tD, rvD] = get_family_orbit('gto', struct('orientDeg', 0));
rv0  = phase_state(tD, rvD);  rv0 = rv0(0.10);
z12  = P.Z8(:, 1, 1, 2);
[tt, yy] = pumpkyn.cr3bp.tfMinProp(z12(8), [rv0(:); 1; z12(1:7)], ndp.Tnd, ndp.cnd, P.meta.muStar);
tg   = linspace(0, z12(8), N+1);
Y    = interp1(tt, yy, tg, 'spline');                         % [N+1 x 14]
lamV = Y(:, 11:13).';
U    = [-lamV ./ max(vecnorm(lamV, 2, 1), eps); ones(1, N+1)];
seedTop = struct('X', Y(:, 1:7).', 'U', U, 'tf', z12(8)*(12/15)^0.6);
o = base;  o.resume = true;  o.cells = [1 1];  o.redo = true;  o.seedTop = seedTop;
P2 = thrust_ladder_library(outMat, o);
tf15b = P2.TF(1, 1, 1);
fprintf('  redo from the 12 N seed: tf(15 N) = %.5f ND, OK = %d, ATT = %d\n', ...
    tf15b, P2.OK(1, 1, 1), P2.ATT(1, 1));
assert(P2.OK(1, 1, 1), 'the re-solved 15 N rung must be OK');
assert(tf15b < 0.29, 'the reverse warm start must land on the faster basin (< 0.29 ND), got %.5f', tf15b);
assert(P2.ATT(1, 1) == 1, 'redo must reset the attempt counter (got %d)', P2.ATT(1, 1));
Q = load(outMat, 'TF', 'OK');
assert(abs(Q.TF(1, 1, 1) - tf15b) < 1e-12 && Q.OK(1, 1, 1), 'the sheet on disk must carry the re-solve');

rmdir(outDir, 's');
ok = true;
fprintf('test_ladder_seedtop: PASS\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Through the MCP tool `run_matlab_file` on `orbit_transfer/costate_common/tests/test_ladder_seedtop.m`.
Expected: check (a) fails with `test:noThrow` (the engine accepts the bad struct because it does not read `seedTop` yet).

- [ ] **Step 3: Implement the options**

In the header's Inputs block, after the `.maxCpuSec` entry, add:

```matlab
%              .seedTop   struct .X [7 x N+1] [r;v;m], .U [4 x N+1]
%                         [alpha;thr], .tf (ND): an EXTERNAL warm start
%                         for the TOP rung of every cell this call
%                         solves (a campaign repair tool's reverse warm
%                         start or phase continuation); [] = the cold
%                         tf0 path, byte-identical to before    [[]]
%              .redo      true = re-solve every cell in .cells even if
%                         it already holds OK rungs: its ATT is reset to
%                         0 and its TF/FLYKM/ACCDZ/RES/WALL/OK/Z8 slices
%                         are cleared before the ladder runs. Refuses
%                         without .cells (thrust_ladder_library:redoCells).
%                         The caller owns the backup of what it clears
%                                                                 [false]
```

Replace the todo block (`todo = d('cells', []);` through the resume `end`) with:

```matlab
todo = d('cells', []);
redo = d('redo', false);
if redo && isempty(todo)
    error('thrust_ladder_library:redoCells', ...
          'opts.redo requires opts.cells (the cells to re-solve)');
end
if isempty(todo)
    [gA, gD] = meshgrid(1:nA, 1:nD);
    todo = [gD(:), gA(:)];
end
seedTop = d('seedTop', []);
if ~isempty(seedTop)
    okSeed = isstruct(seedTop) && all(isfield(seedTop, {'X', 'U', 'tf'})) && ...
             isequal(size(seedTop.X), [7 N+1]) && isequal(size(seedTop.U), [4 N+1]) && ...
             isscalar(seedTop.tf) && seedTop.tf > 0;
    if ~okSeed
        error('thrust_ladder_library:seedTop', ...
              'opts.seedTop must carry .X [7 x %d], .U [4 x %d] and a positive .tf', N+1, N+1);
    end
end

nR = numel(rungs);
TF = nan(nD,nA,nR); FLYKM = TF; ACCDZ = TF; RES = TF; WALL = TF;
OK = false(nD,nA,nR);  Z8 = nan(8,nD,nA,nR);
% ATT counts attempts per cell, written BEFORE the cell is solved so that a
% cell which HANGS (and is killed externally) still leaves a record. Without
% it a hung cell has no verified rung, resume would not skip it, and a
% batched driver would retry it forever.
ATT = zeros(nD,nA);
maxAtt = d('maxAtt', 2);
% RESUME: a prior run's products are data -- reload them and skip cells that
% already produced verified rungs (opts.resume, default true). A redo call
% keeps every cell it names, resets its attempts and clears its slices.
if d('resume', true) && isfile(outMat)
    Q = load(outMat);
    if isequal(size(Q.OK), [nD nA nR])
        TF=Q.TF; FLYKM=Q.FLYKM; ACCDZ=Q.ACCDZ; RES=Q.RES; WALL=Q.WALL;
        OK=Q.OK; Z8=Q.Z8;
        if isfield(Q,'ATT'), ATT = Q.ATT; end
        if redo
            for kq = 1:size(todo,1)
                iq = todo(kq,1);  jq = todo(kq,2);
                ATT(iq,jq) = 0;
                TF(iq,jq,:) = NaN;  FLYKM(iq,jq,:) = NaN;  ACCDZ(iq,jq,:) = NaN;
                RES(iq,jq,:) = NaN; WALL(iq,jq,:) = NaN;   OK(iq,jq,:) = false;
                Z8(:,iq,jq,:) = NaN;
            end
            nSkip = 0;
        else
            keep = false(size(todo,1),1);
            for kq = 1:size(todo,1)
                iq = todo(kq,1);  jq = todo(kq,2);
                keep(kq) = ~any(OK(iq,jq,:)) && ATT(iq,jq) < maxAtt;
            end
            nSkip = nnz(~keep);  todo = todo(keep,:);
        end
    else
        nSkip = 0;
    end
else
    nSkip = 0;
end
```

Note the original `nR = ...` through `maxAtt = ...` lines move INTO this block unchanged (they are repeated above so the block is complete); delete the originals so nothing is defined twice.

Replace the per-cell seed line `seedX = [];  seedU = [];  seedTf = tf0Top;  Tprev = [];` with:

```matlab
    seedX = [];  seedU = [];  seedTf = tf0Top;  Tprev = [];
    if ~isempty(seedTop)
        seedX = seedTop.X;  seedU = seedTop.U;  seedTf = seedTop.tf;
    end
```

Add to the THRUST LADDER opening log line's arguments nothing; add one line after it:

```matlab
if ~isempty(seedTop), lg('  top rung seeded externally (opts.seedTop, tf %.5f)', seedTop.tf); end
if redo, lg('  redo: %d cell(s) cleared and re-solved', size(todo,1)); end
```

- [ ] **Step 4: Run the test to verify it passes**

`run_matlab_file` on the test. Expected: four checks, `test_ladder_seedtop: PASS`, about 6 to 10 minutes (three top-rung solves and two 12 N solves). Record the printed t_f values in the task report.

- [ ] **Step 5: Run the existing engine regressions**

`run_matlab_file` on `orbit_transfer/costate_common/tests/test_mintime_earth_floor.m` and `orbit_transfer/DRO_tulip/direct/certify/tests/test_minenergy_objective.m` (the bitwise default-path regression). Expected: both PASS unchanged. Run `orbit_transfer/costate_common/tests/test_folder_rules.m` and confirm the only failure is the pre-existing `%#ok` in `test_certify_schema.m` (not ours).

- [ ] **Step 6: README paragraph**

In `orbit_transfer/costate_common/README.md`, under the `thrust_ladder_library` entry, add:

```
`opts.seedTop` (struct .X [7 x N+1], .U [4 x N+1], .tf) seeds the TOP rung of
every cell the call solves from an external direct solution (a repair tool's
reverse warm start or phase continuation); `opts.redo` (with `opts.cells`)
re-solves cells that already hold OK rungs, resetting their attempts and
clearing their slices first -- the caller keeps the backup. Both default off
and leave every other call byte-identical (tests/test_ladder_seedtop).
```

- [ ] **Step 7: Commit** (controller)

```bash
git add orbit_transfer/costate_common/thrust_ladder_library.m orbit_transfer/costate_common/tests/test_ladder_seedtop.m orbit_transfer/costate_common/README.md
git commit -m "costate_common: thrust_ladder_library opts.seedTop (external top-rung warm start) and opts.redo (re-solve named cells, attempts reset, slices cleared) -- default off, byte-identical otherwise; test on the GTO->DRO demo cell measures the reverse warm start landing on the faster 15 N basin"
```

---

### Task 2: `seed_top_from_entry` and `earth_clearance_sweep` (+ the cached demo-sheet fixture)

**Files:**
- Create: `orbit_transfer/GTO_DRO/lib/seed_top_from_entry.m`
- Create: `orbit_transfer/GTO_DRO/lib/earth_clearance_sweep.m`
- Create: `orbit_transfer/GTO_DRO/tests/fixture_demo_sheet.m`
- Test: `orbit_transfer/GTO_DRO/tests/test_seed_top_from_entry.m`, `orbit_transfer/GTO_DRO/tests/test_earth_clearance_sweep.m`
- Modify: `orbit_transfer/GTO_DRO/tests/run_tests.m` (add the two names)

**Interfaces:**
- Consumes: Task 1's `seedTop` contract; sheet `.mat` variables `TF, FLYKM, ACCDZ, RES, WALL, OK, Z8, ATT, rungs, sD, sA, meta` (meta has `depFamily, depParams, arrFamily, arrParams, ispS, m0kg, N, floorEarthKm, muStar, lStar, tStar`); `costate_common/phase_state`, `nd_propulsion`, `fly_transfer`, `get_family_orbit`; `GTO_DRO/lib/earth_clearance`.
- Produces:
  - `seedTop = seed_top_from_entry(sheetMat, iD, iA, kr, thrustTargetN)`: flies entry `(iD, iA, kr)`'s `z8` with `pumpkyn.cr3bp.tfMinProp`, resamples on `meta.N + 1` uniform nodes, returns `struct('X' [7 x N+1], 'U' [4 x N+1] with thr = 1, 'tf' = z8(8)*(rungs(kr)/thrustTargetN)^0.6, 'from' struct(iD, iA, kr, thrustN, tf_nd))`. Refuses a cell that is not OK (`seed_top_from_entry:notOK`).
  - `S = earth_clearance_sweep(sheetMat, opts)`: for every OK entry, `fly_transfer` then `earth_clearance` against `opts.floorKm` (default `meta.floorEarthKm`); writes sidecar `<sheet>_earth.mat` with `EARTHKM [nD x nA x nR]` (NaN where not OK), `floorKm`, `done [nD x nA x nR]`, resumable (`opts.batchSec`, `opts.maxEntries`). With `opts.writeback = true` and a complete census, copies the sheet to `<sheet>.bak_earth` (only if no backup exists yet), sets `OK(...) = false` for entries below the floor, appends `meta.earthSweep = struct('date', 'floorKm', 'nDemoted', 'demoted' [k x 3])` and saves. Returns `S` with `.done, .nDone, .nTodo, .nBelow, .minAltKm, .sideMat, .demoted`.
  - `fixture_demo_sheet()` returns the path of `GTO_DRO/results/test_fixture/demo_sheet.mat`, a 1 x 1 sheet (demo cell, rungs `[15 12]`, floor 200 km) built through `thrust_ladder_library` on first use and cached; `fixture_demo_sheet(true)` forces a rebuild.

- [ ] **Step 1: Write the fixture**

```matlab
function sheetMat = fixture_demo_sheet(rebuild)
% FIXTURE_DEMO_SHEET  A cached 1 x 1 GTO -> DRO sheet for the stage 2 tests.
%
%   The stage 1 demo cell (orient 0, DRO tau 1, departure phase 0.10,
%   arrival phase 0, Earth floor 200 km) at rungs [15 12], solved through
%   thrust_ladder_library exactly as the fleet would solve it, saved under
%   results/test_fixture/ (git-ignored) and reused by every test that
%   needs a real entry. The 15 N cold start lands on the 0.3040 ND basin
%   (refuted by the conjugate test, stage 1) and the 12 N rung on 0.3070 ND.
%
% INPUTS:
%   rebuild  - true forces a fresh solve (default false)          [scalar]
%
% OUTPUTS:
%   sheetMat - path of the fixture sheet .mat                     [char]
%
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.

if nargin < 1, rebuild = false; end
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
fixDir   = fullfile(fileparts(here), 'results', 'test_fixture');
if ~isfolder(fixDir), mkdir(fixDir); end
sheetMat = fullfile(fixDir, 'demo_sheet.mat');
if isfile(sheetMat) && ~rebuild
    Q = load(sheetMat, 'OK', 'rungs');
    if isequal(Q.rungs, [15 12]) && all(Q.OK(:)), return, end
end
if isfile(sheetMat), delete(sheetMat); end
opts = struct('rungs', [15 12], 'ispS', 1710, 'm0kg', 150, 'N', 400, ...
    'floorKm', 500, 'floorEarthKm', 200, 'maxIter', 6000, 'maxCpuSec', 600, ...
    'gateKm', 100, 'accTol', 1e-6, 'nD', 1, 'nA', 1, 'sD0', 0.10, 'sA0', 0, ...
    'depFamily', 'gto', 'depParams', struct('orientDeg', 0), ...
    'arrFamily', 'dro', 'arrParams', struct('tau', 1), 'tauDRO', 0, ...
    'thrLock', true, 'tf0', 0.30, 'resume', false, ...
    'logFile', fullfile(fixDir, 'demo_sheet.log'));
P = thrust_ladder_library(sheetMat, opts);
assert(all(P.OK(:)), 'fixture_demo_sheet: the demo cell did not solve at both rungs');
end
```

- [ ] **Step 2: Write the failing tests**

```matlab
function ok = test_seed_top_from_entry()
% TEST_SEED_TOP_FROM_ENTRY  Flown 12 N entry -> engine seedTop struct.
% INPUTS: none.  OUTPUTS: ok - logical.
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(here);
sheetMat = fixture_demo_sheet();
Q = load(sheetMat, 'Z8', 'rungs', 'meta', 'OK');
N = Q.meta.N;

s = seed_top_from_entry(sheetMat, 1, 1, 2, 15);
assert(isequal(size(s.X), [7 N+1]) && isequal(size(s.U), [4 N+1]), 'seedTop shape');
assert(all(s.U(4, :) == 1), 'thrust-locked seed must carry thr = 1');
assert(all(abs(vecnorm(s.U(1:3, :), 2, 1) - 1) < 1e-12), 'direction columns must be unit');
assert(abs(s.X(7, 1) - 1) < 1e-12, 'mass fraction starts at 1');
assert(s.X(7, end) < 1 && s.X(7, end) > 0.3, 'mass fraction ends in (0.3, 1)');
tf12 = Q.Z8(8, 1, 1, 2);
assert(abs(s.tf - tf12*(12/15)^0.6) < 1e-12, 'tf scaled by the engine rule (Tprev/TN)^0.6');
assert(s.from.iD == 1 && s.from.iA == 1 && s.from.kr == 2 && s.from.thrustN == 12, 'provenance');
% the resampled departure state is the sheet's own departure state
[tD, rvD] = get_family_orbit(Q.meta.depFamily, Q.meta.depParams);
dep = phase_state(tD, rvD);  rv0 = dep(0.10);
assert(max(abs(s.X(1:6, 1) - rv0(:))) < 1e-9, 'seed starts at the cell departure state');
% a non-OK cell refuses
Qbad = Q;  Qbad.OK(1, 1, 2) = false;
badMat = fullfile(tempdir, 'seed_top_bad.mat');
S = load(sheetMat);  S.OK(1, 1, 2) = false;  save(badMat, '-struct', 'S');
try
    seed_top_from_entry(badMat, 1, 1, 2, 15);
    error('test:noThrow', 'a non-OK entry must refuse');
catch err
    assert(strcmp(err.identifier, 'seed_top_from_entry:notOK'), 'wrong error: %s', err.identifier);
end
delete(badMat);
ok = true;
fprintf('test_seed_top_from_entry: PASS\n');
end
```

```matlab
function ok = test_earth_clearance_sweep()
% TEST_EARTH_CLEARANCE_SWEEP  Sidecar census, writeback and the demotion path.
% INPUTS: none.  OUTPUTS: ok - logical.
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(here);
src = fixture_demo_sheet();
wd  = fullfile(tempdir, 'test_earth_sweep');
if isfolder(wd), rmdir(wd, 's'); end
mkdir(wd);
sheetMat = fullfile(wd, 'demo_sheet.mat');
copyfile(src, sheetMat);

% census at the campaign floor: both entries clear it (the engine solved
% with the floor and the flown arc is within 1e-4 km of the direct one)
S = earth_clearance_sweep(sheetMat);
assert(S.done && S.nDone == 2 && S.nTodo == 0, 'census must cover both OK entries');
assert(isfile(S.sideMat), 'sidecar written');
E = load(S.sideMat, 'EARTHKM', 'floorKm');
assert(E.floorKm == 200, 'floor defaults to meta.floorEarthKm');
assert(all(E.EARTHKM(1, 1, :) > 200), 'both flown arcs clear 200 km');
assert(S.nBelow == 0, 'nothing below the floor');
fprintf('  min Earth altitude over the two entries: %.0f km\n', S.minAltKm);

% resume: a second call does no flights (done mask) and reports the same
S2 = earth_clearance_sweep(sheetMat);
assert(S2.nDone == 2 && abs(S2.minAltKm - S.minAltKm) < 1e-9, 'resume reproduces the census');

% writeback at the campaign floor: nothing demoted, no backup needed, meta stamped
S3 = earth_clearance_sweep(sheetMat, struct('writeback', true));
Q = load(sheetMat, 'OK', 'meta');
assert(all(Q.OK(:)) && isfield(Q.meta, 'earthSweep') && Q.meta.earthSweep.nDemoted == 0, 'clean writeback');

% impossible floor: both demoted, backup kept, demoted list names them
delete(S.sideMat);
S4 = earth_clearance_sweep(sheetMat, struct('floorKm', 1e6, 'writeback', true));
assert(S4.nBelow == 2, 'both entries below an impossible floor');
Q = load(sheetMat, 'OK', 'meta');
assert(~any(Q.OK(:)), 'demoted entries are OK = false');
assert(isfile([sheetMat(1:end-4) '.bak_earth']), 'pre-sweep backup kept');
assert(isequal(sortrows(Q.meta.earthSweep.demoted), [1 1 1; 1 1 2]), 'demoted list');
B = load([sheetMat(1:end-4) '.bak_earth'], 'OK');
assert(all(B.OK(:)), 'the backup holds the original OK');

rmdir(wd, 's');
ok = true;
fprintf('test_earth_clearance_sweep: PASS\n');
end
```

- [ ] **Step 3: Run both to verify they fail**

`run_matlab_file` on each. Expected: `Unrecognized function` for `seed_top_from_entry` / `earth_clearance_sweep` (the fixture itself builds on the first call, 4 to 6 minutes, and is then cached).

- [ ] **Step 4: Implement `seed_top_from_entry`**

```matlab
function seedTop = seed_top_from_entry(sheetMat, iD, iA, kr, thrustTargetN)
% SEED_TOP_FROM_ENTRY  Fly a sheet entry and resample it as an engine seed.
%
%   Flies entry (iD, iA, kr) of a thrust_ladder_library sheet with
%   pumpkyn.cr3bp.tfMinProp, resamples the flight on the engine's N + 1
%   uniform nodes and returns the struct thrust_ladder_library's
%   opts.seedTop expects. The seed is thrust-locked (thr = 1) and its t_f is
%   scaled to the target rung by the engine's own rule (Tprev/TN)^0.6, so a
%   12 N entry seeds a 15 N re-solve (reverse warm start, basin repair) and
%   a neighbour cell's 15 N entry seeds a 15 N solve (phase continuation).
%
% INPUTS:
%   sheetMat      - sheet .mat written by thrust_ladder_library     [char]
%   iD, iA, kr    - departure phase, arrival phase, rung indices    [scalars]
%   thrustTargetN - rung the seed will be used at, N                [scalar]
%
% OUTPUTS:
%   seedTop       - .X [7 x N+1] [r; v; m], .U [4 x N+1] [alpha; thr],
%                   .tf scalar ND, .from struct(iD, iA, kr, thrustN, tf_nd)
%
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(here, 'tests'));
    sheetMat = fixture_demo_sheet();  iD = 1;  iA = 1;  kr = 2;  thrustTargetN = 15;
end
Q = load(sheetMat, 'Z8', 'OK', 'rungs', 'sD', 'meta');
if ~Q.OK(iD, iA, kr)
    error('seed_top_from_entry:notOK', 'entry (%d,%d,rung %g N) of %s is not OK', ...
          iD, iA, Q.rungs(kr), sheetMat);
end
meta   = Q.meta;
z8     = Q.Z8(:, iD, iA, kr);
TN     = Q.rungs(kr);
ndp    = nd_propulsion(TN, meta.ispS, meta.m0kg, meta.lStar, meta.tStar);
[tD, rvD] = get_family_orbit(meta.depFamily, meta.depParams);
dep    = phase_state(tD, rvD);
rv0    = dep(Q.sD(iD));
[tt, yy] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(:); 1; z8(1:7)], ndp.Tnd, ndp.cnd, meta.muStar);
N      = meta.N;
tg     = linspace(0, z8(8), N+1);
Y      = interp1(tt, yy, tg, 'spline');                       % [N+1 x 14]
lamV   = Y(:, 11:13).';
alpha  = -lamV ./ max(vecnorm(lamV, 2, 1), eps);
seedTop = struct('X', Y(:, 1:7).', 'U', [alpha; ones(1, N+1)], ...
                 'tf', z8(8)*(TN/thrustTargetN)^0.6, ...
                 'from', struct('iD', iD, 'iA', iA, 'kr', kr, 'thrustN', TN, 'tf_nd', z8(8)));
if nargin == 0
    fprintf('seed_top_from_entry: %g N entry tf %.5f -> seed tf %.5f for %g N, m_f/m0 %.4f\n', ...
            TN, z8(8), seedTop.tf, thrustTargetN, seedTop.X(7, end));
end
end
```

- [ ] **Step 5: Implement `earth_clearance_sweep`**

```matlab
function S = earth_clearance_sweep(sheetMat, opts)
% EARTH_CLEARANCE_SWEEP  Flown Earth altitude of every OK entry in a sheet.
%
%   The direct solve enforced the Earth floor at its nodes and midpoints,
%   and the engine's G1b gate flew each accepted z8 for the arrival miss,
%   but nothing recorded the flown arc's minimum Earth altitude. This sweep
%   closes that gap on the sheet (campaign side, no engine change): every
%   OK entry is flown with fly_transfer, earth_clearance measures its
%   minimum altitude, the result goes to a resumable sidecar, and an
%   optional writeback demotes entries below the floor to OK = false after
%   backing the sheet up once. Audit and packaging read the sidecar.
%
% INPUTS:
%   sheetMat - sheet .mat written by thrust_ladder_library        [char]
%   opts     - (optional) .floorKm [meta.floorEarthKm, or 200 if the sheet
%              has none], .sideMat [<sheet>_earth.mat], .batchSec [inf],
%              .maxEntries [inf], .writeback [false], .logFile ['']
%
% OUTPUTS:
%   S        - .done, .nDone, .nTodo, .nBelow, .minAltKm, .sideMat,
%              .demoted [k x 3] (iD, iA, kr) when written back, .floorKm
%
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(here, 'tests'));
    sheetMat = fixture_demo_sheet();  opts = struct();
end
if nargin < 2, opts = struct(); end
Q = load(sheetMat);
meta = Q.meta;
floorDef = 200;
if isfield(meta, 'floorEarthKm') && ~isempty(meta.floorEarthKm), floorDef = meta.floorEarthKm; end
floorKm    = fdef(opts, 'floorKm', floorDef);
sideMat    = fdef(opts, 'sideMat', [sheetMat(1:end-4) '_earth.mat']);
batchSec   = fdef(opts, 'batchSec', inf);
maxEntries = fdef(opts, 'maxEntries', inf);
writeback  = fdef(opts, 'writeback', false);
logFile    = fdef(opts, 'logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

[nD, nA, nR] = size(Q.OK);
EARTHKM = nan(nD, nA, nR);  done = false(nD, nA, nR);
if isfile(sideMat)
    E = load(sideMat);
    if isequal(size(E.EARTHKM), [nD nA nR]) && E.floorKm == floorKm
        EARTHKM = E.EARTHKM;  done = E.done;
    end
end

[tD, rvD] = get_family_orbit(meta.depFamily, meta.depParams);
[tA, rvA] = get_family_orbit(meta.arrFamily, meta.arrParams);
dep = phase_state(tD, rvD);  arr = phase_state(tA, rvA);
pEarth = struct('muStar', meta.muStar, 'lStar', meta.lStar);

tAll = tic;  nThis = 0;
for iD = 1:nD
    for iA = 1:nA
        for kr = 1:nR
            if ~Q.OK(iD, iA, kr) || done(iD, iA, kr), continue, end
            if toc(tAll) > batchSec || nThis >= maxEntries, break, end
            TN  = Q.rungs(kr);
            ndp = nd_propulsion(TN, meta.ispS, meta.m0kg, meta.lStar, meta.tStar);
            phys = struct('Tnd', ndp.Tnd, 'cnd', ndp.cnd, 'muStar', meta.muStar, ...
                          'lStar', meta.lStar, 'tStar', meta.tStar, 'm0kg', meta.m0kg);
            rv0 = dep(Q.sD(iD));  rvf = arr(Q.sA(iA));
            F = fly_transfer(Q.Z8(:, iD, iA, kr), rv0(:), rvf(:), phys);
            hMin = earth_clearance(F.t(:), F.Y, floorKm, pEarth);
            EARTHKM(iD, iA, kr) = hMin;  done(iD, iA, kr) = true;
            nThis = nThis + 1;
            lg('  (%2d,%2d) T=%5.1f  min Earth alt %.1f km %s', iD, iA, TN, hMin, pick(hMin >= floorKm, '', '< FLOOR'));
            save(sideMat, 'EARTHKM', 'done', 'floorKm');
        end
    end
end

okMask  = Q.OK;
todo    = okMask & ~done;
below   = okMask & done & (EARTHKM < floorKm);
S = struct('done', ~any(todo(:)), 'nDone', nnz(okMask & done), 'nTodo', nnz(todo), ...
           'nBelow', nnz(below), 'minAltKm', min(EARTHKM(okMask & done), [], 'all'), ...
           'sideMat', sideMat, 'demoted', zeros(0, 3), 'floorKm', floorKm);
if isempty(S.minAltKm), S.minAltKm = NaN; end
lg('earth sweep: %d done, %d to do, %d below %.0f km, min %.1f km', S.nDone, S.nTodo, S.nBelow, floorKm, S.minAltKm);

if writeback
    assert(S.done, 'earth_clearance_sweep:incomplete', 'writeback needs a complete census (%d to do)', S.nTodo);
    bak = [sheetMat(1:end-4) '.bak_earth'];
    if ~isfile(bak), copyfile(sheetMat, bak); end
    [dI, dJ, dK] = ind2sub(size(below), find(below));
    S.demoted = [dI, dJ, dK];
    Q.OK(below) = false;
    Q.meta.earthSweep = struct('date', datestr(now, 'yyyy-mm-dd HH:MM'), 'floorKm', floorKm, ...
                               'nDemoted', nnz(below), 'demoted', S.demoted);
    save(sheetMat, '-struct', 'Q');
    lg('earth sweep writeback: %d entr%s demoted, backup %s', nnz(below), pick(nnz(below) == 1, 'y', 'ies'), bak);
end
if nargin == 0, disp(S), end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end

function v = pick(c, a, b)
% PICK  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end

function logmsg(f, s)
% LOGMSG  Append-mode file log, stdout when f is empty.  INPUTS: f; s.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); end
end
```

`fly_transfer` returns `.t` and `.Y` (rows = samples, columns 1:3 position); `earth_clearance` takes `(tau [M x 1], y [M x >= 6], floorKm, p)`.

- [ ] **Step 6: Run both tests to verify they pass**

Expected: `test_seed_top_from_entry: PASS` (seconds after the fixture exists), `test_earth_clearance_sweep: PASS` (four flights, under a minute). Add `'test_seed_top_from_entry', 'test_earth_clearance_sweep'` to the `names` cell in `tests/run_tests.m`.

- [ ] **Step 7: Commit** (controller)

```bash
git add orbit_transfer/GTO_DRO/lib/seed_top_from_entry.m orbit_transfer/GTO_DRO/lib/earth_clearance_sweep.m orbit_transfer/GTO_DRO/tests/fixture_demo_sheet.m orbit_transfer/GTO_DRO/tests/test_seed_top_from_entry.m orbit_transfer/GTO_DRO/tests/test_earth_clearance_sweep.m orbit_transfer/GTO_DRO/tests/run_tests.m
git commit -m "GTO_DRO: seed_top_from_entry (flown sheet entry -> engine seedTop, tf scaled by the engine rule) and earth_clearance_sweep (flown Earth altitude per OK entry, resumable sidecar, writeback demotes below the floor with a one-time backup); cached demo-sheet fixture for the stage 2 tests"
```

---

### Task 3: The fleet front door `run_gto_dro_catalog` and its batched driver

**Files:**
- Create: `orbit_transfer/GTO_DRO/catalog/run_gto_dro_catalog.m`
- Create: `orbit_transfer/GTO_DRO/catalog/run_gto_dro_batched.sh` (executable)
- Test: `orbit_transfer/GTO_DRO/tests/test_run_gto_dro_catalog.m`
- Modify: `orbit_transfer/GTO_DRO/tests/run_tests.m`

**Interfaces:**
- Consumes: `thrust_ladder_library(outMat, opts)` as stage 1's `gto_dro_entry` calls it (same option names), `setup_gto_dro_paths` (one level up).
- Produces: `plan = run_gto_dro_catalog(sheetSel, maxCellsIn, batchSecIn, rungsIn, optsIn)`; `optsIn.planOnly = true` returns, without solving, a struct array `plan(k)` with `.tag, .orientDeg, .tau, .outMat, .logFile, .optsA, .optsB` (the exact engine options of pass A and pass B) for the selected sheets; `optsIn.tf0` overrides the top-rung seed (default 0.30; the pilot re-bases it); `optsIn.catDir` overrides the sheet folder (tests). Sheet tags `gto_dro_o%03d_tau%g`; sheet index order: orientation outer `[0 90 180 270]`, tau inner `[0.5 1 2 3]`, so sheet 2 is `(0 deg, tau 1)`, the pilot. Completion markers `CATALOG ALL SHEETS COMPLETE` / `SELECTED SHEETS COMPLETE (...)` in `combined_progress.txt`.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_run_gto_dro_catalog()
% TEST_RUN_GTO_DRO_CATALOG  planOnly contract and the on-disk rung guard.
% INPUTS: none.  OUTPUTS: ok - logical.
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(fullfile(fileparts(here), 'catalog'));

plan = run_gto_dro_catalog([], [], [], [], struct('planOnly', true));
assert(numel(plan) == 16, '16 sheets');
tags = {plan.tag};
assert(numel(unique(tags)) == 16, 'distinct tags');
assert(strcmp(plan(2).tag, 'gto_dro_o000_tau1') && plan(2).orientDeg == 0 && plan(2).tau == 1, 'sheet 2 is the pilot');
assert(strcmp(plan(5).tag, 'gto_dro_o090_tau0.5'), 'tau 0.5 tag spelling');
for k = 1:16
    oA = plan(k).optsA;
    assert(strcmp(oA.depFamily, 'gto') && strcmp(oA.arrFamily, 'dro'), 'families');
    assert(oA.depParams.orientDeg == plan(k).orientDeg, 'orientation');
    assert(oA.arrParams.tau == plan(k).tau, 'DRO period in arrParams');
    assert(oA.tauDRO == plan(k).orientDeg, 'sheet key rule: tauDRO = orientDeg');
    assert(oA.floorEarthKm == 200 && oA.floorKm == 500, 'floors');
    assert(oA.thrLock && abs(oA.tf0 - 0.30) < 1e-12, 'warm recipe');
    assert(oA.ispS == 1710 && oA.m0kg == 150 && oA.N == 400 && oA.maxIter == 6000 && oA.maxCpuSec == 600, 'regime + solver');
    assert(oA.nD == 12 && oA.nA == 6 && oA.sD0 == 0 && oA.sA0 == 0, 'grid');
    assert(isequal(oA.rungs, [15 12 10 7 5]), 'default rungs');
    assert(~isfield(plan(k).optsB, 'tf0') || isempty(plan(k).optsB.tf0), 'pass B is the cold tf0');
    assert(isempty(plan(k).optsB.cells), 'pass B cells are scoped at run time, empty in the plan');
end
plan9 = run_gto_dro_catalog(2, [], [], [15 12 10 7 5 3 2 1.5 1], struct('planOnly', true, 'tf0', 0.27));
assert(numel(plan9) == 1 && isequal(plan9.optsA.rungs, [15 12 10 7 5 3 2 1.5 1]) && abs(plan9.optsA.tf0 - 0.27) < 1e-12, 'overrides');

% on-disk rung guard: a fake sheet with rungs [15 12] refuses a 5-rung call
wd = fullfile(tempdir, 'test_run_gto_dro_catalog');
if isfolder(wd), rmdir(wd, 's'); end
mkdir(wd);
rungs = [15 12];
save(fullfile(wd, 'gto_dro_o000_tau1.mat'), 'rungs');
try
    run_gto_dro_catalog(2, 0, 1, [], struct('catDir', wd));
    error('test:noThrow', 'rung mismatch did not refuse');
catch err
    assert(strcmp(err.identifier, 'run_gto_dro_catalog:rungMismatch'), 'wrong error: %s', err.identifier);
end
% a matching call with maxCells = 0 runs the engine on zero cells and returns
run_gto_dro_catalog(2, 0, 1, [15 12], struct('catDir', wd));
assert(isfile(fullfile(wd, 'combined_progress.txt')), 'driver log written');
rmdir(wd, 's');
ok = true;
fprintf('test_run_gto_dro_catalog: PASS\n');
end
```

- [ ] **Step 2: Run it to verify it fails**

Expected: `Unrecognized function 'run_gto_dro_catalog'`.

- [ ] **Step 3: Write the front door**

Clone `GTO_tulip/catalog/run_gto_catalog.m` into `GTO_DRO/catalog/run_gto_dro_catalog.m` and make these changes (everything else, including the two-pass multistart, the fail-loud try/catch blocks, the open-cell sweep and the completion markers, stays verbatim):

```matlab
function plan = run_gto_dro_catalog(sheetSel, maxCellsIn, batchSecIn, rungsIn, optsIn)
%% Purpose:
%
%   FRONT DOOR for the GTO -> DRO costate catalog (stage 2 of the GTO -> DRO
%   min-time campaign, spec docs/superpowers/specs/2026-11-09-gto-dro-
%   mintime-design.md section 10): the shared three-gate ladder engine
%   (direct solve -> ms_tfmin refinement -> pumpkyn tfMin acceptance)
%   walked down the thrust ladder over (GTO departure orientation x DRO
%   period) sheets. Clone of GTO_tulip/catalog/run_gto_catalog.m with the
%   DRO arrival, the Earth-altitude floor (stage 1 measured a floored-free
%   solve passing 1722 km from Earth's centre) and a planOnly mode. The
%   two-point top-rung multistart (pass A tf0 warm, pass B cold mop-up,
%   first-success-wins) and the SHEET KEY RULE (opts.tauDRO = orientDeg)
%   are inherited unchanged; see that file's header for their rationale.
%
%   WHAT THE FLEET DOES NOT DO (stage 1 findings, GTO_DRO/README.md): it
%   does not select the top-rung basin (the floored cold start can land on
%   an extremal the conjugate test refutes) and it does not reach the
%   perigee column at 15 N. Both are repaired AFTER the fleet by
%   lib/repair_gto_dro_cells (reverse warm start, phase continuation)
%   through the engine's seedTop/redo options, on the sheets this front
%   door writes.
%
%% Inputs:
%
%  sheetSel                 [1 x k]                 Index into the 16-sheet
%                                                   list (default 1:16);
%                                                   orientation outer
%                                                   [0 90 180 270], tau
%                                                   inner [0.5 1 2 3], so
%                                                   sheet 2 = (0 deg, tau 1)
%                                                   is the pilot
%  maxCellsIn               double                  Fresh cells per sheet
%                                                   this call [inf]
%  batchSecIn               double                  Clean-exit wall budget
%                                                   (s) shared by the
%                                                   selected sheets [inf]
%  rungsIn                  [1 x r]                 Thrust ladder override,
%                                                   N high to low
%                                                   [15 12 10 7 5]
%  optsIn                   struct (optional)       .planOnly [false]
%                                                   return the plan without
%                                                   solving; .tf0 [0.30]
%                                                   pass A top-rung seed;
%                                                   .catDir [results/
%                                                   catalog next to this
%                                                   folder's parent]
%
%% Outputs:
%
%  plan                     [k x 1] struct          .tag .orientDeg .tau
%                                                   .outMat .logFile .optsA
%                                                   .optsB (pass B .cells
%                                                   empty: scoped at run
%                                                   time); files under
%                                                   catDir: <tag>.mat +
%                                                   <tag>_progress.txt,
%                                                   combined_progress.txt
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------
```

Argument handling and the parameter block:

```matlab
if ~exist('sheetSel','var') || isempty(sheetSel), sheetSel = 1:16; end
maxCells = inf;  batchSec = inf;
if exist('maxCellsIn','var') && ~isempty(maxCellsIn), maxCells = maxCellsIn; end
if exist('batchSecIn','var') && ~isempty(batchSecIn), batchSec = batchSecIn; end
if ~exist('rungsIn','var'), rungsIn = []; end
if ~exist('optsIn','var') || isempty(optsIn), optsIn = struct(); end
planOnly = fdef(optsIn, 'planOnly', false);

%% ===================== ADJUSTABLE PARAMETERS ============================

%% Departure: GTO pseudo-family, selected by orientation (sheet axis 1):
      orients = [0 90 180 270];               % deg

%% Arrival: DRO periods (sheet axis 2; ND, tau 0.5/1/2/3 = 2.22/4.43/8.87/13.30 d):
         taus = [0.5 1 2 3];

%% Phasing grid per sheet (departure x arrival; uniform in time fraction):
           nD = 12;
           nA = 6;
          sD0 = 0;
          sA0 = 0;

%% Thrust ladder (N, high to low); rungsIn overrides, and must match any
%% sheet already on disk (guard below):
        rungs = [15 12 10 7 5];
if ~isempty(rungsIn), rungs = rungsIn; end
         ispS = 1710;                % specific impulse (s)
         m0kg = 150;                 % initial mass (kg)

%% Solver / gates (the stage 1 recipe, GTO_DRO/README.md):
            N = 400;                 % collocation intervals
      floorKm = 500;                 % lunar altitude floor (km)
 floorEarthKm = 200;                 % Earth altitude floor (km), nodes + midpoints
      maxIter = 6000;                % IPOPT iteration cap
    maxCpuSec = 600;                 % IPOPT CPU cap per rung (s)
       gateKm = 100;                 % flown-arrival gate (km)
       accTol = 1e-6;                % tfMin acceptance
      thrLock = true;                % pin thr == 1 (classical min-time form)
     tf0Fleet = fdef(optsIn, 'tf0', 0.30);   % pass A top-rung seed (ND); the pilot re-bases it

%% ======================= END ADJUSTABLE PARAMETERS ======================
```

Paths and sheet list:

```matlab
here = fileparts(mfilename('fullpath'));                 % GTO_DRO/catalog
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
catDir = fdef(optsIn, 'catDir', fullfile(fileparts(here), 'results', 'catalog'));
if ~planOnly && ~isfolder(catDir), mkdir(catDir); end

sheets = {};
for od = orients
    for ta = taus
        sheets{end+1} = struct('orientDeg', od, 'tau', ta, ...
            'tag', sprintf('gto_dro_o%03d_tau%g', od, ta));
    end
end
nSheets = numel(sheets);
assert(nSheets == 16, 'run_gto_dro_catalog:sheetCount', 'expected 16 sheets, got %d', nSheets);
```

The per-sheet `optsCommon` (replacing the tulip one):

```matlab
    optsCommon = struct('rungs',rungs, 'ispS',ispS, 'm0kg',m0kg, ...
        'N',N, 'floorKm',floorKm, 'floorEarthKm',floorEarthKm, ...
        'maxIter',maxIter, 'maxCpuSec',maxCpuSec, ...
        'gateKm',gateKm, 'nD',nD, 'nA',nA, 'sD0',sD0, 'sA0',sA0, ...
        'depFamily','gto', 'depParams',struct('orientDeg',s.orientDeg), ...
        'arrFamily','dro', 'arrParams',struct('tau',s.tau), ...
        'tauDRO',s.orientDeg, ...          % SHEET KEY RULE: orientDeg, never the DRO period
        'accTol',accTol, 'maxCells',maxCells, ...
        'thrLock',thrLock, 'logFile',logFile);
```

Insert the planOnly branch right after `optsCommon` is built and BEFORE the rung guard (so planOnly never touches disk):

```matlab
    optsA = optsCommon;  optsA.tf0 = tf0Fleet;
    optsB = optsCommon;  optsB.cells = zeros(0,2);
    plan(end+1,1) = struct('tag', s.tag, 'orientDeg', s.orientDeg, 'tau', s.tau, ...
        'outMat', outMat, 'logFile', logFile, 'optsA', optsA, 'optsB', optsB);
    if planOnly, continue, end
```

with `plan = struct('tag',{},'orientDeg',{},'tau',{},'outMat',{},'logFile',{},'optsA',{},'optsB',{});` initialised before the loop and the `fprintf('=== ...')` banner moved below the planOnly `continue`. The error identifiers become `run_gto_dro_catalog:rungMismatch`, `run_gto_dro_catalog:sheetSel`, `run_gto_dro_catalog:sDgrid`. Pass A sets `optsA.batchSec`; pass B reads the mop-up set from disk exactly as the tulip front door does. After the loop, when `planOnly`, return before the completion-marker block. Add the `fdef` subfunction at the end:

```matlab
% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
```

- [ ] **Step 4: Write the batched driver**

`orbit_transfer/GTO_DRO/catalog/run_gto_dro_batched.sh`, a copy of `run_gto_catalog`'s driver with these substitutions: header says GTO -> DRO; `MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab`; `CATDIR="$HERE/../results/catalog"`; `MLOG="$HERE/../results/gtodrocat_matlab.log"`; `LOCK="$HERE/../results/.gtodrocat_driver.lock"`; `count_lines` globs `gto_dro_*_progress.txt`; the MATLAB call is

```bash
  "$MATLAB" -batch "here=pwd; cd('$PUMPKYN'); startup(); cd(here); \
                  addpath('$HERE'); addpath('$HERE/..'); \
                  run_gto_dro_catalog($SHEETSEL, $CELLS, $BSEC$RUNGS_ARG$TF0_ARG);" \
          >> "$MLOG" 2>&1 &
```

with a fifth optional argument `TF0` (`TF0_ARG=", [], struct('tf0', $TF0)"` when set; when `RUNGS` is empty and `TF0` is set, `RUNGS_ARG=", []"` so the argument positions hold). Usage lines:

```
# Usage:   ./run_gto_dro_batched.sh [sheetSel] [cells_per_batch] [batch_seconds] [rungs] [tf0]
# Pilot (sheet 2 = orient 0, tau 1):      ./run_gto_dro_batched.sh 2 6 1800
# Deep probe (6 cells, 9 rungs; a SEPARATE sheet file is made by build_gto_dro_catalog's
#   'deep' stage, not by this driver)
# Fleet:                                  ./run_gto_dro_batched.sh "[1 3:16]" 6 1800 "" 0.30
```

`chmod +x` it. Note the batch budget: a GTO -> DRO cell at 15 N costs 100 to 520 s per rung attempt; 6 cells per 1800 s batch with the OS kill at 1940 s.

- [ ] **Step 5: Run the test**

Expected: PASS in seconds (planOnly builds no sheets; the guard case writes one fake `.mat`; the zero-cell call writes the log and returns).

- [ ] **Step 6: Commit** (controller)

```bash
git add orbit_transfer/GTO_DRO/catalog/run_gto_dro_catalog.m orbit_transfer/GTO_DRO/catalog/run_gto_dro_batched.sh orbit_transfer/GTO_DRO/tests/test_run_gto_dro_catalog.m orbit_transfer/GTO_DRO/tests/run_tests.m
git commit -m "GTO_DRO: run_gto_dro_catalog fleet front door (clone of the GTO->tulip driver with the DRO arrival, Earth floor 200 km, sheet key = orientDeg, planOnly mode, tf0 override) and run_gto_dro_batched.sh OS-timeout driver; test pins the plan contract and the on-disk rung guard"
```

---

### Task 4: `repair_gto_dro_cells` (phase continuation and basin repair through seedTop/redo)

**Files:**
- Create: `orbit_transfer/GTO_DRO/lib/repair_gto_dro_cells.m`
- Test: `orbit_transfer/GTO_DRO/tests/test_repair_gto_dro_cells.m`
- Modify: `orbit_transfer/GTO_DRO/tests/run_tests.m`

**Interfaces:**
- Consumes: Task 1 (`seedTop`, `redo`), Task 2 (`seed_top_from_entry`), the sheet file contract, `seed_from_z8` + `ms_tfmin(..., struct('conjTest', true, 'wallSec', 120))` exactly as `gto_dro_entry` runs its conjugate re-solve (`info.conj.pass`, reproduced only when `info.converged && max(abs(z - z8)) < 1e-6`).
- Produces: `R = repair_gto_dro_cells(sheetMat, mode, opts)` with `mode` in `{'unsolved', 'basin'}`:
  - `'unsolved'`: every cell with no OK rung and `ATT >= maxAtt` (the engine retired it) gets ONE re-solve seeded from the nearest cell in the departure-phase direction (same `iA`; search `iD +- 1, +- 2, ...` with wraparound, first with an OK top rung) through `thrust_ladder_library(sheetMat, engineOpts + cells + redo + seedTop)`. Nothing to back up (the cell was empty).
  - `'basin'`: every cell whose top rung is OK runs the conjugate re-solve on its top-rung entry (`opts.conjOnly = true` stops here and returns the verdicts); a cell refuted (`conj.pass == 0`, reproduced) whose second rung is OK is re-solved with `seed_top_from_entry(sheetMat, iD, iA, 2, rungs(1))`. Before the engine call the cell's slices (TF, FLYKM, ACCDZ, RES, WALL, OK, Z8, ATT) are copied; after it, if the new top rung is not OK or its `t_f` is not smaller than the old one, the copies are restored into the sheet and the attempt is recorded as `'restored'`.
  - Both modes append rows to the sidecar `<sheet>_repairs.mat` (`rows` struct array: `.mode .iD .iA .seedFrom .tfBefore .tfAfter .okAfter .outcome .conjBefore .conjAfter .date .wallS`) and never touch a cell twice in the same mode (`opts.force = true` allows it). `opts.engineOpts` is required: the sheet's own engine options minus `rungs/nD/nA/sD0/sA0/cells/redo/seedTop/tf0` (the tool fills those from the sheet), i.e. `plan(k).optsA` from Task 3 with `tf0` ignored; `opts.maxCells` caps the cells repaired this call; `opts.logFile`.
  - Returns `R` with `.mode, .nCandidates, .nRepaired, .nRestored, .nFailed, .nConjRefuted, .conjPass [nD x nA] (NaN not run / 0 / 1), .rows, .sideMat`.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_repair_gto_dro_cells()
% TEST_REPAIR_GTO_DRO_CELLS  Basin repair on the demo cell, restore on a
% failed repair, and the unsolved-mode neighbour search.
% INPUTS: none.  OUTPUTS: ok - logical.
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(here);  addpath(fullfile(fileparts(here), 'catalog'));
src = fixture_demo_sheet();
wd  = fullfile(tempdir, 'test_repair_gto_dro');
if isfolder(wd), rmdir(wd, 's'); end
mkdir(wd);
sheetMat = fullfile(wd, 'demo_sheet.mat');
copyfile(src, sheetMat);
plan = run_gto_dro_catalog(2, [], [], [15 12], struct('planOnly', true));
eo = plan.optsA;  eo.logFile = fullfile(wd, 'repair.log');

% 1. conjOnly: the fixture's 15 N entry is the refuted basin (stage 1)
R0 = repair_gto_dro_cells(sheetMat, 'basin', struct('engineOpts', eo, 'conjOnly', true));
assert(R0.conjPass(1, 1) == 0, 'the 15 N cold-start entry must be refuted (got %g)', R0.conjPass(1, 1));
assert(R0.nConjRefuted == 1 && R0.nRepaired == 0, 'conjOnly repairs nothing');
Q = load(sheetMat, 'TF', 'OK');
tfOld = Q.TF(1, 1, 1);

% 2. basin repair: reverse warm start from the 12 N entry -> faster basin
R1 = repair_gto_dro_cells(sheetMat, 'basin', struct('engineOpts', eo));
assert(R1.nRepaired == 1 && R1.nRestored == 0, 'one cell repaired');
Q = load(sheetMat, 'TF', 'OK', 'ATT');
assert(Q.OK(1, 1, 1) && Q.TF(1, 1, 1) < 0.29 && Q.TF(1, 1, 1) < tfOld, 'repaired 15 N on the faster basin');
assert(Q.OK(1, 1, 2), 'the 12 N rung re-walked and OK');
S = load(R1.sideMat, 'rows');
assert(numel(S.rows) == 1 && strcmp(S.rows(1).outcome, 'repaired') && S.rows(1).conjBefore == 0, 'sidecar row');
fprintf('  basin repair: tf 15 N %.5f -> %.5f, conj after = %g\n', tfOld, Q.TF(1, 1, 1), S.rows(1).conjAfter);

% 3. a second basin call does nothing (cell already repaired)
R2 = repair_gto_dro_cells(sheetMat, 'basin', struct('engineOpts', eo));
assert(R2.nCandidates == 0 || R2.nRepaired == 0, 'no second repair without force');

% 4. restore on failure: force a repair whose seed is garbage (tf 1e-4 ND
% cannot reach the DRO; the engine's top rung fails) -> originals restored
Qb = load(sheetMat);
R3 = repair_gto_dro_cells(sheetMat, 'basin', struct('engineOpts', eo, 'force', true, ...
        'seedOverride', struct('X', Qb.Z8(1:7, 1, 1, 1)*zeros(1, 401) + [Qb.Z8(1:6,1,1,2); 1]*ones(1, 401), ...
                               'U', [zeros(3, 401); ones(1, 401)], 'tf', 1e-4)));
assert(R3.nRestored == 1, 'a failed repair must restore');
Qa = load(sheetMat);
assert(isequal(Qa.OK, Qb.OK) && max(abs(Qa.TF(:) - Qb.TF(:))) == 0 && isequaln(Qa.Z8, Qb.Z8), 'restored slices');

% 5. unsolved mode on a 2 x 1 sheet: cell 2 empty and retired, seeded from cell 1
S2 = load(sheetMat);
S2.sD = [S2.sD, 0.12];
f = {'TF', 'FLYKM', 'ACCDZ', 'RES', 'WALL'};
for kf = 1:numel(f), S2.(f{kf}) = cat(1, S2.(f{kf}), nan(1, 1, 2)); end
S2.OK  = cat(1, S2.OK, false(1, 1, 2));
S2.Z8  = cat(2, S2.Z8, nan(8, 1, 1, 2));
S2.ATT = [S2.ATT; 2];
sheet2 = fullfile(wd, 'two_cells.mat');
save(sheet2, '-struct', 'S2');
eo2 = eo;  eo2.logFile = fullfile(wd, 'repair2.log');
R4 = repair_gto_dro_cells(sheet2, 'unsolved', struct('engineOpts', eo2));
Q2 = load(sheet2, 'OK', 'TF');
assert(R4.nCandidates == 1 && Q2.OK(2, 1, 1), 'phase continuation solved the neighbour cell at 15 N');
S4 = load(R4.sideMat, 'rows');
assert(S4.rows(end).seedFrom.iD == 1, 'seeded from cell 1');
fprintf('  phase continuation: cell 2 (sD 0.12) tf 15 N = %.5f\n', Q2.TF(2, 1, 1));

rmdir(wd, 's');
ok = true;
fprintf('test_repair_gto_dro_cells: PASS\n');
end
```

Note for the implementer: check 5 hands the engine a sheet whose `sD` vector is not the uniform grid the engine rebuilds from `sD0/nD`. The repair tool therefore passes the engine `nD = numel(sD)`, `sD0 = sD(1)` only when the sheet's `sD` IS uniform (`max(abs(diff(sD) - 1/nD)) < 1e-12`); otherwise it refuses with `repair_gto_dro_cells:grid`. Make check 5's second phase `0.10 + 1/12` is NOT what the test uses (0.12), so the test must instead build the two-cell sheet with `sD = [0.10, 0.10 + 1/12]`... no: the engine recomputes `sD = mod(sD0 + (0:nD-1)/nD, 1)` with `nD = 2`, which gives `[0.10, 0.60]`. **Ruling for the implementer:** the engine owns the grid, so the tool passes `nD, nA, sD0, sA0` from the sheet's `meta` (`meta.sD0`, `meta.sA0`) and `nD = numel(sD)`, `nA = numel(sA)`, and check 5 must set `S2.sD = mod(0.10 + (0:1)/2, 1)` and `S2.meta.sD0 = 0.10` (already true), so the second cell is at phase 0.60 (apogee side; stage 1 did not measure it, so the assertion on `Q2.OK(2,1,1)` is a genuine measurement; if it fails the implementer reports the measured engine log and the test keeps the assertion on `R4.nCandidates == 1` and the sidecar row only). Update the test accordingly before running it.

- [ ] **Step 2: Run it to verify it fails**

Expected: `Unrecognized function 'repair_gto_dro_cells'`.

- [ ] **Step 3: Implement**

```matlab
function R = repair_gto_dro_cells(sheetMat, mode, opts)
% REPAIR_GTO_DRO_CELLS  Re-solve cells of a GTO -> DRO sheet from external seeds.
%
%   Two repairs the fleet cannot do for itself (stage 1 findings,
%   GTO_DRO/README.md), both through thrust_ladder_library's seedTop/redo:
%     'unsolved'  a cell the engine retired without any OK rung (the perigee
%                 column at 15 N burns its attempts at the iteration cap) is
%                 re-solved once, seeded from the nearest solved cell in the
%                 departure-phase direction (phase continuation).
%     'basin'     every OK top-rung entry gets the conjugate re-solve
%                 (seed_from_z8 -> ms_tfmin conjTest, the stage 1 S4
%                 junction test); a refuted entry whose second rung is OK is
%                 re-solved from that second-rung entry (reverse warm start:
%                 stage 1 measured 15 N <- 12 N landing on the faster basin).
%                 The cell's slices are backed up and RESTORED when the
%                 repair fails or is not faster.
%   Every attempt is a row in the sidecar <sheet>_repairs.mat.
%
% INPUTS:
%   sheetMat - sheet .mat written by thrust_ladder_library           [char]
%   mode     - 'unsolved' | 'basin'                                   [char]
%   opts     - .engineOpts (required: the sheet's engine options, e.g.
%              run_gto_dro_catalog planOnly optsA; rungs/nD/nA/sD0/sA0/
%              cells/redo/seedTop/tf0 are overwritten from the sheet),
%              .conjOnly [false] (basin: verdicts only), .force [false]
%              (repair cells the sidecar already lists), .seedOverride []
%              (a seedTop struct used for every repair; tests),
%              .maxCells [inf], .logFile ['']
%
% OUTPUTS:
%   R        - .mode .nCandidates .nRepaired .nRestored .nFailed
%              .nConjRefuted .conjPass [nD x nA] (NaN/0/1) .rows .sideMat
%
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(here, 'tests'), fullfile(here, 'catalog'));
    sheetMat = fixture_demo_sheet();  mode = 'basin';
    plan = run_gto_dro_catalog(2, [], [], [15 12], struct('planOnly', true));
    opts = struct('engineOpts', plan.optsA, 'conjOnly', true);
end
assert(any(strcmp(mode, {'unsolved', 'basin'})), 'repair_gto_dro_cells:mode', 'mode must be unsolved or basin');
assert(isfield(opts, 'engineOpts'), 'repair_gto_dro_cells:engineOpts', 'opts.engineOpts is required');
conjOnly = fdef(opts, 'conjOnly', false);
force    = fdef(opts, 'force', false);
seedOvr  = fdef(opts, 'seedOverride', []);
maxCells = fdef(opts, 'maxCells', inf);
logFile  = fdef(opts, 'logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));
sideMat  = [sheetMat(1:end-4) '_repairs.mat'];

Q = load(sheetMat);
meta = Q.meta;  rungs = Q.rungs;  [nD, nA, nR] = size(Q.OK);
rows = struct('mode', {}, 'iD', {}, 'iA', {}, 'seedFrom', {}, 'tfBefore', {}, 'tfAfter', {}, ...
              'okAfter', {}, 'outcome', {}, 'conjBefore', {}, 'conjAfter', {}, 'date', {}, 'wallS', {});
if isfile(sideMat), S = load(sideMat, 'rows'); rows = S.rows; end
already = false(nD, nA);
for kr = 1:numel(rows)
    if strcmp(rows(kr).mode, mode), already(rows(kr).iD, rows(kr).iA) = true; end
end

% engine options: the sheet owns the grid, the rungs and the key
eo = opts.engineOpts;
for f = {'cells', 'redo', 'seedTop', 'tf0', 'maxCells', 'batchSec'}
    if isfield(eo, f{1}), eo = rmfield(eo, f{1}); end
end
eo.rungs = rungs;  eo.nD = nD;  eo.nA = nA;  eo.sD0 = meta.sD0;  eo.sA0 = meta.sA0;
eo.resume = true;  eo.logFile = logFile;
sDgrid = mod(meta.sD0 + (0:nD-1)/nD, 1);
assert(max(abs(sDgrid - Q.sD(:).')) < 1e-12, 'repair_gto_dro_cells:grid', ...
       'the sheet departure grid is not the engine grid (sD0 %g, nD %d)', meta.sD0, nD);

ndp = nd_propulsion(rungs(1), meta.ispS, meta.m0kg, meta.lStar, meta.tStar);
[tD, rvD] = get_family_orbit(meta.depFamily, meta.depParams);
[tA, rvA] = get_family_orbit(meta.arrFamily, meta.arrParams);
dep = phase_state(tD, rvD);  arr = phase_state(tA, rvA);

conjPass = nan(nD, nA);
cand = zeros(0, 2);  seedOf = {};
switch mode
    case 'unsolved'
        maxAtt = fdef(eo, 'maxAtt', 2);
        for iD = 1:nD
            for iA = 1:nA
                if any(Q.OK(iD, iA, :)) || Q.ATT(iD, iA) < maxAtt, continue, end
                if already(iD, iA) && ~force, continue, end
                nb = [];
                for step = 1:floor(nD/2)
                    for sgn = [-1 1]
                        jD = mod(iD - 1 + sgn*step, nD) + 1;
                        if Q.OK(jD, iA, 1), nb = jD; break, end
                    end
                    if ~isempty(nb), break, end
                end
                if isempty(nb)
                    lg('  (%2d,%2d) unsolved: no solved neighbour in the departure direction -- skipped', iD, iA);
                    continue
                end
                cand(end+1, :) = [iD iA];
                seedOf{end+1} = struct('iD', nb, 'iA', iA, 'kr', 1);
            end
        end
    case 'basin'
        for iD = 1:nD
            for iA = 1:nA
                if ~Q.OK(iD, iA, 1), continue, end
                conjPass(iD, iA) = conj_verdict(Q.Z8(:, iD, iA, 1), dep(Q.sD(iD)), arr(Q.sA(iA)), ndp, meta, lg, iD, iA);
                if conjPass(iD, iA) == 0 && nR >= 2 && Q.OK(iD, iA, 2) && (~already(iD, iA) || force)
                    cand(end+1, :) = [iD iA];
                    seedOf{end+1} = struct('iD', iD, 'iA', iA, 'kr', 2);
                end
            end
        end
end
R = struct('mode', mode, 'nCandidates', size(cand, 1), 'nRepaired', 0, 'nRestored', 0, 'nFailed', 0, ...
           'nConjRefuted', nnz(conjPass == 0), 'conjPass', conjPass, 'rows', rows, 'sideMat', sideMat);
lg('repair %s: %d candidate cell(s), %d refuted top rungs', mode, R.nCandidates, R.nConjRefuted);
if conjOnly || isempty(cand), R.rows = rows; return, end

for kc = 1:min(size(cand, 1), maxCells)
    iD = cand(kc, 1);  iA = cand(kc, 2);  sf = seedOf{kc};
    t0 = tic;
    bk = backup_cell(Q, iD, iA);
    if isempty(seedOvr)
        seedTop = seed_top_from_entry(sheetMat, sf.iD, sf.iA, sf.kr, rungs(1));
    else
        seedTop = seedOvr;
    end
    row = struct('mode', mode, 'iD', iD, 'iA', iA, 'seedFrom', sf, 'tfBefore', bk.TF(1), 'tfAfter', NaN, ...
                 'okAfter', false, 'outcome', '', 'conjBefore', conjPass(iD, iA), 'conjAfter', NaN, ...
                 'date', datestr(now, 'yyyy-mm-dd HH:MM'), 'wallS', NaN);
    e = eo;  e.cells = [iD iA];  e.redo = true;  e.seedTop = seedTop;
    try
        thrust_ladder_library(sheetMat, e);
        Qn = load(sheetMat);
        tfNew = Qn.TF(iD, iA, 1);  okNew = Qn.OK(iD, iA, 1);
        better = okNew && (isnan(row.tfBefore) || tfNew < row.tfBefore);
        if better
            row.outcome = 'repaired';  row.tfAfter = tfNew;  row.okAfter = true;
            if strcmp(mode, 'basin')
                row.conjAfter = conj_verdict(Qn.Z8(:, iD, iA, 1), dep(Q.sD(iD)), arr(Q.sA(iA)), ndp, meta, lg, iD, iA);
            end
            Q = Qn;  R.nRepaired = R.nRepaired + 1;
        else
            Qn = restore_cell(Qn, bk, iD, iA);
            save(sheetMat, '-struct', 'Qn');  Q = Qn;
            row.outcome = pick(okNew, 'restored (not faster)', 'restored (repair failed)');
            row.tfAfter = tfNew;  row.okAfter = okNew;
            R.nRestored = R.nRestored + 1;
        end
    catch ME
        Qn = load(sheetMat);  Qn = restore_cell(Qn, bk, iD, iA);
        save(sheetMat, '-struct', 'Qn');  Q = Qn;
        row.outcome = ['error: ' ME.message];  R.nFailed = R.nFailed + 1;  R.nRestored = R.nRestored + 1;
    end
    row.wallS = toc(t0);
    rows(end+1) = row;
    save(sideMat, 'rows');
    lg('  (%2d,%2d) %s: %s  tf %.5f -> %.5f  (%.0f s)', iD, iA, mode, row.outcome, row.tfBefore, row.tfAfter, row.wallS);
end
R.rows = rows;
if nargin == 0, disp(R), end
end

% ------------------------------------------------------------------------
function cp = conj_verdict(z8, rv0, rvf, ndp, meta, lg, iD, iA)
% CONJ_VERDICT  The stage 1 S4 junction test on one entry: NaN when the
% re-solve does not reproduce the entry.  INPUTS: z8 [8x1]; rv0, rvf [6x1];
% ndp; meta; lg; iD; iA.  OUTPUTS: cp - 1 / 0 / NaN.
cp = NaN;
try
    seed = seed_from_z8(z8, rv0(:), 24, ndp.Tnd, ndp.cnd, meta.muStar);
    [z, info] = ms_tfmin(rv0(:).', rvf(:).', seed, ndp.Tnd, ndp.cnd, meta.muStar, ...
                         struct('conjTest', true, 'wallSec', 120));
    if info.converged && max(abs(z - z8)) < 1e-6
        cp = double(info.conj.pass);
    else
        lg('  (%2d,%2d) conjugate re-solve did not reproduce the entry (|dz| %.1e) -- verdict NaN', iD, iA, max(abs(z - z8)));
    end
catch ME
    lg('  (%2d,%2d) conjugate re-solve error: %s -- verdict NaN', iD, iA, ME.message);
end
end

function bk = backup_cell(Q, iD, iA)
% BACKUP_CELL  Copy one cell's slices.  INPUTS: Q; iD; iA.  OUTPUTS: bk struct.
bk = struct('TF', squeeze(Q.TF(iD, iA, :)), 'FLYKM', squeeze(Q.FLYKM(iD, iA, :)), ...
            'ACCDZ', squeeze(Q.ACCDZ(iD, iA, :)), 'RES', squeeze(Q.RES(iD, iA, :)), ...
            'WALL', squeeze(Q.WALL(iD, iA, :)), 'OK', squeeze(Q.OK(iD, iA, :)), ...
            'Z8', squeeze(Q.Z8(:, iD, iA, :)), 'ATT', Q.ATT(iD, iA));
end

function Q = restore_cell(Q, bk, iD, iA)
% RESTORE_CELL  Put a backed-up cell back.  INPUTS: Q; bk; iD; iA.  OUTPUTS: Q.
nR = numel(bk.TF);
Q.TF(iD, iA, :) = reshape(bk.TF, 1, 1, nR);     Q.FLYKM(iD, iA, :) = reshape(bk.FLYKM, 1, 1, nR);
Q.ACCDZ(iD, iA, :) = reshape(bk.ACCDZ, 1, 1, nR); Q.RES(iD, iA, :) = reshape(bk.RES, 1, 1, nR);
Q.WALL(iD, iA, :) = reshape(bk.WALL, 1, 1, nR);  Q.OK(iD, iA, :) = reshape(bk.OK, 1, 1, nR);
Q.Z8(:, iD, iA, :) = reshape(bk.Z8, 8, 1, 1, nR); Q.ATT(iD, iA) = bk.ATT;
end

function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end

function v = pick(c, a, b)
% PICK  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end

function logmsg(f, s)
% LOGMSG  Append-mode file log, stdout when f is empty.  INPUTS: f; s.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); end
end
```

`squeeze` on a `[1 x 1 x nR]` slice gives `[nR x 1]`; on `Z8(:, iD, iA, :)` it gives `[8 x nR]` (for `nR = 1` it gives `[8 x 1]`, and the `reshape(..., 8, 1, 1, nR)` restore is still correct). `ms_tfmin` takes row endpoints (`rv0(1:6)` rows in the engine); the `.'` above makes rows from the column `phase_state` returns.

- [ ] **Step 4: Run the test**

Expected: PASS in about 10 to 15 minutes (one conjugate re-solve, one 15 + 12 N repair, one failed repair that restores, one phase-continuation solve plus its 12 N rung). Record every printed t_f in the report. If check 5's apogee-side cell does not solve, report the engine log tail and keep the sidecar/candidate assertions (see the note under Step 1).

- [ ] **Step 5: Commit** (controller)

```bash
git add orbit_transfer/GTO_DRO/lib/repair_gto_dro_cells.m orbit_transfer/GTO_DRO/tests/test_repair_gto_dro_cells.m orbit_transfer/GTO_DRO/tests/run_tests.m
git commit -m "GTO_DRO: repair_gto_dro_cells -- 'unsolved' re-solves retired cells by phase continuation from the nearest solved departure-phase neighbour, 'basin' runs the S4 junction test on every top-rung entry and re-solves refuted ones from their own 12 N entry (reverse warm start), backup + restore when a repair fails or is not faster, every attempt in a sidecar; test measures the demo cell's basin repair"
```

---

### Task 5: The post-fleet chain `build_gto_dro_catalog`, the audit, the vendored pickers

**Files:**
- Create: `orbit_transfer/GTO_DRO/catalog/build_gto_dro_catalog.m`
- Create: `orbit_transfer/GTO_DRO/catalog/audit_gto_dro_entries.m`
- Copy: `orbit_transfer/GTO_tulip/catalog/costate_catalog_pick.m` and `costate_lib_describe.m` into `orbit_transfer/GTO_DRO/catalog/` unchanged except the Revision History line (`copied for the GTO -> DRO catalog, 11/10/2026`)
- Test: `orbit_transfer/GTO_DRO/tests/test_build_gto_dro_catalog.m`
- Modify: `orbit_transfer/GTO_DRO/tests/run_tests.m`, `orbit_transfer/GTO_DRO/README.md` (stage 2 section: files, the chain's stages, how to run the pilot and the fleet)

**Interfaces:**
- Consumes: Tasks 2 to 4; `build_costate_catalog_family(catDir, outMat, spec)` with the GTO `complete_dep_params` patch copied verbatim from `build_gto_catalog.m` (it asserts `meta.tauDRO == orientDeg`); `conj_catalog_pass(catMat, opts)` and `gates_catalog_pass(catMat, opts)` (sidecars `_conjprog.mat` / `_gatesprog.mat`, `writeback` after a complete census); `costate_catalog_pick(cat_, tauDep, arrKey, depPhaseDays, arrPhaseDays, thrustN, warnFlag)` where for this catalog `tauDep` is the orientation in degrees and `arrKey` is the DRO period in ND.
- Produces: `chainBlockers = build_gto_dro_catalog(runIn)` with `runIn` a struct of stage switches `census, earth, repairUnsolved, repairBasin, package, conj, gates, audit, pictures, deliverable` (all true by default; `runIn.catDir`, `runIn.outDir` overrides for tests; `runIn.sheetGlob` default `'gto_dro_o*_tau*.mat'`, which excludes `_deep`, `_earth`, `_repairs` and `.bak_*` by the explicit name list below, never by glob alone). Also `D = deep_probe(cells, rungs)` as a separate stage function inside the same file? No: the deep probe is a front-door call with a separate sheet file, run by Task 6 by hand; the chain only reports it if present.
- Files: `results/costate_catalog_gto_dro.mat` (variable `costate_catalog_gto_dro`), `results/audit_gto_dro.mat`, `results/pictures/*.png`, `results/deliverable/costate_catalog_gto_dro_<yyyymmdd>.zip`, `results/chain_<yyyymmdd_HHMM>.log`.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_build_gto_dro_catalog()
% TEST_BUILD_GTO_DRO_CATALOG  The chain on a one-sheet, one-cell catalog
% folder built from the demo fixture: census -> earth -> package -> conj ->
% gates -> audit -> pictures -> deliverable, every stage fenced, blockers
% returned; the repair stages are switched off here (Task 4 tests them).
% INPUTS: none.  OUTPUTS: ok - logical.
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(here);  addpath(fullfile(fileparts(here), 'catalog'));
src = fixture_demo_sheet();
wd  = fullfile(tempdir, 'test_build_gto_dro');
if isfolder(wd), rmdir(wd, 's'); end
mkdir(fullfile(wd, 'catalog'));
copyfile(src, fullfile(wd, 'catalog', 'gto_dro_o000_tau1.mat'));

run_ = struct('census', true, 'earth', true, 'repairUnsolved', false, 'repairBasin', false, ...
              'package', true, 'conj', true, 'gates', true, 'audit', true, 'pictures', true, ...
              'deliverable', true, 'catDir', fullfile(wd, 'catalog'), 'outDir', wd, ...
              'expectSheets', 1);
blockers = build_gto_dro_catalog(run_);
catMat = fullfile(wd, 'costate_catalog_gto_dro.mat');
assert(isfile(catMat), 'package written');
C = load(catMat);  cat_ = C.costate_catalog_gto_dro;
assert(numel(cat_.sheets) == 1 && cat_.sheets(1).tauDRO == 0, 'sheet key = orientDeg');
assert(cat_.sheets(1).dep_params.orientDeg == 0 && isfield(cat_.sheets(1).dep_params, 'sma_km'), 'dep_params completed');
assert(cat_.sheets(1).arr_params.tau == 1, 'arrival DRO period kept');
assert(isfield(cat_.sheets(1), 'conj_pass') && isfield(cat_.sheets(1), 'gate_min_lamv'), 'sweeps written back');
assert(isfile(fullfile(wd, 'audit_gto_dro.mat')), 'audit written');
A = load(fullfile(wd, 'audit_gto_dro.mat'));
assert(all([A.results.passFly]) && all([A.results.passAccept]) && all([A.results.passEarth]), 'audit gates');
assert(~isempty(dir(fullfile(wd, 'pictures', '*.png'))), 'pictures');
z = dir(fullfile(wd, 'deliverable', '*.zip'));
% the 15 N entry of the fixture is the refuted basin, so the conj sweep
% reports one FAIL and the chain must carry it as a blocker and refuse the zip
assert(any(contains(blockers, 'conj')), 'refuted entry must be a blocker');
assert(isempty(z), 'deliverable refused while a blocker stands');
fprintf('  blockers: %s\n', strjoin(blockers, ' | '));

% pickers resolve the sheet by (orientation, DRO period)
[tf_nd, z8, info] = costate_catalog_pick(cat_, 0, 1, 0.10*cat_.sheets(1).derive_period_days_dep, 0, 12, false);
assert(isfinite(tf_nd) && numel(z8) == 8, 'picker returns the 12 N entry');
D = costate_lib_describe(cat_);
assert(isstruct(D), 'describe runs');

rmdir(wd, 's');
ok = true;
fprintf('test_build_gto_dro_catalog: PASS\n');
end
```

Note: `derive_period_days_dep` is a placeholder name for whatever the packaged sheet exposes as the departure period in days; the implementer reads `build_costate_catalog_family`'s output (`cat_.sheets(1)` fields, `cat_.derive`) and `costate_catalog_pick`'s header to pass the departure phase in the picker's units, and replaces that expression with the real one, stating in the report which field it used. The picker's phase convention is DAYS (its header), so the departure phase 0.10 of the GTO's Kepler period must be converted with the period the catalog stores.

- [ ] **Step 2: Run it to verify it fails**

Expected: `Unrecognized function 'build_gto_dro_catalog'`.

- [ ] **Step 3: Write the chain script**

```matlab
function chainBlockers = build_gto_dro_catalog(runIn)
%% Purpose:
%
%   THE CHAIN for the GTO -> DRO costate catalog AFTER the fleet
%   (run_gto_dro_catalog / run_gto_dro_batched.sh write the sheets; this
%   script turns sheets into a shippable catalog), build_70mN_library
%   style: parameter block first, one switch per stage, every stage
%   fenced, a blocker list the deliverable stage refuses on.
%
%     1. census           sheets found, cells solved per rung, retired cells
%     2. earth            earth_clearance_sweep per sheet, writeback (demotes)
%     3. repairUnsolved   repair_gto_dro_cells 'unsolved' (phase continuation)
%     4. repairBasin      repair_gto_dro_cells 'basin' (S4 on every top rung,
%                         reverse warm start on refuted ones)
%     5. package          build_costate_catalog_family (sheet key = orientDeg)
%     6. conj             conj_catalog_pass census -> writeback
%     7. gates            gates_catalog_pass census -> writeback
%     8. audit            audit_gto_dro_entries (flown, accepted, Earth floor)
%     9. pictures         coverage per sheet, t_f and m_f ranges per rung
%    10. deliverable      zip of catalog + pickers + README (refused on blockers)
%
%   Stages 2 to 4 modify the SHEETS (each keeps its own backup/sidecar);
%   stage 5 rebuilds the package from the sheets, and stages 6 to 7 run on
%   the package, so re-running the chain after a repair re-packages and
%   re-sweeps (the sweep sidecars are deleted when the package is rebuilt,
%   since their entry indices would be stale).
%
%% Inputs:
%
%  runIn                    struct (optional)       Stage switches
%                                                   .census .earth
%                                                   .repairUnsolved
%                                                   .repairBasin .package
%                                                   .conj .gates .audit
%                                                   .pictures .deliverable
%                                                   (all true); .catDir
%                                                   [../results/catalog],
%                                                   .outDir [../results],
%                                                   .expectSheets [16],
%                                                   .floorEarthKm [200],
%                                                   .maxRepair [inf]
%
%% Outputs:
%
%  chainBlockers            cell                    Blocker strings (empty
%                                                   = clean); printed too
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, runIn = struct(); end
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(here);

%% ========================================================================
%  0. PARAMETERS AND STAGE SWITCHES
%% ========================================================================
run_ = struct('census', true, 'earth', true, 'repairUnsolved', true, 'repairBasin', true, ...
              'package', true, 'conj', true, 'gates', true, 'audit', true, 'pictures', true, ...
              'deliverable', true);
for f = fieldnames(run_).'
    if isfield(runIn, f{1}), run_.(f{1}) = runIn.(f{1}); end
end
      catDir = fdef(runIn, 'catDir', fullfile(fileparts(here), 'results', 'catalog'));
      outDir = fdef(runIn, 'outDir', fullfile(fileparts(here), 'results'));
expectSheets = fdef(runIn, 'expectSheets', 16);
floorEarthKm = fdef(runIn, 'floorEarthKm', 200);
   maxRepair = fdef(runIn, 'maxRepair', inf);
     orients = [0 90 180 270];
        taus = [0.5 1 2 3];
if ~isfolder(outDir), mkdir(outDir); end
logFile = fullfile(outDir, ['chain_' datestr(now, 'yyyymmdd_HHMM') '.log']);
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));
files = struct('catalog', fullfile(outDir, 'costate_catalog_gto_dro.mat'), ...
               'audit', fullfile(outDir, 'audit_gto_dro.mat'), ...
               'pictures', fullfile(outDir, 'pictures'), ...
               'deliverable', fullfile(outDir, 'deliverable'));
blockers = {};

% the sheets, by their official names (never a bare glob: _deep, _earth,
% _repairs and .bak_* must not enter)
sheetMats = {};
for od = orients
    for ta = taus
        p = fullfile(catDir, sprintf('gto_dro_o%03d_tau%g.mat', od, ta));
        if isfile(p), sheetMats{end+1} = p; end
    end
end

%% ========================================================================
%  1. CENSUS
%% ========================================================================
if run_.census
    try
        nS = numel(sheetMats);
        lg('1. census: %d of %d sheets on disk', nS, expectSheets);
        if nS < expectSheets, blockers{end+1} = sprintf('census: %d of %d sheets', nS, expectSheets); end
        for ks = 1:nS
            Q = load(sheetMats{ks}, 'OK', 'ATT', 'rungs');
            [nD, nA, nR] = size(Q.OK);
            solved = squeeze(sum(sum(Q.OK, 1), 2)).';
            retired = nnz(~any(Q.OK, 3) & Q.ATT >= 2);
            [~, nm] = fileparts(sheetMats{ks});
            lg('   %-22s cells %3d  per rung [%s]  retired %2d', nm, nD*nA, num2str(solved), retired);
        end
    catch ME
        lg('1. census THREW: %s', ME.message);  blockers{end+1} = ['census threw: ' ME.message];
    end
end

%% ========================================================================
%  2. EARTH CLEARANCE SWEEP (sheets; writeback demotes)
%% ========================================================================
if run_.earth
    for ks = 1:numel(sheetMats)
        try
            S = earth_clearance_sweep(sheetMats{ks}, struct('floorKm', floorEarthKm, 'logFile', logFile));
            if S.done
                S = earth_clearance_sweep(sheetMats{ks}, struct('floorKm', floorEarthKm, 'writeback', true, 'logFile', logFile));
                [~, nm] = fileparts(sheetMats{ks});
                lg('2. earth %-22s min %.0f km, %d demoted', nm, S.minAltKm, S.nBelow);
                if S.nBelow > 0, blockers{end+1} = sprintf('earth: %s demoted %d entries', nm, S.nBelow); end
            else
                blockers{end+1} = sprintf('earth: %s incomplete (%d to do)', sheetMats{ks}, S.nTodo);
            end
        catch ME
            lg('2. earth THREW on %s: %s', sheetMats{ks}, ME.message);
            blockers{end+1} = ['earth threw: ' ME.message];
        end
    end
end

%% ========================================================================
%  3./4. REPAIRS (sheets; sidecars)
%% ========================================================================
plan = run_gto_dro_catalog([], [], [], [], struct('planOnly', true));
for mode = {'unsolved', 'basin'}
    sw = pick(strcmp(mode{1}, 'unsolved'), run_.repairUnsolved, run_.repairBasin);
    if ~sw, continue, end
    for ks = 1:numel(sheetMats)
        try
            [~, nm] = fileparts(sheetMats{ks});
            kp = find(strcmp({plan.tag}, nm), 1);
            assert(~isempty(kp), 'build_gto_dro_catalog:plan', 'no plan row for %s', nm);
            R = repair_gto_dro_cells(sheetMats{ks}, mode{1}, ...
                    struct('engineOpts', plan(kp).optsA, 'maxCells', maxRepair, 'logFile', logFile));
            lg('3/4. repair %-8s %-22s candidates %d repaired %d restored %d failed %d refuted %d', ...
               mode{1}, nm, R.nCandidates, R.nRepaired, R.nRestored, R.nFailed, R.nConjRefuted);
            if R.nFailed > 0, blockers{end+1} = sprintf('repair %s: %s %d failed', mode{1}, nm, R.nFailed); end
        catch ME
            lg('3/4. repair %s THREW on %s: %s', mode{1}, sheetMats{ks}, ME.message);
            blockers{end+1} = ['repair threw: ' ME.message];
        end
    end
end

%% ========================================================================
%  5. PACKAGE (sheet key = orientDeg; dep_params completed on disk)
%% ========================================================================
if run_.package
    try
        complete_dep_params(catDir, sheetMats);
        spec = struct( ...
          'glob', 'gto_dro_o*_tau*.mat', ...
          'name', 'costate_catalog_gto_dro', ...
          'description', ['Minimum-time low-thrust GTO->DRO costate CATALOG: one (GTO orientation ', ...
            'x DRO period) phasing sheet per pair, 12 x 6 phasing, thrust rungs walked from 15 N. ', ...
            'COMPACT format: canonical ND quantities only; derived values via .derive. Each entry ', ...
            'z8 = [lambda_r(3); lambda_v(3); lambda_m; tf] is a converged PMP solution accepted by ', ...
            'pumpkyn.cr3bp.tfMin, solved with a 200 km Earth-altitude floor and swept for the ', ...
            'flown Earth clearance. Departure is the GTO PSEUDO-family (algebraic Kepler locus, ', ...
            'see get_family_orbit), 350 x 35786 km, orientDeg = argument of perigee from the ', ...
            'Earth->Moon line, prograde; sheets(k).tauDRO/.tau_dep hold orientDeg IN DEGREES (the ', ...
            'GTO sheet key rule), the DRO period is arr_params.tau (ND). conj_pass marks entries ', ...
            'the conjugate test refutes: consult it, a refuted entry is an extremal, not a minimum.'], ...
          'provenance', ['GTO -> DRO catalog run (run_gto_dro_catalog), 2026-11; pilot sheet ', ...
            '(0 deg, tau 1) then the fleet; two-point top-rung multistart, first-success-wins; ', ...
            'Earth sweep, phase-continuation and reverse-warm-start repairs recorded per sheet ', ...
            'in _earth.mat / _repairs.mat sidecars. Casey/Koblick 2026-11.'], ...
          'depReconstruction', ['gto: [tau,rv] = get_family_orbit(''gto'', dep_params) -> algebraic ', ...
            'Kepler locus (NO propagation); state at phase f: phase_state(tau, rv)(f); dro: ', ...
            '[tau,rv] = get_family_orbit(''dro'', arr_params) -> getDRO seed, cont_np, prop; ', ...
            'state at phase f: phase_state(tau, rv)(f)']);
        % the sweep sidecars index entries of the OLD package: delete them
        for sc = {'_conjprog.mat', '_gatesprog.mat'}
            p = [files.catalog(1:end-4) sc{1}];
            if isfile(p), delete(p); end
        end
        % the glob above must not pick up the deep probe or backups: the
        % packager globs, so move nothing -- assert the extras are not matched
        extra = dir(fullfile(catDir, 'gto_dro_o*_tau*.mat'));
        extra = extra(~ismember(fullfile({extra.folder}, {extra.name}), sheetMats));
        assert(isempty(extra), 'build_gto_dro_catalog:glob', ...
               'the package glob would pick up %s -- rename it', strjoin({extra.name}, ', '));
        cat_ = build_costate_catalog_family(catDir, files.catalog, spec);
        lg('5. package: %d sheets, %d entries -> %s', numel(cat_.sheets), count_entries(cat_), files.catalog);
    catch ME
        lg('5. package THREW: %s', ME.message);  blockers{end+1} = ['package threw: ' ME.message];
    end
end

%% ========================================================================
%  6./7. CONJUGATE AND HYPOTHESIS SWEEPS (package; sidecars; writeback)
%% ========================================================================
if run_.conj && isfile(files.catalog)
    try
        Sc = conj_catalog_pass(files.catalog, struct('logFile', logFile));
        if Sc.done
            Sc = conj_catalog_pass(files.catalog, struct('logFile', logFile, 'writeback', true));
            lg('6. conj: %d pass, %d FAIL, %d not run', Sc.nPass, Sc.nFail, Sc.nNotrun);
            if Sc.nFail > 0, blockers{end+1} = sprintf('conj: %d refuted entries (consult conj_pass; repairBasin re-solves refuted TOP rungs)', Sc.nFail); end
        else
            blockers{end+1} = sprintf('conj incomplete: %d to do', Sc.nTodo);
        end
    catch ME
        lg('6. conj THREW: %s', ME.message);  blockers{end+1} = ['conj threw: ' ME.message];
    end
end
if run_.gates && isfile(files.catalog)
    try
        Sg = gates_catalog_pass(files.catalog, struct('logFile', logFile));
        if Sg.done
            Sg = gates_catalog_pass(files.catalog, struct('logFile', logFile, 'writeback', true));
            lg('7. gates: %d done, H2 fail %d, H3 fail %d, abnormal %d', Sg.nDone, Sg.nH2fail, Sg.nH3fail, Sg.nAbnormal);
            if Sg.nH2fail + Sg.nH3fail + Sg.nAbnormal > 0
                blockers{end+1} = sprintf('gates: %d entries fail a hypothesis gate', Sg.nH2fail + Sg.nH3fail + Sg.nAbnormal);
            end
        else
            blockers{end+1} = sprintf('gates incomplete: %d to do', Sg.nTodo);
        end
    catch ME
        lg('7. gates THREW: %s', ME.message);  blockers{end+1} = ['gates threw: ' ME.message];
    end
end

%% ========================================================================
%  8. AUDIT (3-entry replay with the Earth floor as a gate)
%% ========================================================================
if run_.audit
    try
        results = audit_gto_dro_entries(catDir, sheetMats, floorEarthKm);
        save(files.audit, 'results');
        lg('8. audit: fly %d/3, accept %d/3, earth %d/3', nnz([results.passFly]), nnz([results.passAccept]), nnz([results.passEarth]));
        if ~all([results.passFly] & [results.passAccept] & [results.passEarth])
            blockers{end+1} = 'audit: a replayed entry failed a gate';
        end
    catch ME
        lg('8. audit THREW: %s', ME.message);  blockers{end+1} = ['audit threw: ' ME.message];
    end
end

%% ========================================================================
%  9. PICTURES (coverage per sheet and rung; t_f and m_f ranges)
%% ========================================================================
if run_.pictures
    try
        if ~isfolder(files.pictures), mkdir(files.pictures); end
        figsBefore = findall(0, 'Type', 'figure');
        set(0, 'DefaultFigureVisible', 'off');
        catalog_pictures(sheetMats, files.pictures);
        close(setdiff(findall(0, 'Type', 'figure'), figsBefore));
        set(0, 'DefaultFigureVisible', 'on');
        lg('9. pictures -> %s', files.pictures);
    catch ME
        set(0, 'DefaultFigureVisible', 'on');
        lg('9. pictures THREW: %s', ME.message);  blockers{end+1} = ['pictures threw: ' ME.message];
    end
end

%% ========================================================================
%  10. DELIVERABLE (refused on any blocker but a pictures failure)
%% ========================================================================
shipBlockers = blockers(~contains(blockers, 'pictures threw'));
if run_.deliverable
    if ~isempty(shipBlockers)
        lg('10. deliverable REFUSED -- %d blocker(s):', numel(shipBlockers));
        for kb = 1:numel(shipBlockers), lg('     - %s', shipBlockers{kb}); end
    else
        if ~isfolder(files.deliverable), mkdir(files.deliverable); end
        stamp = datestr(now, 'yyyymmdd');
        stage = fullfile(files.deliverable, ['costate_catalog_gto_dro_' stamp]);
        if isfolder(stage), rmdir(stage, 's'); end
        mkdir(stage);
        copyfile(files.catalog, stage);
        copyfile(fullfile(here, 'costate_catalog_pick.m'), stage);
        copyfile(fullfile(here, 'costate_lib_describe.m'), stage);
        copyfile(fullfile(fileparts(here), 'README.md'), fullfile(stage, 'README.md'));
        zipPath = [stage '.zip'];
        zip(zipPath, stage);
        lg('10. deliverable: %s', zipPath);
    end
end

%% ========================================================================
%  BLOCKER SUMMARY
%% ========================================================================
if isempty(blockers)
    lg('CHAIN CLEAN: every stage that ran, passed.');
else
    lg('CHAIN COMPLETE WITH %d BLOCKER(S):', numel(blockers));
    for kb = 1:numel(blockers), lg('  - %s', blockers{kb}); end
    lg('Nothing ships until these clear; the measurements above are still valid.');
end
chainBlockers = blockers;
end

% ------------------------------------------------------------------------
function complete_dep_params(catDir, sheetMats)
% COMPLETE_DEP_PARAMS  Patch each sheet's meta.depParams to the full
% {sma_km, ecc, orientDeg} triple, idempotently and atomically (copy of
% GTO_tulip/catalog/build_gto_catalog.m's subfunction, by explicit file
% list instead of a glob).  INPUTS: catDir char; sheetMats cell.  OUTPUTS: none.
for kf = 1:numel(sheetMats)
    fpath = sheetMats{kf};
    S = load(fpath);
    dp = S.meta.depParams;
    if isfield(dp, 'sma_km') && isfield(dp, 'ecc'), continue, end
    assert(strcmpi(S.meta.depFamily, 'gto'), 'complete_dep_params:family', ...
        'sheet %s: depFamily is ''%s'', expected ''gto''', fpath, S.meta.depFamily);
    [~, ~, info] = get_family_orbit('gto', struct('orientDeg', dp.orientDeg));
    S.meta.depParams = struct('sma_km', info.sma_km, 'ecc', info.ecc, 'orientDeg', info.orientDeg);
    assert(abs(info.orientDeg - S.meta.tauDRO) < 1e-12, 'complete_dep_params:key', ...
        'sheet %s: orientDeg (%.6g) ~= tauDRO key (%.6g)', fpath, info.orientDeg, S.meta.tauDRO);
    tmpPath = [tempname(fileparts(fpath)) '.mat'];
    save(tmpPath, '-struct', 'S');
    Schk = load(tmpPath, 'meta');
    okTmp = isfield(Schk.meta.depParams, 'sma_km') && isfield(Schk.meta.depParams, 'ecc');
    if ~okTmp
        if isfile(tmpPath), delete(tmpPath); end
        error('complete_dep_params:tmpVerify', 'sheet %s: temp file failed verification', fpath);
    end
    movefile(tmpPath, fpath);
    fprintf('complete_dep_params: patched %s -> sma_km=%.6g ecc=%.6g orientDeg=%g\n', ...
        fpath, info.sma_km, info.ecc, info.orientDeg);
end
% unused argument kept for the call-site symmetry with the tulip chain
if isempty(catDir), return, end
end

function n = count_entries(cat_)
% COUNT_ENTRIES  Accepted entries across sheets.  INPUTS: cat_.  OUTPUTS: n.
n = 0;
for ks = 1:numel(cat_.sheets), n = n + nnz(cat_.sheets(ks).has_solution); end
end

function catalog_pictures(sheetMats, picDir)
% CATALOG_PICTURES  One coverage figure per sheet (OK per rung on the
% 12 x 6 torus) and one summary figure (t_f and m_f ranges per rung over
% every sheet).  INPUTS: sheetMats cell; picDir char.  OUTPUTS: none (PNGs).
nS = numel(sheetMats);
rungAll = [];  tfAll = [];  mfAll = [];
for ks = 1:nS
    Q = load(sheetMats{ks}, 'OK', 'TF', 'rungs', 'meta', 'sD', 'sA');
    [nD, nA, nR] = size(Q.OK);
    [~, nm] = fileparts(sheetMats{ks});
    fig = figure('Color', 'w', 'Position', [100 100 1100 320]);
    tl = tiledlayout(fig, 1, nR, 'TileSpacing', 'compact', 'Padding', 'compact');
    for kr = 1:nR
        ax = nexttile(tl);
        imagesc(ax, Q.sA, Q.sD, double(Q.OK(:, :, kr)));
        set(ax, 'YDir', 'normal');  clim(ax, [0 1]);  colormap(ax, [0.92 0.92 0.92; 0.1 0.5 0.2]);
        title(ax, sprintf('%g N: %d/%d', Q.rungs(kr), nnz(Q.OK(:, :, kr)), nD*nA), 'FontSize', 11);
        xlabel(ax, 'arrival phase');  if kr == 1, ylabel(ax, 'departure phase'); end
    end
    title(tl, strrep(nm, '_', ' '), 'FontSize', 13);
    exportgraphics(fig, fullfile(picDir, [nm '_coverage.png']), 'Resolution', 120);
    close(fig);
    for kr = 1:nR
        m = Q.OK(:, :, kr);
        tf = Q.TF(:, :, kr);  tf = tf(m);
        ndp = nd_propulsion(Q.rungs(kr), Q.meta.ispS, Q.meta.m0kg, Q.meta.lStar, Q.meta.tStar);
        mf = 1 - ndp.Tnd*tf/ndp.cnd;
        rungAll = [rungAll; Q.rungs(kr)*ones(numel(tf), 1)];
        tfAll   = [tfAll; tf(:)*Q.meta.tStar/86400];
        mfAll   = [mfAll; mf(:)*Q.meta.m0kg];
    end
end
if isempty(tfAll), return, end
fig = figure('Color', 'w', 'Position', [100 100 900 380]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(tl);
scatter(ax1, rungAll, tfAll, 14, 'filled', 'MarkerFaceAlpha', 0.35, 'DisplayName', 'accepted entries');
set(ax1, 'XScale', 'log');  grid(ax1, 'on');
xlabel(ax1, 'thrust (N)');  ylabel(ax1, 't_f (days)');  title(ax1, 'transfer time per rung');
legend(ax1, 'Location', 'northeast');
ax2 = nexttile(tl);
scatter(ax2, rungAll, mfAll, 14, 'filled', 'MarkerFaceAlpha', 0.35, 'DisplayName', 'accepted entries');
set(ax2, 'XScale', 'log');  grid(ax2, 'on');
xlabel(ax2, 'thrust (N)');  ylabel(ax2, 'm_f (kg)');  title(ax2, 'final mass per rung');
legend(ax2, 'Location', 'southeast');
exportgraphics(fig, fullfile(picDir, 'catalog_ranges.png'), 'Resolution', 120);
close(fig);
end

function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end

function v = pick(c, a, b)
% PICK  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end

function logmsg(f, s)
% LOGMSG  Append-mode file log AND stdout.  INPUTS: f; s.
fprintf('%s\n', s);
if ~isempty(f), fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); end
end
```

Implementer notes: `has_solution` is the coverage authority in the packaged schema (deliverable 9 README); confirm the field name on `cat_.sheets(1)` and adjust `count_entries` if the v2 packager names it differently (state which in the report). `clim` exists in R2026a (`caxis` is its old name).

- [ ] **Step 4: Write the audit**

`audit_gto_dro_entries(catDirIn, sheetMats, floorEarthKm)` is a copy of `GTO_tulip/catalog/audit_gto_entries.m` with: the sheet list taken from `sheetMats` (the chain's explicit list; when called with no arguments it builds the 16 official names under `../results/catalog` and uses `floorEarthKm = 200`); `pick_three` unchanged (low rung = 5 or 7 N; two orientations; two rung classes); the Earth altitude read from the sheet's `_earth.mat` sidecar when present (`EARTHKM(iD, iA, kr)`) and otherwise computed as the tulip audit does, and GATED: `passEarth = minEarthAltKm >= floorEarthKm`, asserted with `audit_gto_dro_entries:earthFloor`; results rows gain `.passEarth`. Header Purpose says it closes the certifier's missing propagated Earth clearance for the audited entries (the sweep closes it for all).

- [ ] **Step 5: Copy the pickers**

```bash
cp orbit_transfer/GTO_tulip/catalog/costate_catalog_pick.m orbit_transfer/GTO_DRO/catalog/
cp orbit_transfer/GTO_tulip/catalog/costate_lib_describe.m orbit_transfer/GTO_DRO/catalog/
```

Edit only the Revision History line of each to add `copied for the GTO -> DRO catalog 11/10/2026`. Read `costate_catalog_pick.m`'s header to confirm it selects sheets by nearest `(tauDRO, arrKey)` and that `arrKey` compares against the sheet's arrival key field; the DRO catalog's arrival key is the DRO period (ND), so check which sheet field the picker compares `arrKey` to (`Np` for tulips; for a DRO arrival the packager stores `arr_params.tau` and possibly `tau_arr`). If the picker cannot resolve a DRO sheet without a code change, make the minimal change in the COPY only, and say so in the report.

- [ ] **Step 6: Run the test**

Expected: PASS in 5 to 10 minutes (the package, one conjugate re-solve per entry, the gates, the audit's three replays on the two entries, two PNGs). The 15 N refuted entry must appear as a `conj` blocker and the zip must be refused.

- [ ] **Step 7: README stage 2 section**

Append to `orbit_transfer/GTO_DRO/README.md`:

```
## Stage 2: the catalog (built 2026-11-10; pilot and fleet status below)

| file | role |
|---|---|
| `catalog/run_gto_dro_catalog.m` | fleet front door: 16 sheets (orientDeg {0,90,180,270} x DRO tau {0.5,1,2,3}), 12 x 6 phasing, rungs [15 12 10 7 5], Earth floor 200 km, two-point top-rung multistart; `planOnly` returns the engine options without solving |
| `catalog/run_gto_dro_batched.sh` | nohup'd OS-timeout driver (matlab-campaign discipline) |
| `lib/earth_clearance_sweep.m` | flown Earth altitude of every accepted entry, sidecar + writeback (demotes below the floor) |
| `lib/seed_top_from_entry.m` | a flown entry as the engine's top-rung seed (reverse warm start, phase continuation) |
| `lib/repair_gto_dro_cells.m` | 'unsolved' (phase continuation from the nearest solved departure phase) and 'basin' (S4 junction test on every top rung; refuted ones re-solved from their 12 N entry; backup + restore) |
| `catalog/build_gto_dro_catalog.m` | the post-fleet chain: census, earth, repairs, package, conj, gates, audit, pictures, deliverable (refused on blockers) |
| `catalog/audit_gto_dro_entries.m` | 3-entry replay: flown miss, tfMin acceptance, Earth floor (gated) |
| `catalog/costate_catalog_pick.m`, `costate_lib_describe.m` | vendored pickers (deliverable convention) |

Run the pilot (sheet 2), then the chain on it, then the fleet:

    cd orbit_transfer/GTO_DRO/catalog
    nohup ./run_gto_dro_batched.sh 2 6 1800 > /dev/null 2>&1 &
    # MATLAB, after the pilot's completion marker:
    build_gto_dro_catalog(struct('expectSheets', 1, 'deliverable', false))
    nohup ./run_gto_dro_batched.sh "[1 3:16]" 6 1800 "" <tf0 from the pilot> > /dev/null 2>&1 &

Engine options added for this stage (costate_common, default off): `seedTop`, `redo`.
```

- [ ] **Step 8: Commit** (controller)

```bash
git add orbit_transfer/GTO_DRO/catalog/build_gto_dro_catalog.m orbit_transfer/GTO_DRO/catalog/audit_gto_dro_entries.m orbit_transfer/GTO_DRO/catalog/costate_catalog_pick.m orbit_transfer/GTO_DRO/catalog/costate_lib_describe.m orbit_transfer/GTO_DRO/tests/test_build_gto_dro_catalog.m orbit_transfer/GTO_DRO/tests/run_tests.m orbit_transfer/GTO_DRO/README.md
git commit -m "GTO_DRO: build_gto_dro_catalog post-fleet chain (census, Earth sweep, repairs, package with the orientDeg sheet key, conjugate + hypothesis sweeps, audit with the Earth floor gated, pictures, deliverable refused on blockers), audit_gto_dro_entries, vendored pickers; test runs the chain on the one-cell fixture and checks the refuted 15 N entry blocks the zip"
```

---

### Task 6: The pilot (controller task: hours of background MATLAB, decisions recorded)

**Files:**
- Modify: `orbit_transfer/GTO_DRO/README.md` ("Pilot findings (stage 2, <date>)")
- Results (not committed): `results/catalog/gto_dro_o000_tau1.mat`, `gto_dro_o000_tau1_deep.mat`, sidecars, `results/pilot_summary.txt`

**Interfaces:** Tasks 3 to 5 as built. This task is executed by the controller (the session that owns the shell), not a subagent, because it waits on hours of background compute.

- [ ] **Step 1: Launch the pilot sheet in the background**

```bash
cd /Users/msc/Desktop/optimal_control/orbit_transfer/GTO_DRO/catalog && nohup ./run_gto_dro_batched.sh 2 6 1800 > /dev/null 2>&1 &
```

Expected cost: 72 cells x (15 N 100 to 520 s + four lower rungs ~60 to 120 s each) ~ 8 to 12 h plus pass B mop-ups on failed cells. Monitor with a stall-aware watcher on `results/catalog/gto_dro_o000_tau1_progress.txt` (log age > 1500 s with no line growth = STALLED; the driver's own `[pass N timed out]` lines mean a solver hang was killed).

- [ ] **Step 2: Launch the deep probe in the shared session once the pilot sheet has at least 12 OK cells at 5 N**

Six cells with 9 rungs in a SEPARATE sheet file, through the front door's engine options (never through the fleet file, whose rungs are fixed):

```matlab
cd(fullfile(getenv('HOME'), 'Desktop/optimal_control/orbit_transfer/GTO_DRO/catalog'))
plan = run_gto_dro_catalog(2, [], [], [15 12 10 7 5 3 2 1.5 1], struct('planOnly', true));
o = plan.optsA;  o.cells = [2 1; 2 4; 5 1; 5 4; 8 1; 8 4];  o.maxCells = 6;  o.resume = true;
o.logFile = fullfile(fileparts(plan.outMat), 'gto_dro_o000_tau1_deep_progress.txt');
P = thrust_ladder_library(fullfile(fileparts(plan.outMat), 'gto_dro_o000_tau1_deep.mat'), o);
squeeze(sum(sum(P.OK, 1), 2)).'       % OK count per rung: where does the ladder stop?
```

Record which rung the ladder stops at per cell and the engine's log line for the first failed rung (defect, iteration cap, flown miss, acceptance).

- [ ] **Step 3: When the pilot's completion marker appears, run the chain on it**

```matlab
blockers = build_gto_dro_catalog(struct('expectSheets', 1, 'deliverable', false));
```

Then compute and record in `results/pilot_summary.txt` and the README:
- coverage per rung on the pilot sheet before and after repairs; the perigee column (iD = 1) outcome at 15 N: how many of its 6 cells the engine retired, how many `'unsolved'` repaired by phase continuation and their t_f;
- the top-rung basin census: `conjPass` from the `'basin'` stage (fraction refuted at 15 N), how many repairs landed on a faster basin, the t_f before/after, and whether the repaired entries pass S4 (`conjAfter`);
- tf0 calibration: median and range of accepted 15 N `t_f` over the sheet; the ruling `tf0Fleet = round(median, 2)` if it differs from 0.30 by more than 0.02, else keep 0.30;
- the Earth sweep: minimum flown altitude over the sheet, demotions (expected 0);
- the deep probe: lowest rung reached per cell;
- the 180-degree question is NOT answered by the pilot; the fleet order puts sheets `[9 10 11 12]` (180 deg) FIRST so the weak corner is measured early: launch the fleet as `"[9:12 1 3:8 13:16]"`.

- [ ] **Step 4: README "Pilot findings" and commit** (controller)

Numbers only as measured; each bullet names the stage that measured it. Commit the README with message `GTO_DRO: pilot findings (sheet 0 deg / tau 1): coverage per rung, perigee column by phase continuation, top-rung basin census and repairs, tf0 calibration, Earth sweep, deep probe`.

---

### Task 7: Astra review before the fleet (controller task)

**Files:**
- Create: `orbit_transfer/GTO_DRO/reviews/astra_stage2_<date>.md` (the review verbatim) and `reviews/astra_stage2_<date>_dispositions.md` (accepted / changed / refused, each with the reason)

- [ ] **Step 1: Assemble the review package** into the scratchpad: stage 1 README findings, the pilot findings, the stage 2 plan (this file) and the three new library functions' headers (`repair_gto_dro_cells`, `earth_clearance_sweep`, the engine's new option block). Keep the package under 15 KB of prose plus the headers; Astra hangs above that in crush, so use the raw API runner `~/.claude/skills/opencode-headless/references/astra_raw.sh`.

- [ ] **Step 2: Ask** for: (1) holes in the basin-repair logic (first-success-wins plus reverse warm start: can a repaired 15 N entry be a worse basin that happens to be faster? what does the S4 verdict after repair prove?), (2) the phase-continuation seed for the perigee column (is a neighbour at phase 1/12 a valid seed for phase 0 when the perigee cell burns 95% of the mass?), (3) the Earth sweep as the closure of the certifier gap (is flown-arc sampling at `fly_transfer`'s resolution enough near perigee?), (4) anything in the fleet recipe that the pilot did not measure.

- [ ] **Step 3: Disposition** every finding in the dispositions file; apply the accepted ones as fix commits before the fleet launches; refused ones carry the reason. Commit both files.

---

### Task 8: The fleet, the chain, the ship (controller task; waits on ~40 h of background compute)

- [ ] **Step 1: Launch the fleet** with the pilot's tf0 and the 180-degree sheets first:

```bash
cd /Users/msc/Desktop/optimal_control/orbit_transfer/GTO_DRO/catalog && nohup ./run_gto_dro_batched.sh "[9:12 1 3:8 13:16]" 6 1800 "" <tf0Fleet> > /dev/null 2>&1 &
```

Monitor as in Task 6; the driver stops on `CATALOG ALL SHEETS COMPLETE` or three dry passes. The shared MATLAB session stays free for the mission-simulation work.

- [ ] **Step 2: Run the chain** when the marker appears: `blockers = build_gto_dro_catalog()`; iterate the repair stages (`maxRepair` in batches) until the blocker list holds only what is measured and accepted (refuted entries that the repair could not improve are shipped with `conj_pass = 0`, the catalog's standing convention since the 0.5 N extension; the deliverable stage then needs `runIn.allowConjBlockers = true`, which the implementer of Task 5 adds as a switch that drops `conj:` blockers from `shipBlockers` and records the count in the log).

- [ ] **Step 3: Record**: README "Fleet findings" (coverage by orientation and DRO period, the 180-degree corner, perigee column, basin census over 16 sheets, deep probe), `orbit_transfer/README.md` row, memory note; copy the zip to `~/Documents/orbit_transfer_bundles/`. Commit and push.

---

## Self-review (done while writing)

- Spec coverage: section 10's pilot -> calibration -> fleet -> deep rungs -> conjugate -> gates -> package -> audit -> pictures -> zip -> batched driver -> pickers are Tasks 3, 6, 8, 6 (deep probe as a 6-cell probe, see Deviations), 5, 5, 5, 5, 5, 5, 3, 5. Section 2's sheet key rule: Task 3 and 5 (asserted). Section 11's risks: tf0 (Task 6 calibration), 180 deg (fleet order, Task 6/8), 3 to 1 N wall (deep probe), Earth constraint (engine floor from stage 1 + Task 2 sweep), 180-degree weakness unmeasured until the fleet (stated). Stage 1's open items: basin selection (Task 4 'basin'), perigee by continuation (Task 4 'unsolved'), tf0 (Task 6), propagated Earth clearance (Task 2 sweep, Task 5 audit gate).
- Deviations from the spec, for the record: (a) stage 2 files live under `GTO_DRO/catalog/` (the spec's folder list was stage 1's); (b) the deep rungs are a six-cell probe in a separate sheet file rather than a fleet pass, because the engine's rung vector is fixed per sheet file and a fleet-wide 9-rung ladder would charge up to 600 s per cell at the wall; the fleet rung vector stays [15 12 10 7 5] as the GTO -> tulip catalog shipped; (c) two engine options (seedTop, redo) are added, default off, on the stage 1 precedent; (d) an Earth sweep stage is inserted between the fleet and the package.
- Placeholder scan: one deliberate implementer-resolved expression in Task 5's test (`derive_period_days_dep`, flagged in the step text with the rule for replacing it) and one in Task 5 step 5 (the picker's DRO arrival key, flagged with the rule). Both are reads of existing code the implementer must do; the plan states what to do in each outcome.
- Type consistency: `seedTop` struct `{X, U, tf, from}` (Tasks 1, 2, 4); `repair_gto_dro_cells(sheetMat, mode, opts)` returning `R` with the listed fields (Tasks 4, 5); `earth_clearance_sweep(sheetMat, opts)` returning `S` (Tasks 2, 5); `run_gto_dro_catalog(..., optsIn)` returning `plan` with `.optsA/.optsB` (Tasks 3, 4, 5, 6); sheet tags `gto_dro_o%03d_tau%g` everywhere.
- Review Focus: each line names its pinning test (Tasks 4, 1, 3, 2, 3).
