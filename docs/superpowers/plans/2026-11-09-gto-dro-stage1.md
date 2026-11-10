# GTO -> DRO minimum-time, stage 1 (single-transfer study) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One GTO -> DRO minimum-time low-thrust transfer in the Earth-Moon CR3BP, solved through the existing shared ladder engine from user-set parameters, independently verified, gated (PASS/FAIL/N.A.) and drawn, by a study script Mike runs with no arguments.

**Architecture:** A new campaign folder `orbit_transfer/GTO_DRO/` holding a thin single-cell driver (`lib/gto_dro_entry.m`, a copy of the GTO -> tulip catalog's `gto_entry.m` with the arrival family switched to the DRO), one new geometric instrument (`lib/earth_clearance.m`), and the study script `gto_dro_study.m` (SDA-style `%% knobs` block, numbered sections, stable gate IDs, light-theme figures, `.mat` with run metadata). No solver code is written and no `costate_common` file is changed.

**Tech Stack:** MATLAB R2026a (shared session through the MCP tools), CasADi 3.7.0 + IPOPT (inside `casadi_mintime_dro`), pumpkyn + pumpkynPie (`get_family_orbit`, `tfMin`, `tfMinProp`, `fromPCI`, `toPCI`, `lagrangePts`, `jacobi`, `minDeltaV`), the costate_common engine (`thrust_ladder_library`, `ms_tfmin`, `seed_from_z8`, `fly_transfer`, `pmp_pointwise_checks`, `mintime_hypothesis_gates`, `lift_margin`, `conj_spectrum`, `nd_propulsion`, `phase_state`), oclib (`oc.*`).

**Spec:** `docs/superpowers/specs/2026-11-09-gto-dro-mintime-design.md`

## Global Constraints

- Catalog regime only: defaults m0 = 150 kg, Isp = 1710 s, top rung 15 N; `thrustN` outside [1, 15] N is refused by name (`gto_dro_entry:thrust`).
- Standard rungs `[15 12 10 7 5 3 2 1.5 1]`; the ladder walks from 15 N DOWN to `thrustN`; higher rungs are stepping stones, the returned entry is the bottom rung. No deep rung is ever cold-solved.
- Warm recipe defaults: `thrLock = true`, `tf0 = 0.30` (top rung only), `maxCpuSec = 600`, `maxIter = 6000`, `N = 400`, `floorKm = 500`, `gateKm = 100`, `accTol = 1e-6`.
- Pipeline unchanged: direct solve -> fly -> `ms_tfmin` -> pumpkyn `tfMin` acceptance. The study script runs the acceptance again in the open (X1).
- No file under `orbit_transfer/costate_common/`, `oclib/`, or any other campaign folder is modified. No DRO_tulip or GTO_tulip folder goes on the path (`plot_transfer_3d` and `certify_root` live there and are NOT used).
- MATLAB house style: header quartet (`%% Purpose:` / `%% Inputs:` / `%% Outputs:` / `%% Revision History:` with `%  M. Casey                                                   (c) 11/09/2026` and `%  Copyright Coorbital Inc.`, then `%% ------------------------ Begin Code Sequence ---------------------------`), `if nargin == 0` self-demo in every function, no `i`/`j` as variables, no `norm` (use `sqrt(sum(.^2))` or `max(abs())`), no `%#ok`, no emojis, `=` aligned in parameter blocks, ASCII only in README.
- Study-script discipline: numbered sections, every condition computed IN the script with a stable ID and PASS / FAIL / NOT APPLICABLE / UNRESOLVED, one tolerance block, a final verdict table and an assertion on the code gates.
- Figures: `figure('Theme', 'light', 'Color', 'w')`, `tiledlayout`, every call with an explicit axes handle, `DisplayName` on every line, fonts >= 12 pt, saved PNG at 150 dpi; the shared session defaults axes to dark, so `set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')`.
- MATLAB is run ONLY through the MCP tools (`mcp__matlab__run_matlab_file`, `evaluate_matlab_code`, `run_matlab_test_file`, `check_matlab_code`) in the shared R2026a session. Never `matlab -batch` from the shell.
- Results are git-ignored (`orbit_transfer/GTO_DRO/.gitignore` with `results/`); `*.mat` is ignored repo-wide already.
- Subagents never commit. The controller commits after each clean review; commit messages end `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`; git through `/Library/Developer/CommandLineTools/usr/bin/git`; never stage the unrelated `missiles/` or `papers/` deletions or `.DS_Store` changes.
- Ruling recorded here (spec section 7, item 6): the direct solution's control is not saved by `thrust_ladder_library`, so the "second flight" of the direct control is NOT re-flown; the script reports the engine's own G1b flown miss for that rung from `P.FLYKM` beside the `fly_transfer` flight of the accepted z8. Cost if wrong: one number is a read-back instead of a re-measurement; the z8 flight (the one every gate uses) is still independent of the solve.

## Review Focus

1. A `thrustN` that is not a standard rung (e.g. 4 N): the ladder must be `[15 12 10 7 5 4]`, the returned entry at 4 N, and `P.rungs(end) == 4`. Tested in Task 2 (`test_gto_dro_entry`, plan-only mode).
2. `study_opts` overrides must REACH the solve (the advertised-but-ignored-option bug): the test sets `thrustN = 12` through `study_opts` and checks the printed ladder ends at 12 N. Tested in Task 5.
3. The Earth-clearance gate must be able to FAIL: an arc that crosses the Earth's radius must report FAIL, and the study rerun with `tol.earthKm = 1e6` must print `N7a ... FAIL`. Tested in Tasks 1 and 5.
4. `out.reuse` (loading a saved result instead of solving) must refuse a file whose knobs (orientDeg, tau, phases, thrust, m0, Isp) differ from the current ones, so a stale result is never gated as if it were fresh. Tested in Task 5.
5. `tfMin` acceptance without an fsolve licence: X1 must print NOT CHECKED (not PASS) when `tfMin` throws, and the fallback `solverType = 'lsqnonlin'` knob must be honoured. Tested in Task 4 (unit-level, by forcing the knob).

---

### Task 1: Scaffold, paths and the Earth-clearance instrument

**Files:**
- Create: `orbit_transfer/GTO_DRO/.gitignore`
- Create: `orbit_transfer/GTO_DRO/setup_gto_dro_paths.m`
- Create: `orbit_transfer/GTO_DRO/lib/earth_clearance.m`
- Create: `orbit_transfer/GTO_DRO/tests/test_earth_clearance.m`
- Create: `orbit_transfer/GTO_DRO/tests/run_tests.m`

**Interfaces:**
- Produces: `setup_gto_dro_paths()` (no args; adds GTO_DRO, GTO_DRO/lib, costate_common, oclib parent, casadi; pumpkynPie startup if `pumpkyn.cr3bp.tfMin` is absent).
- Produces: `[hMinKm, tAtMin, pass, out] = earth_clearance(tau, y, floorKm, p)`; `p.muStar`, `p.lStar` required; optional `p.center` [1 x 3] ND (default Earth `[-muStar 0 0]`), `p.radiusKm` (default 6378.137); `out.hKm` [M x 1] altitude history, `out.center`, `out.radiusKm`.
- Produces: `ok = run_tests()` running `{'test_earth_clearance'}` for now; later tasks append names.

- [ ] **Step 1: Write `.gitignore`**

```
results/
```

- [ ] **Step 2: Write `setup_gto_dro_paths.m`**

```matlab
function setup_gto_dro_paths()
%% Purpose:
%
%   GTO_DRO campaign paths: this dir + lib/ + costate_common (family
%   provider get_family_orbit, the shared ladder engine
%   thrust_ladder_library / casadi_mintime_dro / certify_dro_mintime, the
%   verification instruments fly_transfer / pmp_pointwise_checks /
%   mintime_hypothesis_gates / lift_margin / conj_spectrum) + the top-level
%   oclib/ dir (parent of +oc; certify_dro_mintime -> dro_residual calls
%   oc.local_residual) + casadi + pumpkyn/pumpkynPie through pumpkynPie's
%   own startup if they are not already on the path. No DRO_tulip or
%   GTO_tulip folder is added: the engine left DRO_tulip on 2026-09-16 and
%   GTO_tulip/catalog's helpers are tulip-keyed.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none
%
%% Revision History:
%  M. Casey                                                   (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
ot   = fileparts(here);                        % orbit_transfer/
addpath(here, fullfile(here, 'lib'), fullfile(ot, 'costate_common'), ...
        fullfile(fileparts(ot), 'oclib'), ...
        fullfile(getenv('HOME'), 'casadi-3.7.0'));
if isempty(which('pumpkyn.cr3bp.tfMin'))
    pp = '/Users/msc/Desktop/proj7/external/pumpkynPie';
    od = cd(pp); startup(); cd(od);
end
end
```

- [ ] **Step 3: Write the failing test `tests/test_earth_clearance.m`**

```matlab
function ok = test_earth_clearance()
%% Purpose:
%
%   earth_clearance on three arcs: (1) a circular 6678.137 km Earth orbit
%   converted with pumpkyn fromPCI reads 300 km altitude to 1e-6 km and
%   PASSES a 200 km floor; (2) the same orbit with a 1000 km floor FAILS
%   (the gate can fail); (3) a straight radial arc that crosses the
%   Earth's radius reports a negative minimum altitude and FAILS; (4) the
%   Moon-centre mode on the tau = 1 DRO reproduces get_family_orbit's
%   perilune to 1 km.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every check passed
%
%% Revision History:
%  M. Casey                                                   (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
muStar = 0.012150585609624;  lStar = 389703.264829278;  tStar = 382981.289129055;
p = struct('muStar', muStar, 'lStar', lStar);

% (1) circular Earth orbit, radius 6678.137 km, 64 samples, converted at t = 0
rC   = 6678.137;
muE  = 6.67384e-20*(1 - muStar)*(5.9736E24 + 7.35E22);
vC   = sqrt(muE/rC);
th   = linspace(0, 2*pi, 64).';
rvKm = [rC*cos(th), rC*sin(th), zeros(64,1), -vC*sin(th), vC*cos(th), zeros(64,1)];
rvND = zeros(64, 6);
for kk = 1:64
    rvND(kk,:) = pumpkyn.cr3bp.fromPCI(0, rvKm(kk,:), muStar, tStar, lStar, 1);
end
tau = linspace(0, 1, 64).';
[hMin, tAt, pass1] = earth_clearance(tau, rvND, 200, p);
ok = check(ok, abs(hMin - 300) < 1e-6, sprintf('circular 6678.137 km orbit reads %.9f km altitude (300 expected)', hMin));
ok = check(ok, pass1, '300 km orbit passes a 200 km floor');
ok = check(ok, tAt >= 0 && tAt <= 1, 'time of minimum lies in the arc');

% (2) the gate can fail
[~, ~, pass2] = earth_clearance(tau, rvND, 1000, p);
ok = check(ok, ~pass2, '300 km orbit FAILS a 1000 km floor');

% (3) a radial arc through the Earth's centre
xs  = linspace(0.05, -0.05, 21).';                       % ND, crosses x = -muStar
arc = [xs - muStar, zeros(21, 5)];
[hMin3, ~, pass3] = earth_clearance(linspace(0, 1, 21).', arc, 200, p);
ok = check(ok, hMin3 < 0 && ~pass3, sprintf('radial arc through the centre reads %.1f km and FAILS', hMin3));

% (4) Moon-centre mode on the tau = 1 DRO vs get_family_orbit's own geometry
[tA, rvA] = get_family_orbit('dro', struct('tau', 1));
pM = struct('muStar', muStar, 'lStar', lStar, 'center', [1 - muStar, 0, 0], 'radiusKm', 1737.4);
[hMoon, ~, passM, outM] = earth_clearance(tA, rvA, 500, pM);
dMoon = lStar*sqrt(sum((rvA(:,1:3) - [1 - muStar, 0, 0]).^2, 2));
ok = check(ok, abs(hMoon - (min(dMoon) - 1737.4)) < 1, sprintf('Moon mode perilune altitude %.1f km matches the direct computation', hMoon));
ok = check(ok, passM, 'tau = 1 DRO clears a 500 km lunar floor');
ok = check(ok, numel(outM.hKm) == size(rvA, 1), 'altitude history has one entry per sample');

fprintf('test_earth_clearance: %s\n', tern(ok, 'PASS', 'FAIL'));
end

function ok = check(ok, c, msg)
% CHECK  Print one line and fold the verdict.  INPUTS: ok; c; msg.  OUTPUTS: ok.
if c, s = 'ok  '; else, s = 'FAIL'; end
fprintf('  %s %s\n', s, msg);
ok = ok && c;
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end
```

- [ ] **Step 4: Run the test to see it fail**

Use `mcp__matlab__evaluate_matlab_code` with
```matlab
cd('/Users/msc/Desktop/optimal_control/orbit_transfer/GTO_DRO/tests'); test_earth_clearance
```
Expected: error `Unrecognized function or variable 'earth_clearance'`.

- [ ] **Step 5: Write `lib/earth_clearance.m`**

```matlab
function [hMinKm, tAtMin, pass, out] = earth_clearance(tau, y, floorKm, p)
%% Purpose:
%
%   Altitude of a flown arc above a primary, in km, and whether it stays
%   above a floor. Pure geometry on the caller's samples: no propagation,
%   no interpolation between samples (the number is the SAMPLE minimum,
%   and says so). Default primary is the EARTH at [-muStar 0 0] with the
%   equatorial radius 6378.137 km -- the check the shared ladder engine
%   lacks for an Earth departure (it enforces a LUNAR floor only). With
%   p.center = [1-muStar 0 0], p.radiusKm = 1737.4 it is the lunar
%   altitude, so a study script measures both clearances with ONE
%   instrument.
%
%   Self-demo (nargin == 0): a circular 6678.137 km Earth orbit converted
%   with pumpkyn fromPCI at t = 0 reads 300 km.
%
%% Inputs:
%
%  tau                      [M x 1]                 sample times, ND
%
%  y                        [M x >= 6]              rotating-frame
%                                                   barycentric ND states,
%                                                   one row per time;
%                                                   columns 1:3 are used
%
%  floorKm                  double                  altitude floor, km
%
%  p                        struct                  .muStar, .lStar
%                                                   (required); .center
%                                                   [1 x 3] ND
%                                                   [-muStar 0 0]; .radiusKm
%                                                   [6378.137]
%
%% Outputs:
%
%  hMinKm                   double                  minimum altitude over the
%                                                   samples, km (negative
%                                                   when the arc crosses the
%                                                   surface)
%
%  tAtMin                   double                  tau at that sample
%
%  pass                     logical                 hMinKm >= floorKm
%
%  out                      struct                  .hKm [M x 1] altitude
%                                                   history, .center,
%                                                   .radiusKm, .iMin
%
%% Revision History:
%  M. Casey                                                   (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    muStar = 0.012150585609624;  lStar = 389703.264829278;  tStar = 382981.289129055;
    muE = 6.67384e-20*(1 - muStar)*(5.9736E24 + 7.35E22);
    rC  = 6678.137;  vC = sqrt(muE/rC);  th = linspace(0, 2*pi, 64).';
    rvKm = [rC*cos(th), rC*sin(th), zeros(64,1), -vC*sin(th), vC*cos(th), zeros(64,1)];
    y = zeros(64, 6);
    for kk = 1:64
        y(kk,:) = pumpkyn.cr3bp.fromPCI(0, rvKm(kk,:), muStar, tStar, lStar, 1);
    end
    tau = linspace(0, 1, 64).';  floorKm = 200;
    p = struct('muStar', muStar, 'lStar', lStar);
end
if nargin < 4 || ~isfield(p, 'muStar') || ~isfield(p, 'lStar')
    error('earth_clearance:params', 'p.muStar and p.lStar are required');
end
center   = fieldd(p, 'center',   [-p.muStar, 0, 0]);
radiusKm = fieldd(p, 'radiusKm', 6378.137);
tau = tau(:);
if size(y, 2) < 6 || size(y, 1) ~= numel(tau)
    error('earth_clearance:shape', 'y must be [numel(tau) x >= 6]; got [%d x %d] for %d times', ...
          size(y, 1), size(y, 2), numel(tau));
end

dKm  = p.lStar*sqrt(sum((y(:, 1:3) - center).^2, 2));     % distance to the centre, km
hKm  = dKm - radiusKm;
[hMinKm, iMin] = min(hKm);
tAtMin = tau(iMin);
pass   = hMinKm >= floorKm;
out    = struct('hKm', hKm, 'center', center, 'radiusKm', radiusKm, 'iMin', iMin);

if nargin == 0 || nargout == 0
    fprintf('earth_clearance: min altitude %.6f km at tau = %.4f (floor %.0f km) -> %s\n', ...
            hMinKm, tAtMin, floorKm, tern(pass, 'PASS', 'FAIL'));
end
end

function v = fieldd(s, f, d_)
% FIELDD  s.(f) if present else d_.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end
```

- [ ] **Step 6: Run the test to see it pass**

Same MCP call as step 4. Expected: eight `ok` lines and `test_earth_clearance: PASS`. Also run `mcp__matlab__check_matlab_code` on `lib/earth_clearance.m`; fix any analyzer message that is not a style-only note.

- [ ] **Step 7: Write `tests/run_tests.m`**

```matlab
function ok = run_tests()
%% Purpose:
%
%   Every test in this folder, in dependency order, with ONE logical
%   verdict and a throw on failure so a scripted run returns a real exit
%   code. The study-script test solves one transfer at 15 N and takes
%   several minutes; the others are seconds.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every test passed (the
%                                                   function also THROWS
%                                                   on failure)
%
%% Revision History:
%  M. Casey                                                   (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(here);

names = {'test_earth_clearance'};
res   = false(size(names));
secs  = zeros(size(names));
for kk = 1:numel(names)
    fprintf('\n######## %s\n', names{kk});
    t0 = tic;
    try
        res(kk) = feval(names{kk});
    catch err
        fprintf('  ERROR %s\n', err.message);
        res(kk) = false;
    end
    secs(kk) = toc(t0);
end

fprintf('\n%-28s %-6s %8s\n', 'test', 'result', 'seconds');
for kk = 1:numel(names)
    if res(kk), s = 'PASS'; else, s = 'FAIL'; end
    fprintf('%-28s %-6s %8.1f\n', names{kk}, s, secs(kk));
end
ok = all(res);
fprintf('\nGTO_DRO: %d of %d passed\n', nnz(res), numel(res));
if ~ok
    error('GTO_DRO:run_tests:fail', '%d test(s) failed', nnz(~res));
end
end
```

- [ ] **Step 8: Run `run_tests` through MCP; expected `GTO_DRO: 1 of 1 passed`.**

- [ ] **Step 9: Report** (the controller commits): files created, test output lines.

---

### Task 2: The single-cell driver `lib/gto_dro_entry.m`

**Files:**
- Create: `orbit_transfer/GTO_DRO/lib/gto_dro_entry.m`
- Create: `orbit_transfer/GTO_DRO/tests/test_gto_dro_entry.m`
- Modify: `orbit_transfer/GTO_DRO/tests/run_tests.m` (append the test name)

**Interfaces:**
- Consumes: `thrust_ladder_library(outMat, opts)` with `depFamily 'gto'`, `arrFamily 'dro'`; `nd_propulsion(thrustN, ispS, m0kg, lStar, tStar)` -> `.cnd .Tnd`; `get_family_orbit`, `phase_state`, `seed_from_z8`, `ms_tfmin`.
- Produces: `[entry, gates, P] = gto_dro_entry(orientDeg, tauDRO, depFrac, arrFrac, thrustN, opts)` exactly as spec section 5, plus `opts.planOnly` [false] (return the planned ladder in `entry.rungsWalked` and everything else empty, no solve) and `opts.outMat`, `opts.keepScratch`, `opts.gtoParams`, `opts.K` [24].

- [ ] **Step 1: Write the failing test `tests/test_gto_dro_entry.m`** (fast checks only; the slow solve is exercised by Task 5's study test)

```matlab
function ok = test_gto_dro_entry()
%% Purpose:
%
%   gto_dro_entry's PLANNING and REFUSALS, without a solve: (1) thrustN
%   outside [1, 15] is refused by name; (2) the ladder plan is
%   [15] for 15 N, [15 12 10 7] for 7 N, all nine rungs for 1 N and
%   [15 12 10 7 5 4] for 4 N (a non-standard rung is appended, never
%   cold-solved); (3) the GTO parameter passthrough and the DRO period
%   land in the engine options the driver builds (planOnly returns them in
%   entry.engineOpts). The solve itself is exercised end to end by
%   test_gto_dro_study.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every check passed
%
%% Revision History:
%  M. Casey                                                   (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
po = struct('planOnly', true);

% (1) refusals
for bad = [0.5, 16, -1]
    threw = false;
    try
        gto_dro_entry(0, 1, 0, 0, bad, po);
    catch err
        threw = strcmp(err.identifier, 'gto_dro_entry:thrust');
    end
    ok = check(ok, threw, sprintf('thrustN = %g refused with gto_dro_entry:thrust', bad));
end

% (2) ladder plans
plans = {15, [15];  7, [15 12 10 7];  1, [15 12 10 7 5 3 2 1.5 1];  4, [15 12 10 7 5 4]};
for kk = 1:size(plans, 1)
    e = gto_dro_entry(0, 1, 0, 0, plans{kk, 1}, po);
    ok = check(ok, isequal(e.rungsWalked, plans{kk, 2}), ...
               sprintf('thrustN = %g plans [%s]', plans{kk, 1}, num2str(e.rungsWalked)));
end

% (3) passthrough: GTO params, DRO period, phases, recipe reach the engine options
e = gto_dro_entry(90, 2, 0.25, 0.5, 10, struct('planOnly', true, ...
        'gtoParams', struct('sma_km', 25000, 'ecc', 0.7), 'tf0', 0.4, 'm0kg', 120, 'ispS', 2000));
eo = e.engineOpts;
ok = check(ok, strcmp(eo.depFamily, 'gto') && strcmp(eo.arrFamily, 'dro'), 'families gto -> dro');
ok = check(ok, eo.depParams.orientDeg == 90 && eo.depParams.sma_km == 25000 && abs(eo.depParams.ecc - 0.7) < 1e-15, 'GTO params reach depParams');
ok = check(ok, eo.arrParams.tau == 2, 'DRO period reaches arrParams.tau');
ok = check(ok, eo.sD0 == 0.25 && eo.sA0 == 0.5 && eo.nD == 1 && eo.nA == 1, 'phases reach a 1x1 grid');
ok = check(ok, eo.tf0 == 0.4 && eo.m0kg == 120 && eo.ispS == 2000 && eo.thrLock, 'recipe, mass and Isp reach the engine');
ok = check(ok, eo.tauDRO == 90, 'sheet key carries orientDeg (the GTO catalog rule)');
ok = check(ok, isequal(eo.rungs, [15 12 10]), 'planned rungs are the engine rungs');

fprintf('test_gto_dro_entry: %s\n', tern(ok, 'PASS', 'FAIL'));
end

function ok = check(ok, c, msg)
% CHECK  Print one line and fold the verdict.  INPUTS: ok; c; msg.  OUTPUTS: ok.
if c, s = 'ok  '; else, s = 'FAIL'; end
fprintf('  %s %s\n', s, msg);
ok = ok && c;
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end
```

- [ ] **Step 2: Run it (MCP) and see it fail with `Unrecognized function 'gto_dro_entry'`.**

- [ ] **Step 3: Write `lib/gto_dro_entry.m`**

```matlab
function [entry, gates, P] = gto_dro_entry(orientDeg, tauDRO, depFrac, arrFrac, thrustN, opts)
%% Purpose:
%
%   THE SINGLE-ENTRY DRIVER for GTO -> DRO minimum-time transfers: solve
%   ONE (departure phase x arrival phase) cell at a requested thrust level
%   through the SAME shared engine every costate catalog uses
%   (costate_common/thrust_ladder_library, family-agnostic via
%   get_family_orbit) -- never a bespoke solver. A copy of
%   GTO_tulip/catalog/gto_entry.m (2026-08-26) with the arrival family
%   switched from the tulip to a DRO of period tauDRO; the sheet-writing
%   path is dropped (stage 2's chain script owns sheets).
%
%   LADDER DIRECTION: solves from 15 N DOWN to thrustN; higher standard
%   rungs are warm-start stepping stones (never returned) and the entry
%   AT thrustN is returned. A deep rung is never cold-solved: the pilot
%   of the GTO -> tulip campaign measured a cold deep start as exactly
%   the non-convergent case.
%
%   WARM RECIPE (inherited from run_gto_catalog's pass A, 2026-08-26
%   diagnostic A; a DRO-specific recalibration is stage 2's pilot):
%   thrLock = true, tf0 = 0.30 ND (seeds the TOP rung only),
%   maxCpuSec = 600, maxIter = 6000. The engine's cold defaults (free
%   throttle, tf0 = 4.0) measured 67/72 GTO cells unreachable; this driver
%   never falls back to them.
%
%   Self-demo (nargin == 0): orient 0 deg, tau 1, PERIGEE departure
%   (phase 0), arrival phase 0, 15 N -- one solve. Perigee, not apogee:
%   the GTO campaign measured an apogee cold start at 15 N landing in a
%   bad local minimum (11.95e6 km flown miss).
%
%% Inputs:
%
%  orientDeg                double                  GTO locus orientation,
%                                                   deg: Earth->Moon line to
%                                                   the perigee direction,
%                                                   rotation sense
%
%  tauDRO                   double                  DRO period, ND (0.5, 1,
%                                                   2, 3 are the catalog
%                                                   periods; 1 = 4.43 d)
%
%  depFrac                  double                  departure phase: TIME
%                                                   fraction since perigee
%                                                   over one Kepler period
%                                                   (0 perigee, 0.5 apogee)
%
%  arrFrac                  double                  arrival phase: time
%                                                   fraction along the DRO
%                                                   from its +x crossing
%
%  thrustN                  double                  thrust, N, in [1, 15]:
%                                                   the BOTTOM (returned)
%                                                   rung
%
%  opts                     struct (optional)       .m0kg [150], .ispS [1710],
%                                                   .thrLock [true], .tf0
%                                                   [0.30], .maxCpuSec [600],
%                                                   .maxIter [6000], .N [400],
%                                                   .floorKm [500] lunar,
%                                                   .gateKm [100], .accTol
%                                                   [1e-6], .gtoParams
%                                                   (struct, .sma_km .ecc;
%                                                   default = get_family_
%                                                   orbit's GTO), .conjTest
%                                                   [false] re-solve ms_tfmin
%                                                   at the accepted z8 with
%                                                   the conjugate test, .K
%                                                   [24] segments for it,
%                                                   .outMat scratch path
%                                                   [results/scratch/...],
%                                                   .keepScratch [false],
%                                                   .planOnly [false] return
%                                                   the plan, no solve
%
%% Outputs:
%
%  entry                    struct                  .ok, .z8 [8 x 1], .tf_nd,
%                                                   .tf_d, .mf_kg, .dV_kms,
%                                                   .thrustN, .rungsWalked
%                                                   [1 x nR], .tfPerRung,
%                                                   .okPerRung, .rv0, .rvf
%                                                   [1 x 6], .depFrac,
%                                                   .arrFrac, .orientDeg,
%                                                   .tauDRO, .m0kg, .ispS,
%                                                   .engineOpts (the options
%                                                   handed to the engine)
%
%  gates                    struct                  .msNormR, .flownKm,
%                                                   .accDz, .wallS, .conjPass
%                                                   (NaN unless .conjTest
%                                                   reproduced z8),
%                                                   .ratioScreen [1 x nR-1],
%                                                   .ratioPass
%
%  P                        struct                  thrust_ladder_library's
%                                                   output for the cell
%                                                   (every rung); empty when
%                                                   .planOnly
%
%% Revision History:
%  M. Casey                                                   (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    orientDeg = 0;  tauDRO = 1;  depFrac = 0;  arrFrac = 0;  thrustN = 15;
end
if nargin < 6, opts = struct(); end
m0kg       = fdef(opts, 'm0kg',       150);
ispS       = fdef(opts, 'ispS',       1710);
thrLock    = fdef(opts, 'thrLock',    true);    % warm recipe
tf0Warm    = fdef(opts, 'tf0',        0.30);    % warm recipe (top rung only)
maxCpuSec  = fdef(opts, 'maxCpuSec',  600);     % warm recipe
maxIterOpt = fdef(opts, 'maxIter',    6000);    % warm recipe (engine default 3000 fell short)
N          = fdef(opts, 'N',          400);
floorKm    = fdef(opts, 'floorKm',    500);     % LUNAR floor (the engine's only one)
gateKm     = fdef(opts, 'gateKm',     100);
accTol     = fdef(opts, 'accTol',     1e-6);
gtoParams  = fdef(opts, 'gtoParams',  struct());
conjTest   = fdef(opts, 'conjTest',   false);
Kseg       = fdef(opts, 'K',          24);
keepScratch = fdef(opts, 'keepScratch', false);
planOnly   = fdef(opts, 'planOnly',   false);

if ~(isscalar(thrustN) && isreal(thrustN) && thrustN >= 1 && thrustN <= 15)
    error('gto_dro_entry:thrust', ...
          ['thrustN = %g N is outside [1, 15]: above 15 N there is no standard rung to anchor on; ' ...
           'below 1 N is the many-rev regime that ships junction states, not bare z8 (spec section 9)'], thrustN);
end

here   = fileparts(fileparts(mfilename('fullpath')));      % GTO_DRO/
run(fullfile(here, 'setup_gto_dro_paths.m'));

%% Ladder from 15 N DOWN to thrustN (stepping stones above, thrustN last)
rungsStd = [15 12 10 7 5 3 2 1.5 1];
rungs    = [rungsStd(rungsStd > thrustN), thrustN];
krEntry  = numel(rungs);

muStar = 0.012150585609624;
lStar  = 389703.264829278;
tStar  = 382981.289129055;

depParams = struct('orientDeg', orientDeg);
for f = {'sma_km', 'ecc'}
    if isfield(gtoParams, f{1}), depParams.(f{1}) = gtoParams.(f{1}); end
end

scratchDir = fullfile(here, 'results', 'scratch');
outMat  = fdef(opts, 'outMat', fullfile(scratchDir, sprintf('gto_dro_o%03d_tau%g_%gN.mat', orientDeg, tauDRO, thrustN)));
logFile = [outMat(1:end-4) '.log'];

optsEngine = struct( ...
    'rungs',     rungs, ...
    'ispS',      ispS, ...
    'm0kg',      m0kg, ...
    'N',         N, ...
    'floorKm',   floorKm, ...
    'maxIter',   maxIterOpt, ...
    'gateKm',    gateKm, ...
    'accTol',    accTol, ...
    'nD',        1, ...
    'nA',        1, ...
    'sD0',       depFrac, ...
    'sA0',       arrFrac, ...
    'muStar',    muStar, 'lStar', lStar, 'tStar', tStar, ...
    'depFamily', 'gto', ...
    'depParams', depParams, ...
    'arrFamily', 'dro', ...
    'arrParams', struct('tau', tauDRO), ...
    'tauDRO',    orientDeg, ...        % SHEET KEY RULE (GTO catalog): the key is the orientation
    'thrLock',   thrLock, 'tf0', tf0Warm, 'maxCpuSec', maxCpuSec, ...
    'resume',    false, ...
    'logFile',   logFile);

entry = struct('ok', false, 'z8', nan(8, 1), 'tf_nd', NaN, 'tf_d', NaN, 'mf_kg', NaN, 'dV_kms', NaN, ...
    'thrustN', thrustN, 'rungsWalked', rungs, 'tfPerRung', nan(1, numel(rungs)), 'okPerRung', false(1, numel(rungs)), ...
    'rv0', nan(1, 6), 'rvf', nan(1, 6), 'depFrac', depFrac, 'arrFrac', arrFrac, ...
    'orientDeg', orientDeg, 'tauDRO', tauDRO, 'm0kg', m0kg, 'ispS', ispS, 'engineOpts', optsEngine);
gates = struct('msNormR', NaN, 'flownKm', NaN, 'accDz', NaN, 'wallS', NaN, 'conjPass', NaN, ...
    'ratioScreen', nan(1, max(numel(rungs) - 1, 0)), 'ratioPass', true);
P = struct();
if planOnly, return; end

%% Endpoints by THE shared endpoint rule (costate_common/phase_state)
[tD, rvD] = get_family_orbit('gto', depParams);
[tA, rvA] = get_family_orbit('dro', struct('tau', tauDRO));
depAt = phase_state(tD, rvD);   arrAt = phase_state(tA, rvA);
entry.rv0 = depAt(depFrac).';
entry.rvf = arrAt(arrFrac).';

%% Solve
if ~isfolder(fileparts(outMat)), mkdir(fileparts(outMat)); end
if isfile(outMat),  delete(outMat);  end
if isfile(logFile), delete(logFile); end
thrust_ladder_library(outMat, optsEngine);
P = load(outMat);
assert(isequal(P.rungs(krEntry), thrustN), 'gto_dro_entry:rungIndex', ...
    'P.rungs(%d) = %.4g does not match the requested thrustN = %.4g', krEntry, P.rungs(krEntry), thrustN);

entry.ok        = logical(P.OK(1, 1, krEntry));
entry.z8        = P.Z8(:, 1, 1, krEntry);
entry.tf_nd     = P.TF(1, 1, krEntry);
entry.tf_d      = entry.tf_nd*tStar/86400;
entry.tfPerRung = squeeze(P.TF(1, 1, :)).';
entry.okPerRung = logical(squeeze(P.OK(1, 1, :)).');
gates.msNormR   = P.RES(1, 1, krEntry);
gates.flownKm   = P.FLYKM(1, 1, krEntry);
gates.accDz     = P.ACCDZ(1, 1, krEntry);
gates.wallS     = P.WALL(1, 1, krEntry);

% rung-to-rung screen on STEPS only (the top rung's guess has no branch behind it)
if numel(rungs) > 1
    gates.ratioScreen = entry.tfPerRung(2:end)./entry.tfPerRung(1:end-1);
    gates.ratioPass   = all(gates.ratioScreen > 1/3 & gates.ratioScreen < 3 | ~entry.okPerRung(2:end));
end

if entry.ok
    ndp = nd_propulsion(thrustN, ispS, m0kg, lStar, tStar);
    cnd = ndp.cnd;   Tnd = ndp.Tnd;
    mfFrac       = 1 - Tnd*entry.tf_nd/cnd;                 % all-burn: m(t) = 1 - T t / c
    entry.mf_kg  = m0kg*mfFrac;
    entry.dV_kms = cnd*log(1/mfFrac)*lStar/tStar;

    if conjTest
        seed = seed_from_z8(entry.z8, entry.rv0(1:6), Kseg, Tnd, cnd, muStar);
        [z, info] = ms_tfmin(entry.rv0(1:6), entry.rvf(1:6), seed, Tnd, cnd, muStar, ...
                             struct('conjTest', true, 'wallSec', 120));
        dz = max(abs(z - entry.z8));
        if info.converged && dz < 1e-6
            gates.conjPass = double(info.conj.pass);
        else
            fprintf('gto_dro_entry: conjTest re-solve did not reproduce z8 (|dz| = %.1e) -- conjPass left NaN\n', dz);
        end
    end
end

if ~keepScratch
    if isfile(outMat),  delete(outMat);  end
    if isfile(logFile), delete(logFile); end
end

if nargin == 0 || nargout == 0
    fprintf('gto_dro_entry: orient %g deg, DRO tau %g, phases (%.3f, %.3f), %g N -> ok = %d\n', ...
            orientDeg, tauDRO, depFrac, arrFrac, thrustN, entry.ok);
    fprintf('   t_f = %.4f ND = %.3f d, m_f = %.2f kg, dV = %.3f km/s; ms |R| %.1e, flown %.3f km, accept |dz| %.1e\n', ...
            entry.tf_nd, entry.tf_d, entry.mf_kg, entry.dV_kms, gates.msNormR, gates.flownKm, gates.accDz);
end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
```

- [ ] **Step 4: Run the test (MCP); expected every line `ok` and `test_gto_dro_entry: PASS`.** Run `check_matlab_code` on the new file.

- [ ] **Step 5: Append `'test_gto_dro_entry'` to `names` in `tests/run_tests.m`; run `run_tests` (MCP): `GTO_DRO: 2 of 2 passed`.**

- [ ] **Step 6: Smoke the real solve once** (this is the first GTO -> DRO transfer ever solved here; it is also the number stage 2's pilot calibrates against). MCP `evaluate_matlab_code`:
```matlab
cd('/Users/msc/Desktop/optimal_control/orbit_transfer/GTO_DRO'); setup_gto_dro_paths;
t0 = tic; [e, g, P] = gto_dro_entry(0, 1, 0, 0, 15, struct('keepScratch', true)); toc(t0)
disp(e); disp(g)
```
Expected: `e.ok == true`, `g.flownKm < 100`, `g.accDz < 1e-6`, `g.msNormR < 1e-6`, t_f of order 0.2 to 0.8 ND. Record t_f, m_f, dV, wall time and t_f / 0.30 in the report. If `ok` is false, read the scratch log at `results/scratch/gto_dro_o000_tau1_15N.log`, try `tf0` in {0.20, 0.45, 0.60} (knob, not code), and report which converged; DONE_WITH_CONCERNS if none does (the controller rules on the recipe).

- [ ] **Step 7: Report**: test output, smoke numbers, any recipe finding.

---

### Task 3: The study script, sections 0 to 7 (knobs, orbits, endpoints, solve, flight, clearances)

**Files:**
- Create: `orbit_transfer/GTO_DRO/gto_dro_study.m`

**Interfaces:**
- Consumes: `gto_dro_entry` (Task 2), `earth_clearance` (Task 1), `nd_propulsion`, `get_family_orbit`, `phase_state`, `fly_transfer(z8, rv0, rvf, phys)` -> `.t .Y .flyKm .flyVms .tfDays .dvKms .finalMassKg .propellantKg`, `pumpkyn.cr3bp.{jacobi, toPCI, minDeltaV}`.
- Produces: workspace variables Task 4 extends: `tol`, `knobs`, `entry`, `gates`, `P`, `rv0`, `rvf`, `flight`, `Tnd`, `cnd`, `muStar`, `lStar`, `tStar`, `verdicts` (cell array `{id, status, value, threshold}`), `hEarth`, `hMoon`, local functions `fieldd`, `pass3`, `addv`.
- The `study_opts` hook: any field of `study_opts` named like a knob (`thrustN`, `m0kg`, `ispS`, `gto`, `dro`, `mesh`, `recipe`, `tol`, `out`, `selfCheck`) overrides the default; struct knobs are merged field by field.

- [ ] **Step 1: Write the script, sections 0 to 7, with a provisional save so it runs end to end**

```matlab
%% GTO_DRO_STUDY  One minimum-time low-thrust transfer: GTO -> DRO
%   Edit the knobs, press Run. Every step is done here rather than behind a
%   front door so the machinery is visible; the SOLVE goes through the one
%   shared engine (lib/gto_dro_entry -> costate_common/thrust_ladder_library:
%   direct solve -> fly -> ms_tfmin -> pumpkyn tfMin acceptance) and is never
%   reimplemented. Five instruments stay in the library because a second copy
%   here would be a second, unverified implementation: the pointwise PMP
%   checks (N2, N4-N6), the hypothesis gates (S1-S3, X2), lift_margin (S3),
%   conj_spectrum (S4) and fly_transfer (the ONE flight every number below
%   reads).
%
%   The GTO "orbit" is get_family_orbit's algebraic Kepler locus frozen at
%   epoch zero in the rotating frame: a SET of departure states keyed by
%   orientation and time-since-perigee fraction, not a trajectory. The DRO
%   is pumpkynPie's getDRO seed refined by cont_np and propagated one period.
%
%     0  knobs                             (one source; study_opts overrides)
%     1  paths and constants
%     2  the DEPARTURE locus (GTO)         V1 geometry
%     3  the ARRIVAL orbit (DRO)           V2 closure
%     4  endpoints and the energy bound
%     5  solve down the ladder             N0 ladder, V3 ratio screen
%     6  the flight                        N3 arrival
%     7  clearances                        N7a Earth, N7b Moon
%     8  necessary conditions              N1, N2, N4, N5, N6
%     9  sufficiency hypotheses            S1-S4, X2
%    10  acceptance by pumpkyn             X1
%    11  figures
%    12  save and the verdict table
%   IDs: N necessary (fail one and this is not an extremal), S sufficiency
%   hypotheses (Bonnard-Caillau-Trelat; results, never thrown), V validity
%   of the inputs and instruments, X cross-checks against pumpkyn.
%
%   study_opts: a struct in the workspace before running overrides any knob
%   (e.g. study_opts.thrustN = 7; study_opts.out.visible = 'off'). Used by
%   tests/test_gto_dro_study and by anyone scripting several runs.
%
%  M. Casey                                                   (c) 11/09/2026
%  Copyright Coorbital Inc.

clearvars -except study_opts; clc
if ~exist('study_opts', 'var'), study_opts = struct(); end
here = fileparts(mfilename('fullpath'));

%% ========================================================================
%  0. KNOBS -- every adjustable parameter, with units; one source.
%% ========================================================================
thrustN = 15;                        % returned rung [N]; the ladder walks 15 N down to it
m0kg    = 150;                       % initial mass [kg]
ispS    = 1710;                      % specific impulse [s]
gto     = struct( ...
    'orientDeg', 0, ...              % Earth->Moon line to perigee [deg], rotating sense
    'sma_km',    24446, ...          % semi-major axis [km] (350 x 35786 km GTO)
    'ecc',       0.724781150290436, ...
    'depFrac',   0);                 % departure phase: TIME fraction since perigee (0.5 = apogee)
dro     = struct( ...
    'tau',       1, ...              % DRO period [ND]; 0.5 1 2 3 are the catalog periods (1 = 4.43 d)
    'arrFrac',   0);                 % arrival phase: time fraction along the DRO from its +x crossing
mesh    = struct('N', 400);          % Hermite-Simpson intervals per rung
recipe  = struct( ...
    'thrLock',   true, ...           % pin the throttle to 1 (classical min-time)
    'tf0',       0.30, ...           % top-rung t_f guess [ND] (tuned for tulips; printed against the result)
    'maxCpuSec', 600, ...            % IPOPT wall cap per rung [s]
    'K',         24);                % multiple-shooting segments for the conjugate re-solve
tol     = struct( ...
    'gtoKm',     1, ...              % V1 perigee/apogee radius vs a(1 -+ e) [km]
    'closure',   1e-7, ...           % V2 DRO periodicity [ND]
    'ratioLo',   1/3, 'ratioHi', 3, ...   % V3 t_f ratio between consecutive rungs (steps only)
    'msNormR',   1e-8, ...           % N1 multiple-shooting residual, inf-norm
    'H',         1e-6, ...           % N2 |H| along the arc
    'km',        100, 'ms', 10, ...  % N3 flown arrival position [km] and velocity [m/s]
    'lamm',      1e-6, ...           % N4 |lambda_m(t_f)|
    'adj',       1e-7, ...           % N5 adjoint equations, relative
    'gap',       1e-12, 'throttle', 1e-10, ...   % N6 minimum-principle gap, applied throttle
    'earthKm',   200, ...            % N7a Earth altitude floor [km] (below the 350 km perigee by design)
    'moonKm',    500, ...            % N7b lunar altitude floor [km] (the engine's own floorKm)
    'lift',      1e-6, 'liftMargin', 10, ...     % S3 lift backward error, Eckart-Young margin
    'accDz',     1e-6, ...           % X1 pumpkyn tfMin must reproduce z8 to this
    'field',     1e-10);             % X2 pumpkyn field vs independent field, relative
out     = struct( ...
    'dir',       '', ...             % results folder ('' = results/study_o<orient>_tau<tau>_<thrust>N)
    'visible',   'on', ...           % figure visibility
    'conjTest',  true, ...           % run the conjugate re-solve (S4); false = NOT APPLICABLE
    'reuse',     '', ...             % path of a saved gto_dro_study.mat to gate instead of solving
    'solverType', 'fsolve');         % X1: pumpkyn tfMin solver ('fsolve' | 'lsqnonlin')
selfCheck = struct('strict', false); % true: any FAIL or UNRESOLVED gate throws (tests); false: code gates only

% study_opts overrides (struct knobs merge field by field)
for f = fieldnames(study_opts).'
    name = f{1};
    if exist(name, 'var') && isstruct(study_opts.(name)) && isstruct(eval(name))
        cur = eval(name);
        for g = fieldnames(study_opts.(name)).'
            cur.(g{1}) = study_opts.(name).(g{1});
        end
        eval([name ' = cur;']);
    else
        eval([name ' = study_opts.(name);']);
    end
end
if isempty(out.dir)
    out.dir = fullfile(here, 'results', sprintf('study_o%03d_tau%g_%gN', gto.orientDeg, dro.tau, thrustN));
end
if ~isfolder(out.dir), mkdir(out.dir); end
knobs = struct('thrustN', thrustN, 'm0kg', m0kg, 'ispS', ispS, 'gto', gto, 'dro', dro, ...
               'mesh', mesh, 'recipe', recipe, 'tol', tol, 'out', out);
verdicts = cell(0, 4);                                       % {id, status, value, threshold}
tStudy = tic;

%% ========================================================================
%  1. PATHS AND CONSTANTS
%% ========================================================================
fprintf('1. PATHS AND CONSTANTS\n');
run(fullfile(here, 'setup_gto_dro_paths.m'));
muStar = 0.012150585609624;                 % Earth-Moon mass ratio
lStar  = 389703.264829278;                  % km
tStar  = 382981.289129055;                  % s
rEarthKm = 6378.137;  rMoonKm = 1737.4;
day    = @(tnd) tnd*tStar/86400;
ndp    = nd_propulsion(thrustN, ispS, m0kg, lStar, tStar);
Tnd    = ndp.Tnd;   cnd = ndp.cnd;
fprintf('   mu* = %.12f, L* = %.3f km, T* = %.3f s (1 TU = %.4f d)\n', muStar, lStar, tStar, day(1));
fprintf('   thrust %.3g N on %g kg at Isp %g s -> T/m0 = %.3e ND, c = %.4f ND (%.3f km/s)\n', ...
        thrustN, m0kg, ispS, Tnd, cnd, cnd*lStar/tStar);

%% ========================================================================
%  2. THE DEPARTURE LOCUS (GTO)
%% ========================================================================
fprintf('\n2. DEPARTURE: GTO locus, orientation %g deg\n', gto.orientDeg);
[tD, rvD, infoD] = get_family_orbit('gto', struct('orientDeg', gto.orientDeg, 'sma_km', gto.sma_km, 'ecc', gto.ecc));
rDkm   = lStar*sqrt(sum((rvD(:, 1:3) - [-muStar, 0, 0]).^2, 2));
rPeri  = min(rDkm);  rApo = max(rDkm);
TkepD  = day(tD(end));
fprintf('   %d locus points; perigee altitude %.1f km, apogee altitude %.1f km, Kepler period %.4f d\n', ...
        size(rvD, 1), rPeri - rEarthKm, rApo - rEarthKm, TkepD);
fprintf('   (a SET of departure states frozen at epoch zero in the rotating frame, not a trajectory)\n');
v1val = max(abs(rPeri - gto.sma_km*(1 - gto.ecc)), abs(rApo - gto.sma_km*(1 + gto.ecc)));
verdicts = addv(verdicts, 'V1 GTO geometry', pass3(v1val < tol.gtoKm), v1val, tol.gtoKm);

%% ========================================================================
%  3. THE ARRIVAL ORBIT (DRO)
%% ========================================================================
fprintf('\n3. ARRIVAL: DRO, period tau = %g ND\n', dro.tau);
[tA, rvA, infoA] = get_family_orbit('dro', struct('tau', dro.tau));
closure = max(abs(rvA(end, :) - rvA(1, :)));
CJ      = pumpkyn.cr3bp.jacobi(rvA(1, :), muStar);
dAkm    = lStar*sqrt(sum((rvA(:, 1:3) - [1 - muStar, 0, 0]).^2, 2));
fprintf('   period %.4f d, Jacobi C = %.5f, perilune altitude %.0f km, apolune altitude %.0f km, closure %.1e ND\n', ...
        day(infoA.periodND), CJ, min(dAkm) - rMoonKm, max(dAkm) - rMoonKm, closure);
verdicts = addv(verdicts, 'V2 DRO closure', pass3(closure < tol.closure), closure, tol.closure);

%% ========================================================================
%  4. ENDPOINTS AND THE ENERGY BOUND
%% ========================================================================
fprintf('\n4. ENDPOINTS (costate_common/phase_state, spline rule)\n');
depAt = phase_state(tD, rvD);   arrAt = phase_state(tA, rvA);
rv0 = depAt(gto.depFrac).';     rvf = arrAt(dro.arrFrac).';
rv0km = pumpkyn.cr3bp.toPCI(0, rv0, muStar, tStar, lStar, 1);     % Earth-centred, km and km/s
rvfkm = pumpkyn.cr3bp.toPCI(0, rvf, muStar, tStar, lStar, 2);     % Moon-centred
sepKm = lStar*sqrt(sum((rvf(1:3) - rv0(1:3)).^2));
dVlb  = pumpkyn.cr3bp.minDeltaV(rv0, rvf, muStar)*lStar/tStar;
fprintf('   depart  (phase %.3f) ND [%s]\n', gto.depFrac, sprintf('%+.6f ', rv0));
fprintf('           Earth-centred: r = %.1f km, v = %.4f km/s\n', sqrt(sum(rv0km(1:3).^2)), sqrt(sum(rv0km(4:6).^2)));
fprintf('   arrive  (phase %.3f) ND [%s]\n', dro.arrFrac, sprintf('%+.6f ', rvf));
fprintf('           Moon-centred:  r = %.1f km, v = %.4f km/s\n', sqrt(sum(rvfkm(1:3).^2)), sqrt(sum(rvfkm(4:6).^2)));
fprintf('   straight-line separation %.0f km; pumpkyn minDeltaV energy bound %.3f km/s\n', sepKm, dVlb);

%% ========================================================================
%  5. SOLVE DOWN THE LADDER (lib/gto_dro_entry -> the shared engine)
%% ========================================================================
fprintf('\n5. SOLVE: 15 N down to %g N\n', thrustN);
entryOpts = struct('m0kg', m0kg, 'ispS', ispS, 'thrLock', recipe.thrLock, 'tf0', recipe.tf0, ...
    'maxCpuSec', recipe.maxCpuSec, 'N', mesh.N, 'K', recipe.K, 'conjTest', out.conjTest, ...
    'gtoParams', struct('sma_km', gto.sma_km, 'ecc', gto.ecc), ...
    'outMat', fullfile(out.dir, 'ladder_scratch.mat'), 'keepScratch', true);
if isempty(out.reuse)
    tSolve = tic;
    [entry, gates, P] = gto_dro_entry(gto.orientDeg, dro.tau, gto.depFrac, dro.arrFrac, thrustN, entryOpts);
    solveWall = toc(tSolve);
    fprintf('   solved in %.0f s wall\n', solveWall);
else
    S = load(out.reuse, 'entry', 'gates', 'P', 'knobs');
    same = S.knobs.thrustN == thrustN && S.knobs.m0kg == m0kg && S.knobs.ispS == ispS && ...
           S.knobs.gto.orientDeg == gto.orientDeg && S.knobs.gto.depFrac == gto.depFrac && ...
           S.knobs.gto.sma_km == gto.sma_km && S.knobs.gto.ecc == gto.ecc && ...
           S.knobs.dro.tau == dro.tau && S.knobs.dro.arrFrac == dro.arrFrac;
    assert(same, 'gto_dro_study:reuse', ...
           '%s was solved with different knobs (thrust, mass, Isp, GTO or DRO parameters); refusing to gate it as this case', out.reuse);
    entry = S.entry;  gates = S.gates;  P = S.P;  solveWall = NaN;
    fprintf('   REUSED %s (no solve)\n', out.reuse);
end
fprintf('   %6s %9s %8s %8s %8s %9s %9s %9s %6s\n', 'N', 't_f ND', 't_f d', 'm_f kg', 'dV km/s', 'flown km', 'ms |R|', 'acc |dz|', 'ok');
for kr = 1:numel(P.rungs)
    tfk = P.TF(1, 1, kr);  mfk = m0kg*(1 - nd_propulsion(P.rungs(kr), ispS, m0kg, lStar, tStar).Tnd*tfk/cnd);
    dVk = cnd*log(m0kg/mfk)*lStar/tStar;
    fprintf('   %6.3g %9.4f %8.3f %8.2f %8.3f %9.3f %9.1e %9.1e %6d\n', P.rungs(kr), tfk, day(tfk), mfk, dVk, ...
            P.FLYKM(1, 1, kr), P.RES(1, 1, kr), P.ACCDZ(1, 1, kr), P.OK(1, 1, kr));
end
fprintf('   top-rung t_f / tf0 = %.2f (tf0 = %.2f ND was tuned for tulips)\n', entry.tfPerRung(1)/recipe.tf0, recipe.tf0);
verdicts = addv(verdicts, 'N0 ladder', pass3(all(entry.okPerRung)), nnz(entry.okPerRung), numel(entry.okPerRung));
if numel(P.rungs) > 1
    rs = gates.ratioScreen;
    fprintf('   step ratios t_f(k+1)/t_f(k): %s\n', sprintf('%.3f ', rs));
    verdicts = addv(verdicts, 'V3 ratio screen', pass3(all(rs > tol.ratioLo & rs < tol.ratioHi)), max(abs(log(rs))), log(tol.ratioHi));
else
    verdicts = addv(verdicts, 'V3 ratio screen', 'NOT APPLICABLE', NaN, NaN);
end
assert(entry.ok, 'gto_dro_study:noRoot', 'the bottom rung did not produce an accepted entry; nothing below can be gated');

%% ========================================================================
%  6. THE FLIGHT -- pumpkyn tfMinProp from the accepted z8 (ONE flight; every
%     number in sections 7-11 reads it). The engine's own G1b flown miss of
%     the DIRECT solution at this rung is printed beside it (P.FLYKM): the
%     direct control is not saved by the engine, so it is read, not re-flown.
%% ========================================================================
fprintf('\n6. FLIGHT of the accepted costates\n');
phys   = struct('Tnd', Tnd, 'cnd', cnd, 'muStar', muStar, 'lStar', lStar, 'tStar', tStar, 'm0kg', m0kg);
flight = fly_transfer(entry.z8, rv0(1:6).', rvf(1:6).', phys);
fprintf('   t_f = %.4f ND = %.3f d; m_f = %.2f kg (%.2f kg propellant); dV = %.3f km/s (bound %.3f)\n', ...
        entry.z8(8), flight.tfDays, flight.finalMassKg, flight.propellantKg, flight.dvKms, dVlb);
fprintf('   arrival miss %.4f km, %.4f m/s (z8 flight); engine G1b flown miss of the direct solution %.4f km\n', ...
        flight.flyKm, flight.flyVms, gates.flownKm);
verdicts = addv(verdicts, 'N3 arrival', pass3(flight.flyKm < tol.km && flight.flyVms < tol.ms), flight.flyKm, tol.km);

%% ========================================================================
%  7. CLEARANCES on the flown arc (lib/earth_clearance, one instrument twice)
%% ========================================================================
fprintf('\n7. CLEARANCES\n');
[hEarth, tE, pE, oE] = earth_clearance(flight.t, flight.Y, tol.earthKm, struct('muStar', muStar, 'lStar', lStar));
[hMoon,  tM, pM, oM] = earth_clearance(flight.t, flight.Y, tol.moonKm, ...
                           struct('muStar', muStar, 'lStar', lStar, 'center', [1 - muStar, 0, 0], 'radiusKm', rMoonKm));
fprintf('   N7a Earth altitude min %9.1f km at t/t_f = %.3f  (floor %4.0f)  %s\n', hEarth, tE/flight.t(end), tol.earthKm, pass3(pE));
fprintf('   N7b Moon  altitude min %9.1f km at t/t_f = %.3f  (floor %4.0f)  %s\n', hMoon,  tM/flight.t(end), tol.moonKm,  pass3(pM));
verdicts = addv(verdicts, 'N7a Earth floor', pass3(pE), hEarth, tol.earthKm);
verdicts = addv(verdicts, 'N7b Moon floor',  pass3(pM), hMoon,  tol.moonKm);

%% ======================= sections 8-12 follow (Task 4) ===================
save(fullfile(out.dir, 'gto_dro_study.mat'), 'entry', 'gates', 'P', 'knobs', 'verdicts', 'flight');
fprintf('\n   (provisional) saved %s after %.0f s\n', fullfile(out.dir, 'gto_dro_study.mat'), toc(tStudy));

%% ------------------------------------------------------------------------
function v = fieldd(s, f, d_)
% FIELDD  s.(f) if present else d_.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = d_; end
end

function s = pass3(c)
% PASS3  PASS / FAIL from a logical.  INPUTS: c.  OUTPUTS: s.
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function V = addv(V, id, status, value, threshold)
% ADDV  Append one verdict row {id, status, value, threshold}.  INPUTS: V; id; status; value; threshold.  OUTPUTS: V.
V(end+1, :) = {id, status, value, threshold};
end
```

- [ ] **Step 2: Run the script (MCP `run_matlab_file`) at the demo point.** Expected: sections 1 to 7 print, `V1`, `V2`, `N0`, `N3`, `N7a`, `N7b` show PASS (V3 NOT APPLICABLE at 15 N), the `.mat` is saved. Takes a few minutes (one 15 N solve plus the conjugate re-solve). Record every printed number in the report, especially `top-rung t_f / tf0`, the Earth minimum altitude and where it occurs.

- [ ] **Step 3: Run `check_matlab_code` on the script; fix real messages. Report.**

---

### Task 4: The study script, sections 8 to 12 (conditions, hypotheses, acceptance, figures, verdict table)

**Files:**
- Modify: `orbit_transfer/GTO_DRO/gto_dro_study.m` (replace the provisional block from `%% ======================= sections 8-12 follow` through the provisional `save`/`fprintf`, keeping the local functions)

**Interfaces:**
- Consumes: `pmp_pointwise_checks(t, Y, Tnd, cnd, mu)` -> `.Hmax .lamMf .adjErr .nSample .fdAgree .fullGap .fieldGap .gapMin .gapMax .throttleErr .throttleAccErr .throttleMassErr .minQmt`; `mintime_hypothesis_gates(z8, rv0, Tmax, c, muStar, opts)` -> `.minLamV .tMinLamV .minQmt .tMinQmt .dimS .nullResid .nullResidRel .Hresid .sv .svRatio .C .fieldErr .adjErrRef .h6Ok .h6Margin .h6LamM0 .h6Threshold .h6Clearance`; `lift_margin(C_tight, C_loose, lam0, struct('marginMin', ., 'liftTol', .))` -> `.certified .sigma6 .errEst .margin`; `conj_spectrum(z8, rv0, Tmax, c, muStar, struct('K', ., 'nSub', .))` -> `.clear .testable .nZero .nUnresolved .nInterior .multiplicity .t`; `pumpkyn.cr3bp.tfMin(rv0, rvf, z8, Tmax, c, muStar, solverType)`; `pumpkyn.cr3bp.lagrangePts(mu)` [5 x 3].

- [ ] **Step 1: Replace the provisional block with sections 8 to 12**

```matlab
%% ========================================================================
%  8. NECESSARY CONDITIONS (the production instrument: pmp_pointwise_checks,
%     the same call certify_root makes on every certified catalog entry)
%% ========================================================================
fprintf('\n8. NECESSARY CONDITIONS      (value / threshold)\n');
PW = pmp_pointwise_checks(flight.t, flight.Y, Tnd, cnd, muStar);
n1 = gates.msNormR < tol.msNormR;
fprintf('   N1 BVP residual   |R|_inf   %9.2e / %-9.0e  %s   (ms_tfmin, the engine''s own)\n', gates.msNormR, tol.msNormR, pass3(n1));
fprintf('   N2 Hamiltonian    max|H|    %9.2e / %-9.0e  %s\n', PW.Hmax, tol.H, pass3(PW.Hmax < tol.H));
fprintf('   N4 transversality |lam_m(t_f)| %6.2e / %-9.0e  %s\n', PW.lamMf, tol.lamm, pass3(PW.lamMf < tol.lamm));
fprintf('   N5 adjoint eqns   rel err   %9.2e / %-9.0e  %s   (%d samples; FD steps agree to %.1e)\n', ...
        PW.adjErr, tol.adj, pass3(PW.adjErr < tol.adj), PW.nSample, PW.fdAgree);
n6gap = max(PW.fullGap, PW.fieldGap);
n6 = n6gap <= tol.gap && PW.throttleErr < tol.throttle && PW.minQmt >= 0;
fprintf('   N6 min principle  |gap|     %9.2e / %-9.0e  %s   (signed %+.1e..%+.1e; throttle 1 to %.1e acc / %.1e mass; min Q = %.3f)\n', ...
        n6gap, tol.gap, pass3(n6), PW.gapMin, PW.gapMax, PW.throttleAccErr, PW.throttleMassErr, PW.minQmt);
for fv = {{'N1', gates.msNormR}, {'N2', PW.Hmax}, {'N4', PW.lamMf}, {'N5', PW.adjErr}, {'N6 gap', n6gap}, {'N6 throttle', PW.throttleErr}}
    assert(isscalar(fv{1}{2}) && isreal(fv{1}{2}) && isfinite(fv{1}{2}) && fv{1}{2} >= 0, 'gto_dro_study:malformed', ...
           '%s is not a real finite non-negative scalar (%g): the gate cannot be judged', fv{1}{1}, fv{1}{2});
end
verdicts = addv(verdicts, 'N1 BVP residual',   pass3(n1),                 gates.msNormR, tol.msNormR);
verdicts = addv(verdicts, 'N2 Hamiltonian',    pass3(PW.Hmax < tol.H),    PW.Hmax,       tol.H);
verdicts = addv(verdicts, 'N4 transversality', pass3(PW.lamMf < tol.lamm), PW.lamMf,     tol.lamm);
verdicts = addv(verdicts, 'N5 adjoint eqs',    pass3(PW.adjErr < tol.adj), PW.adjErr,    tol.adj);
verdicts = addv(verdicts, 'N6 min principle',  pass3(n6),                 n6gap,         tol.gap);

%% ========================================================================
%  9. SUFFICIENCY HYPOTHESES (Bonnard-Caillau-Trelat) and the physics check
%% ========================================================================
fprintf('\n9. SUFFICIENCY HYPOTHESES\n');
G  = mintime_hypothesis_gates(entry.z8, rv0(1:6).', Tnd, cnd, muStar, struct('keepC', true, 'relTol', 1e-12));
GL = mintime_hypothesis_gates(entry.z8, rv0(1:6).', Tnd, cnd, muStar, struct('keepC', true, 'relTol', 1e-9));
LM = lift_margin(G.C, GL.C, entry.z8(1:7), struct('marginMin', tol.liftMargin, 'liftTol', tol.lift));
s1 = G.minLamV > 0;   s2 = G.minQmt > 0;
s3 = G.dimS == 1 && G.nullResidRel < tol.lift && G.Hresid < tol.H && LM.certified;
fprintf('   S1 Legendre      min|lam_v| = %.4e at t/t_f = %.3f   %s\n', G.minLamV, G.tMinLamV/entry.z8(8), pass3(s1));
fprintf('   S2 strict bang   min Q = %.4e at t/t_f = %.3f        %s\n', G.minQmt, G.tMinQmt/entry.z8(8), pass3(s2));
fprintf('   S3 normality     dim S = %d (1 = no abnormal lift); lift backward error %.1e / %.0e; |lam.f + 1| %.1e; Eckart-Young margin %.0fx / %gx   %s\n', ...
        G.dimS, G.nullResidRel, tol.lift, G.Hresid, LM.margin, tol.liftMargin, pass3(s3));
% S4: the conjugate re-solve's verdict (gto_dro_entry, conjTest) AND the dense spectrum scan
if ~out.conjTest
    s4 = 'NOT APPLICABLE';  nScan = 0;
    fprintf('   S4 conjugate     NOT APPLICABLE (out.conjTest = false)\n');
else
    CS = conj_spectrum(entry.z8, rv0(1:6).', Tnd, cnd, muStar, struct('K', recipe.K, 'nSub', 8));
    nScan = numel(CS.t);
    if isnan(gates.conjPass) || ~CS.testable || CS.nUnresolved > 0
        s4 = 'UNRESOLVED';
    elseif gates.conjPass == 1 && CS.clear && CS.nZero == 0 && CS.multiplicity == 0 && CS.nInterior == 0
        s4 = 'PASS';
    else
        s4 = 'FAIL';
    end
    fprintf('   S4 conjugate     junction test pass = %g; dense scan (%d samples): %d zero, %d unresolved, %d interior sign change(s), clear = %d   %s\n', ...
            gates.conjPass, nScan, CS.nZero, CS.nUnresolved, CS.nInterior, CS.clear, s4);
end
x2 = G.fieldErr < tol.field && G.adjErrRef < tol.field;
fprintf('   X2 physics       pumpkyn vs independent field: state rows %.1e, adjoint rows %.1e / %-6.0e %s\n', G.fieldErr, G.adjErrRef, tol.field, pass3(x2));
verdicts = addv(verdicts, 'S1 strong Legendre', pass3(s1), G.minLamV, 0);
verdicts = addv(verdicts, 'S2 Q_mt > 0',        pass3(s2), G.minQmt,  0);
verdicts = addv(verdicts, 'S3 dim S = 1',       pass3(s3), G.dimS,    1);
verdicts = addv(verdicts, 'S4 conjugate',       s4,        gates.conjPass, 1);
verdicts = addv(verdicts, 'X2 field agreement', pass3(x2), max(G.fieldErr, G.adjErrRef), tol.field);

%% ========================================================================
% 10. ACCEPTANCE BY PUMPKYN -- tfMin called at the accepted z8, in the open
%% ========================================================================
fprintf('\n10. ACCEPTANCE: pumpkyn.cr3bp.tfMin at the accepted z8 (%s)\n', out.solverType);
x1 = 'NOT CHECKED';  dzAcc = NaN;
try
    sol   = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), entry.z8, Tnd, cnd, muStar, out.solverType);
    dzAcc = max(abs(sol(:) - entry.z8));
    x1    = pass3(dzAcc < tol.accDz);
    fprintf('   X1 tfMin accept  |dz| %9.2e / %-9.0e  %s   (engine recorded %.1e)\n', dzAcc, tol.accDz, x1, gates.accDz);
catch err
    fprintf('   X1 tfMin accept  NOT CHECKED: tfMin threw "%s" (set out.solverType = ''lsqnonlin'' if fsolve is unlicensed)\n', err.message);
end
verdicts = addv(verdicts, 'X1 tfMin acceptance', x1, dzAcc, tol.accDz);

%% ========================================================================
% 11. FIGURES (light theme, explicit axes; saved at 150 dpi)
%% ========================================================================
fprintf('\n11. FIGURES -> %s\n', out.dir);
FS = 13;
mFrac = flight.Y(:, 7);
Lp = pumpkyn.cr3bp.lagrangePts(muStar);
[sx, sy, sz] = sphere(24);
hf1 = figure('Theme', 'light', 'Color', 'w', 'Position', [60 60 1400 650], 'Visible', out.visible);
tl1 = tiledlayout(hf1, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for view3 = [true false]
    ax = nexttile(tl1);  hold(ax, 'on');
    plot3(ax, rvD(:, 1), rvD(:, 2), rvD(:, 3), '-', 'Color', [0.1 0.3 0.8], 'LineWidth', 1.6, 'DisplayName', sprintf('GTO locus, %g deg', gto.orientDeg));
    plot3(ax, rvA(:, 1), rvA(:, 2), rvA(:, 3), '-', 'Color', [0.1 0.6 0.2], 'LineWidth', 1.6, 'DisplayName', sprintf('DRO tau = %g', dro.tau));
    scatter3(ax, flight.Y(:, 1), flight.Y(:, 2), flight.Y(:, 3), 9, mFrac, 'filled', 'DisplayName', 'transfer (colour = m/m_0)');
    surf(ax, -muStar + rEarthKm/lStar*sx, rEarthKm/lStar*sy, rEarthKm/lStar*sz, 'FaceColor', [0.3 0.5 0.9], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    surf(ax, 1 - muStar + rMoonKm/lStar*sx, rMoonKm/lStar*sy, rMoonKm/lStar*sz, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot3(ax, Lp(1:2, 1), Lp(1:2, 2), Lp(1:2, 3), '^k', 'MarkerFaceColor', 'k', 'MarkerSize', 7, 'DisplayName', 'L1, L2');
    plot3(ax, rv0(1), rv0(2), rv0(3), 'o', 'Color', [0.1 0.3 0.8], 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'departure');
    plot3(ax, rvf(1), rvf(2), rvf(3), 's', 'Color', [0.1 0.6 0.2], 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'arrival');
    axis(ax, 'equal'); grid(ax, 'on');
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', 'FontSize', FS);
    xlabel(ax, 'x [ND]'); ylabel(ax, 'y [ND]');
    if view3, zlabel(ax, 'z [ND]'); view(ax, [-35 22]); legend(ax, 'Location', 'northwest', 'FontSize', FS - 2);
    else,     view(ax, 2); end
    cb = colorbar(ax); cb.Label.String = 'm / m_0'; cb.Color = 'k';
end
title(tl1, sprintf('GTO (%g deg, phase %.2f) -> DRO (tau %g, phase %.2f) at %g N: t_f = %.3f d, dV = %.3f km/s', ...
      gto.orientDeg, gto.depFrac, dro.tau, dro.arrFrac, thrustN, flight.tfDays, flight.dvKms), 'FontSize', FS + 1);
print(hf1, fullfile(out.dir, 'transfer_3d.png'), '-dpng', '-r150');

td = day(flight.t);
lamV = sqrt(sum(flight.Y(:, 11:13).^2, 2));
alphaDir = -flight.Y(:, 11:13)./lamV;
Qmt = lamV./mFrac + flight.Y(:, 14)/cnd;
Hres = nan(size(td));
for kk = 1:numel(td)
    % tfMinEoM's second output is the Hamiltonian; it takes one augmented state.
    % If it wants a column (check its header), pass flight.Y(kk, :).' instead.
    [~, Hres(kk)] = pumpkyn.cr3bp.tfMinEoM(flight.t(kk), flight.Y(kk, :), Tnd, cnd, muStar);
end
hf2 = figure('Theme', 'light', 'Color', 'w', 'Position', [80 40 1000 1100], 'Visible', out.visible);
tl2 = tiledlayout(hf2, 5, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl2); hold(ax, 'on');
plot(ax, td, alphaDir(:, 1), 'LineWidth', 1.4, 'DisplayName', '\alpha_x');
plot(ax, td, alphaDir(:, 2), 'LineWidth', 1.4, 'DisplayName', '\alpha_y');
plot(ax, td, alphaDir(:, 3), 'LineWidth', 1.4, 'DisplayName', '\alpha_z');
ylabel(ax, 'thrust direction'); legend(ax, 'Location', 'eastoutside'); grid(ax, 'on');
ax = nexttile(tl2); hold(ax, 'on');
plot(ax, td, lamV, 'LineWidth', 1.4, 'DisplayName', '|\lambda_v|');
plot(ax, td, Qmt, 'LineWidth', 1.4, 'DisplayName', 'Q_{mt}');
ylabel(ax, 'switching quantities'); legend(ax, 'Location', 'eastoutside'); grid(ax, 'on');
ax = nexttile(tl2);
plot(ax, td, m0kg*mFrac, 'LineWidth', 1.4, 'DisplayName', 'mass'); ylabel(ax, 'mass [kg]'); grid(ax, 'on'); legend(ax, 'Location', 'eastoutside');
ax = nexttile(tl2);
semilogy(ax, td, max(abs(Hres), 1e-18), 'LineWidth', 1.4, 'DisplayName', '|H|'); ylabel(ax, 'Hamiltonian'); grid(ax, 'on'); legend(ax, 'Location', 'eastoutside');
ax = nexttile(tl2); hold(ax, 'on');
semilogy(ax, td, oE.hKm, 'LineWidth', 1.4, 'DisplayName', 'Earth altitude');
semilogy(ax, td, oM.hKm, 'LineWidth', 1.4, 'DisplayName', 'Moon altitude');
yline(ax, tol.earthKm, '--', 'HandleVisibility', 'off'); yline(ax, tol.moonKm, ':', 'HandleVisibility', 'off');
ylabel(ax, 'altitude [km]'); xlabel(ax, 'time [d]'); legend(ax, 'Location', 'eastoutside'); grid(ax, 'on');
for ax = findall(hf2, 'Type', 'axes').'
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', FS);
end
title(tl2, sprintf('GTO -> DRO at %g N: the flight of the accepted costates', thrustN), 'FontSize', FS + 1);
print(hf2, fullfile(out.dir, 'timeseries.png'), '-dpng', '-r150');

hf3 = figure('Theme', 'light', 'Color', 'w', 'Position', [100 100 900 420], 'Visible', out.visible);
tl3 = tiledlayout(hf3, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
okR = entry.okPerRung;
dVr = nan(size(P.rungs));
for kr = 1:numel(P.rungs)
    mfk = m0kg*(1 - nd_propulsion(P.rungs(kr), ispS, m0kg, lStar, tStar).Tnd*P.TF(1, 1, kr)/cnd);
    dVr(kr) = cnd*log(m0kg/mfk)*lStar/tStar;
end
ax = nexttile(tl3); semilogx(ax, P.rungs(okR), day(entry.tfPerRung(okR)), 'o-', 'LineWidth', 1.6, 'MarkerFaceColor', 'w', 'DisplayName', 't_f');
xlabel(ax, 'thrust [N]'); ylabel(ax, 't_f [d]'); grid(ax, 'on'); set(ax, 'XDir', 'reverse', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', FS);
ax = nexttile(tl3); semilogx(ax, P.rungs(okR), dVr(okR), 's-', 'LineWidth', 1.6, 'MarkerFaceColor', 'w', 'DisplayName', 'dV');
xlabel(ax, 'thrust [N]'); ylabel(ax, '\DeltaV [km/s]'); grid(ax, 'on'); set(ax, 'XDir', 'reverse', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', FS);
title(tl3, 'the ladder walked (accepted rungs)', 'FontSize', FS + 1);
print(hf3, fullfile(out.dir, 'ladder.png'), '-dpng', '-r150');
fprintf('   transfer_3d.png, timeseries.png, ladder.png written\n');

%% ========================================================================
% 12. SAVE AND THE VERDICT TABLE
%% ========================================================================
[~, gitOC] = system('/Library/Developer/CommandLineTools/usr/bin/git -C /Users/msc/Desktop/optimal_control rev-parse --short HEAD');
[~, gitPk] = system('/Library/Developer/CommandLineTools/usr/bin/git -C /Users/msc/Desktop/proj7/external/pumpkyn rev-parse --short HEAD');
[~, gitPP] = system('/Library/Developer/CommandLineTools/usr/bin/git -C /Users/msc/Desktop/proj7/external/pumpkynPie rev-parse --short HEAD');
[~, hostName] = system('hostname');
run_meta = struct('date', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), 'matlab', version, 'host', strtrim(hostName), ...
    'optimal_control', strtrim(gitOC), 'pumpkyn', strtrim(gitPk), 'pumpkynPie', strtrim(gitPP), ...
    'solveWallS', solveWall, 'studyWallS', toc(tStudy), 'reused', out.reuse);
save(fullfile(out.dir, 'gto_dro_study.mat'), 'entry', 'gates', 'P', 'knobs', 'verdicts', 'flight', 'PW', 'G', 'run_meta');

fprintf('\n== SUMMARY  GTO (%g deg, phase %.3f) -> DRO (tau %g, phase %.3f), %g N, %g kg, Isp %g s\n', ...
        gto.orientDeg, gto.depFrac, dro.tau, dro.arrFrac, thrustN, m0kg, ispS);
fprintf('   t_f = %.4f ND = %.3f d; m_f = %.2f kg; dV = %.3f km/s (energy bound %.3f)\n', ...
        entry.z8(8), flight.tfDays, flight.finalMassKg, flight.dvKms, dVlb);
fprintf('   %-22s %-15s %12s %12s\n', 'gate', 'status', 'value', 'threshold');
for kk = 1:size(verdicts, 1)
    fprintf('   %-22s %-15s %12.4g %12.4g\n', verdicts{kk, :});
end
isCode  = ~cellfun(@isempty, regexp(verdicts(:, 1), '^(V|N|X)', 'once'));
codeOK  = all(strcmp(verdicts(isCode, 2), 'PASS') | strcmp(verdicts(isCode, 2), 'NOT APPLICABLE'));
allOK   = all(strcmp(verdicts(:, 2), 'PASS') | strcmp(verdicts(:, 2), 'NOT APPLICABLE'));
failedIds = verdicts(~(strcmp(verdicts(:, 2), 'PASS') | strcmp(verdicts(:, 2), 'NOT APPLICABLE')), 1);
fprintf('   code gates (V, N, X): %s; all gates: %s; %.0f s\n', pass3(codeOK), pass3(allOK), toc(tStudy));
if allOK
    fprintf(['   VERDICT: a certified local minimum-time extremal that pumpkyn''s tfMin accepts unchanged;\n' ...
             '            the sufficiency hypotheses hold on the sampled arc (numerical evidence, not a certificate).\n']);
elseif codeOK
    fprintf('   VERDICT: an extremal; sufficiency not claimed (%s).\n', strjoin(failedIds.', ', '));
else
    fprintf('   VERDICT: the root or an instrument is in doubt (%s); no claim.\n', strjoin(failedIds.', ', '));
end
assert(codeOK, 'gto_dro_study:brokenRoot', 'code gates did not pass: %s', strjoin(failedIds.', ', '));
if selfCheck.strict
    assert(allOK, 'gto_dro_study:notClaimed', 'strict self-check: %s', strjoin(failedIds.', ', '));
end
```

- [ ] **Step 2: Run the script (MCP `run_matlab_file`)** at the demo point. Expected: sections 8 to 12 print; N1, N2, N4, N5, N6, X1, X2 PASS; S1 to S3 PASS; S4 PASS or UNRESOLVED (record which); the summary table; three PNGs in `results/study_o000_tau1_15N/`. If X1 throws for a licence reason, rerun with `study_opts.out.solverType = 'lsqnonlin'` and record both outcomes.

- [ ] **Step 3: Read the three PNGs (Read tool) and fix layout problems** (legend over data, unreadable fonts, dark axes). Confirm the 3-D panel shows the Earth-side GTO and the Moon-side DRO at the right scale.

- [ ] **Step 4: `check_matlab_code` on the script; fix real messages. Report** every printed number (the report is the first record of a GTO -> DRO transfer).

---

### Task 5: Headless test of the study, README, repo README row, memory

**Files:**
- Create: `orbit_transfer/GTO_DRO/tests/test_gto_dro_study.m`
- Modify: `orbit_transfer/GTO_DRO/tests/run_tests.m` (append the test name)
- Create: `orbit_transfer/GTO_DRO/README.md`
- Modify: `orbit_transfer/README.md` (one row/line pointing at GTO_DRO, in the campaign list where GTO_tulip and orbit_network appear)

**Interfaces:**
- Consumes: the `study_opts` hook (`thrustN`, `out.dir`, `out.visible`, `out.reuse`, `tol.earthKm`, `selfCheck.strict`), the printed table format `   <id padded 22> <status padded 15>`.

- [ ] **Step 1: Write `tests/test_gto_dro_study.m`**

```matlab
function ok = test_gto_dro_study()
%% Purpose:
%
%   The study script headless, three ways: (1) a full run at thrustN = 12
%   through study_opts (so the override is proven to REACH the solve: the
%   printed ladder must end at 12 N and the summary line must say 12 N)
%   with strict self-check, every code gate PASS, three PNGs and the .mat
%   written; (2) a REUSE run on that .mat with tol.earthKm = 1e6 and
%   strict off must print N7a FAIL and throw gto_dro_study:brokenRoot (the
%   gate can fail, and the script fails when it should); (3) a REUSE run
%   with a different DRO period must be refused (gto_dro_study:reuse). The
%   script parses clean (checkcode, parse errors only). Scripts run in the
%   base workspace, so they are driven with evalin.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every check passed
%
%% Revision History:
%  M. Casey                                                   (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok   = true;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
run(fullfile(root, 'setup_gto_dro_paths.m'));
file   = fullfile(root, 'gto_dro_study.m');
outDir = fullfile(tempdir, 'test_gto_dro_study');
if ~isfolder(outDir), mkdir(outDir); end
for f = {'transfer_3d.png', 'timeseries.png', 'ladder.png', 'gto_dro_study.mat'}
    if isfile(fullfile(outDir, f{1})), delete(fullfile(outDir, f{1})); end
end

msgs = checkcode(file, '-id');
bad  = msgs(arrayfun(@(m) ~isempty(regexp(m.id, '^(SYNER|PARSE)', 'once')) || contains(m.message, 'Parse error'), msgs));
ok = check(ok, isempty(bad), sprintf('gto_dro_study parses clean (%d analyzer messages, %d parse errors)', numel(msgs), numel(bad)));

% (1) full run at 12 N, strict
assignin('base', 'study_opts', struct('thrustN', 12, 'out', struct('dir', outDir, 'visible', 'off'), 'selfCheck', struct('strict', true)));
t0 = tic;  ran = true;
try
    txt = evalc(sprintf('evalin(''base'', ''run(''''%s'''')'')', file));
catch err
    txt = err.message;  ran = false;
end
evalin('base', 'clear study_opts');
ok = check(ok, ran, sprintf('full run at 12 N completed (%.0f s)', toc(t0)));
if ~ran, fprintf('      %s\n', strrep(txt, newline, [newline '      '])); return; end
ok = check(ok, ~isempty(regexp(txt, '== SUMMARY.*12 N', 'once')), 'summary line reports 12 N (study_opts reached the solve)');
ok = check(ok, ~isempty(regexp(txt, '\n\s+12\s+[\d.]+\s+[\d.]+.*\n\s+top-rung', 'once')), 'ladder table''s last row is the 12 N rung');
for id = {'V1 GTO geometry', 'V2 DRO closure', 'V3 ratio screen', 'N0 ladder', 'N1 BVP residual', 'N2 Hamiltonian', ...
          'N3 arrival', 'N4 transversality', 'N5 adjoint eqs', 'N6 min principle', 'N7a Earth floor', 'N7b Moon floor', ...
          'X1 tfMin acceptance', 'X2 field agreement', 'S1 strong Legendre', 'S2 Q_mt > 0', 'S3 dim S = 1', 'S4 conjugate'}
    pat = ['\n\s+' regexptranslate('escape', id{1}) '\s+(PASS|NOT APPLICABLE)'];
    ok = check(ok, ~isempty(regexp(txt, pat, 'once')), sprintf('gate %s is PASS or NOT APPLICABLE', id{1}));
end
for f = {'transfer_3d.png', 'timeseries.png', 'ladder.png', 'gto_dro_study.mat'}
    ok = check(ok, isfile(fullfile(outDir, f{1})), sprintf('%s written', f{1}));
end

% (2) reuse with an impossible Earth floor: N7a must FAIL and the script must throw brokenRoot
matFile = fullfile(outDir, 'gto_dro_study.mat');
assignin('base', 'study_opts', struct('thrustN', 12, 'out', struct('dir', outDir, 'visible', 'off', 'reuse', matFile), ...
                                       'tol', struct('earthKm', 1e6), 'selfCheck', struct('strict', false)));
threw = '';  txt2 = '';
try
    txt2 = evalc(sprintf('evalin(''base'', ''run(''''%s'''')'')', file));
catch err
    threw = err.identifier;  txt2 = err.message;
end
evalin('base', 'clear study_opts');
ok = check(ok, strcmp(threw, 'gto_dro_study:brokenRoot'), sprintf('impossible Earth floor throws gto_dro_study:brokenRoot (got ''%s'')', threw));
ok = check(ok, contains(txt2, 'N7a Earth floor') || contains(txt2, 'N7a'), 'the failure names N7a');

% (3) reuse with a different DRO period is refused
assignin('base', 'study_opts', struct('thrustN', 12, 'dro', struct('tau', 2), 'out', struct('dir', outDir, 'visible', 'off', 'reuse', matFile)));
threw = '';
try
    evalc(sprintf('evalin(''base'', ''run(''''%s'''')'')', file));
catch err
    threw = err.identifier;
end
evalin('base', 'clear study_opts');
ok = check(ok, strcmp(threw, 'gto_dro_study:reuse'), sprintf('reuse with tau = 2 refused with gto_dro_study:reuse (got ''%s'')', threw));

fprintf('test_gto_dro_study: %s\n', tern(ok, 'PASS', 'FAIL'));
end

function ok = check(ok, c, msg)
% CHECK  Print one line and fold the verdict.  INPUTS: ok; c; msg.  OUTPUTS: ok.
if c, s = 'ok  '; else, s = 'FAIL'; end
fprintf('  %s %s\n', s, msg);
ok = ok && c;
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end
```

Note for the implementer: check (2) relies on the N7a failure surfacing in the thrown message; the script's assert message lists the failed IDs, so `contains(txt2, 'N7a')` holds. If the 12 N ladder (two rungs, 15 then 12) exposes a V3 ratio problem or an S4 UNRESOLVED, the strict run fails: report the printed table verbatim as DONE_WITH_CONCERNS rather than loosening a tolerance.

- [ ] **Step 2: Append `'test_gto_dro_study'` to `names` in `run_tests.m`; run `run_tests` (MCP). Expected `GTO_DRO: 3 of 3 passed`.** Record wall time.

- [ ] **Step 3: Write `README.md`**

```markdown
# GTO_DRO -- GTO -> DRO minimum-time low-thrust transfers

The Earth -> DRO leg the orbital network lacks (see `orbit_network/`),
solved in the CATALOG REGIME (150 kg, Isp 1710 s, thrust ladder 15 N down)
through the shared ladder engine in `costate_common/` -- the same
pipeline every costate catalog used: direct solve -> fly the control ->
`ms_tfmin` multiple shooting -> pumpkyn `tfMin` acceptance unchanged. No
solver code lives here. Design: `docs/superpowers/specs/2026-11-09-gto-dro-mintime-design.md`.

Stage 1 (this folder today): ONE transfer with its parameters exposed.
Stage 2 (planned): the catalog chain script over 16 sheets. Stage 3
(planned): the GTO catalogs enter the orbital network.

## Run

```matlab
cd orbit_transfer/GTO_DRO
gto_dro_study              % one transfer: edit the knobs block, press Run (~minutes)
cd tests; run_tests        % clearance instrument, driver planning, the study headless
```

Paths come from `setup_gto_dro_paths` (this dir, lib/, costate_common,
oclib, casadi, pumpkyn + pumpkynPie). Outputs go to `results/<run>/`
(git-ignored): three PNGs and `gto_dro_study.mat` with `run_meta`.

## Knobs (section 0 of the study script)

| knob | default | meaning |
|---|---|---|
| `thrustN` | 15 | returned rung [N]; the ladder walks 15 N down to it; [1, 15] only |
| `m0kg`, `ispS` | 150, 1710 | spacecraft |
| `gto.orientDeg` | 0 | Earth->Moon line to perigee [deg], rotating sense |
| `gto.sma_km`, `gto.ecc` | 24446, 0.7248 | the 350 x 35786 km GTO |
| `gto.depFrac` | 0 | departure phase, TIME fraction since perigee (0.5 = apogee; an apogee cold start at 15 N found a bad local minimum on the tulip campaign) |
| `dro.tau` | 1 | DRO period [ND]; 0.5, 1, 2, 3 are the catalog periods (1 = 4.43 d) |
| `dro.arrFrac` | 0 | arrival phase along the DRO from its +x crossing |
| `recipe.tf0` | 0.30 | top-rung t_f guess [ND], inherited from the tulip campaign; the script prints t_f / tf0 |
| `out.reuse` | '' | gate a saved `gto_dro_study.mat` instead of solving (knobs must match) |
| `out.conjTest` | true | run the conjugate re-solve (S4); false = NOT APPLICABLE |

`study_opts` in the workspace overrides any knob (the tests use it).

## What the study script checks

| id | condition |
|---|---|
| V1 | GTO perigee and apogee radii match a(1 -+ e) to 1 km |
| V2 | DRO closes to 1e-7 ND after cont_np |
| V3 | t_f ratio between consecutive rungs in [1/3, 3] (steps only) |
| N0 | every walked rung accepted |
| N1 | multiple-shooting residual < 1e-8 |
| N2 | max |H| along the flown arc < 1e-6 (autonomous, free time) |
| N3 | flown arrival within 100 km and 10 m/s |
| N4 | |lambda_m(t_f)| < 1e-6 |
| N5 | adjoint equations hold (finite differences of H) |
| N6 | minimum principle for the applied control, throttle 1, Q_mt >= 0 |
| N7a | Earth altitude >= 200 km on the flown arc (the check the engine lacks) |
| N7b | lunar altitude >= 500 km on the flown arc |
| S1-S3 | strong Legendre, strict bang, dim S = 1 (hypothesis gates + lift margin) |
| S4 | no conjugate point: junction test and dense spectrum scan; UNRESOLVED when either cannot decide |
| X1 | pumpkyn tfMin at the accepted z8 returns within 1e-6 (NOT CHECKED if it throws) |
| X2 | pumpkyn's field agrees with an independent field row by row |

V, N and X are code gates: the script throws when one fails. S1-S4 are
results and are reported either way (strict mode throws on them too).

## What it does not claim

- No new solver and no change to `costate_common`; a violated Earth floor
  is REPORTED, not constrained (a stage 2 item).
- No free arrival or departure phase: rendezvous at the chosen states.
- No thrust below 1 N (many-rev regime: junction states, Sundman mesh).
- No sheet writing; the catalog is stage 2.
- A full PASS table says "a certified local minimum-time extremal that
  pumpkyn's tfMin accepts unchanged", never global optimality.
- The GTO is an algebraic Kepler locus frozen at epoch zero in the
  rotating frame: a set of departure states keyed by orientation and
  time-since-perigee, not a trajectory.
```

- [ ] **Step 4: Add the GTO_DRO row to `orbit_transfer/README.md`** next to the GTO_tulip / orbit_network entries (read the file first; match its table or list format; one line: "`GTO_DRO/` -- GTO -> DRO min-time, catalog regime, through the shared ladder engine; stage 1 = single-transfer study script (`gto_dro_study.m`); stages 2-3 planned (catalog, network)").

- [ ] **Step 5: Report** (the controller commits, then writes the memory note `gto-dro-stage1-<date>.md` with the first GTO -> DRO numbers, the t_f / tf0 ratio, the Earth minimum altitude, S4's status, and the stage 2 open items, and adds the MEMORY.md line).

---

## Self-review notes

- Spec coverage: section 3 (folder) Tasks 1 to 5; section 5 (`gto_dro_entry`) Task 2 (plus `planOnly`, `engineOpts` added for testability); section 6 (`earth_clearance`) Task 1; section 7 (sections 0 to 12) Tasks 3 and 4 (the direct-control re-flight is replaced by the engine's G1b read-back, ruling recorded in Global Constraints); section 8 (tests) Tasks 1, 2, 5 (the spec's slow driver test is folded into the study test, which solves at 12 N and so exercises a two-rung ladder, covering V3 as well); section 9 refusals: `thrust` refusal Task 2, no sheet writing Task 2, README Task 5; section 12 acceptance: Task 4 step 2 and Task 5 step 2.
- Type consistency: `gto_dro_entry` returns `entry.z8` [8 x 1], `entry.rv0`/`rvf` [1 x 6] rows; the study passes `rv0(1:6).'` (column) to `fly_transfer`, `mintime_hypothesis_gates`, `conj_spectrum` and rows to `tfMin`/`toPCI`/`minDeltaV` (pumpkyn takes rows). `gates.conjPass` is a double 1/0/NaN. `verdicts` is `{id, status, value, threshold}` throughout and the test regexes match the `%-22s %-15s` print.
- The `study_opts` merge uses `eval` on knob names; that is confined to section 0 and the names are the script's own variables. An alternative without eval would need a knob struct everywhere; kept for readability of the knobs block.
