# Orbit Network (topic 2, data-only) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `orbit_network_study.m`, its lib helpers, tests, figures and the LaTeX note: the low-thrust orbital network over the 56 certified min-time sheets, with the (thrust, time-cap) accessibility plane, L1-halo-to-tulip relay closure, and an off-by-default mirrored-edge option.

**Architecture:** Every computation with logic of its own lives in `orbit_network/lib/` as a function with an `if nargin == 0` demo and a hand-checked test (`edge_tables`, `accessibility_plane`, `relay_closure`, `mirror_state`, `symmetry_residual`); the study script `orbit_network_study.m` wires them in numbered sections, prints every gate PASS/FAIL, writes four PNGs and a `results/network_summary.txt` that the note cites. Topic 1 (`stability_cost_study.m`, `load_catalog_sheets`, `orbit_node_table`, `lib/partial_f_test`) is already built and committed (d5b87d5); this plan consumes it.

**Tech Stack:** MATLAB R2026a in the shared session (MCP tools), `digraph`/`centrality`/`conncomp`, pumpkyn (`pumpkyn.cr3bp.prop`, `.monodromy`, `.jacobi`), costate_common `get_family_orbit`, pdflatex.

**Spec:** `docs/superpowers/specs/2026-11-05-orbit-network-design.md` (sections 7-11). Plan refinement of the spec's layout: two more lib files, `edge_tables.m` and `relay_closure.m`, plus `symmetry_residual.m` beside `mirror_state.m`, because each is a testable unit the spec describes inside the script's sections 2, 5 and 6.

## Global Constraints

- Data-only: no new solves, no changes under `DRO_tulip/`, `HALO_tulip/`, `DPO_tulip/`, `HALO_HALO/`, `costate_common/`.
- Headline numbers never use mirrored edges; `useMirror` defaults to `false` and every number from that section prints with the word APPROXIMATE.
- Relay-versus-direct COST is never claimed from data; only relay CLOSURE is computed.
- MATLAB house style: header quartet (Purpose / Inputs / Outputs / Revision History with `M. Casey (c) 11/05/2026` and `Copyright Coorbital Inc.`), `if nargin == 0` demo on every function, no `i`/`j` as variables, no `norm`, no `%#ok`, no emojis, ASCII-only READMEs and comments. No `M. Casey / D. Koblick` signature lines in docs. Note uses `\author{Michael Casey}` only.
- Tests are function files returning a logical `ok`, printed `PASS`/`FAIL` per check through the same `check`/`verdict` helpers as `tests/test_lib.m`; run from the shared session by `cd` into `tests/` and calling by name. `tests/run_tests.m` lists them in dependency order.
- Paths: scripts and tests add `here`, `here/lib`, `../costate_common`; pumpkyn/pumpkynPie come from the session path. No campaign folder on the path.
- Figures: set axes `Color`, `XColor`, `YColor` explicitly (the shared session defaults axes to dark); node labels and titles with `'Interpreter', 'none'` where a key contains an underscore; PNG at 150 dpi via `exportgraphics`.
- Outputs go to `results/` (git-ignored). Commit messages end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`; implementers never commit, the controller commits after a clean review; git via `/Library/Developer/CommandLineTools/usr/bin/git`; never stage the unrelated `missiles/` deletions.
- Clean LaTeX aux files after every compile: `rm -f *.aux *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz`.
- Constants: muStar 0.012150585609624, tStar_s 382981.289129055 (read from `T(1).constants`, never retyped); days = t_nd * tStar_s / 86400.

## Review Focus

1. A rung on which a sheet has NO closed cell: `tf_min_days` and `tf_med_days` must be NaN, `closed_frac` 0, `n_certified` 0, the accessibility counts must not count it, the graph figure must drop that edge rather than plot a NaN width. Pinned in Task 2 (edge_tables synthetic sheet with an empty rung) and Task 3 (plane counts with NaN).
2. A time-cap grid that lies entirely below the smallest closed t_f: every count is 0, nothing errors. Pinned in Task 3.
3. A node whose stored initial state already sits on the symmetry plane (y = vx = vz = 0): `mirror_state(x0) == x0`, the residual is ~0 at t = 0 AND the search must also find the mirror image near t = T without wrapping errors; and a deliberately perturbed state (y shifted by 1e-3) must FAIL the gate with a residual above 1e-4, so the gate can fail. Pinned in Task 1.
4. A relay route with exactly one closed leg at a rung: `closed` false, `t_total_days` NaN, and the G4 cross-check against `edge_tables` agrees. Pinned in Task 4.
5. The `useMirror = true` path is never exercised by default: the pipeline test runs the script a second time with `study_opts.useMirror = true` and requires G5 PASS and the APPROXIMATE label in the output. Pinned in Task 7.
6. Node keys contain dots and underscores (`halo_L2_tau2.2`): `digraph` node names must round-trip exactly (`findnode` by key) and figure labels must not TeX-interpret the underscore. Pinned in Task 5 (findnode check in the script's G2) and Task 7 (Interpreter none, verified by reading the PNG).

---

### Task 1: `lib/mirror_state.m` and `lib/symmetry_residual.m`

**Files:**
- Create: `orbit_transfer/orbit_network/lib/mirror_state.m`
- Create: `orbit_transfer/orbit_network/lib/symmetry_residual.m`
- Modify: `orbit_transfer/orbit_network/tests/test_lib.m` (append a section before the final `fprintf('test_lib: ...')`)

**Interfaces:**
- Consumes: `pumpkyn.cr3bp.prop(tauVec, x0row, mu)` -> `[tau, X]` with `X` [N x 6] when `x0row` is [1 x 6]; `get_family_orbit('dro', struct('tau', 1))` -> `[tgrid, rv, info]`, initial state `rv(1, :)`, period `info.periodND`.
- Produces: `R = mirror_state(x)` for `x` [6 x n] -> [6 x n], `(x, -y, z, -vx, vy, -vz)`; refusal `orbit_network:mirror_state:shape`. `[res, tStar] = symmetry_residual(x0, T, mu, opts)` for `x0` [6 x 1], returns the minimum over t in [0, T] of `max(abs(phi(x0, t) - mirror_state(x0)))` and the minimizing time; `opts.nGrid` [400].

- [ ] **Step 1: Append the failing tests to `tests/test_lib.m`** (insert immediately before the line `fprintf('test_lib: %s\n', verdict(ok));`):

```matlab
%% mirror_state: the hand vector, a matrix, and the refusal
x  = [1; 2; 3; 4; 5; 6];
ok = check(ok, isequal(mirror_state(x), [1; -2; 3; -4; 5; -6]), 'mirror_state flips y, vx, vz');
ok = check(ok, isequal(mirror_state([x, 2*x]), [[1; -2; 3; -4; 5; -6], 2*[1; -2; 3; -4; 5; -6]]), 'mirror_state acts column-wise on [6 x n]');
ok = check(ok, isequal(mirror_state(mirror_state(x)), x), 'mirror_state is an involution');
try
    mirror_state(x(1:5));
    ok = check(ok, false, 'mirror_state :shape refusal');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:mirror_state:shape'), 'mirror_state :shape refusal on a 5-vector');
end

%% symmetry_residual: a DRO is x-axis symmetric; a perturbed state is not
addpath(fullfile(fileparts(fileparts(here)), 'costate_common'));
[~, ~, rv, info] = evalc('get_family_orbit(''dro'', struct(''tau'', 1))');
x0 = rv(1, :).';  Tp = info.periodND;  mu = 0.012150585609624;
[res, tS] = symmetry_residual(x0, Tp, mu);
ok = check(ok, res < 1e-8, sprintf('DRO tau 1: mirror image lies on the orbit (residual %.1e at t = %.4f of T = %.4f)', res, tS, Tp));
ok = check(ok, tS >= 0 && tS <= Tp, 'minimizing time inside [0, T]');
xp = x0;  xp(2) = xp(2) + 1e-3;                            % off the symmetric orbit
resP = symmetry_residual(xp, Tp, mu);
ok = check(ok, resP > 1e-4, sprintf('perturbed state is NOT symmetric (residual %.1e): the gate can fail', resP));
```

- [ ] **Step 2: Run the test to verify it fails**

In the shared MATLAB session:
```matlab
cd('/Users/msc/Desktop/optimal_control/orbit_transfer/orbit_network/tests'); ok = test_lib();
```
Expected: error `Unrecognized function or variable 'mirror_state'`.

- [ ] **Step 3: Write `lib/mirror_state.m`**

```matlab
function R = mirror_state(x)
%% Purpose:
%
%   The x-axis (x-z plane) symmetry of the CR3BP applied to a state:
%   (x, y, z, vx, vy, vz) -> (x, -y, z, -vx, vy, -vz). With time reversed,
%   the image of a trajectory is a trajectory; the image of a symmetric
%   periodic orbit is the same orbit. Used to build REVERSE transfer edges
%   as an approximation (mass flow is not reversed by the symmetry).
%
%% Inputs:
%
%  x                        [6 x n]                 States as columns
%
%% Outputs:
%
%  R                        [6 x n]                 Mirrored states
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a DRO initial state (x-axis crossing) is its own mirror image
     x0 = [0.85; 0; 0; 0; 0.5; 0];
     R  = mirror_state(x0);
     fprintf('x0 = %s\nR  = %s\n', mat2str(x0.'), mat2str(R.'));
     return
end

if size(x, 1) ~= 6
    error('orbit_network:mirror_state:shape', 'x must be [6 x n] (got %d x %d)', size(x, 1), size(x, 2));
end
R = [1; -1; 1; -1; 1; -1].*x;
end
```

- [ ] **Step 4: Write `lib/symmetry_residual.m`**

```matlab
function [res, tStar] = symmetry_residual(x0, T, mu, opts)
%% Purpose:
%
%   How far the mirror image of a state is from the orbit through it:
%   min over t in [0, T] of max|phi(x0, t) - mirror_state(x0)|. Zero (to
%   integrator accuracy) when the periodic orbit is x-axis symmetric,
%   which is what makes a mirrored transfer edge land on the same orbit.
%   A coarse grid over one period locates the minimum, fminbnd refines it
%   on the bracketing interval.
%
%% Inputs:
%
%  x0                       [6 x 1]                 State on the orbit
%
%  T                        double                  Period (ND)
%
%  mu                       double                  Mass ratio
%
%  opts                     struct (optional)       .nGrid [400] coarse
%                                                   grid points over [0, T]
%
%% Outputs:
%
%  res                      double                  The minimum residual
%                                                   (max-abs over the six
%                                                   components, ND)
%
%  tStar                    double                  The minimizing time
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the L1 halo of pumpkyn's own prop demo
     x0 = [9.2294247990951594E-1; 0; 2.1592634764719765E-1; 4.2453615136710976E-13; 1.2786666755130519E-1; -2.3343675990673871E-12];
     [res, tStar] = symmetry_residual(x0, 1.8036821222727220, 1.215058560962404E-2);
     fprintf('residual %.2e at t = %.4f\n', res, tStar);
     return
end

if nargin < 4 || isempty(opts), opts = struct(); end
if isfield(opts, 'nGrid') && ~isempty(opts.nGrid), nGrid = opts.nGrid; else, nGrid = 400; end
x0  = x0(:);
Rx  = mirror_state(x0);
tg  = linspace(0, T, nGrid);
[~, X] = pumpkyn.cr3bp.prop(tg, x0.', mu);                 % [nGrid x 6]
d   = max(abs(X(:, 1:6).' - Rx), [], 1);                  % [1 x nGrid]
[~, k] = min(d);
lo  = tg(max(k - 1, 1));  hi = tg(min(k + 1, nGrid));
if hi > lo
    f = @(t) resid_at(t, x0, Rx, mu);
    [tStar, res] = fminbnd(f, lo, hi, optimset('TolX', 1e-12));
    if d(k) < res, res = d(k); tStar = tg(k); end           % never worse than the grid
else
    res = d(k);  tStar = tg(k);
end
end

% ------------------------------------------------------------------------
function r = resid_at(t, x0, Rx, mu)
%% Purpose:
%
%   The residual at one time (t = 0 needs no propagation).
%
if t <= 0
    r = max(abs(x0 - Rx));
else
    [~, X] = pumpkyn.cr3bp.prop([0 t], x0.', mu);
    r = max(abs(X(end, 1:6).' - Rx));
end
end
```

- [ ] **Step 5: Run the test to verify it passes**

```matlab
cd('/Users/msc/Desktop/optimal_control/orbit_transfer/orbit_network/tests'); ok = test_lib();
```
Expected: every line PASS, `test_lib: PASS`. If the DRO residual is above 1e-8, raise `nGrid` in the test's call to 1000 before loosening anything else; the DRO closes to 1e-10 in the catalogs.

- [ ] **Step 6: Code Analyzer on both new files** (MCP `check_matlab_code`): no warnings.

- [ ] **Step 7: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/lib/mirror_state.m orbit_transfer/orbit_network/lib/symmetry_residual.m orbit_transfer/orbit_network/tests/test_lib.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: mirror_state (x-axis symmetry) and symmetry_residual (mirror image on the orbit: grid + fminbnd over one period); tests incl. a deliberately perturbed state that fails

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 2: `lib/edge_tables.m`

**Files:**
- Create: `orbit_transfer/orbit_network/lib/edge_tables.m`
- Modify: `orbit_transfer/orbit_network/tests/test_lib.m` (append before the final fprintf)

**Interfaces:**
- Consumes: the sheet table `T` from `load_catalog_sheets([])` (fields `dep_key`, `arr_key`, `dep_family`, `arr_family`, `rungs_N` [1 x nR], `has_solution`/`tf_nd`/`conj_pass` [nD x nA x nR], `constants.tStar_s`).
- Produces: `E = edge_tables(T)`, struct [numel(T) x 1] with `dep_key`, `arr_key`, `dep_family`, `arr_family`, `rungs_N` [1 x nR], `closed_frac` [1 x nR] (closed cells / (nD*nA)), `tf_min_days`, `tf_med_days` [1 x nR] (NaN where nothing closed), `n_certified` [1 x nR] (`has_solution & conj_pass`), `F_floor_N` (smallest rung with closed_frac > 0; NaN if none), `F_floor_str` (`'<= 0.5'` when the last rung closes, else `sprintf('%g', F_floor_N)`).

- [ ] **Step 1: Append the failing test** (before the final fprintf of `test_lib.m`):

```matlab
%% edge_tables: a synthetic 2-sheet table with a hand-built grid
tStar = 382981.289129055;
S = struct('dep_key', 'a', 'arr_key', 'b', 'dep_family', 'dro', 'arr_family', 'tulip', 'rungs_N', [2 1 0.5], ...
           'has_solution', false(2, 2, 3), 'tf_nd', nan(2, 2, 3), 'conj_pass', false(2, 2, 3), ...
           'constants', struct('muStar', 0.012150585609624, 'lStar_km', 389703.264829278, 'tStar_s', tStar));
S.has_solution(:, :, 1) = [true true; true false];  S.tf_nd(:, :, 1) = [1 2; 3 NaN];  S.conj_pass(:, :, 1) = [true false; true false];
S.has_solution(:, :, 2) = [true false; false false]; S.tf_nd(:, :, 2) = [4 NaN; NaN NaN]; S.conj_pass(:, :, 2) = [true false; false false];
S2 = S;  S2.arr_key = 'c';  S2.has_solution(:) = false;  S2.has_solution(1, 1, 3) = true;  S2.tf_nd(1, 1, 3) = 5;
E = edge_tables([S; S2]);
ok = check(ok, numel(E) == 2 && strcmp(E(2).arr_key, 'c'), 'edge_tables: one row per sheet, keys carried');
ok = check(ok, isequal(E(1).closed_frac, [3/4 1/4 0]), 'edge_tables: closed_frac per rung');
ok = check(ok, abs(E(1).tf_min_days(1) - 1*tStar/86400) < 1e-12 && abs(E(1).tf_med_days(1) - 2*tStar/86400) < 1e-12, 'edge_tables: min and median over CLOSED cells only, in days');
ok = check(ok, isnan(E(1).tf_min_days(3)) && isnan(E(1).tf_med_days(3)), 'edge_tables: NaN on a rung with nothing closed');
ok = check(ok, isequal(E(1).n_certified, [2 1 0]), 'edge_tables: n_certified = has_solution & conj_pass');
ok = check(ok, E(1).F_floor_N == 1 && strcmp(E(1).F_floor_str, '1'), 'edge_tables: floor is the deepest rung that closes (1 N here)');
ok = check(ok, E(2).F_floor_N == 0.5 && strcmp(E(2).F_floor_str, '<= 0.5'), 'edge_tables: a sheet closing at the last rung reports the bracket <= 0.5');
```

- [ ] **Step 2: Run to verify it fails** (`ok = test_lib();` -> `Unrecognized function or variable 'edge_tables'`).

- [ ] **Step 3: Write `lib/edge_tables.m`**

```matlab
function E = edge_tables(T)
%% Purpose:
%
%   Per-edge, per-rung attributes of the network built from a sheet
%   table: the closed fraction of the phasing grid, the min and median
%   time of flight over the CLOSED cells (NaN when none), the certified
%   count (closed and conjugate-pass), and the thrust floor -- the deepest
%   rung that closes, reported as a bracket when it is the catalog's last
%   rung, because the data cannot say how much deeper the edge survives.
%
%% Inputs:
%
%  T                        struct [n x 1]          Sheet table
%                                                   (load_catalog_sheets)
%
%% Outputs:
%
%  E                        struct [n x 1]          dep_key, arr_key,
%                                                   dep_family,
%                                                   arr_family, rungs_N
%                                                   [1 x nR], closed_frac,
%                                                   tf_min_days,
%                                                   tf_med_days,
%                                                   n_certified [1 x nR],
%                                                   F_floor_N, F_floor_str
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the 56 edges at 1 N and 0.5 N
     here = fileparts(fileparts(mfilename('fullpath')));
     addpath(here);
     E = edge_tables(load_catalog_sheets([]));
     r1 = find(E(1).rungs_N == 1);  r05 = find(E(1).rungs_N == 0.5);
     fprintf('%-18s %-14s %7s %8s %7s %8s %8s\n', 'dep', 'arr', 'frac1N', 'tf1N[d]', 'frac.5', 'tf.5[d]', 'floor');
     for kk = 1:numel(E)
         fprintf('%-18s %-14s %7.3f %8.3f %7.3f %8.3f %8s\n', E(kk).dep_key, E(kk).arr_key, E(kk).closed_frac(r1), E(kk).tf_min_days(r1), E(kk).closed_frac(r05), E(kk).tf_min_days(r05), E(kk).F_floor_str);
     end
     return
end

n    = numel(T);
rows = cell(n, 1);
for kk = 1:n
    s     = T(kk);
    nR    = numel(s.rungs_N);
    nCell = size(s.has_solution, 1)*size(s.has_solution, 2);
    toDay = s.constants.tStar_s/86400;
    e.dep_key     = s.dep_key;
    e.arr_key     = s.arr_key;
    e.dep_family  = s.dep_family;
    e.arr_family  = s.arr_family;
    e.rungs_N     = s.rungs_N;
    e.closed_frac = zeros(1, nR);
    e.tf_min_days = nan(1, nR);
    e.tf_med_days = nan(1, nR);
    e.n_certified = zeros(1, nR);
    for rr = 1:nR
        h  = s.has_solution(:, :, rr);
        tf = s.tf_nd(:, :, rr);
        e.closed_frac(rr) = nnz(h)/nCell;
        e.n_certified(rr) = nnz(h & s.conj_pass(:, :, rr));
        if any(h(:))
            e.tf_min_days(rr) = min(tf(h))*toDay;
            e.tf_med_days(rr) = median(tf(h))*toDay;
        end
    end
    deep = find(e.closed_frac > 0, 1, 'last');            % rungs run shallow -> deep
    if isempty(deep)
        e.F_floor_N = NaN;  e.F_floor_str = 'none';
    elseif deep == nR
        e.F_floor_N = s.rungs_N(nR);  e.F_floor_str = sprintf('<= %g', s.rungs_N(nR));
    else
        e.F_floor_N = s.rungs_N(deep);  e.F_floor_str = sprintf('%g', s.rungs_N(deep));
    end
    rows{kk} = e;
end
E = vertcat(rows{:});
end
```

- [ ] **Step 4: Run to verify it passes** (`ok = test_lib();` all PASS).
- [ ] **Step 5: Code Analyzer clean on `edge_tables.m`.**
- [ ] **Step 6: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/lib/edge_tables.m orbit_transfer/orbit_network/tests/test_lib.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: edge_tables -- per-edge per-rung closed fraction, min/median t_f over closed cells (NaN when none), certified count, thrust floor as a bracket at the last rung; hand-built synthetic test

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 3: `lib/accessibility_plane.m`

**Files:**
- Create: `orbit_transfer/orbit_network/lib/accessibility_plane.m`
- Modify: `orbit_transfer/orbit_network/tests/test_lib.m` (append before the final fprintf)

**Interfaces:**
- Consumes: sheet table `T` (as in Task 2).
- Produces: `P = accessibility_plane(T, caps_days)` with `P.rungs_N` [1 x nR], `P.caps_days` [1 x nC], `P.pairs` [nR x nC] (sheets whose min closed t_f at that rung is <= cap), `P.cells` [nR x nC] (closed cells with t_f <= cap, summed over sheets), `P.closed` {nR x nC} (row vector of sheet indices counted in `pairs`).

- [ ] **Step 1: Append the failing test** (reuses `S`, `S2`, `tStar` from the Task 2 block, which precedes it in the file):

```matlab
%% accessibility_plane: counts on the synthetic table, monotone in the cap, zero below the data
d = tStar/86400;
P = accessibility_plane([S; S2], d*[0.5 1 2.5 10]);
ok = check(ok, isequal(size(P.pairs), [3 4]) && isequal(size(P.cells), [3 4]) && isequal(size(P.closed), [3 4]), 'accessibility_plane: [nR x nC] arrays');
ok = check(ok, isequal(P.pairs(1, :), [0 1 1 1]), 'accessibility_plane: pairs at rung 1 (min t_f = 1 day for sheet 1; sheet 2 closed nothing)');
ok = check(ok, isequal(P.cells(1, :), [0 1 2 3]), 'accessibility_plane: cells at rung 1 count closed cells under the cap');
ok = check(ok, isequal(P.pairs(3, :), [0 0 0 1]) && isequal(P.cells(3, :), [0 0 0 1]), 'accessibility_plane: rung 3 only sheet 2 (t_f = 5 days) under the 10-day cap');
ok = check(ok, isequal(P.closed{1, 4}, 1) && isequal(P.closed{3, 4}, 2) && isempty(P.closed{1, 1}), 'accessibility_plane: closed sheet lists');
ok = check(ok, all(all(diff(P.pairs, 1, 2) >= 0)) && all(all(diff(P.cells, 1, 2) >= 0)), 'accessibility_plane: non-decreasing along the cap axis');
P0 = accessibility_plane([S; S2], d*[0.1 0.2]);
ok = check(ok, ~any(P0.pairs(:)) && ~any(P0.cells(:)), 'accessibility_plane: caps below the data give all zeros, no error');
```

- [ ] **Step 2: Run to verify it fails** (`Unrecognized function or variable 'accessibility_plane'`).

- [ ] **Step 3: Write `lib/accessibility_plane.m`**

```matlab
function P = accessibility_plane(T, caps_days)
%% Purpose:
%
%   The (thrust rung, time cap) accessibility plane of a sheet table: at
%   each rung and cap, how many (departure, arrival) pairs have a closed
%   cell whose time of flight is under the cap, and how many cells. The
%   pair count is the direct-access regime of Braik & Ross with our axes;
%   the gap between the two counts is phase robustness. Rungs are the
%   catalog's, shallow to deep; caps are in days.
%
%% Inputs:
%
%  T                        struct [n x 1]          Sheet table
%                                                   (load_catalog_sheets)
%
%  caps_days                [1 x nC]                Time caps, days,
%                                                   increasing
%
%% Outputs:
%
%  P                        struct                  .rungs_N [1 x nR];
%                                                   .caps_days [1 x nC];
%                                                   .pairs [nR x nC];
%                                                   .cells [nR x nC];
%                                                   .closed {nR x nC}
%                                                   sheet indices counted
%                                                   in .pairs
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the four catalogs on 10 caps between 1 and 4 days
     here = fileparts(fileparts(mfilename('fullpath')));
     addpath(here);
     P = accessibility_plane(load_catalog_sheets([]), linspace(1, 4, 10));
     fprintf('pairs (rows = rungs %s N, cols = caps %s d):\n', mat2str(P.rungs_N), mat2str(P.caps_days, 3));
     disp(P.pairs);
     return
end

caps  = caps_days(:).';
nC    = numel(caps);
rungs = T(1).rungs_N;
nR    = numel(rungs);
P.rungs_N   = rungs;
P.caps_days = caps;
P.pairs     = zeros(nR, nC);
P.cells     = zeros(nR, nC);
P.closed    = cell(nR, nC);
for rr = 1:nR
    for cc = 1:nC
        P.closed{rr, cc} = zeros(1, 0);
    end
end
for kk = 1:numel(T)
    toDay = T(kk).constants.tStar_s/86400;
    for rr = 1:nR
        h  = T(kk).has_solution(:, :, rr);
        if ~any(h(:)), continue; end
        tf = T(kk).tf_nd(:, :, rr)*toDay;
        tc = tf(h);                                        % closed cells only
        for cc = 1:nC
            nUnder = nnz(tc <= caps(cc));
            P.cells(rr, cc) = P.cells(rr, cc) + nUnder;
            if nUnder > 0
                P.pairs(rr, cc)  = P.pairs(rr, cc) + 1;
                P.closed{rr, cc} = [P.closed{rr, cc}, kk];
            end
        end
    end
end
end
```

- [ ] **Step 4: Run to verify it passes.** 
- [ ] **Step 5: Code Analyzer clean.**
- [ ] **Step 6: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/lib/accessibility_plane.m orbit_transfer/orbit_network/tests/test_lib.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: accessibility_plane -- (rung, time cap) pair and cell counts with the closed-sheet lists; synthetic test incl. caps below the data

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 4: `lib/relay_closure.m`

**Files:**
- Create: `orbit_transfer/orbit_network/lib/relay_closure.m`
- Modify: `orbit_transfer/orbit_network/tests/test_lib.m` (append before the final fprintf)

**Interfaces:**
- Consumes: `E` from `edge_tables` (Task 2) and the node table `N` from `orbit_node_table` (fields `key`, `period_nd`), plus `tStar_s`.
- Produces: `R = relay_closure(E, N, tStar_s, src_prefix, relay_prefix, dst_prefix)` -> struct [nRoutes x 1] with `src_key`, `relay_key`, `dst_key`, `e1`, `e2` (indices into `E`), `closed` [1 x nR] logical (both legs have closed_frac > 0), `t_total_days` [1 x nR] (`tf_min` leg 1 + leg 2, NaN unless closed), `t_wait_days` [1 x nR] (`t_total_days` + the relay's period in days). Routes are every (src, relay, dst) with `src_key` starting `src_prefix` ('halo_L1'), `relay_key` starting `relay_prefix` ('halo_L2'), `dst_key` starting `dst_prefix` ('tulip'), for which both edges exist in `E`. Prefix match is `strncmp(key, prefix, numel(prefix))`.

- [ ] **Step 1: Append the failing test**:

```matlab
%% relay_closure: a 3-node toy (a -> b -> c) with one leg open at rung 2
Ea = struct('dep_key', 'L1_a', 'arr_key', 'L2_b', 'closed_frac', [1 0.5 0], 'tf_min_days', [1 2 NaN]);
Eb = struct('dep_key', 'L2_b', 'arr_key', 'tulip_c', 'closed_frac', [1 0 0.5], 'tf_min_days', [3 NaN 4]);
Ex = struct('dep_key', 'L1_a', 'arr_key', 'tulip_c', 'closed_frac', [1 1 1], 'tf_min_days', [1 1 1]);  % a direct edge: not a relay route
Nt = struct('key', {'L1_a'; 'L2_b'; 'tulip_c'}, 'period_nd', {1; 2; 3});
R  = relay_closure([Ea; Eb; Ex], Nt, tStar, 'L1', 'L2', 'tulip');
ok = check(ok, numel(R) == 1 && strcmp(R.relay_key, 'L2_b') && R.e1 == 1 && R.e2 == 2, 'relay_closure: one route a -> b -> c, edges indexed into E');
ok = check(ok, isequal(R.closed, [true false false]), 'relay_closure: closed iff BOTH legs close (rung 2: leg 2 open; rung 3: leg 1 open)');
ok = check(ok, R.t_total_days(1) == 4 && all(isnan(R.t_total_days(2:3))), 'relay_closure: t_total = sum of leg minima, NaN when not closed');
ok = check(ok, abs(R.t_wait_days(1) - (4 + 2*d)) < 1e-12, 'relay_closure: t_wait adds one relay period in days');
```

- [ ] **Step 2: Run to verify it fails.**

- [ ] **Step 3: Write `lib/relay_closure.m`**

```matlab
function R = relay_closure(E, N, tStar_s, src_prefix, relay_prefix, dst_prefix)
%% Purpose:
%
%   Two-leg routes src -> relay -> dst through the edge tables: closed at
%   a rung iff both legs have a closed cell; the route time is the sum of
%   the legs' min-over-phases times, and the phase hand-over is bounded
%   by one relay period (leg 2's minimum is over its departure phases, so
%   the wait at the relay is at most one period). This measures relay
%   CLOSURE only; relay-versus-direct cost is not a question this data
%   can answer (no pair has both).
%
%% Inputs:
%
%  E                        struct [m x 1]          edge_tables output
%                                                   (dep_key, arr_key,
%                                                   closed_frac,
%                                                   tf_min_days)
%
%  N                        struct [n x 1]          Node table (key,
%                                                   period_nd)
%
%  tStar_s                  double                  Time unit, s
%
%  src_prefix               char                    Key prefix of the
%                                                   sources ('halo_L1')
%
%  relay_prefix             char                    Of the relays
%                                                   ('halo_L2')
%
%  dst_prefix               char                    Of the destinations
%                                                   ('tulip')
%
%% Outputs:
%
%  R                        struct [k x 1]          src_key, relay_key,
%                                                   dst_key, e1, e2,
%                                                   closed [1 x nR],
%                                                   t_total_days,
%                                                   t_wait_days [1 x nR]
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the 32 L1 halo -> L2 halo -> tulip routes of the catalogs, at 1 N
     here = fileparts(fileparts(mfilename('fullpath')));
     addpath(here, fullfile(fileparts(here), 'costate_common'));
     T = load_catalog_sheets([]);
     R = relay_closure(edge_tables(T), orbit_node_table(T, struct('verbose', false)), T(1).constants.tStar_s, 'halo_L1', 'halo_L2', 'tulip');
     r1 = find(T(1).rungs_N == 1);
     for kk = 1:numel(R)
         fprintf('%-18s -> %-16s -> %-10s  closed %d  t_total %.3f d  t_wait %.3f d\n', R(kk).src_key, R(kk).relay_key, R(kk).dst_key, R(kk).closed(r1), R(kk).t_total_days(r1), R(kk).t_wait_days(r1));
     end
     return
end

dep  = {E.dep_key};  arr = {E.arr_key};
keys = {N.key};
isS  = @(k, p) strncmp(k, p, numel(p));
rows = {};
for e1 = 1:numel(E)
    if ~(isS(dep{e1}, src_prefix) && isS(arr{e1}, relay_prefix)), continue; end
    for e2 = 1:numel(E)
        if ~(strcmp(dep{e2}, arr{e1}) && isS(arr{e2}, dst_prefix)), continue; end
        r.src_key   = dep{e1};
        r.relay_key = arr{e1};
        r.dst_key   = arr{e2};
        r.e1        = e1;
        r.e2        = e2;
        r.closed    = E(e1).closed_frac > 0 & E(e2).closed_frac > 0;
        tt          = E(e1).tf_min_days + E(e2).tf_min_days;
        tt(~r.closed) = NaN;
        r.t_total_days = tt;
        kN = find(strcmp(keys, r.relay_key), 1);
        if isempty(kN)
            error('orbit_network:relay_closure:node', 'relay %s is not in the node table', r.relay_key);
        end
        r.t_wait_days  = tt + N(kN).period_nd*tStar_s/86400;
        rows{end+1, 1} = r;
    end
end
if isempty(rows)
    R = struct('src_key', {}, 'relay_key', {}, 'dst_key', {}, 'e1', {}, 'e2', {}, 'closed', {}, 't_total_days', {}, 't_wait_days', {});
else
    R = vertcat(rows{:});
end
end
```

(`rows{end+1, 1}` grows in a loop; the Code Analyzer reports it as an info-level note, which is acceptable; do NOT add `%#ok`. If the reviewer objects, count the routes in a first pass and preallocate.)

- [ ] **Step 4: Run to verify it passes.**
- [ ] **Step 5: Code Analyzer: no warnings (info-level growth note tolerated).**
- [ ] **Step 6: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/lib/relay_closure.m orbit_transfer/orbit_network/tests/test_lib.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: relay_closure -- two-leg routes by key prefix, closed iff both legs close, t_total and the one-relay-period wait bound; toy test with a one-leg-open rung

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 5: `orbit_network_study.m`, sections 0-4 (settings, load, graph, plane, centrality) with G1-G3

**Files:**
- Create: `orbit_transfer/orbit_network/orbit_network_study.m`
- Modify: `orbit_transfer/orbit_network/tests/test_pipeline.m` (the `scripts` cell gets a second row)

**Interfaces:**
- Consumes: `load_catalog_sheets([])`, `orbit_node_table(T, struct('verbose', false))`, `edge_tables(T)`, `accessibility_plane(T, caps)`.
- Produces: the base-workspace variables later sections use: `T`, `N`, `E`, `G` (digraph), `src`, `dst` [56 x 1] node indices, `rungs`, `rC` (index of the 1 N rung), `caps`, `P` (with `P.largest` [nR x nC] added), `gates` struct array, `lines` cellstr (the summary file's lines), `outDir`, `visible`, `useMirror`. The script honours `study_opts` (`.outDir`, `.visible`, `.useMirror`) from the base workspace exactly as `stability_cost_study.m` does.

- [ ] **Step 1: Extend the pipeline test** so it expects the new script. In `tests/test_pipeline.m` replace the single-row `scripts = {...};` with:

```matlab
scripts = {'stability_cost_study', {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7'}, {'H1', 'H2'}, {'stability_cost.png'}
           'orbit_network_study',  {'G1', 'G2', 'G3'},                         {},           {}};
```

- [ ] **Step 2: Run `test_pipeline` to verify it fails** (`orbit_network_study ran without error` FAIL: file not found).

- [ ] **Step 3: Write `orbit_network_study.m` with sections 0-4 and the summary block.** The relay, mirror and figure sections arrive in Tasks 6 and 7; leave the markers `%% 5  relay closure`, `%% 6  mirror option`, `%% 7  figures` as comment lines in place so later tasks insert under them.

```matlab
%% ORBIT_NETWORK_STUDY  The low-thrust orbital network over the costate catalogs
%
%   Braik & Ross (arXiv 2605.31543, 2026) organize periodic orbits into a
%   weighted graph whose edges come from reachable-set overlap under
%   heading-change impulses, then read hubs, relays and budget regimes off
%   it. Here the nodes are the 18 pumpkyn orbits of our four min-time
%   catalogs and the edges are the 56 certified sheets. No proxies, no new
%   solves. The edge attribute that matters for an optimal control problem
%   is FEASIBILITY ACROSS THRUST: a thrust floor and a phase-resolved time
%   of flight; their (dV cap, time cap) budget plane becomes our (thrust,
%   time cap) accessibility plane.
%
%   What this data cannot say, said once: no pair has both a direct edge
%   and a two-leg route, so relay-versus-direct cost is not tested here
%   (in a genuine OCP the direct feasible set contains every relay route,
%   so direct never loses). Relay CLOSURE is tested: L1 halo -> L2 halo
%   -> tulip, the only multileg structure the catalogs hold.
%
%     0  settings                          (one place for every choice)
%     1  load                              (G1 periods, G2 edges + 1 N / 0.5 N nesting)
%     2  graph and edge tables             (digraph; per-rung closed fraction, t_f, floor)
%     3  accessibility plane               (G3 monotone in the cap; connectivity overlay)
%     4  centrality                        (strength, betweenness; only L2 halos can relay)
%     5  relay closure                     (G4 routes agree with their legs)
%     6  mirror option                     (G5 every node is x-axis symmetric; APPROXIMATE edges)
%     7  figures                           (results/*.png)
%
%   G gates are code gates: a FAIL means the data or the computation is
%   not what this script assumes, and the summary throws. Every printed
%   number from section 6 carries the word APPROXIMATE and never enters
%   sections 3-5.
%
%   A pipeline test runs this headless: if a struct `study_opts` exists in
%   the base workspace, its .outDir, .visible and .useMirror are honoured.
%
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.

%% 0  settings
clearvars -except study_opts; clc
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'lib'), fullfile(fileparts(here), 'costate_common'));

F_ref     = 1;                    % N, the rung every sheet closes at: centrality and the graph figure use it
F_deep    = 0.5;                  % N, the catalogs' last rung: the graph figure's second panel
nCaps     = 25;                   % time caps from the smallest to the largest closed t_f (days)
useMirror = false;                % mirrored (reverse) edges: APPROXIMATE, off by default
relay     = struct('src', 'halo_L1', 'via', 'halo_L2', 'dst', 'tulip');   % key prefixes of the only multileg structure in the data
outDir    = fullfile(here, 'results');
visible   = 'on';
if exist('study_opts', 'var')
    if isfield(study_opts, 'outDir'),    outDir    = study_opts.outDir;    end
    if isfield(study_opts, 'visible'),   visible   = study_opts.visible;   end
    if isfield(study_opts, 'useMirror'), useMirror = study_opts.useMirror; end
end
if ~exist(outDir, 'dir'), mkdir(outDir); end
gates = struct('id', {}, 'ok', {}, 'msg', {});
lines = {};                                                 % the summary file, one printed fact per line

%% 1  load
T     = load_catalog_sheets([]);
N     = orbit_node_table(T, struct('verbose', false));
nS    = numel(T);  nN = numel(N);
rungs = T(1).rungs_N;  nR = numel(rungs);
rC    = find(rungs == F_ref, 1);  rD = find(rungs == F_deep, 1);
tStar = T(1).constants.tStar_s;
keyN  = {N.key};
src   = cellfun(@(k) find(strcmp(keyN, k), 1), {T.dep_key}).';
dst   = cellfun(@(k) find(strcmp(keyN, k), 1), {T.arr_key}).';
fprintf('\n== 1  load: %d sheets, %d orbits, rungs %s N\n', nS, nN, mat2str(rungs));

ePer = max([abs([N(src).period_nd] - [T.T_dep]), abs([N(dst).period_nd] - [T.T_arr])]);
nNest = 0;  nViolAll = 0;
for kk = 1:nS
    h = T(kk).has_solution;
    nNest    = nNest    + nnz(h(:, :, rD) & ~h(:, :, rC));
    nViolAll = nViolAll + nnz(h(:, :, 2:end) & ~h(:, :, 1:end-1));
end
perRung = sum(cell2mat(arrayfun(@(r) reshape(sum(sum(r.has_solution, 1), 2), 1, []), T, 'UniformOutput', false)), 1);
lines{end+1} = sprintf('closed cells per rung %s: %s', mat2str(rungs), mat2str(perRung));
lines{end+1} = sprintf('whole-ladder non-monotone cells: %d (15 N starting rung closes fewer than 12 N: solver coverage, not gated)', nViolAll);
fprintf('%s\n%s\n', lines{end-1}, lines{end});
gates(end+1) = struct('id', 'G1', 'ok', ePer < 1e-9, 'msg', sprintf('node periods match the catalogs (max |dT| = %.1e)', ePer));
gates(end+1) = struct('id', 'G2', 'ok', nS == 56 && nNest == 0 && all(src > 0) && all(dst > 0), 'msg', sprintf('56 edges on named nodes; every cell closed at %g N is closed at %g N (%d violations)', F_deep, F_ref, nNest));

%% 2  graph and edge tables
E = edge_tables(T);
nodeTab = table(keyN.', {N.family}.', [N.period_nd].', [N.CJ].', [N.sigma].', [N.stable].', ...
                'VariableNames', {'Name', 'family', 'period_nd', 'CJ', 'sigma', 'stable'});
edgeTab = table([src dst], vertcat(E.closed_frac), vertcat(E.tf_min_days), vertcat(E.tf_med_days), vertcat(E.n_certified), [E.F_floor_N].', ...
                'VariableNames', {'EndNodes', 'closed_frac', 'tf_min_days', 'tf_med_days', 'n_certified', 'F_floor_N'});
G = digraph(edgeTab, nodeTab);
fprintf('\n== 2  graph: %d nodes, %d edges\n', numnodes(G), numedges(G));
for rr = [rC rD]
    fprintf('\n  edge table at %g N\n', rungs(rr));
    fprintf('  %-18s %-16s %7s %8s %8s %6s %8s\n', 'departure', 'arrival', 'closed', 'tf_min', 'tf_med', 'cert', 'floor');
    for kk = 1:nS
        fprintf('  %-18s %-16s %7.3f %8.3f %8.3f %6d %8s\n', E(kk).dep_key, E(kk).arr_key, E(kk).closed_frac(rr), E(kk).tf_min_days(rr), E(kk).tf_med_days(rr), E(kk).n_certified(rr), E(kk).F_floor_str);
    end
end
nFloorDeep = nnz(strncmp({E.F_floor_str}, '<=', 2));
lines{end+1} = sprintf('edges whose thrust floor is a bracket at the last rung (<= %g N): %d of %d', F_deep, nFloorDeep, nS);
fprintf('%s\n', lines{end});
% G2 also pins the keys round-tripping through the digraph (dots and underscores in names)
gates(end).ok  = gates(end).ok && all(findnode(G, {T.dep_key}) == src.') && all(findnode(G, {T.arr_key}) == dst.');

%% 3  accessibility plane
tfAll = vertcat(E.tf_min_days);
caps  = linspace(min(tfAll(:), [], 'omitnan'), max(tfAll(:), [], 'omitnan'), nCaps);
P     = accessibility_plane(T, caps);
P.largest = zeros(nR, nCaps);                              % largest weak component at (rung, cap)
for rr = 1:nR
    for cc = 1:nCaps
        idx = P.closed{rr, cc};
        Grc = digraph(src(idx), dst(idx), ones(numel(idx), 1), nN);
        [~, binsize] = conncomp(Grc, 'Type', 'weak');
        P.largest(rr, cc) = max(binsize);
    end
end
fprintf('\n== 3  accessibility plane: %d rungs x %d caps (%.3f .. %.3f days)\n', nR, nCaps, caps(1), caps(end));
fprintf('  pairs closed (rows = rungs, cols = caps):\n');  disp(P.pairs);
fprintf('  largest weak component (of %d nodes):\n', nN);  disp(P.largest);
monoCap = all(all(diff(P.pairs, 1, 2) >= 0)) && all(all(diff(P.cells, 1, 2) >= 0)) && all(all(diff(P.largest, 1, 2) >= 0));
[~, cFull] = max(P.largest(rC, :) == nN);
lines{end+1} = sprintf('at %g N the graph is weakly connected (all %d nodes) from a %.3f-day cap on', F_ref, nN, caps(cFull));
lines{end+1} = sprintf('at %g N, %d of %d pairs close under the largest cap', F_deep, P.pairs(rD, end), nS);
fprintf('%s\n%s\n', lines{end-1}, lines{end});
gates(end+1) = struct('id', 'G3', 'ok', monoCap, 'msg', 'pair, cell and component counts non-decreasing along the cap axis');

%% 4  centrality (at the reference rung, where every edge is closed)
cost  = tfAll(:, rC);
outS  = accumarray(src, 1./cost, [nN 1]);
inS   = accumarray(dst, 1./cost, [nN 1]);
bc    = centrality(G, 'betweenness', 'Cost', cost);
fprintf('\n== 4  centrality at %g N (weights 1 / t_f[days]; betweenness with t_f as cost)\n', F_ref);
fprintf('  %-20s %6s %6s %9s %9s %11s\n', 'orbit', 'in', 'out', 'in-str', 'out-str', 'betweenness');
for kk = 1:nN
    fprintf('  %-20s %6d %6d %9.3f %9.3f %11.3f\n', keyN{kk}, N(kk).in_degree, N(kk).out_degree, inS(kk), outS(kk), bc(kk));
end
isL2 = strncmp(keyN, relay.via, numel(relay.via));
lines{end+1} = sprintf('betweenness is nonzero only on %s nodes (%d of them): they are the only orbits with both in- and out-edges', relay.via, nnz(isL2));
fprintf('%s\n', lines{end});
if any(bc(~isL2) > 0)
    lines{end+1} = 'NOTE: a non-relay node has nonzero betweenness -- the graph has grown beyond the four catalogs; re-read section 4';
    fprintf('%s\n', lines{end});
end

%% 5  relay closure

%% 6  mirror option

%% 7  figures

%% summary
fprintf('\n== SUMMARY\n');
for kk = 1:numel(gates)
    if gates(kk).ok, s = 'PASS'; else, s = 'FAIL'; end
    fprintf('  %s  %s  %s\n', gates(kk).id, s, gates(kk).msg);
end
fid = fopen(fullfile(outDir, 'network_summary.txt'), 'w');
fprintf(fid, 'orbit_network_study summary, %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
for kk = 1:numel(lines), fprintf(fid, '%s\n', lines{kk}); end
for kk = 1:numel(gates)
    if gates(kk).ok, s = 'PASS'; else, s = 'FAIL'; end
    fprintf(fid, '%s %s %s\n', gates(kk).id, s, gates(kk).msg);
end
fclose(fid);
fprintf('summary: %s\n', fullfile(outDir, 'network_summary.txt'));
if ~all([gates.ok])
    error('orbit_network:orbit_network_study:gate', 'a G gate failed; see the summary');
end
```

- [ ] **Step 4: Run the script in the shared session** (`run_matlab_file`) and read the output: G1-G3 PASS, the 1 N edge table has 56 rows, the plane prints, betweenness nonzero only on `halo_L2_*`. Then `cd tests; ok = test_pipeline();` -> all PASS including the three G gates for the new script.

- [ ] **Step 5: Code Analyzer on the script**: no warnings (info-level growth notes on `gates(end+1)`/`lines{end+1}` tolerated, no `%#ok`).

- [ ] **Step 6: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/orbit_network_study.m orbit_transfer/orbit_network/tests/test_pipeline.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: orbit_network_study sections 0-4 -- digraph of 18 orbits / 56 certified edges with per-rung closed fraction, t_f and thrust-floor attributes, the (thrust, time-cap) accessibility plane with the weak-connectivity overlay, strength and betweenness at 1 N; G1-G3; summary file; pipeline test row

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 6: sections 5-6 (relay closure with G4, mirror option with G5)

**Files:**
- Modify: `orbit_transfer/orbit_network/orbit_network_study.m` (fill the `%% 5` and `%% 6` markers)
- Modify: `orbit_transfer/orbit_network/tests/test_pipeline.m` (gate list for the network row becomes `{'G1', 'G2', 'G3', 'G4', 'G5'}`; add the second, mirrored run described in Step 1)

**Interfaces:**
- Consumes: Task 5's workspace (`T`, `N`, `E`, `G`, `src`, `dst`, `P`, `caps`, `rungs`, `rC`, `tStar`, `relay`, `useMirror`, `gates`, `lines`), `relay_closure(E, N, tStar, relay.src, relay.via, relay.dst)`, `symmetry_residual(x0, T, mu)`, `mirror_state`.
- Produces: `R` (routes), `P.multileg` [nR x nCaps] (number of (src, dst) pairs whose best route time is under the cap), `best` struct array over (src, dst) pairs with `t_total_days` [1 x nR] (min over relays, NaN when none closed), `deepest_N` (deepest rung with any closed route), and, when `useMirror`, `Gm` (the digraph with the 56 reversed edges, edge variable `approximate` logical).

- [ ] **Step 1: Extend the pipeline test.** Change the network row's gate list to `{'G1', 'G2', 'G3', 'G4', 'G5'}`. Then, after the `for ss = 1:size(scripts, 1)` loop and before the final `fprintf('test_pipeline: ...')`, add the mirrored run:

```matlab
%% the mirror option, which the default run never exercises
assignin('base', 'study_opts', struct('outDir', outDir, 'visible', 'off', 'useMirror', true));
try
    txt = evalc(sprintf('evalin(''base'', ''run(''''%s'''')'')', fullfile(root, 'orbit_network_study.m')));
    ran = true;
catch err
    txt = err.message;  ran = false;
end
evalin('base', 'clear study_opts');
ok = check(ok, ran, 'orbit_network_study ran with useMirror = true');
ok = check(ok, ran && ~isempty(regexp(txt, '\n\s*G5\s+PASS', 'once')), 'mirrored run: G5 PASS');
ok = check(ok, ran && contains(txt, 'APPROXIMATE'), 'mirrored run labels its numbers APPROXIMATE');
ok = check(ok, ran && contains(txt, '112 edges'), 'mirrored run reports 112 edges (56 + 56 reversed)');
```

- [ ] **Step 2: Run `test_pipeline` to verify the new checks fail** (G4/G5 absent, mirrored run has no `112 edges`).

- [ ] **Step 3: Fill section 5** (replace the `%% 5  relay closure` marker line with):

```matlab
%% 5  relay closure (L1 halo -> L2 halo -> tulip: the only multileg structure in the data)
R  = relay_closure(E, N, tStar, relay.src, relay.via, relay.dst);
fprintf('\n== 5  relay closure: %d routes %s -> %s -> %s\n', numel(R), relay.src, relay.via, relay.dst);
pairKeys = unique(strcat({R.src_key}, ' -> ', {R.dst_key}), 'stable');
best = struct('pair', pairKeys, 't_total_days', [], 't_wait_days', [], 'relay_key', [], 'deepest_N', []);
for pp = 1:numel(pairKeys)
    sel = find(strcmp(strcat({R.src_key}, ' -> ', {R.dst_key}), pairKeys{pp}));
    tt  = vertcat(R(sel).t_total_days);                     % [nRelays x nR]
    tw  = vertcat(R(sel).t_wait_days);
    [tb, ib] = min(tt, [], 1, 'omitnan');                   % best relay per rung
    tb(all(isnan(tt), 1)) = NaN;
    best(pp).t_total_days = tb;
    best(pp).t_wait_days  = arrayfun(@(rr) tw(ib(rr), rr), 1:nR);
    best(pp).relay_key    = arrayfun(@(rr) R(sel(ib(rr))).relay_key, 1:nR, 'UniformOutput', false);
    deep = find(~isnan(tb), 1, 'last');
    if isempty(deep), best(pp).deepest_N = NaN; else, best(pp).deepest_N = rungs(deep); end
    fprintf('  %-32s deepest rung closed %5g N   at %g N: best relay %-16s t_total %.3f d, t_wait %.3f d\n', ...
            pairKeys{pp}, best(pp).deepest_N, F_ref, best(pp).relay_key{rC}, tb(rC), best(pp).t_wait_days(rC));
end
P.multileg = zeros(nR, nCaps);
for rr = 1:nR
    for cc = 1:nCaps
        P.multileg(rr, cc) = nnz(arrayfun(@(b) b.t_total_days(rr) <= caps(cc), best));
    end
end
% G4: a route is closed at a rung iff both of its edges are, read back from E
g4 = true;
for kk = 1:numel(R)
    legs = E(R(kk).e1).closed_frac > 0 & E(R(kk).e2).closed_frac > 0;
    g4 = g4 && isequal(R(kk).closed, legs) && all(isnan(R(kk).t_total_days(~legs))) && all(~isnan(R(kk).t_total_days(legs)));
end
lines{end+1} = sprintf('%d relay routes over %d (src, dst) pairs; pairs reaching %g N by relay: %d of %d', numel(R), numel(best), F_deep, nnz([best.deepest_N] == F_deep), numel(best));
fprintf('%s\n', lines{end});
gates(end+1) = struct('id', 'G4', 'ok', g4, 'msg', sprintf('all %d routes closed iff both legs closed, times NaN exactly when open', numel(R)));
```

- [ ] **Step 4: Fill section 6** (replace the `%% 6  mirror option` marker line with):

```matlab
%% 6  mirror option: reverse edges by the x-axis symmetry (APPROXIMATE: mass flow is not reversed)
mu   = T(1).constants.muStar;
resS = zeros(nN, 1);
for kk = 1:nN
    resS(kk) = symmetry_residual(N(kk).x0, N(kk).period_nd, mu);
end
fprintf('\n== 6  mirror option: x-axis symmetry residual per node (max %.1e)\n', max(resS));
gates(end+1) = struct('id', 'G5', 'ok', all(resS < 1e-8), 'msg', sprintf('every node is x-axis symmetric: mirror image on the orbit to %.1e (bar 1e-8)', max(resS)));
if useMirror
    edgeM = G.Edges;
    edgeM.EndNodes = fliplr(edgeM.EndNodes);
    edgeAll = [G.Edges; edgeM];
    edgeAll.approximate = [false(nS, 1); true(nS, 1)];
    Gm = digraph(edgeAll, G.Nodes);
    A1 = adjacency(Gm) > 0;  A1(logical(eye(nN))) = false;
    A2 = (A1*A1) > 0;  A2(logical(eye(nN))) = false;
    nTwoLeg = nnz(A2 & ~A1);
    fprintf('  APPROXIMATE: %d edges (%d + %d reversed); mutually reachable pairs %d; ordered pairs reachable in exactly two legs but not one: %d\n', numedges(Gm), nS, nS, nS, nTwoLeg);
    lines{end+1} = sprintf('APPROXIMATE (mirror on): %d edges; two-leg-only ordered pairs %d; never used in sections 3-5', numedges(Gm), nTwoLeg);
else
    fprintf('  mirror edges off (useMirror = false); G5 is still checked so the option is always available\n');
end
```

- [ ] **Step 5: Run the script, then `test_pipeline`.** Expected: G4 and G5 PASS; the mirrored run prints `112 edges` and `APPROXIMATE`; all pipeline checks PASS. If G5 fails on a node, print `resS` with keys and check that node's `x0` is the catalog's stored initial state (it is `rv(1, :)` of `get_family_orbit`), then raise `nGrid` to 1000 before anything else; a halo residual of order 1e-9 is normal, 1e-3 is a bug.

- [ ] **Step 6: Code Analyzer clean (same tolerance as Task 5).**
- [ ] **Step 7: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/orbit_network_study.m orbit_transfer/orbit_network/tests/test_pipeline.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: sections 5-6 -- relay closure L1 halo -> L2 halo -> tulip (best relay per rung, deepest rung reached, one-period wait bound, multileg overlay for the plane; G4) and the mirror option (G5 x-axis symmetry of every node; reversed edges APPROXIMATE, off by default, exercised by the pipeline test)

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 7: section 7, the four figures

**Files:**
- Modify: `orbit_transfer/orbit_network/orbit_network_study.m` (fill the `%% 7  figures` marker)
- Modify: `orbit_transfer/orbit_network/tests/test_pipeline.m` (the network row's figure list becomes `{'network_graph.png', 'accessibility_plane.png', 'edge_survival.png', 'relay_closure.png'}`)

**Interfaces:**
- Consumes: `G`, `E`, `P` (`pairs`, `largest`, `multileg`, `caps_days`), `best`, `rungs`, `rC`, `rD`, `nN`, `keyN`, `outDir`, `visible`, `useMirror`, `Gm` when present.
- Produces: the four PNGs in `outDir`, and `lines` entries naming them.

- [ ] **Step 1: Extend the pipeline test's figure list** for the network row as above; run it to see the four `wrote ...png` checks FAIL.

- [ ] **Step 2: Fill section 7** (replace the `%% 7  figures` marker line with):

```matlab
%% 7  figures
axStyle = @(ax) set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
famOf   = @(k) regexprep(k, '_tau.*|_Np.*', '');          % dro, dpo, halo_L1, halo_L2, tulip
depFam  = cellfun(famOf, {E.dep_key}, 'UniformOutput', false);
famList = {'dro', 'dpo', 'halo_L1', 'halo_L2'};

% F1: the graph at the reference rung and at the deepest rung; width = closed fraction, color = t_f
fig = figure('Color', 'w', 'Position', [60 60 1400 620], 'Visible', visible);
for pp = 1:2
    rr = [rC rD];  rr = rr(pp);
    cfAll = vertcat(E.closed_frac);  tfMin = vertcat(E.tf_min_days);
    keep = find(cfAll(:, rr) > 0);
    Gr = digraph(src(keep), dst(keep), ones(numel(keep), 1), nN);
    ax = subplot(1, 2, pp);  axStyle(ax);
    h  = plot(ax, Gr, 'Layout', 'layered', 'Direction', 'right', 'NodeLabel', keyN, 'Interpreter', 'none', ...
              'MarkerSize', 6, 'NodeColor', [0.1 0.3 0.7], 'ArrowSize', 8);
    h.LineWidth = 0.5 + 4*cfAll(keep, rr);
    h.EdgeCData = tfMin(keep, rr);
    colormap(ax, 'parula');  cb = colorbar(ax);  cb.Label.String = 'min t_f over phases [days]';
    title(ax, sprintf('%g N: %d edges, width = closed fraction of the 6 x 6 grid', rungs(rr), numel(keep)), 'Color', 'k', 'Interpreter', 'none');
end
f1 = fullfile(outDir, 'network_graph.png');  exportgraphics(fig, f1, 'Resolution', 150);  if strcmp(visible, 'off'), close(fig); end

% F2: the accessibility plane with the connectivity and multileg-closure contours
fig = figure('Color', 'w', 'Position', [60 60 900 560], 'Visible', visible);
ax  = axes(fig);  axStyle(ax);  hold(ax, 'on');
imagesc(ax, caps, 1:nR, P.pairs);  axis(ax, 'tight');  set(ax, 'YDir', 'normal', 'YTick', 1:nR, 'YTickLabel', arrayfun(@(r) sprintf('%g', r), rungs, 'UniformOutput', false));
colormap(ax, 'hot');  cb = colorbar(ax);  cb.Label.String = sprintf('pairs closed (of %d)', nS);
contour(ax, caps, 1:nR, P.largest, [nN nN], 'c-', 'LineWidth', 2);
if max(P.multileg(:)) > 0, contour(ax, caps, 1:nR, P.multileg, [0.5 0.5], 'g--', 'LineWidth', 2); end
xlabel(ax, 'time cap [days]', 'Color', 'k');  ylabel(ax, 'thrust rung [N]', 'Color', 'k');
title(ax, 'accessibility: pairs closed; cyan = all 18 nodes weakly connected; green dashed = first relay route closed', 'Color', 'k');
f2 = fullfile(outDir, 'accessibility_plane.png');  exportgraphics(fig, f2, 'Resolution', 150);  if strcmp(visible, 'off'), close(fig); end

% F3: edge survival down the ladder, by departure family
fig = figure('Color', 'w', 'Position', [60 60 900 480], 'Visible', visible);
ax  = axes(fig);  axStyle(ax);  hold(ax, 'on');  grid(ax, 'on');
cf  = vertcat(E.closed_frac);
for ff = 1:numel(famList)
    sel = strcmp(depFam, famList{ff});
    plot(ax, 1:nR, mean(cf(sel, :), 1), '-o', 'LineWidth', 1.8, 'DisplayName', sprintf('%s departures (%d edges)', famList{ff}, nnz(sel)));
end
set(ax, 'XTick', 1:nR, 'XTickLabel', arrayfun(@(r) sprintf('%g', r), rungs, 'UniformOutput', false));
xlabel(ax, 'thrust rung [N] (shallow to deep)', 'Color', 'k');  ylabel(ax, 'mean closed fraction of the phasing grid', 'Color', 'k');
legend(ax, 'Location', 'southwest', 'Interpreter', 'none', 'TextColor', 'k', 'Color', 'w');
title(ax, 'edge survival down the ladder (the 15 N dip is ladder coverage, not physics)', 'Color', 'k');
f3 = fullfile(outDir, 'edge_survival.png');  exportgraphics(fig, f3, 'Resolution', 150);  if strcmp(visible, 'off'), close(fig); end

% F4: relay closure, best route per (L1 halo, tulip) pair against thrust
fig = figure('Color', 'w', 'Position', [60 60 900 520], 'Visible', visible);
ax  = axes(fig);  axStyle(ax);  hold(ax, 'on');  grid(ax, 'on');
cols = parula(numel(best));                                % NOT lines(): that name is the script's summary cellstr
for pp = 1:numel(best)
    plot(ax, 1:nR, best(pp).t_total_days, '-o',  'Color', cols(pp, :), 'LineWidth', 1.6, 'DisplayName', [best(pp).pair ' (t_total)']);
    plot(ax, 1:nR, best(pp).t_wait_days,  '--',  'Color', cols(pp, :), 'LineWidth', 1.0, 'HandleVisibility', 'off');
end
set(ax, 'XTick', 1:nR, 'XTickLabel', arrayfun(@(r) sprintf('%g', r), rungs, 'UniformOutput', false));
xlabel(ax, 'thrust rung [N]', 'Color', 'k');  ylabel(ax, 'two-leg time [days]; dashed = + one relay period', 'Color', 'k');
legend(ax, 'Location', 'northwest', 'Interpreter', 'none', 'TextColor', 'k', 'Color', 'w', 'FontSize', 8);
title(ax, 'relay closure L1 halo -> L2 halo -> tulip: best relay per pair and rung', 'Color', 'k', 'Interpreter', 'none');
f4 = fullfile(outDir, 'relay_closure.png');  exportgraphics(fig, f4, 'Resolution', 150);  if strcmp(visible, 'off'), close(fig); end
lines{end+1} = sprintf('figures: %s, %s, %s, %s', f1, f2, f3, f4);
fprintf('%s\n', lines{end});
```

The relay figure's colors come from `parula`, not `lines`: `lines` is the script's summary cellstr, and calling it as a function there would index the cell and error.

- [ ] **Step 3: Run the script with figures visible, then read each PNG** with the Read tool: F1 node labels show underscores literally (no subscripts), no NaN-width edges at 0.5 N; F2 has the cyan contour; F3 four lines with the 15 N dip; F4 eight solid lines with dashed companions. Fix what the eye catches before running the test.

- [ ] **Step 4: Run `test_pipeline` and `run_tests`**: all PASS.
- [ ] **Step 5: Code Analyzer clean.**
- [ ] **Step 6: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/orbit_network_study.m orbit_transfer/orbit_network/tests/test_pipeline.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: section 7 figures -- the graph at 1 N and 0.5 N (width = closed fraction, color = t_f), the accessibility plane with connectivity and multileg contours, edge survival by departure family, relay closure per (L1 halo, tulip) pair; pipeline test checks all four PNGs

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 8: README and the folder's runner registration

**Files:**
- Modify: `orbit_transfer/orbit_network/README.md` (replace the paragraph beginning `Topic 2, the network itself` and extend the file table)
- Modify: `orbit_transfer/README.md` (the `**Across the catalogs (2026-11-05)` paragraph: drop "and, next," so it reads as built)

**Interfaces:** none (documentation).

- [ ] **Step 1: In `orbit_network/README.md`** add rows to the `What is here` table:

```markdown
| `orbit_network_study.m` | topic 2: digraph of 18 orbits and 56 certified edges with per-rung closed fraction, min/median t_f and thrust floor; (thrust, time-cap) accessibility plane with weak-connectivity and relay-closure overlays; strength and betweenness at 1 N; L1 halo -> L2 halo -> tulip relay closure; x-axis mirror option (APPROXIMATE, off by default); G1-G5 gates; four figures and `results/network_summary.txt` |
| `lib/edge_tables.m`, `lib/accessibility_plane.m`, `lib/relay_closure.m`, `lib/mirror_state.m`, `lib/symmetry_residual.m` | the network's computations, each with a demo and a hand-checked test in `tests/test_lib.m` |
| `doc/orbit_network_note.tex` | the note (Task 9) |
```

and replace the `Topic 2 ... is not yet built.` paragraph with:

```markdown
Both scripts print every gate PASS/FAIL and throw on a failed code gate;
hypothesis verdicts (H1-H2) are results either way. `results/network_summary.txt`
holds the numbers the note cites.
```

Add under `## Run`:

```matlab
orbit_network_study         % topic 2: the network, the accessibility plane, relay closure
```

- [ ] **Step 2: In `orbit_transfer/README.md`** change `and, next, the data-only low-thrust network` to `and the data-only low-thrust network`.

- [ ] **Step 3: Check both files are ASCII** (`LC_ALL=C grep -nP '[^\x00-\x7F]' orbit_transfer/orbit_network/README.md` prints nothing; the top-level README already contains arrows and is left as it is).

- [ ] **Step 4: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/README.md orbit_transfer/README.md
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: README for topic 2 (what the network script computes, what it does not claim); orbit_transfer README wording

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 9: the note, `doc/orbit_network_note.tex`

**Files:**
- Create: `orbit_transfer/orbit_network/doc/orbit_network_note.tex`
- Create: `orbit_transfer/orbit_network/doc/figures/` (copies of the five PNGs: `stability_cost.png`, `network_graph.png`, `accessibility_plane.png`, `edge_survival.png`, `relay_closure.png`, copied from `results/` after a fresh run of both scripts)
- Produce: `orbit_transfer/orbit_network/doc/orbit_network_note.pdf` (PDFs are git-ignored by the repo's `*.pdf` rule; the `.tex` and the figures are committed)

**Interfaces:**
- Consumes: `results/network_summary.txt` and the printed output of `stability_cost_study.m` (the H1/H2 lines and the two regression blocks) for every number quoted. No number in the note is typed from memory; each comes from one of those two sources, run fresh in this task.

- [ ] **Step 1: Run both scripts fresh** in the shared session (`run_matlab_file` on each), copy the five PNGs into `doc/figures/`, and keep the two outputs open while writing.

- [ ] **Step 2: Write the note** with this skeleton (fill each section's prose from the sources named; the lemma and the fixed sentences are given verbatim):

```latex
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb,amsthm,graphicx,booktabs,hyperref}
\newtheorem{lemma}{Lemma}
\title{A Low-Thrust Orbital Network over the Costate Catalogs:\\ what Braik and Ross's framing says about certified min-time transfers}
\author{Michael Casey}
\date{November 2026}
\begin{document}
\maketitle

\section{Motivation}
% One paragraph: Braik and Ross (arXiv 2605.31543) build a weighted graph of 13 planar Earth-Moon
% periodic orbits on one Jacobi level from reachable-set overlap under heading-change impulses
% (turn cost 2 v sin(|delta|/2), symmetric half-budget split), and read hubs (C32), relays, and three
% budget regimes (direct access, connectedness, multileg closure) off it. Their stability reading:
% unstable orbits share a chaotic sea and connect cheaply; stable ones sit in tori and are hard.
% One sentence on why we care: the same question for our DRO -> tulip transfers, whose thrust wall at
% 0.09 -> 0.067 N we have measured.

\section{Their problem and ours}
% Table: planar vs spatial; one Jacobi level vs Jacobi gaps up to 0.33; energy-preserving impulses vs
% mass-varying continuous thrust; proxy overlap cost vs certified min-time extremals (multiple shooting,
% conjugate test, hypothesis gates); phase as future work vs 6 x 6 phasing grids per sheet.
\begin{lemma}[direct never loses]
Let $\mathcal{F}_{AB}$ be the set of admissible controls steering orbit $A$ to orbit $B$ in the
mass-varying low-thrust model, and let $J$ be any objective that is additive along a trajectory
(time of flight, propellant). For any relay orbit $C$, every concatenation of an admissible
$A \to C$ transfer, a coast on $C$, and an admissible $C \to B$ transfer is itself admissible for
$A \to B$. Hence $\inf_{\mathcal{F}_{AB}} J \le \inf J(\text{relay through } C)$.
\end{lemma}
\begin{proof}
Concatenation of admissible controls with a zero-thrust coast is admissible (the state is continuous and
the thrust bound holds piecewise), and $J$ adds over the pieces. The relay route is one element of
$\mathcal{F}_{AB}$; the infimum over the whole set is no larger.
\end{proof}
% One sentence: so a network over an OCP is about FEASIBILITY ACROSS THRUST (floors, phase robustness),
% never about relay savings; the paper's multileg-cheaper headline is a property of its restricted
% maneuver model.

\section{The data}
% The four catalogs, the sheet table (56 rows: 16 + 16 + 16 + 8), rungs 15 -> 0.5 N, the 6 x 6 grids,
% certification. The node table (18 orbits) as a booktabs table: key, period, C_J, |lambda|max, sigma,
% stable -- copied from stability_cost_study section 2.
% The closure fact: nested between 1 N and 0.5 N, non-monotone along the ladder at 15 N (quote the
% per-rung counts line from network_summary.txt).

\section{The graph}
% Nodes, edges, attributes (closed fraction, min/median t_f, certified count, thrust floor as a bracket).
% Figure: network_graph.png. Quote: how many edges have floor <= 0.5 N (summary line).
% Centrality: the table from section 4 of the script, and the sentence that only L2 halos can have
% nonzero betweenness because they are the only nodes with both in- and out-edges.

\section{The accessibility plane}
% Figure: accessibility_plane.png. The three regimes in our axes: direct access (pairs closed),
% connectedness (cyan contour: the cap from which all 18 nodes are one weak component at 1 N, quoted),
% multileg closure (green: first relay route). Figure: edge_survival.png, with the 15 N dip explained.

\section{Relay closure}
% Figure: relay_closure.png. The 32 routes, the best relay per (L1 halo, tulip) pair at 1 N, the deepest
% rung each pair reaches by relay (summary line), the one-period wait bound. State again that
% relay-versus-direct cost is untestable here and why (no direct L1 halo -> tulip edge exists).

\section{Stability versus cost}
% Figure: stability_cost.png. The two responses, the nested regression, the H1/H2 verdict lines verbatim,
% the group means, the sigma_arr confound (only the 8 halo-halo sheets), and the regime caveat
% (1 N and 0.5 N sit far above the 70 mN regime of the DRO wall).

\section{What the legs campaign would add}
% DRO -> DPO and DRO -> L2 halo legs down toward 70 mN would make DRO -> tulip the first pair with both a
% direct edge and a relay route; the question is relay CLOSURE below the 67 mN direct wall, which the
% lemma does not settle because the wall is a basin/topology wall, not a cost. Which relay orbit the
% present data points to (read off section 5's best relays and section 4's strengths).

\bibliographystyle{plain}
\begin{thebibliography}{9}
\bibitem{braik2026} S.~Braik and S.~D.~Ross, ``Orbital Networks in the Three-Body Problem,'' arXiv:2605.31543, 2026.
\end{thebibliography}
\end{document}
```

Write the prose under every `%` comment from the named source; delete the comments as each is filled. Six pages is the target; the figures are `\includegraphics[width=\textwidth]{figures/<name>.png}`.

- [ ] **Step 3: Compile and clean**

```bash
cd /Users/msc/Desktop/optimal_control/orbit_transfer/orbit_network/doc
/Library/TeX/texbin/pdflatex -interaction=nonstopmode orbit_network_note.tex > /dev/null && /Library/TeX/texbin/pdflatex -interaction=nonstopmode orbit_network_note.tex > /dev/null
rm -f *.aux *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz
ls -la orbit_network_note.pdf
```
Expected: a PDF of 5 to 7 pages, no `?` references (grep the second pdflatex run's output for `undefined` before deleting the log if anything looks off).

- [ ] **Step 4: Read the PDF** (Read tool, pages 1-7) and check: every figure present, every quoted number matches the summary file or the study output, no comment stubs left, `\author{Michael Casey}` only.

- [ ] **Step 5: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/doc/orbit_network_note.tex orbit_transfer/orbit_network/doc/figures
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: the note (doc/orbit_network_note.tex, ~6 pp) -- Braik-Ross framing against a mass-varying min-time OCP, the direct-never-loses lemma, the sheet and node tables, the graph and its thrust-floor attributes, the accessibility plane and its three regimes, relay closure L1 halo -> L2 halo -> tulip, the stability-versus-cost verdict, and what the legs campaign would add

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

## Self-review (done while writing)

- Spec coverage: section 7.0-7.7 -> Tasks 5, 6, 7; 7.8 lib -> Tasks 1-4 (plus `edge_tables`, `relay_closure`, `symmetry_residual` as refinements); 7.9 tests -> every task extends `test_lib` or `test_pipeline`, the runner already lists both; 7.10 note -> Task 9; 7.11 out of scope respected (no solves, no API, mirror never in headline numbers).
- Type consistency: `E` fields (`closed_frac`, `tf_min_days`, `tf_med_days`, `n_certified`, `F_floor_N`, `F_floor_str`) are the same in Tasks 2, 4, 5, 6, 7; `P` fields (`pairs`, `cells`, `closed`, `largest`, `multileg`, `caps_days`) in Tasks 3, 5, 6, 7; `R`/`best` fields in Tasks 4, 6, 7; `relay_closure(E, N, tStar_s, src, via, dst)` signature identical in Tasks 4 and 6.
- Review Focus items each have an owning task test: 1 -> Task 2 and Task 3; 2 -> Task 3; 3 -> Task 1; 4 -> Task 4; 5 -> Task 6; 6 -> Task 5 (findnode in G2) and Task 7 (visual check).
- Known hazard written into the plan: the name `lines` is both the script's summary cellstr and a MATLAB colormap function; Task 7 uses `parula` for the relay figure's colors for that reason.
