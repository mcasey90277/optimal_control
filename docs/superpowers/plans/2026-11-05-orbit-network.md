# Orbit Network (topic 2, data-only) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `orbit_network_study.m`, its lib helpers, tests, figures and the LaTeX note: the low-thrust orbital network over the 56 accepted min-time sheets, with the (thrust capability, time cap) coverage and envelope planes, the phase-compatible L1-halo-to-tulip relay closure, and an off-by-default hypothetical reverse graph.

**Architecture:** Every computation with logic of its own lives in `orbit_network/lib/` as a function with an `if nargin == 0` demo and a hand-checked test (`edge_tables`, `network_graph`, `accessibility_plane`, `relay_closure`, `mirror_state`, `symmetry_residual`); the study script `orbit_network_study.m` wires them in numbered sections, prints every gate PASS/FAIL, writes four PNGs and a `results/network_summary.txt` that the note cites. Topic 1 (`stability_cost_study.m`, `load_catalog_sheets` with its `accepted` mask, `orbit_node_table`, `lib/partial_f_test`) is already built; this plan consumes it.

**Tech Stack:** MATLAB R2026a in the shared session (MCP tools), `digraph`/`conncomp`, pumpkyn (`pumpkyn.cr3bp.prop`, `.monodromy`, `.jacobi`), costate_common `get_family_orbit`, pdflatex.

**Spec:** `docs/superpowers/specs/2026-11-05-orbit-network-design.md` (revised 2026-11-05 after the Astra review; sections 7-11). Review record: `orbit_transfer/orbit_network/reviews/`.

## Global Constraints

- Data-only: no new solves, no changes under `DRO_tulip/`, `HALO_tulip/`, `DPO_tulip/`, `HALO_HALO/`, `costate_common/`.
- Every count, minimum, plane and route uses the sheet table's `accepted` mask (`has_solution & conj_pass == 1 & gate_min_lamv > 0 & gate_min_qmt > 0 & gate_dimS == 1`), never `has_solution` alone; solved counts are reported beside accepted counts.
- The hypothetical reverse graph: `hypothetical` defaults to `false`; every printed number from that section is prefixed HYPOTHETICAL; those edges carry `n_accepted = 0`; nothing from that section enters sections 3-5 or the figures.
- Relay-versus-direct COST is never claimed; relay CLOSURE and the phase-compatible assembled time are computed.
- Wording: "lowest accepted rung" (never floor, bracket, threshold); "catalog closure down the ladder" (never survival of physical accessibility); "diagnostics" (never regimes); "coverage" (never robustness).
- MATLAB `digraph` sorts edges: every edge attribute is read from `G.Edges` (or a filtered copy of it), never from a parallel array in sheet order.
- MATLAB house style: header quartet (Purpose / Inputs / Outputs / Revision History with `M. Casey (c) 11/05/2026` and `Copyright Coorbital Inc.`), `if nargin == 0` demo on every function, no `i`/`j`, no `norm`, no `%#ok`, no emojis, ASCII-only READMEs and comments. No `M. Casey / D. Koblick` signature lines in docs. Note uses `\author{Michael Casey}` only.
- Tests are function files returning a logical `ok`, printed `PASS`/`FAIL` per check through the same `check`/`verdict` helpers as `tests/test_lib.m`; run from the shared session by `cd` into `tests/` and calling by name.
- Paths: scripts and tests add `here`, `here/lib`, `../costate_common`; pumpkyn/pumpkynPie come from the session path. No campaign folder on the path.
- Figures: set axes `Color`, `XColor`, `YColor` explicitly (the shared session defaults axes to dark); labels and titles with `'Interpreter', 'none'` where a key contains an underscore; PNG at 150 dpi via `exportgraphics`.
- Outputs go to `results/` (git-ignored). Commit messages end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`; implementers never commit, the controller commits after a clean review; git via `/Library/Developer/CommandLineTools/usr/bin/git`; never stage the unrelated `missiles/` deletions.
- Clean LaTeX aux files after every compile: `rm -f *.aux *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz`.
- Constants: muStar 0.012150585609624, tStar_s 382981.289129055 (read from `T(1).constants`, never retyped); days = t_nd * tStar_s / 86400.

## Review Focus

1. A rung on which a sheet has NO accepted cell: `tf_min_days`/`tf_med_days` NaN, `closed_frac` 0, the plane counts skip it, the graph figure drops that edge. Pinned in Task 2 and Task 4.
2. The fastest SOLVED cell on a rung fails certification: it must not set the edge's minimum. Pinned in Task 2.
3. Edge order: the graph's edge table is sorted by MATLAB; an attribute looked up by (dep, arr) must be that edge's own. Pinned in Task 3 with scrambled input.
4. A relay pair whose smallest sum of leg minima is NOT its smallest phase-compatible assembled time: the assembled minimizer must be chosen, and t_route >= t_sum. Pinned in Task 5.
5. The hypothetical graph path is never exercised by default: the pipeline test runs the script a second time with `study_opts.hypothetical = true` and requires G5 PASS, the HYPOTHETICAL prefix and the reversed-edge count. Pinned in Task 8.
6. A symmetry check at a stored initial state is vacuous (y = vx = vz = 0): G5 uses an off-plane point and an independent periodic-closure check, and a non-periodic state must fail it. Pinned in Task 1 and Task 7.
7. No time cap reaches full direct access at some rung: the script prints "no sampled cap" instead of index 1. Pinned in Task 6.

---

### Task 1: `lib/mirror_state.m` and `lib/symmetry_residual.m`

**Files:**
- Create: `orbit_transfer/orbit_network/lib/mirror_state.m`
- Create: `orbit_transfer/orbit_network/lib/symmetry_residual.m`
- Modify: `orbit_transfer/orbit_network/tests/test_lib.m` (append a section before the final `fprintf('test_lib: ...')`)

**Interfaces:**
- Consumes: `pumpkyn.cr3bp.prop(tauVec, x0row, mu)` -> `[tau, X]` with `X` [N x 6] when `x0row` is [1 x 6]; `get_family_orbit('dro', struct('tau', 1))` -> `[tgrid, rv, info]`, initial state `rv(1, :)`, period `info.periodND`.
- Produces: `R = mirror_state(x)` for `x` [6 x n] -> [6 x n], `(x, -y, z, -vx, vy, -vz)`; refusal `orbit_network:mirror_state:shape`. `[res, tStar] = symmetry_residual(x1, T, mu, opts)` for `x1` [6 x 1]: the minimum over t in [0, T] of `max(abs(phi(x1, t) - mirror_state(x1)))` and the minimizing time; `opts.nGrid` [400].

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

%% symmetry_residual: at an OFF-PLANE point a symmetric orbit maps onto itself; a non-periodic state does not
addpath(fullfile(fileparts(fileparts(here)), 'costate_common'));
[~, ~, rv, info] = evalc('get_family_orbit(''dro'', struct(''tau'', 1))');
x0 = rv(1, :).';  Tp = info.periodND;  mu = 0.012150585609624;
ok = check(ok, max(abs(x0([2 4 6]))) < 1e-12, 'the stored initial state sits on the symmetry plane (so a check there would be vacuous)');
[~, X] = pumpkyn.cr3bp.prop([0 0.37*Tp], x0.', mu);
x1 = X(end, 1:6).';
ok = check(ok, max(abs(x1([2 4 6]))) > 1e-3, 'phi(x0, 0.37 T) is off the plane');
[res, tS] = symmetry_residual(x1, Tp, mu);
ok = check(ok, res < 1e-8, sprintf('DRO tau 1: mirror image of the off-plane point lies on the orbit (residual %.1e at t = %.4f of T = %.4f)', res, tS, Tp));
ok = check(ok, tS >= 0 && tS <= Tp, 'minimizing time inside [0, T]');
xn = x1;  xn(5) = xn(5)*1.02;                              % 2% faster: not on the periodic orbit
resN = symmetry_residual(xn, Tp, mu);
ok = check(ok, resN > 1e-6, sprintf('a non-periodic state is NOT mapped onto an orbit (residual %.1e): the gate can fail', resN));
```

- [ ] **Step 2: Run the test to verify it fails**

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
%   (x, y, z, vx, vy, vz) -> (x, -y, z, -vx, vy, -vz). With time reversed
%   the image of a ballistic trajectory is a ballistic trajectory, and the
%   image of a symmetric periodic orbit is the same orbit. For a THRUSTING
%   arc the thrust acceleration is reflected by diag(1, -1, 1) and the
%   mass runs backward (a mass-gain model), so a mirrored transfer is a
%   geometric object only -- see the spec's "hypothetical reverse graph".
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
   %Demo: a state on the symmetry plane is its own image; an off-plane one is not
     x0 = [0.85; 0; 0; 0; 0.5; 0];  x1 = [0.8; 0.1; 0.02; -0.1; 0.4; 0.01];
     fprintf('on-plane : %s -> %s\n', mat2str(x0.'), mat2str(mirror_state(x0).'));
     fprintf('off-plane: %s -> %s\n', mat2str(x1.'), mat2str(mirror_state(x1).'));
     R = mirror_state(x1);
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
function [res, tStar] = symmetry_residual(x1, T, mu, opts)
%% Purpose:
%
%   How far the mirror image of a state is from the orbit through it:
%   min over t in [0, T] of max|phi(x1, t) - mirror_state(x1)|. Zero (to
%   integrator accuracy) when the periodic orbit through x1 is x-axis
%   symmetric. The check is meaningful only for an OFF-PLANE x1: a state
%   with y = vx = vz = 0 is its own image and the residual is zero at
%   t = 0 whatever the orbit does. A coarse grid over one period locates
%   the minimum, fminbnd refines it on the bracketing interval. This does
%   not test periodicity; the caller checks closure separately.
%
%% Inputs:
%
%  x1                       [6 x 1]                 An off-plane state on
%                                                   the orbit
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
   %Demo: pumpkyn's own L1 halo, from a point 0.37 T along it
     x0 = [9.2294247990951594E-1; 0; 2.1592634764719765E-1; 4.2453615136710976E-13; 1.2786666755130519E-1; -2.3343675990673871E-12];
     T0 = 1.8036821222727220;  mu0 = 1.215058560962404E-2;
     [~, X] = pumpkyn.cr3bp.prop([0 0.37*T0], x0.', mu0);
     [res, tStar] = symmetry_residual(X(end, 1:6).', T0, mu0);
     fprintf('residual %.2e at t = %.4f\n', res, tStar);
     return
end

if nargin < 4 || isempty(opts), opts = struct(); end
if isfield(opts, 'nGrid') && ~isempty(opts.nGrid), nGrid = opts.nGrid; else, nGrid = 400; end
x1  = x1(:);
Rx  = mirror_state(x1);
tg  = linspace(0, T, nGrid);
[~, X] = pumpkyn.cr3bp.prop(tg, x1.', mu);                 % [nGrid x 6]
d   = max(abs(X(:, 1:6).' - Rx), [], 1);                  % [1 x nGrid]
[~, k] = min(d);
lo  = tg(max(k - 1, 1));  hi = tg(min(k + 1, nGrid));
f   = @(t) resid_at(t, x1, Rx, mu);
[tStar, res] = fminbnd(f, lo, hi, optimset('TolX', 1e-12));
if d(k) < res, res = d(k); tStar = tg(k); end               % never worse than the grid
end

% ------------------------------------------------------------------------
function r = resid_at(t, x1, Rx, mu)
%% Purpose:
%
%   The residual at one time (t = 0 needs no propagation).
%
if t <= 0
    r = max(abs(x1 - Rx));
else
    [~, X] = pumpkyn.cr3bp.prop([0 t], x1.', mu);
    r = max(abs(X(end, 1:6).' - Rx));
end
end
```

- [ ] **Step 5: Run the test to verify it passes** (`ok = test_lib();` every line PASS). If the DRO residual is above 1e-8, raise `nGrid` to 1000 in the test's call before loosening anything else.
- [ ] **Step 6: Code Analyzer on both new files** (MCP `check_matlab_code`): no warnings.
- [ ] **Step 7: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/lib/mirror_state.m orbit_transfer/orbit_network/lib/symmetry_residual.m orbit_transfer/orbit_network/tests/test_lib.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: mirror_state (x-axis symmetry) and symmetry_residual (mirror image of an OFF-PLANE orbit point lies on the orbit: grid + fminbnd over one period); tests incl. the vacuous-on-plane fact and a non-periodic state that fails

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 2: `lib/edge_tables.m`

**Files:**
- Create: `orbit_transfer/orbit_network/lib/edge_tables.m`
- Modify: `orbit_transfer/orbit_network/tests/test_lib.m` (append before the final fprintf)

**Interfaces:**
- Consumes: the sheet table `T` from `load_catalog_sheets([])` (fields `dep_key`, `arr_key`, `dep_family`, `arr_family`, `rungs_N` [1 x nR], `has_solution`/`accepted`/`tf_nd` [nD x nA x nR], `constants.tStar_s`).
- Produces: `E = edge_tables(T)`, struct [numel(T) x 1] with `sheet_index` (position in `T`), `dep_key`, `arr_key`, `dep_family`, `arr_family`, `rungs_N` [1 x nR], `n_solved` [1 x nR] (has_solution count), `n_accepted` [1 x nR], `closed_frac` [1 x nR] (= n_accepted / (nD*nA)), `tf_min_days`, `tf_med_days` [1 x nR] over ACCEPTED cells (NaN when none), `F_low_N` (the smallest rung with n_accepted > 0; NaN if none), `F_low_str` (`'0.5 (lowest tested)'` when it is the last rung, else `sprintf('%g', F_low_N)`, `'none'` if none).

- [ ] **Step 1: Append the failing test** (before the final fprintf of `test_lib.m`):

```matlab
%% edge_tables: synthetic 2-sheet table; the fastest SOLVED cell is refuted and must not be selected
tStar = 382981.289129055;  d = tStar/86400;
S = struct('dep_key', 'a', 'arr_key', 'b', 'dep_family', 'dro', 'arr_family', 'tulip', 'rungs_N', [2 1 0.5], ...
           'has_solution', false(2, 2, 3), 'accepted', false(2, 2, 3), 'tf_nd', nan(2, 2, 3), ...
           'sD_frac', [0 0.5], 'sA_frac', [0 0.5], ...
           'constants', struct('muStar', 0.012150585609624, 'lStar_km', 389703.264829278, 'tStar_s', tStar));
S.has_solution(:, :, 1) = [true true; true false];  S.accepted(:, :, 1) = [false true; true false];  S.tf_nd(:, :, 1) = [1 2; 3 NaN];
S.has_solution(:, :, 2) = [true false; false false]; S.accepted(:, :, 2) = [true false; false false]; S.tf_nd(:, :, 2) = [4 NaN; NaN NaN];
S2 = S;  S2.arr_key = 'c';  S2.has_solution(:) = false;  S2.accepted(:) = false;  S2.has_solution(1, 1, 3) = true;  S2.accepted(1, 1, 3) = true;  S2.tf_nd(1, 1, 3) = 5;
E = edge_tables([S; S2]);
ok = check(ok, numel(E) == 2 && E(2).sheet_index == 2 && strcmp(E(2).arr_key, 'c'), 'edge_tables: one row per sheet, sheet_index and keys carried');
ok = check(ok, isequal(E(1).n_solved, [3 1 0]) && isequal(E(1).n_accepted, [2 1 0]) && isequal(E(1).closed_frac, [2/4 1/4 0]), 'edge_tables: solved and accepted counts, closed_frac over accepted');
ok = check(ok, abs(E(1).tf_min_days(1) - 2*d) < 1e-12, 'edge_tables: the refuted 1-day cell is NOT the minimum; min over ACCEPTED is 2 days');
ok = check(ok, abs(E(1).tf_med_days(1) - 2.5*d) < 1e-12, 'edge_tables: median over accepted cells (2, 3 days -> 2.5)');
ok = check(ok, isnan(E(1).tf_min_days(3)) && isnan(E(1).tf_med_days(3)), 'edge_tables: NaN on a rung with nothing accepted');
ok = check(ok, E(1).F_low_N == 1 && strcmp(E(1).F_low_str, '1'), 'edge_tables: lowest accepted rung 1 N');
ok = check(ok, E(2).F_low_N == 0.5 && strcmp(E(2).F_low_str, '0.5 (lowest tested)'), 'edge_tables: a sheet accepted at the last rung says so, not a bracket');
```

- [ ] **Step 2: Run to verify it fails** (`Unrecognized function or variable 'edge_tables'`).

- [ ] **Step 3: Write `lib/edge_tables.m`**

```matlab
function E = edge_tables(T)
%% Purpose:
%
%   Per-edge, per-rung attributes of the network built from a sheet
%   table, over ACCEPTED cells (certified witnesses): solved and accepted
%   counts, the accepted fraction of the phasing grid, the min and median
%   time of flight over accepted cells (NaN when none), and the lowest
%   accepted rung. An unaccepted rung is unresolved, not infeasible, so
%   the lowest accepted rung is reported as what it is and never as a
%   floor or a bracket.
%
%% Inputs:
%
%  T                        struct [n x 1]          Sheet table
%                                                   (load_catalog_sheets)
%
%% Outputs:
%
%  E                        struct [n x 1]          sheet_index, dep_key,
%                                                   arr_key, dep_family,
%                                                   arr_family, rungs_N
%                                                   [1 x nR], n_solved,
%                                                   n_accepted,
%                                                   closed_frac,
%                                                   tf_min_days,
%                                                   tf_med_days [1 x nR],
%                                                   F_low_N, F_low_str
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
     fprintf('%-18s %-14s %7s %8s %7s %8s %20s\n', 'dep', 'arr', 'frac1N', 'tf1N[d]', 'frac.5', 'tf.5[d]', 'lowest accepted');
     for kk = 1:numel(E)
         fprintf('%-18s %-14s %7.3f %8.3f %7.3f %8.3f %20s\n', E(kk).dep_key, E(kk).arr_key, E(kk).closed_frac(r1), E(kk).tf_min_days(r1), E(kk).closed_frac(r05), E(kk).tf_min_days(r05), E(kk).F_low_str);
     end
     return
end

n    = numel(T);
rows = cell(n, 1);
for kk = 1:n
    s     = T(kk);
    nR    = numel(s.rungs_N);
    nCell = size(s.accepted, 1)*size(s.accepted, 2);
    toDay = s.constants.tStar_s/86400;
    e.sheet_index = kk;
    e.dep_key     = s.dep_key;
    e.arr_key     = s.arr_key;
    e.dep_family  = s.dep_family;
    e.arr_family  = s.arr_family;
    e.rungs_N     = s.rungs_N;
    e.n_solved    = zeros(1, nR);
    e.n_accepted  = zeros(1, nR);
    e.closed_frac = zeros(1, nR);
    e.tf_min_days = nan(1, nR);
    e.tf_med_days = nan(1, nR);
    for rr = 1:nR
        a  = s.accepted(:, :, rr);
        tf = s.tf_nd(:, :, rr);
        e.n_solved(rr)    = nnz(s.has_solution(:, :, rr));
        e.n_accepted(rr)  = nnz(a);
        e.closed_frac(rr) = nnz(a)/nCell;
        if any(a(:))
            e.tf_min_days(rr) = min(tf(a))*toDay;
            e.tf_med_days(rr) = median(tf(a))*toDay;
        end
    end
    low = find(e.n_accepted > 0, 1, 'last');              % rungs run shallow -> deep
    if isempty(low)
        e.F_low_N = NaN;  e.F_low_str = 'none';
    elseif low == nR
        e.F_low_N = s.rungs_N(nR);  e.F_low_str = sprintf('%g (lowest tested)', s.rungs_N(nR));
    else
        e.F_low_N = s.rungs_N(low);  e.F_low_str = sprintf('%g', s.rungs_N(low));
    end
    rows{kk} = e;
end
E = vertcat(rows{:});
end
```

- [ ] **Step 4: Run to verify it passes.**
- [ ] **Step 5: Code Analyzer clean.**
- [ ] **Step 6: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/lib/edge_tables.m orbit_transfer/orbit_network/tests/test_lib.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: edge_tables -- per-edge per-rung solved/accepted counts, min/median t_f over ACCEPTED cells (a refuted fastest cell is not selected), lowest accepted rung stated as such; synthetic test

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 3: `lib/network_graph.m`

**Files:**
- Create: `orbit_transfer/orbit_network/lib/network_graph.m`
- Modify: `orbit_transfer/orbit_network/tests/test_lib.m` (append before the final fprintf)

**Interfaces:**
- Consumes: `E` from Task 2; node table `N` (fields `key`, `family`, `period_nd`, `CJ`, `sigma`, `stable`).
- Produces: `G = network_graph(E, N)`, a `digraph` whose `G.Nodes` has variables `Name` (the keys), `family`, `period_nd`, `CJ`, `sigma`, `stable`, and whose `G.Edges` has `EndNodes` (names), `sheet_index`, `closed_frac`, `tf_min_days`, `tf_med_days`, `n_accepted`, `n_solved` (each [nEdges x nR]), `F_low_N`, `hypothetical` (false). Refusal `orbit_network:network_graph:node` if an edge key is not in `N`.

- [ ] **Step 1: Append the failing test** (uses `E` from the Task 2 block and `d`):

```matlab
%% network_graph: scrambled input order; attributes looked up by (dep, arr) are the right ones
Nt = struct('key', {'b'; 'c'; 'a'}, 'family', {'tulip'; 'tulip'; 'dro'}, 'period_nd', {1; 2; 3}, 'CJ', {3; 3; 3}, 'sigma', {0; 0; 0}, 'stable', {true; true; true});
G  = network_graph([E(2); E(1)], Nt);                       % scrambled: sheet 2 first
ok = check(ok, numnodes(G) == 3 && numedges(G) == 2 && isequal(sort(G.Nodes.Name), {'a'; 'b'; 'c'}), 'network_graph: 3 named nodes, 2 edges');
ia = findedge(G, 'a', 'b');  ic = findedge(G, 'a', 'c');
ok = check(ok, G.Edges.sheet_index(ia) == 1 && G.Edges.sheet_index(ic) == 2, 'network_graph: sheet_index follows the edge, not the input order');
ok = check(ok, abs(G.Edges.tf_min_days(ia, 1) - 2*d) < 1e-12 && isnan(G.Edges.tf_min_days(ic, 1)) && abs(G.Edges.tf_min_days(ic, 3) - 5*d) < 1e-12, 'network_graph: per-rung t_f rows belong to their edges');
ok = check(ok, G.Edges.F_low_N(ia) == 1 && G.Edges.F_low_N(ic) == 0.5 && ~any(G.Edges.hypothetical), 'network_graph: scalar attributes follow too; hypothetical false');
ok = check(ok, strcmp(G.Nodes.family{findnode(G, 'a')}, 'dro') && G.Nodes.period_nd(findnode(G, 'c')) == 2, 'network_graph: node variables by name');
try
    network_graph(E(1), Nt(1:2));
    ok = check(ok, false, 'network_graph :node refusal');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:network_graph:node'), 'network_graph :node refusal when an edge key is missing');
end
```

- [ ] **Step 2: Run to verify it fails.**

- [ ] **Step 3: Write `lib/network_graph.m`**

```matlab
function G = network_graph(E, N)
%% Purpose:
%
%   The directed graph of the catalogs: nodes are orbits (named by key),
%   edges are sheets. MATLAB sorts a digraph's edges by (source, target),
%   so every attribute lives IN the edge table and is read back from
%   G.Edges; sheet_index maps an edge to its row of the sheet table.
%
%% Inputs:
%
%  E                        struct [m x 1]          edge_tables output
%
%  N                        struct [n x 1]          Node table (key,
%                                                   family, period_nd, CJ,
%                                                   sigma, stable)
%
%% Outputs:
%
%  G                        digraph                 Nodes: Name, family,
%                                                   period_nd, CJ, sigma,
%                                                   stable. Edges:
%                                                   EndNodes (names),
%                                                   sheet_index,
%                                                   closed_frac,
%                                                   tf_min_days,
%                                                   tf_med_days,
%                                                   n_accepted, n_solved
%                                                   [m x nR], F_low_N,
%                                                   hypothetical (false)
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the catalogs' graph, degrees per node
     here = fileparts(fileparts(mfilename('fullpath')));
     addpath(here, fullfile(fileparts(here), 'costate_common'));
     T = load_catalog_sheets([]);
     G = network_graph(edge_tables(T), orbit_node_table(T, struct('verbose', false)));
     disp([G.Nodes(:, 'Name'), table(indegree(G), outdegree(G), 'VariableNames', {'in', 'out'})]);
     return
end

keys = {N.key};
dep  = {E.dep_key};  arr = {E.arr_key};
miss = [dep(~ismember(dep, keys)), arr(~ismember(arr, keys))];
if ~isempty(miss)
    error('orbit_network:network_graph:node', 'edge key(s) not in the node table: %s', strjoin(unique(miss), ', '));
end
nodeTab = table(keys(:), {N.family}.', [N.period_nd].', [N.CJ].', [N.sigma].', logical([N.stable]).', ...
                'VariableNames', {'Name', 'family', 'period_nd', 'CJ', 'sigma', 'stable'});
edgeTab = table([dep(:), arr(:)], [E.sheet_index].', vertcat(E.closed_frac), vertcat(E.tf_min_days), vertcat(E.tf_med_days), ...
                vertcat(E.n_accepted), vertcat(E.n_solved), [E.F_low_N].', false(numel(E), 1), ...
                'VariableNames', {'EndNodes', 'sheet_index', 'closed_frac', 'tf_min_days', 'tf_med_days', 'n_accepted', 'n_solved', 'F_low_N', 'hypothetical'});
G = digraph(edgeTab, nodeTab);
end
```

- [ ] **Step 4: Run to verify it passes.**
- [ ] **Step 5: Code Analyzer clean.**
- [ ] **Step 6: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/lib/network_graph.m orbit_transfer/orbit_network/tests/test_lib.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: network_graph -- the digraph with every attribute in the edge table (MATLAB sorts edges) and a sheet_index back-reference; scrambled-order test

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 4: `lib/accessibility_plane.m`

**Files:**
- Create: `orbit_transfer/orbit_network/lib/accessibility_plane.m`
- Modify: `orbit_transfer/orbit_network/tests/test_lib.m` (append before the final fprintf)

**Interfaces:**
- Consumes: sheet table `T` (`accepted`, `tf_nd`, `rungs_N`, `constants.tStar_s`).
- Produces: `P = accessibility_plane(T, caps_days)` with `P.rungs_N` [1 x nR], `P.caps_days` [1 x nC], raw coverage `P.pairs`, `P.cells` [nR x nC] (accepted cells AT that rung with t_f <= cap; pairs = sheets with at least one), `P.closed` {nR x nC} (sheet indices in `pairs`), and the known-feasible envelope `P.pairs_env`, `P.cells_env` [nR x nC] (a trajectory accepted at rung s with F_s <= F_r is feasible under the bound F_r, so row r unions rows s >= r: per cell, the min of t_f over accepted rungs s >= r, counted if <= cap).

- [ ] **Step 1: Append the failing test** (reuses `S`, `S2`, `d`):

```matlab
%% accessibility_plane: raw coverage vs the known-feasible envelope; monotone in the cap; zeros below the data
P = accessibility_plane([S; S2], d*[0.5 2 2.5 10]);
ok = check(ok, isequal(size(P.pairs), [3 4]) && isequal(size(P.pairs_env), [3 4]) && isequal(size(P.closed), [3 4]), 'accessibility_plane: [nR x nC] arrays');
ok = check(ok, isequal(P.pairs(1, :), [0 1 1 1]) && isequal(P.cells(1, :), [0 1 1 2]), 'accessibility_plane: raw at rung 1 counts ACCEPTED cells (2 and 3 days; the refuted 1-day cell is not there)');
ok = check(ok, isequal(P.pairs(3, :), [0 0 0 1]) && isequal(P.cells(3, :), [0 0 0 1]), 'accessibility_plane: raw at rung 3 only sheet 2 (5 days)');
ok = check(ok, isequal(P.closed{1, 4}, 1) && isequal(P.closed{3, 4}, 2) && isempty(P.closed{1, 1}), 'accessibility_plane: closed sheet lists');
ok = check(ok, isequal(P.pairs_env(1, :), [0 1 1 2]), 'accessibility_plane: the envelope at the 2 N bound includes sheet 2''s 0.5 N witness (feasible under a higher bound)');
ok = check(ok, isequal(P.cells_env(1, :), [0 1 1 4]) && isequal(P.cells_env(3, :), P.cells(3, :)), 'accessibility_plane: envelope cells take the min over deeper rungs (2, 3, 4 d of sheet 1 and 5 d of sheet 2 under the 10-day cap); the deepest row equals the raw row');
ok = check(ok, all(all(diff(P.pairs, 1, 2) >= 0)) && all(all(diff(P.cells_env, 1, 2) >= 0)) && all(all(diff(P.pairs_env, 1, 1) <= 0)), 'accessibility_plane: non-decreasing in the cap; envelope non-increasing toward deeper rungs');
P0 = accessibility_plane([S; S2], d*[0.1 0.2]);
ok = check(ok, ~any(P0.pairs(:)) && ~any(P0.cells_env(:)), 'accessibility_plane: caps below the data give all zeros, no error');
```

- [ ] **Step 2: Run to verify it fails.**

- [ ] **Step 3: Write `lib/accessibility_plane.m`**

```matlab
function P = accessibility_plane(T, caps_days)
%% Purpose:
%
%   The (thrust bound, time cap) plane of a sheet table, two ways. RAW
%   COVERAGE: at each rung and cap, how many pairs have an ACCEPTED cell
%   solved at that rung with time of flight under the cap, and how many
%   cells -- what the catalogs hold at that thrust. KNOWN-FEASIBLE
%   ENVELOPE: a trajectory accepted at a lower thrust bound is feasible
%   under a higher one (throttleable thrust), so the envelope at a rung
%   unions every rung at or below its thrust. Thrust here is a capability
%   axis for one spacecraft, not a consumable budget. Rungs are the
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
%                                                   .pairs, .cells [nR x
%                                                   nC] raw; .closed {nR
%                                                   x nC} sheet indices
%                                                   in .pairs; .pairs_env,
%                                                   .cells_env [nR x nC]
%                                                   envelope
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
     fprintf('raw pairs (rows = rungs %s N, cols = caps %s d):\n', mat2str(P.rungs_N), mat2str(P.caps_days, 3));
     disp(P.pairs);
     fprintf('envelope pairs:\n');  disp(P.pairs_env);
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
P.closed    = repmat({zeros(1, 0)}, nR, nC);
P.pairs_env = zeros(nR, nC);
P.cells_env = zeros(nR, nC);
for kk = 1:numel(T)
    toDay = T(kk).constants.tStar_s/86400;
    tfA   = T(kk).tf_nd*toDay;
    tfA(~T(kk).accepted) = Inf;                            % unaccepted cells never count
    for rr = 1:nR
        tc = tfA(:, :, rr);  tc = tc(isfinite(tc));        % raw: this rung's accepted cells
        te = min(tfA(:, :, rr:end), [], 3);                % envelope: best over this and every deeper rung, per cell
        te = te(isfinite(te));
        for cc = 1:nC
            nRaw = nnz(tc <= caps(cc));
            nEnv = nnz(te <= caps(cc));
            P.cells(rr, cc)     = P.cells(rr, cc) + nRaw;
            P.cells_env(rr, cc) = P.cells_env(rr, cc) + nEnv;
            if nRaw > 0
                P.pairs(rr, cc)  = P.pairs(rr, cc) + 1;
                P.closed{rr, cc} = [P.closed{rr, cc}, kk];
            end
            if nEnv > 0
                P.pairs_env(rr, cc) = P.pairs_env(rr, cc) + 1;
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
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: accessibility_plane -- raw catalog coverage and the known-feasible envelope over (thrust bound, time cap), accepted cells only; synthetic test incl. a deeper-rung witness under a shallower bound and caps below the data

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 5: `lib/relay_closure.m` (phase-compatible assembled time)

**Files:**
- Create: `orbit_transfer/orbit_network/lib/relay_closure.m`
- Modify: `orbit_transfer/orbit_network/tests/test_lib.m` (append before the final fprintf)

**Interfaces:**
- Consumes: sheet table `T` (`dep_key`, `arr_key`, `accepted`, `tf_nd`, `sD_frac`, `sA_frac`, `rungs_N`, `constants.tStar_s`), node table `N` (`key`, `period_nd`).
- Produces: `R = relay_closure(T, N, src_prefix, relay_prefix, dst_prefix)` -> struct [nRoutes x 1] with `src_key`, `relay_key`, `dst_key`, `s1`, `s2` (indices into `T`), `closed` [1 x nR] (both legs have an accepted cell), `t_sum_days` [1 x nR] (min t1 + min t2, the optimistic lower bound; NaN unless closed), `t_route_days` [1 x nR] (min over accepted cell pairs of t1(a,b) + T_via * mod(sD2(c) - sA1(b), 1) + t2(c,d); NaN unless closed), `wait_days` [1 x nR] (the wait inside that minimizer), `cells` [4 x nR] ([a; b; c; d] of the minimizer, 0 when open). Routes are every (src, relay, dst) with keys matching the prefixes (`strncmp`) for which both sheets exist. Refusals: `orbit_network:relay_closure:node` (relay key not in `N`), `orbit_network:relay_closure:grid` (leg 1's `sA_frac` and leg 2's `sD_frac` are not the same set to 1e-12: the phase origins would not be comparable).

- [ ] **Step 1: Append the failing test**:

```matlab
%% relay_closure: a toy a -> b -> c where the smallest SUM of leg minima is not the smallest ASSEMBLED time
mk = @(dk, ak, A, tf) struct('dep_key', dk, 'arr_key', ak, 'dep_family', 'x', 'arr_family', 'x', 'rungs_N', [2 1], ...
                           'has_solution', A, 'accepted', A, 'tf_nd', tf, 'sD_frac', [0 0.5], 'sA_frac', [0 0.5], ...
                           'constants', struct('muStar', 0.012150585609624, 'lStar_km', 389703.264829278, 'tStar_s', tStar));
% leg 1 (a -> b): at rung 1 two accepted cells: arrive at phase 0 in 3 days, or at phase 0.5 in 3.2 days; rung 2 open
A1 = false(2, 2, 2);  A1(1, 1, 1) = true;  A1(1, 2, 1) = true;
t1 = nan(2, 2, 2);    t1(1, 1, 1) = 3/d;   t1(1, 2, 1) = 3.2/d;
% leg 2 (b -> c): at rung 1 one accepted cell departing at phase 0.5 in 2 days; at rung 2 one cell departing at phase 0 in 2.5 days
A2 = false(2, 2, 2);  A2(2, 1, 1) = true;  A2(1, 1, 2) = true;
t2 = nan(2, 2, 2);    t2(2, 1, 1) = 2/d;   t2(1, 1, 2) = 2.5/d;
Tt = [mk('L1_a', 'L2_b', A1, t1); mk('L2_b', 'tulip_c', A2, t2); mk('L1_a', 'tulip_c', A1, t1)];   % the third is a direct edge, not a route
Nt = struct('key', {'L1_a'; 'L2_b'; 'tulip_c'}, 'period_nd', {1; 4/d; 3});                        % relay period = 4 days
R  = relay_closure(Tt, Nt, 'L1', 'L2', 'tulip');
ok = check(ok, numel(R) == 1 && strcmp(R.relay_key, 'L2_b') && R.s1 == 1 && R.s2 == 2, 'relay_closure: one route a -> b -> c, sheets indexed into T');
ok = check(ok, isequal(R.closed, [true false]), 'relay_closure: closed iff BOTH legs have an accepted cell (rung 2: leg 1 open)');
ok = check(ok, abs(R.t_sum_days(1) - 5) < 1e-9, 'relay_closure: optimistic sum = 3 + 2 days');
% assembled: arrive phase 0 (3 d) then wait to phase 0.5 = 2 d, +2 d = 7 d; or arrive phase 0.5 (3.2 d), no wait, +2 d = 5.2 d -> choose 5.2
ok = check(ok, abs(R.t_route_days(1) - 5.2) < 1e-9 && abs(R.wait_days(1)) < 1e-9 && isequal(R.cells(:, 1), [1; 2; 2; 1]), 'relay_closure: the assembled minimizer is the phase-compatible 5.2-day pairing, not the 5-day sum');
ok = check(ok, R.t_route_days(1) >= R.t_sum_days(1) && isnan(R.t_route_days(2)) && isnan(R.t_sum_days(2)) && all(R.cells(:, 2) == 0), 'relay_closure: t_route >= t_sum; NaN and zero cells when open');
Tbad = Tt;  Tbad(2).sD_frac = [0.1 0.6];
try
    relay_closure(Tbad, Nt, 'L1', 'L2', 'tulip');
    ok = check(ok, false, 'relay_closure :grid refusal');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:relay_closure:grid'), 'relay_closure :grid refusal when the two legs'' phase grids on the relay differ');
end
```

- [ ] **Step 2: Run to verify it fails.**

- [ ] **Step 3: Write `lib/relay_closure.m`**

```matlab
function R = relay_closure(T, N, src_prefix, relay_prefix, dst_prefix)
%% Purpose:
%
%   Two-leg routes src -> relay -> dst assembled from catalog sheets
%   without new solves. A route is closed at a rung iff both legs have an
%   accepted cell. Its time is PHASE-COMPATIBLE: leg 1 arrives on the
%   relay at phase sA1(b), leg 2 departs at phase sD2(c), and the coast
%   between them on the relay orbit takes T_relay * mod(sD2(c) - sA1(b),
%   1); the route time is the minimum over accepted cell pairs of
%   t1(a, b) + wait + t2(c, d). The sum of the two legs' independent
%   minima is kept as the optimistic lower bound. Phases of both legs
%   must be measured from the same relay-orbit origin (the same
%   get_family_orbit state), which the grid check enforces. Mass: a leg
%   solved at m0 is flyable at the previous leg's arrival mass by scaling
%   thrust and mass by q = m1/m0 <= 1 (same path, bound respected); the
%   note states this. Relay-versus-direct cost is not a question for
%   this function.
%
%% Inputs:
%
%  T                        struct [n x 1]          Sheet table
%                                                   (load_catalog_sheets)
%
%  N                        struct [m x 1]          Node table (key,
%                                                   period_nd)
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
%                                                   dst_key, s1, s2,
%                                                   closed [1 x nR],
%                                                   t_sum_days,
%                                                   t_route_days,
%                                                   wait_days [1 x nR],
%                                                   cells [4 x nR]
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
     R = relay_closure(T, orbit_node_table(T, struct('verbose', false)), 'halo_L1', 'halo_L2', 'tulip');
     r1 = find(T(1).rungs_N == 1);
     for kk = 1:numel(R)
         fprintf('%-18s -> %-16s -> %-10s  closed %d  t_sum %.3f d  t_route %.3f d (wait %.3f d)\n', R(kk).src_key, R(kk).relay_key, R(kk).dst_key, R(kk).closed(r1), R(kk).t_sum_days(r1), R(kk).t_route_days(r1), R(kk).wait_days(r1));
     end
     return
end

dep  = {T.dep_key};  arr = {T.arr_key};
keys = {N.key};
isS  = @(k, p) strncmp(k, p, numel(p));
pairs = zeros(0, 2);
for s1 = 1:numel(T)
    if ~(isS(dep{s1}, src_prefix) && isS(arr{s1}, relay_prefix)), continue; end
    for s2 = 1:numel(T)
        if strcmp(dep{s2}, arr{s1}) && isS(arr{s2}, dst_prefix)
            pairs(end+1, :) = [s1 s2];
        end
    end
end
nRoute = size(pairs, 1);
rows   = cell(nRoute, 1);
for qq = 1:nRoute
    s1 = pairs(qq, 1);  s2 = pairs(qq, 2);
    L1 = T(s1);  L2 = T(s2);
    nR = numel(L1.rungs_N);
    kN = find(strcmp(keys, L1.arr_key), 1);
    if isempty(kN)
        error('orbit_network:relay_closure:node', 'relay %s is not in the node table', L1.arr_key);
    end
    if numel(L1.sA_frac) ~= numel(L2.sD_frac) || max(abs(sort(L1.sA_frac(:)) - sort(L2.sD_frac(:)))) > 1e-12
        error('orbit_network:relay_closure:grid', 'route %s -> %s -> %s: leg 1 arrival phases %s and leg 2 departure phases %s differ', ...
              L1.dep_key, L1.arr_key, L2.arr_key, mat2str(L1.sA_frac, 4), mat2str(L2.sD_frac, 4));
    end
    toDay = L1.constants.tStar_s/86400;
    Trel  = N(kN).period_nd*toDay;
    W     = Trel*mod(L2.sD_frac(:).' - L1.sA_frac(:), 1);  % [nA1 x nD2] wait from arrival phase b to departure phase c
    r.src_key   = L1.dep_key;
    r.relay_key = L1.arr_key;
    r.dst_key   = L2.arr_key;
    r.s1        = s1;
    r.s2        = s2;
    r.closed       = false(1, nR);
    r.t_sum_days   = nan(1, nR);
    r.t_route_days = nan(1, nR);
    r.wait_days    = nan(1, nR);
    r.cells        = zeros(4, nR);
    for rr = 1:nR
        a1 = L1.accepted(:, :, rr);  t1 = L1.tf_nd(:, :, rr)*toDay;
        a2 = L2.accepted(:, :, rr);  t2 = L2.tf_nd(:, :, rr)*toDay;
        if ~any(a1(:)) || ~any(a2(:)), continue; end
        r.closed(rr)     = true;
        r.t_sum_days(rr) = min(t1(a1)) + min(t2(a2));
        best = Inf;
        [aa, bb] = find(a1);  [cc, dd] = find(a2);
        for p1 = 1:numel(aa)
            for p2 = 1:numel(cc)
                tot = t1(aa(p1), bb(p1)) + W(bb(p1), cc(p2)) + t2(cc(p2), dd(p2));
                if tot < best
                    best = tot;
                    r.t_route_days(rr) = tot;
                    r.wait_days(rr)    = W(bb(p1), cc(p2));
                    r.cells(:, rr)     = [aa(p1); bb(p1); cc(p2); dd(p2)];
                end
            end
        end
    end
    rows{qq} = r;
end
if nRoute == 0
    R = struct('src_key', {}, 'relay_key', {}, 'dst_key', {}, 's1', {}, 's2', {}, 'closed', {}, 't_sum_days', {}, 't_route_days', {}, 'wait_days', {}, 'cells', {});
else
    R = vertcat(rows{:});
end
end
```

(`pairs(end+1, :)` grows in a loop; the Code Analyzer notes it at info level; do NOT add `%#ok`.)

- [ ] **Step 4: Run to verify it passes.**
- [ ] **Step 5: Code Analyzer: no warnings.**
- [ ] **Step 6: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/lib/relay_closure.m orbit_transfer/orbit_network/tests/test_lib.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: relay_closure -- two-leg routes with the PHASE-COMPATIBLE assembled time (leg 1 arrival phase, coast on the relay, leg 2 departure phase) over accepted cell pairs, the optimistic sum kept as a labelled lower bound, grid-origin refusal; toy test where the smallest sum is not the smallest assembled time

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 6: `orbit_network_study.m`, sections 0-4 with G1-G3

**Files:**
- Create: `orbit_transfer/orbit_network/orbit_network_study.m`
- Modify: `orbit_transfer/orbit_network/tests/test_pipeline.m` (the `scripts` cell gets a second row)

**Interfaces:**
- Consumes: `load_catalog_sheets([])`, `orbit_node_table(T, struct('verbose', false))`, `edge_tables(T)`, `network_graph(E, N)`, `accessibility_plane(T, caps)`.
- Produces (base-workspace variables later sections use): `T`, `N`, `E`, `G`, `rungs`, `nR`, `rC` (reference rung index), `rD` (deepest), `tStar`, `caps`, `nCaps`, `P`, `relay` (prefix struct), `hypothetical`, `gates`, `lines`, `outDir`, `visible`, `nS`, `nN`, `keyN`. The script honours `study_opts` (`.outDir`, `.visible`, `.hypothetical`) from the base workspace exactly as `stability_cost_study.m` does.

- [ ] **Step 1: Extend the pipeline test.** In `tests/test_pipeline.m` replace the single-row `scripts = {...};` with:

```matlab
scripts = {'stability_cost_study', {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7'}, {'H1', 'H2'}, {'stability_cost.png'}
           'orbit_network_study',  {'G1', 'G2', 'G3'},                         {},           {}};
```

- [ ] **Step 2: Run `test_pipeline` to verify it fails** (the network script does not exist).

- [ ] **Step 3: Write `orbit_network_study.m` with sections 0-4 and the summary block.** Leave the lines `%% 5  relay closure`, `%% 6  hypothetical reverse graph`, `%% 7  figures` as markers for Tasks 7 and 8.

```matlab
%% ORBIT_NETWORK_STUDY  The low-thrust orbital network over the costate catalogs
%
%   Braik & Ross (arXiv 2605.31543, 2026) organize periodic orbits into a
%   weighted graph whose edges come from reachable-set overlap under
%   heading-change impulses, then read hubs, relays and budget regimes off
%   it. Here the nodes are the 18 pumpkyn orbits of our four min-time
%   catalogs and the edges are the 56 sheets, every number over ACCEPTED
%   cells (certified witnesses). No proxies, no new solves. The edge
%   attributes are catalog coverage across thrust: the lowest accepted
%   rung and a phase-resolved time of flight; their (dV cap, time cap)
%   budget plane becomes our (thrust capability, time cap) plane, shown
%   as raw coverage and as the known-feasible envelope.
%
%   What this data cannot say, said once: an unaccepted cell or rung is
%   unresolved, not infeasible; no pair has both a direct edge and a
%   two-leg route, so relay-versus-direct cost is not tested (the lemma in
%   the note: the unrestricted infimum is never above the relay infimum);
%   relay CLOSURE is tested, with the phase-compatible assembled time, on
%   the only multileg structure the catalogs hold, L1 halo -> L2 halo ->
%   tulip.
%
%     0  settings                          (one place for every choice)
%     1  load                              (G1 periods; G2 edges, keys, 1 N / 0.5 N nesting)
%     2  graph and edge tables             (digraph; per-rung accepted fraction, t_f, lowest accepted rung)
%     3  planes                            (raw coverage and envelope; G3 monotone in the cap)
%     4  strengths                         (in/out strength at 1 N; a topology caveat)
%     5  relay closure                     (phase-compatible assembled times; G4)
%     6  hypothetical reverse graph        (G5 symmetry at an off-plane point; HYPOTHETICAL edges)
%     7  figures                           (results/*.png)
%
%   G gates are code gates: a FAIL of G1-G4 means the data or the
%   computation is not what this script assumes, and the summary throws;
%   G5 throws only when the hypothetical graph is on. Every printed number
%   from section 6 carries the word HYPOTHETICAL and never enters sections
%   3-5 or the figures.
%
%   A pipeline test runs this headless: if a struct `study_opts` exists in
%   the base workspace, its .outDir, .visible and .hypothetical are
%   honoured.
%
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.

%% 0  settings
clearvars -except study_opts; clc
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'lib'), fullfile(fileparts(here), 'costate_common'));

F_ref        = 1;                 % N, the rung every sheet has an accepted cell at: strengths and the graph figure use it
F_deep       = 0.5;               % N, the catalogs' deepest rung
nCaps        = 25;                % time caps, from the smallest accepted cell time to the largest (or the largest relay time)
hypothetical = false;             % the reverse graph by the x-axis symmetry: geometry only, off by default
relay        = struct('src', 'halo_L1', 'via', 'halo_L2', 'dst', 'tulip');   % key prefixes of the only multileg structure in the data
outDir       = fullfile(here, 'results');
visible      = 'on';
if exist('study_opts', 'var')
    if isfield(study_opts, 'outDir'),       outDir       = study_opts.outDir;       end
    if isfield(study_opts, 'visible'),      visible      = study_opts.visible;      end
    if isfield(study_opts, 'hypothetical'), hypothetical = study_opts.hypothetical; end
end
if ~exist(outDir, 'dir'), mkdir(outDir); end
gates = struct('id', {}, 'ok', {}, 'msg', {}, 'throws', {});
lines = {};                                                 % the summary file, one printed fact per line

%% 1  load
T     = load_catalog_sheets([]);
N     = orbit_node_table(T, struct('verbose', false));
nS    = numel(T);  nN = numel(N);
rungs = T(1).rungs_N;  nR = numel(rungs);
rC    = find(rungs == F_ref, 1);  rD = find(rungs == F_deep, 1);
tStar = T(1).constants.tStar_s;
keyN  = {N.key};
fprintf('\n== 1  load: %d sheets, %d orbits, rungs %s N\n', nS, nN, mat2str(rungs));
ePer = max([abs(cellfun(@(k) N(strcmp(keyN, k)).period_nd, {T.dep_key}) - [T.T_dep]), ...
            abs(cellfun(@(k) N(strcmp(keyN, k)).period_nd, {T.arr_key}) - [T.T_arr])]);
nNest = 0;  nViolAll = 0;  nAccX = 0;
for kk = 1:nS
    a = T(kk).accepted;  h = T(kk).has_solution;
    nNest    = nNest    + nnz(h(:, :, rD) & ~h(:, :, rC));   % SOLVED sets nest (the ladder descends from 1 N)
    nAccX    = nAccX    + nnz(a(:, :, rD) & ~a(:, :, rC));   % accepted at 0.5 N but refuted/gated at 1 N: a fact, not a gate
    nViolAll = nViolAll + nnz(h(:, :, 2:end) & ~h(:, :, 1:end-1));
end
lines{end+1} = sprintf('cells accepted at %g N but not at %g N (refuted or gated at the shallower rung): %d', F_deep, F_ref, nAccX);
fprintf('%s\n', lines{end});
perSolved   = sum(cell2mat(arrayfun(@(r) reshape(sum(sum(r.has_solution, 1), 2), 1, []), T, 'UniformOutput', false)), 1);
perAccepted = sum(cell2mat(arrayfun(@(r) reshape(sum(sum(r.accepted, 1), 2), 1, []), T, 'UniformOutput', false)), 1);
lines{end+1} = sprintf('solved cells per rung %s: %s', mat2str(rungs), mat2str(perSolved));
lines{end+1} = sprintf('accepted cells per rung: %s (certified witnesses; every count below uses these)', mat2str(perAccepted));
lines{end+1} = sprintf('whole-ladder non-monotone solved cells: %d (the 15 N starting rung closes fewer than 12 N: solver coverage, not gated)', nViolAll);
fprintf('%s\n%s\n%s\n', lines{end-2}, lines{end-1}, lines{end});
gates(end+1) = struct('id', 'G1', 'ok', ePer < 1e-9, 'msg', sprintf('node periods match the catalogs (max |dT| = %.1e)', ePer), 'throws', true);

%% 2  graph and edge tables
E = edge_tables(T);
G = network_graph(E, N);
src = findnode(G, {T.dep_key});  dst = findnode(G, {T.arr_key});          % node indices in SHEET order (used only for key checks)
fprintf('\n== 2  graph: %d nodes, %d edges (every attribute read from G.Edges; MATLAB sorts edges)\n', numnodes(G), numedges(G));
for rr = [rC rD]
    fprintf('\n  edge table at %g N (accepted cells)\n', rungs(rr));
    fprintf('  %-18s %-16s %7s %6s %8s %8s %20s\n', 'departure', 'arrival', 'frac', 'solved', 'tf_min', 'tf_med', 'lowest accepted');
    for kk = 1:nS
        fprintf('  %-18s %-16s %7.3f %6d %8.3f %8.3f %20s\n', E(kk).dep_key, E(kk).arr_key, E(kk).closed_frac(rr), E(kk).n_solved(rr), E(kk).tf_min_days(rr), E(kk).tf_med_days(rr), E(kk).F_low_str);
    end
end
nLowest = nnz(contains({E.F_low_str}, 'lowest tested'));
lines{end+1} = sprintf('edges accepted at the lowest tested rung (%g N): %d of %d; edges with at least one accepted cell at %g N: %d of %d', F_deep, nLowest, nS, F_ref, nnz([E.n_accepted(rC)] > 0), nS);
fprintf('%s\n', lines{end});
gates(end+1) = struct('id', 'G2', 'ok', nS == 56 && numedges(G) == 56 && all(src > 0) && all(dst > 0) && nNest == 0, ...
                      'msg', sprintf('56 edges on named nodes (findnode round-trips every key); every cell solved at %g N is solved at %g N (%d violations)', F_deep, F_ref, nNest), 'throws', true);

%% 3  planes: raw coverage and the known-feasible envelope
tfCells = [];
for kk = 1:nS
    tf = T(kk).tf_nd(T(kk).accepted)*tStar/86400;
    tfCells = [tfCells; tf(:)];
end
capLo = min(tfCells);  capHi = max(tfCells);                % section 5 may extend capHi for the relay times
caps  = linspace(capLo, capHi, nCaps);
P     = accessibility_plane(T, caps);
fprintf('\n== 3  planes: %d rungs x %d caps (%.3f .. %.3f days)\n', nR, nCaps, caps(1), caps(end));
fprintf('  raw coverage, pairs (rows = rungs, cols = caps):\n');  disp(P.pairs);
fprintf('  known-feasible envelope, pairs:\n');  disp(P.pairs_env);
monoCap  = all(all(diff(P.pairs, 1, 2) >= 0)) && all(all(diff(P.cells, 1, 2) >= 0)) && all(all(diff(P.pairs_env, 1, 2) >= 0)) && all(all(diff(P.cells_env, 1, 2) >= 0));
monoEnv  = all(all(diff(P.pairs_env, 1, 1) <= 0)) && all(all(diff(P.cells_env, 1, 1) <= 0));
cFull = find(P.pairs(rC, :) == nS, 1, 'first');
if isempty(cFull)
    lines{end+1} = sprintf('at %g N no sampled cap gives all %d direct pairs (max %d of %d under %.3f days)', F_ref, nS, max(P.pairs(rC, :)), nS, caps(end));
else
    lines{end+1} = sprintf('at %g N all %d direct pairs are covered from the first sampled cap of %.3f days', F_ref, nS, caps(cFull));
end
lines{end+1} = sprintf('at %g N, raw coverage %d of %d pairs under the largest cap; envelope at the %g N bound %d of %d', F_deep, P.pairs(rD, end), nS, rungs(1), P.pairs_env(1, end), nS);
fprintf('%s\n%s\n', lines{end-1}, lines{end});
gates(end+1) = struct('id', 'G3', 'ok', monoCap && monoEnv, 'msg', 'raw and envelope counts non-decreasing along the cap; envelope non-increasing toward deeper rungs', 'throws', true);

%% 4  strengths at the reference rung (descriptive; the between-family ranking is the catalog's topology)
tfRef = G.Edges.tf_min_days(:, rC);
eSrc  = findnode(G, G.Edges.EndNodes(:, 1));  eDst = findnode(G, G.Edges.EndNodes(:, 2));
outS  = accumarray(eSrc, 1./tfRef, [nN 1]);
inS   = accumarray(eDst, 1./tfRef, [nN 1]);
fprintf('\n== 4  strengths at %g N (weights 1 / t_f[days]); each tulip has 12 in-edges, each L2 halo 2 in and 4 out: rankings reflect the catalog, not transport\n', F_ref);
hdr = sprintf('  %-20s %4s %4s %9s %9s', 'orbit', 'in', 'out', 'in-str', 'out-str');
fprintf('%s\n', hdr);  lines{end+1} = ['strengths at the reference rung:' hdr];
for kk = 1:nN
    row = sprintf('  %-20s %4d %4d %9.3f %9.3f', G.Nodes.Name{kk}, indegree(G, kk), outdegree(G, kk), inS(kk), outS(kk));
    fprintf('%s\n', row);  lines{end+1} = row;
end

%% 5  relay closure

%% 6  hypothetical reverse graph

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
if any(~[gates.ok] & [gates.throws])
    error('orbit_network:orbit_network_study:gate', 'a G gate failed; see the summary');
end
```

- [ ] **Step 4: Run the script** (`run_matlab_file`) and read the output: G1-G3 PASS, the two edge tables have 56 rows, the planes print, the envelope row at 15 N is at least the raw row. Then `cd tests; ok = test_pipeline();` -> all PASS.
- [ ] **Step 5: Code Analyzer on the script**: no warnings (info-level growth notes on `gates(end+1)`/`lines{end+1}`/`tfCells` tolerated, no `%#ok`).
- [ ] **Step 6: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/orbit_network_study.m orbit_transfer/orbit_network/tests/test_pipeline.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: orbit_network_study sections 0-4 -- digraph of 18 orbits / 56 edges over accepted cells with per-rung coverage, t_f and lowest-accepted-rung attributes (read from G.Edges), the raw-coverage and known-feasible-envelope planes, in/out strengths at 1 N with the topology caveat; G1-G3; summary file; pipeline test row

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 7: sections 5-6 (relay closure with G4, hypothetical reverse graph with G5)

**Files:**
- Modify: `orbit_transfer/orbit_network/orbit_network_study.m` (fill the `%% 5` and `%% 6` markers)
- Modify: `orbit_transfer/orbit_network/tests/test_pipeline.m` (gate list for the network row becomes `{'G1', 'G2', 'G3', 'G4', 'G5'}`; add the second, hypothetical run described in Step 1)

**Interfaces:**
- Consumes: Task 6's workspace, `relay_closure(T, N, relay.src, relay.via, relay.dst)`, `symmetry_residual(x1, T, mu)`, `mirror_state`, `pumpkyn.cr3bp.prop`.
- Produces: `R` (routes), `best` struct array over the 8 (src, dst) pairs with `pair`, `t_route_days`, `t_sum_days` [1 x nR] (min over relays), `relay_key` {1 x nR}, `deepest_N`; `P.reach` [nR x nCaps] = `P.pairs` + the number of (src, dst) pairs with `t_route` under the cap; `usage` [nN x 1] relay-usage counts at the reference rung; when `hypothetical`, `Gh` (the digraph with the 56 reversed edges).

- [ ] **Step 1: Extend the pipeline test.** Change the network row's gate list to `{'G1', 'G2', 'G3', 'G4', 'G5'}`. After the `for ss = 1:size(scripts, 1)` loop and before the final `fprintf('test_pipeline: ...')`, add:

```matlab
%% the hypothetical reverse graph, which the default run never exercises
assignin('base', 'study_opts', struct('outDir', outDir, 'visible', 'off', 'hypothetical', true));
try
    txt = evalc(sprintf('evalin(''base'', ''run(''''%s'''')'')', fullfile(root, 'orbit_network_study.m')));
    ran = true;
catch err
    txt = err.message;  ran = false;
end
evalin('base', 'clear study_opts');
ok = check(ok, ran, 'orbit_network_study ran with hypothetical = true');
ok = check(ok, ran && ~isempty(regexp(txt, '\n\s*G5\s+PASS', 'once')), 'hypothetical run: G5 PASS');
ok = check(ok, ran && contains(txt, 'HYPOTHETICAL'), 'hypothetical run labels its numbers HYPOTHETICAL');
ok = check(ok, ran && contains(txt, '112 edges'), 'hypothetical run reports 112 edges (56 + 56 reversed)');
```

- [ ] **Step 2: Run `test_pipeline` to verify the new checks fail.**

- [ ] **Step 3: Fill section 5** (replace the `%% 5  relay closure` marker with):

```matlab
%% 5  relay closure (L1 halo -> L2 halo -> tulip: the only multileg structure in the data)
R = relay_closure(T, N, relay.src, relay.via, relay.dst);
fprintf('\n== 5  relay closure: %d routes %s -> %s -> %s, phase-compatible assembled times\n', numel(R), relay.src, relay.via, relay.dst);
pairKeys = unique(strcat({R.src_key}, ' -> ', {R.dst_key}), 'stable');
best = struct('pair', pairKeys, 't_route_days', [], 't_sum_days', [], 'relay_key', [], 'deepest_N', []);
for pp = 1:numel(pairKeys)
    sel = find(strcmp(strcat({R.src_key}, ' -> ', {R.dst_key}), pairKeys{pp}));
    tr  = vertcat(R(sel).t_route_days);                     % [nRelays x nR]
    ts  = vertcat(R(sel).t_sum_days);
    [tb, ib] = min(tr, [], 1, 'omitnan');                   % best relay per rung by the ASSEMBLED time
    tb(all(isnan(tr), 1)) = NaN;
    best(pp).t_route_days = tb;
    best(pp).t_sum_days   = min(ts, [], 1, 'omitnan');
    best(pp).relay_key    = arrayfun(@(rr) R(sel(ib(rr))).relay_key, 1:nR, 'UniformOutput', false);
    deep = find(~isnan(tb), 1, 'last');
    if isempty(deep), best(pp).deepest_N = NaN; else, best(pp).deepest_N = rungs(deep); end
    row = sprintf('  %-32s deepest rung closed %5g N | at %g N: best relay %-16s t_route %.3f d (t_sum %.3f d)', ...
                  pairKeys{pp}, best(pp).deepest_N, F_ref, best(pp).relay_key{rC}, tb(rC), best(pp).t_sum_days(rC));
    fprintf('%s\n', row);  lines{end+1} = row;
end
% extend the caps if the relay times exceed the direct ones, then the directed-reachability overlay
tRel = [best.t_route_days];
if max(tRel, [], 'omitnan') > caps(end)
    caps = linspace(capLo, max(tRel, [], 'omitnan'), nCaps);
    P = accessibility_plane(T, caps);
    lines{end+1} = sprintf('caps extended to %.3f days to cover the relay routes', caps(end));
    fprintf('%s\n', lines{end});
end
P.multileg = zeros(nR, nCaps);
for rr = 1:nR
    for cc = 1:nCaps
        P.multileg(rr, cc) = nnz(arrayfun(@(b) b.t_route_days(rr) <= caps(cc), best));
    end
end
P.reach = P.pairs + P.multileg;                            % directed reachable ordered pairs under the cap (direct + two-leg)
usage = zeros(nN, 1);
for pp = 1:numel(best)
    if ~isnan(best(pp).t_route_days(rC))
        kU = find(strcmp(keyN, best(pp).relay_key{rC}), 1);
        usage(kU) = usage(kU) + 1;
    end
end
for kk = find(usage > 0).'
    row = sprintf('  relay usage at %g N: %-16s serves %d of %d (L1 halo, tulip) pairs as the best assembled relay', F_ref, keyN{kk}, usage(kk), numel(best));
    fprintf('%s\n', row);  lines{end+1} = row;
end
lines{end+1} = sprintf('%d relay routes over %d (src, dst) pairs; pairs reaching %g N by relay: %d of %d; directed reachable ordered pairs at %g N under the largest cap: %d (direct %d + two-leg %d) of %d possible', ...
                       numel(R), numel(best), F_deep, nnz([best.deepest_N] == F_deep), numel(best), F_ref, P.reach(rC, end), P.pairs(rC, end), P.multileg(rC, end), nS + numel(best));
fprintf('%s\n', lines{end});
% G4: enumeration and consistency with the sheets (not a physics check)
g4 = numel(R) == 32 && numel(best) == 8;
for kk = 1:numel(R)
    L1 = T(R(kk).s1);  L2 = T(R(kk).s2);
    legs = squeeze(any(any(L1.accepted, 1), 2)).' & squeeze(any(any(L2.accepted, 1), 2)).';
    g4 = g4 && strcmp(L1.arr_key, L2.dep_key) && strncmp(L1.dep_key, relay.src, numel(relay.src)) && strncmp(L2.arr_key, relay.dst, numel(relay.dst)) ...
            && isequal(R(kk).closed, legs) && all(isnan(R(kk).t_route_days(~legs))) && all(~isnan(R(kk).t_route_days(legs))) ...
            && all(R(kk).t_route_days(legs) >= R(kk).t_sum_days(legs) - 1e-9) && isequal(sort(L1.sA_frac), sort(L2.sD_frac));
end
gates(end+1) = struct('id', 'G4', 'ok', g4, 'msg', '32 routes over 8 pairs, legs share the relay key and its phase grid, closed iff both legs accepted, t_route >= t_sum and NaN exactly when open', 'throws', true);
```

- [ ] **Step 4: Fill section 6** (replace the `%% 6  hypothetical reverse graph` marker with):

```matlab
%% 6  hypothetical reverse graph (x-axis symmetry: geometry only; the mirrored arc belongs to a mass-GAIN model)
mu    = T(1).constants.muStar;
resP  = zeros(nN, 1);  resS = zeros(nN, 1);
for kk = 1:nN
    x0 = N(kk).x0;  Tp = N(kk).period_nd;
    [~, X] = pumpkyn.cr3bp.prop([0 0.37*Tp Tp], x0.', mu);
    resP(kk) = max(abs(X(end, 1:6).' - x0));                % periodic closure, independent of the symmetry
    resS(kk) = symmetry_residual(X(2, 1:6).', Tp, mu);      % mirror image of an OFF-PLANE point lies on the orbit
end
fprintf('\n== 6  symmetry: periodic closure max %.1e; off-plane mirror residual max %.1e (every stored x0 is on the plane, so the check is made at phi(x0, 0.37 T))\n', max(resP), max(resS));
gates(end+1) = struct('id', 'G5', 'ok', all(resP < 1e-8) && all(resS < 1e-8), ...
                      'msg', sprintf('every node closes (%.1e) and is x-axis symmetric at an off-plane point (%.1e), bars 1e-8; throws only when the hypothetical graph is on', max(resP), max(resS)), 'throws', hypothetical);
if hypothetical
    Eh = G.Edges;
    Eh.EndNodes = fliplr(Eh.EndNodes);
    Eh.hypothetical = true(height(Eh), 1);
    Eh.n_accepted   = zeros(size(Eh.n_accepted));          % no certificate travels with a mirrored arc
    Gh = digraph([G.Edges; Eh], G.Nodes);
    [~, sccSize] = conncomp(Gh, 'Type', 'strong');
    D  = distances(Gh);  D(logical(eye(nN))) = Inf;
    nReach = nnz(isfinite(D));
    row = sprintf('  HYPOTHETICAL: %d edges (%d + %d reversed); strongly connected components %d (largest %d nodes); ordered node pairs joined by a directed path %d of %d', ...
                  numedges(Gh), nS, nS, numel(sccSize), max(sccSize), nReach, nN*(nN - 1));
    fprintf('%s\n', row);  lines{end+1} = row;
    lines{end+1} = 'HYPOTHETICAL numbers never enter sections 3-5 or the figures';
else
    fprintf('  hypothetical reverse graph off (hypothetical = false); G5 is still reported\n');
end
```

- [ ] **Step 5: Run the script, then `test_pipeline`.** Expected: G4 and G5 PASS; the hypothetical run prints `112 edges` and `HYPOTHETICAL`; all checks PASS. If G5 fails on a node, print `resP` and `resS` with keys; a halo residual of order 1e-9 is normal, 1e-3 is a bug (check that `x0` is `rv(1, :)` of `get_family_orbit`).
- [ ] **Step 6: Code Analyzer clean (same tolerance as Task 6).**
- [ ] **Step 7: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/orbit_network_study.m orbit_transfer/orbit_network/tests/test_pipeline.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: sections 5-6 -- relay closure L1 halo -> L2 halo -> tulip with phase-compatible assembled times (best relay per rung, deepest rung reached, directed-reachability overlay, relay usage; G4) and the hypothetical reverse graph (G5 periodic closure + off-plane symmetry; reversed edges carry no certificate; SCC and reachability counts, HYPOTHETICAL-prefixed; exercised by the pipeline test)

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 8: section 7, the four figures

**Files:**
- Modify: `orbit_transfer/orbit_network/orbit_network_study.m` (fill the `%% 7  figures` marker)
- Modify: `orbit_transfer/orbit_network/tests/test_pipeline.m` (the network row's figure list becomes `{'network_graph.png', 'accessibility_plane.png', 'ladder_closure.png', 'relay_closure.png'}`)

**Interfaces:**
- Consumes: `G`, `E`, `P` (`pairs`, `pairs_env`, `multileg`, `caps_days`), `best`, `rungs`, `rC`, `rD`, `nS`, `outDir`, `visible`.
- Produces: the four PNGs in `outDir`, and a `lines` entry naming them.

- [ ] **Step 1: Extend the pipeline test's figure list** as above; run it to see the four `wrote ...png` checks FAIL.

- [ ] **Step 2: Fill section 7** (replace the `%% 7  figures` marker with):

```matlab
%% 7  figures (never the hypothetical edges)
axStyle = @(ax) set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
rungLab = arrayfun(@(r) sprintf('%g', r), rungs, 'UniformOutput', false);
famOf   = @(k) regexprep(k, '_tau.*|_Np.*', '');          % dro, dpo, halo_L1, halo_L2, tulip
famList = {'dro', 'dpo', 'halo_L1', 'halo_L2'};

% F1: the graph at the reference rung and at the deepest rung; width = accepted fraction, color = t_f; built from the edge table
fig = figure('Color', 'w', 'Position', [60 60 1400 620], 'Visible', visible);
for pp = 1:2
    rr = [rC rD];  rr = rr(pp);
    Er = G.Edges(G.Edges.closed_frac(:, rr) > 0, :);        % accepted edges at this rung, attributes travel with the rows
    Gr = digraph(Er, G.Nodes);
    ax = subplot(1, 2, pp);  axStyle(ax);
    h  = plot(ax, Gr, 'Layout', 'layered', 'Direction', 'right', 'NodeLabel', Gr.Nodes.Name, 'Interpreter', 'none', ...
              'MarkerSize', 6, 'NodeColor', [0.1 0.3 0.7], 'ArrowSize', 8);
    h.LineWidth = 0.5 + 4*Gr.Edges.closed_frac(:, rr);
    h.EdgeCData = Gr.Edges.tf_min_days(:, rr);
    colormap(ax, 'parula');  cb = colorbar(ax);  cb.Label.String = 'min t_f over accepted phases [days]';
    title(ax, sprintf('%g N: %d edges, width = accepted fraction of the 6 x 6 grid', rungs(rr), numedges(Gr)), 'Color', 'k', 'Interpreter', 'none');
end
f1 = fullfile(outDir, 'network_graph.png');  exportgraphics(fig, f1, 'Resolution', 150);  if strcmp(visible, 'off'), close(fig); end

% F2: raw coverage and the envelope, with boolean-mask contours (all direct pairs; first relay route)
fig = figure('Color', 'w', 'Position', [60 60 1400 560], 'Visible', visible);
planes = {P.pairs, 'raw coverage: pairs with an accepted cell AT this thrust under the cap'; P.pairs_env, 'known-feasible envelope: pairs with an accepted cell at this thrust OR LOWER under the cap'};
for pp = 1:2
    ax = subplot(1, 2, pp);  axStyle(ax);  hold(ax, 'on');
    imagesc(ax, caps, 1:nR, planes{pp, 1});  axis(ax, 'tight');  set(ax, 'YDir', 'normal', 'YTick', 1:nR, 'YTickLabel', rungLab);
    colormap(ax, 'hot');  cb = colorbar(ax);  cb.Label.String = sprintf('pairs (of %d)', nS);
    mAll = double(planes{pp, 1} == nS);
    if any(mAll(:)) && ~all(mAll(:)), contour(ax, caps, 1:nR, mAll, [0.5 0.5], 'c-', 'LineWidth', 2); end
    mRel = double(P.multileg >= 1);
    if any(mRel(:)) && ~all(mRel(:)), contour(ax, caps, 1:nR, mRel, [0.5 0.5], 'g--', 'LineWidth', 2); end
    xlabel(ax, 'time cap [days]', 'Color', 'k');  ylabel(ax, 'thrust bound [N]', 'Color', 'k');
    title(ax, {planes{pp, 2}, 'cyan: all 56 direct pairs; green dashed: first relay route (assembled time)'}, 'Color', 'k', 'FontSize', 9);
end
f2 = fullfile(outDir, 'accessibility_plane.png');  exportgraphics(fig, f2, 'Resolution', 150);  if strcmp(visible, 'off'), close(fig); end

% F3: catalog closure down the ladder, by departure family (accepted fraction; the 15 N dip is ladder coverage)
fig = figure('Color', 'w', 'Position', [60 60 900 480], 'Visible', visible);
ax  = axes(fig);  axStyle(ax);  hold(ax, 'on');  grid(ax, 'on');
cf  = vertcat(E.closed_frac);  depFam = cellfun(famOf, {E.dep_key}, 'UniformOutput', false);
for ff = 1:numel(famList)
    sel = strcmp(depFam, famList{ff});
    plot(ax, 1:nR, mean(cf(sel, :), 1), '-o', 'LineWidth', 1.8, 'DisplayName', sprintf('%s departures (%d edges)', famList{ff}, nnz(sel)));
end
set(ax, 'XTick', 1:nR, 'XTickLabel', rungLab);
xlabel(ax, 'thrust rung [N] (shallow to deep)', 'Color', 'k');  ylabel(ax, 'mean accepted fraction of the phasing grid', 'Color', 'k');
legend(ax, 'Location', 'southwest', 'Interpreter', 'none', 'TextColor', 'k', 'Color', 'w');
title(ax, 'catalog closure down the ladder (an unaccepted cell is unresolved, not infeasible; the 15 N dip is ladder coverage)', 'Color', 'k', 'FontSize', 9);
f3 = fullfile(outDir, 'ladder_closure.png');  exportgraphics(fig, f3, 'Resolution', 150);  if strcmp(visible, 'off'), close(fig); end

% F4: relay closure, best route per (L1 halo, tulip) pair: assembled time (solid) and the optimistic sum (dotted)
fig = figure('Color', 'w', 'Position', [60 60 900 520], 'Visible', visible);
ax  = axes(fig);  axStyle(ax);  hold(ax, 'on');  grid(ax, 'on');
cols = parula(numel(best));                                % NOT lines(): that name is the script's summary cellstr
for pp = 1:numel(best)
    plot(ax, 1:nR, best(pp).t_route_days, '-o', 'Color', cols(pp, :), 'LineWidth', 1.6, 'DisplayName', best(pp).pair);
    plot(ax, 1:nR, best(pp).t_sum_days,   ':',  'Color', cols(pp, :), 'LineWidth', 1.0, 'HandleVisibility', 'off');
end
set(ax, 'XTick', 1:nR, 'XTickLabel', rungLab);
xlabel(ax, 'thrust rung [N]', 'Color', 'k');  ylabel(ax, 'two-leg time [days]: solid = phase-compatible assembled, dotted = sum of leg minima', 'Color', 'k');
legend(ax, 'Location', 'northwest', 'Interpreter', 'none', 'TextColor', 'k', 'Color', 'w', 'FontSize', 8);
title(ax, 'relay closure L1 halo -> L2 halo -> tulip: best assembled relay per pair and rung', 'Color', 'k', 'Interpreter', 'none');
f4 = fullfile(outDir, 'relay_closure.png');  exportgraphics(fig, f4, 'Resolution', 150);  if strcmp(visible, 'off'), close(fig); end
lines{end+1} = sprintf('figures: %s, %s, %s, %s', f1, f2, f3, f4);
fprintf('%s\n', lines{end});
```

- [ ] **Step 3: Run the script with figures visible, then read each PNG** with the Read tool: F1 node labels show underscores literally, no NaN-width edges at 0.5 N; F2 two panels, the envelope never below the raw plane, contours present only where the masks change; F3 four lines with the 15 N dip; F4 eight solid lines with dotted companions never above them. Fix what the eye catches before running the test.
- [ ] **Step 4: Run `test_pipeline` and `run_tests`**: all PASS.
- [ ] **Step 5: Code Analyzer clean.**
- [ ] **Step 6: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/orbit_network_study.m orbit_transfer/orbit_network/tests/test_pipeline.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: section 7 figures -- the graph at 1 N and 0.5 N from the edge table (width = accepted fraction, color = t_f), raw-coverage and envelope planes with boolean-mask contours, catalog closure down the ladder by departure family, relay closure (assembled vs sum) per (L1 halo, tulip) pair; pipeline test checks all four PNGs

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 9: README and the orbit_transfer README line

**Files:**
- Modify: `orbit_transfer/orbit_network/README.md`
- Modify: `orbit_transfer/README.md` (the `**Across the catalogs (2026-11-05)` paragraph: drop "and, next," so it reads as built)

- [ ] **Step 1: In `orbit_network/README.md`** add rows to the `What is here` table:

```markdown
| `orbit_network_study.m` | topic 2: digraph of 18 orbits and 56 edges over ACCEPTED cells with per-rung accepted fraction, min/median t_f and the lowest accepted rung; raw-coverage and known-feasible-envelope planes over (thrust bound, time cap) with the directed-reachability overlay; in/out strengths at 1 N with the topology caveat; L1 halo -> L2 halo -> tulip relay closure with PHASE-COMPATIBLE assembled times; the hypothetical reverse graph (off by default, geometry only); G1-G5 gates; four figures and `results/network_summary.txt` |
| `lib/edge_tables.m`, `lib/network_graph.m`, `lib/accessibility_plane.m`, `lib/relay_closure.m`, `lib/mirror_state.m`, `lib/symmetry_residual.m` | the network's computations, each with a demo and a hand-checked test in `tests/test_lib.m` |
| `reviews/` | the GPT-6 Astra review of the spec and plan (2026-11-05) and its dispositions |
| `doc/orbit_network_note.tex` | the note (Task 10) |
```

Replace the paragraph beginning `Topic 2, the network itself` with:

```markdown
Both scripts print every gate PASS/FAIL and throw on a failed code gate;
hypothesis verdicts (H1-H2) are results either way. `results/network_summary.txt`
holds every number the note cites.
```

Add under `## Run`:

```matlab
orbit_network_study         % topic 2: the network, the coverage and envelope planes, relay closure
```

Extend `## What it does not claim` with:

```markdown
- An unaccepted cell or rung is unresolved, never a demonstrated infeasibility;
  the "lowest accepted rung" is not a thrust floor.
- Weak connectivity and betweenness are not reported: only L2 halos have both
  in- and out-edges, so those statistics would restate the catalog's topology.
- The hypothetical reverse graph is geometry: a mirrored thrusting arc belongs
  to a mass-GAIN model and carries no certificate.
```

- [ ] **Step 2: In `orbit_transfer/README.md`** change `and, next, the data-only low-thrust network` to `and the data-only low-thrust network`.
- [ ] **Step 3: Check `orbit_network/README.md` is ASCII** (`LC_ALL=C grep -nP '[^\x00-\x7F]' orbit_transfer/orbit_network/README.md` prints nothing).
- [ ] **Step 4: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/README.md orbit_transfer/README.md
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: README for topic 2 (what the network script computes, what it does not claim); orbit_transfer README wording

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 10: the note, `doc/orbit_network_note.tex`

**Files:**
- Create: `orbit_transfer/orbit_network/doc/orbit_network_note.tex`
- Create: `orbit_transfer/orbit_network/doc/figures/` (copies of the five PNGs: `stability_cost.png`, `network_graph.png`, `accessibility_plane.png`, `ladder_closure.png`, `relay_closure.png`, copied from `results/` after a fresh run of both scripts with a REBUILT node table)
- Produce: `orbit_transfer/orbit_network/doc/orbit_network_note.pdf` (PDFs are git-ignored; the `.tex` and the figures are committed)

**Interfaces:**
- Consumes: `results/network_summary.txt` and the printed output of `stability_cost_study.m` for every number quoted. No number is typed from memory.

- [ ] **Step 1: Rebuild and run fresh.** In the shared session:

```matlab
cd('/Users/msc/Desktop/optimal_control/orbit_transfer/orbit_network');
addpath(pwd, fullfile(pwd, 'lib'), fullfile(fileparts(pwd), 'costate_common'));
orbit_node_table(load_catalog_sheets([]), struct('rebuild', true, 'verbose', false));
```
then `run_matlab_file` on `stability_cost_study.m` and on `orbit_network_study.m`; copy the five PNGs into `doc/figures/`; keep the two outputs and `results/network_summary.txt` open while writing.

- [ ] **Step 2: Write the note** with this skeleton; the lemma, its proof and the fixed sentences are verbatim; the prose under each `%` comment comes from the named source, and the comment is deleted once filled:

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
% budget regimes off it. Their stability reading: unstable orbits share a chaotic sea and connect
% cheaply; stable ones sit in tori and are hard. Why we care: the same question for our DRO -> tulip
% transfers, whose direct min-time branch we have continued down to 0.067 N.

\section{Their problem and ours}
% Table: planar vs spatial; one Jacobi level vs Jacobi gaps from -0.16 to 0.33; energy-preserving
% impulses vs mass-varying continuous thrust; proxy overlap cost vs accepted min-time extremals
% (multiple shooting, conjugate test, hypothesis gates); phase as future work vs 6 x 6 phasing grids.
\begin{lemma}[unrestricted versus relay-constrained optimum]
Let $\mathcal{F}_{AB}$ be the set of feasible trajectories of a specified $A \to B$ optimal control
problem, including the spacecraft mass and every boundary and path constraint. Let
$\mathcal{R}_{ACB} \subseteq \mathcal{F}_{AB}$ be the subset of trajectories that pass through the
relay orbit $C$ with compatible full states at every join and no relay-only operation. Then for the
same cost functional $J$,
\[
\inf_{x \in \mathcal{F}_{AB}} J[x] \;\le\; \inf_{x \in \mathcal{R}_{ACB}} J[x].
\]
\end{lemma}
\begin{proof}
$\mathcal{R}_{ACB}$ is a subset of $\mathcal{F}_{AB}$, and an infimum over a set is no larger than
the infimum over any of its subsets.
\end{proof}
% Fixed sentences (verbatim): "The lemma is a statement about global infima, not about the costs
% returned by particular direct or relay solution procedures. A multileg advantage in Braik and Ross
% is relative to their restricted direct-maneuver class; it is not a comparison against the global
% optimum over a feasible set that already contains those multileg trajectories. A relay witness
% below the termination thrust of our computed direct branch would establish another feasible
% unrestricted transfer; it would diagnose a limitation of that branch or of the search, not passage
% through a proven global feasibility barrier."
% Mass composition (verbatim): "A leg solved from m_0 can be flown from the previous leg's arrival
% mass m_1 by scaling its thrust and mass histories by q = m_1/m_0 <= 1: the acceleration F/m and the
% path are unchanged, the mass equation still holds, and the thrust bound is respected. The assembled
% route is feasible without a new solve; the scaled leg is not thereby optimal or certified for the
% changed initial-mass problem."

\section{The data}
% The four catalogs, the sheet table (56 rows: 16 + 16 + 16 + 8), rungs 15 -> 0.5 N, the 6 x 6 grids,
% the accepted mask (has_solution, conj_pass = 1, the three hypothesis-gate margins), and the counts
% (quote the "solved cells per rung" and "accepted cells per rung" lines of network_summary.txt).
% The node table (18 orbits) as a booktabs table: key, period, C_J, |lambda|max, sigma, stable
% (tolerance 1e-3 stated; "no resolved exponential instability in the computed multipliers").
% The closure fact: nested between 1 N and 0.5 N; the 15 N starting rung closes fewer cells than 12 N.

\section{The graph}
% Nodes, edges, attributes (accepted fraction, min/median t_f, lowest accepted rung: a successful
% rung bounds nothing below it). Figure: network_graph.png. Quote the "edges accepted at the lowest
% tested rung" line. The strength table (from the summary) with the caveat: each tulip has 12
% in-edges, each L2 halo 2 in and 4 out; between-family rankings reflect the catalog's topology.

\section{The planes}
% Figure: accessibility_plane.png. Thrust is a capability axis for one spacecraft, not a consumable
% budget. Raw coverage (what the catalogs hold at that thrust) against the known-feasible envelope
% (every accepted trajectory at that thrust or lower). The diagnostics, named as diagnostics: direct
% availability, and directed reachable ordered pairs with the relay routes (at most 56 + 8 of 306).
% Figure: ladder_closure.png with the 15 N dip explained.

\section{Relay closure}
% Figure: relay_closure.png. The 32 routes, the phase-compatible assembled time
% t_1(a,b) + T_C mod(s_D2(c) - s_A1(b), 1) + t_2(c,d), the optimistic sum beneath it, the best relay
% per (L1 halo, tulip) pair at 1 N and the deepest rung each pair reaches (quote the relay lines and
% the usage lines from the summary). State again that relay-versus-direct cost is untestable here
% and why (no direct L1 halo -> tulip edge exists).

\section{Stability versus cost}
% Figure: stability_cost.png. The two responses over accepted cells, the nested regression, the
% H1/H2 lines, the coefficient t tests and 95% confidence intervals, the group means labelled
% descriptive. Fixed wording (verbatim): "Departure instability was not associated with shorter
% sampled minimum durations or greater catalog retention at the nominal significance threshold in
% this exploratory sheet-level model. The time model's significant joint stability term is
% associated primarily with arrival instability on the halo-to-halo sheets. These results neither
% establish absence of a departure-stability effect nor test the mechanism in the substantially
% lower-thrust regime. There are 14 distinct departure orbits and 8 distinct arrivals; sheets share
% orbit nodes, and the reported OLS p-values are nominal model-based summaries."

\section{The hypothetical reverse graph}
% One paragraph: the x-axis symmetry with time reversal maps a thrusting arc to a thrusting arc of a
% mass-GAIN model (thrust acceleration reflected by diag(1,-1,1), mass running m_f -> m_0); the
% reversed edges are geometry, carry no certificate, and are off by default; G5 and what it checks
% (periodic closure and symmetry at an off-plane point); the HYPOTHETICAL counts if quoted at all.

\section{What the legs campaign would add}
% DRO -> DPO and DRO -> L2 halo legs down toward 70 mN would make DRO -> tulip the first pair with
% both a direct edge and a relay route. A relay witness below the direct branch's termination would
% be a new feasible unrestricted transfer (lemma). The present data supports candidate
% prioritization among the L2 halos (relay usage, strengths), not identification of the best DRO
% relay: the existing L1 -> L2 legs say nothing about the missing DRO -> L2 leg costs.

\begin{thebibliography}{9}
\bibitem{braik2026} S.~Braik and S.~D.~Ross, ``Orbital Networks in the Three-Body Problem,'' arXiv:2605.31543, 2026.
\end{thebibliography}
\end{document}
```

Six pages is the target; figures are `\includegraphics[width=\textwidth]{figures/<name>.png}`.

- [ ] **Step 3: Compile and clean**

```bash
cd /Users/msc/Desktop/optimal_control/orbit_transfer/orbit_network/doc
/Library/TeX/texbin/pdflatex -interaction=nonstopmode orbit_network_note.tex > /dev/null && /Library/TeX/texbin/pdflatex -interaction=nonstopmode orbit_network_note.tex | grep -i "undefined\|error" ; rm -f *.aux *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz
ls -la orbit_network_note.pdf
```
Expected: a PDF of 5 to 7 pages, nothing printed by the grep.

- [ ] **Step 4: Read the PDF** (Read tool, pages 1-7): every figure present, every quoted number matches the summary file or the study output, no comment stubs left, the verbatim sentences present, `\author{Michael Casey}` only.

- [ ] **Step 5: Commit (controller, after review)**

```bash
cd /Users/msc/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add orbit_transfer/orbit_network/doc/orbit_network_note.tex orbit_transfer/orbit_network/doc/figures
/Library/Developer/CommandLineTools/usr/bin/git commit -m "orbit_network: the note (doc/orbit_network_note.tex, ~6 pp) -- Braik-Ross framing against a mass-varying min-time OCP, the set-inclusion lemma and the mass-composition remark, the sheet and node tables over accepted cells, the graph and its lowest-accepted-rung attribute, raw-coverage and envelope planes, phase-compatible relay closure, the stability-versus-cost result in its reviewed wording, the hypothetical reverse graph, and what the legs campaign would add

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

## Self-review (done while writing, after the Astra revision)

- Spec coverage: 7.0-7.7 -> Tasks 6, 7, 8; 8 lib -> Tasks 1-5 (`partial_f_test` already built and revised); 9 tests -> every task extends `test_lib` or `test_pipeline`, the runner already lists both; 10 note -> Task 10; 11 out of scope respected.
- Astra dispositions each have a home: 1, 2 -> Task 5/7/10; 3 -> Task 1/7/10; 4 -> Task 2 and the built loader; 5 -> built code and Task 10; 6 -> Task 2; 7 -> Task 4; 8, 9 -> Task 6/7; 10 -> Task 3; 11A-C -> Task 6/7; 12A-B -> Task 7; 13, 14 -> built code; 15A-D -> Tasks 7, 8, and the tests named in each task.
- Type consistency: `E` fields (`sheet_index`, `n_solved`, `n_accepted`, `closed_frac`, `tf_min_days`, `tf_med_days`, `F_low_N`, `F_low_str`) are the same in Tasks 2, 3, 6, 8; `P` fields (`pairs`, `cells`, `closed`, `pairs_env`, `cells_env`, `multileg`, `reach`, `caps_days`) in Tasks 4, 6, 7, 8; `R` fields (`s1`, `s2`, `closed`, `t_sum_days`, `t_route_days`, `wait_days`, `cells`) and `best` fields in Tasks 5, 7, 8; `relay_closure(T, N, src, via, dst)` identical in Tasks 5 and 7; `network_graph(E, N)` identical in Tasks 3, 6, 8.
- Review Focus items each have an owning test: 1 -> Tasks 2, 4; 2 -> Task 2; 3 -> Task 3; 4 -> Task 5; 5 -> Task 7; 6 -> Tasks 1, 7; 7 -> Task 6 (the `isempty(cFull)` branch).
- Known hazard: the name `lines` is the script's summary cellstr and a MATLAB colormap function; Task 8 uses `parula`.
