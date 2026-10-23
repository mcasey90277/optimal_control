# Drone Route 1 (Direct NLP) Implementation Plan (Stage 1)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Solve the paper's multi-agent kinematic drone problem (arXiv:2610.30388, eq. 19) as a CasADi/IPOPT trapezoidal NLP with hard pairwise clearance (the reference) and with the paper's own log barrier (Route 2's twin), gate every solve, reproduce Table 1 for M = 4..32 with certified numbers, and ship the M = 4 and M = 8 fixtures Route 2 seeds from.

**Architecture:** `optimal_control_examples/drone/route1_direct/`: `solve_drones_direct` (one NLP, modes `hard`/`barrier`, 2D or 3D-ball, vectorised pair constraints at nodes and Hermite midpoints), `verify_drones_direct` (gates G1-G5 with its own `.ok/.why`), `drone_table1` (the Table 1 loop writing a tracked CSV and the fixtures), `plot_drones`/`movie_drones`, a batch launcher for the M = 64/128 stretch. Everything reads the plant and geometry from `drone_dynamics/` (Stage 0, DONE before this plan starts).

**Tech Stack:** MATLAB R2026a, CasADi 3.7.0 (`~/casadi-3.7.0`, resolved from `getenv('HOME')` if `casadi.Opti` is not on the path), IPOPT with MUMPS. No oclib use in this stage.

**Spec:** `docs/superpowers/specs/2026-10-23-drone-multiagent-design.md` (Stage 1 section).

## Global Constraints

- House MATLAB style (header quartet, `if nargin == 0` demo in library functions, no `i`/`j`, no `norm`, no `%#ok`, no emojis); paths relative to `mfilename('fullpath')`; no `orbit_transfer` references.
- Defect constraints are added FIRST to the `Opti` so their multipliers are `opti.lam_g(1 : 3M N)`; never read `opti.dual` (canonicalised orientation, memory `casadi-optidual-sign-bug`).
- The barrier in the NLP is a second copy of `barrier_shape`'s two formulas (CasADi needs `if_else`); `solve_drones_direct` returns the barrier evaluated at the solution nodes and the test pins it to `barrier_terms` at 1e-12.
- Verbatim numbers: `T = 10`, `N = 100` default, bump amplitude `A = 0.05`, IPOPT `tol 1e-10`, `acceptable_tol 1e-8`, `max_iter 3000`; gates G1 `−1e-2·2r`, G3 defect `1e-9`, G4 `1e-3`, G5 `1e-6`, G6 `5e-4` at N = 100 and `5e-5` at N = 400.
- Slow solves (M = 32 mesh gate, M = 64/128) never run inside tests; they run from the batch launcher under the campaign discipline (nohup script file, log to file, `JOB DONE` marker).
- `.mat` fixtures are force-added (`git add -f`); each has a generator (`drone_table1`).
- Subagents never commit; commit messages end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`; do not stage the unrelated `missiles/` deletions.

## Review Focus

1. Straight-line initial guess with no bump: all agents coincide at T/2, the barrier's `sqrt` has an infinite derivative there and `hard` mode starts from a maximally infeasible point. The bump is mandatory; `opts.bump = 0` must be refused by identifier. (Task 1 test.)
2. 2D layouts solved as 3D: `v_z` is pinned by equality, so the z rows of X must come back exactly zero and no agent may dodge in z. (Task 1 test: `max |X(3:3:end,:)| == 0`.)
3. A solve that IPOPT reports as failed must not be verified as a solution: `verify_drones_direct` reads `out.success` first and its `.why` names the IPOPT status. (Task 2 test with `max_iter = 2`.)
4. The mirror gate on a non-symmetric layout (`sphere`) is meaningless; `verify_drones_direct` must skip it with a note rather than fail. (Task 2.)
5. Fixture staleness: `test_table1_fixture` re-derives the gates from the stored X, V and compares the stored J with a trapezoid re-integration, so an edited fixture cannot sit green. (Task 3.)

---

### Task 1: `solve_drones_direct`

**Files:**
- Create: `optimal_control_examples/drone/route1_direct/solve_drones_direct.m`
- Create: `optimal_control_examples/drone/route1_direct/casadi_path.m` (one-line path helper, shared by the folder)
- Test: `optimal_control_examples/drone/route1_direct/tests/test_direct_oracle.m`

**Interfaces:**
- Consumes: `drone_layouts`, `pair_list`, `clearance_terms`, `barrier_terms`, `lower_bound` (Stage 0).
- Produces: `out = solve_drones_direct(x0, xT, geom, opts)` with `opts` fields `.N [100] .T [10] .mode ['hard'|'barrier'] .sense [+1] .bump [0.05] .pairs [pair_list(M)] .midpoints [true] .maxIter [3000] .print [0] .X0 .V0 [warm start, same shapes as X, V]`; `out` fields `.X [3M x N+1] .V [3M x N+1] .tN [1 x N+1] .J .JLB .success .ipoptStatus .iterations .maxDefect .lamDefect [3M x N] .minHnodes .Unodes [1 x N+1] (barrier mode, else []) .mode .N .T .geom .pairs .wall`.

- [ ] **Step 1: Write the failing oracle test**

```matlab
function ok = test_direct_oracle()
%% Purpose:
%
%   THE closed-form oracle for the direct solver: two agents head-on,
%   (-0.4,0) <-> (0.4,0), r = 0.02, T = 10. In centre-of-mass and relative
%   coordinates the cost is integral |cdot|^2 + 1/4 |wdot|^2, c rests, and
%   w = x1 - x2 follows the shortest path from (-0.8,0) to (0.8,0) outside
%   the disk |w| < 0.04 at constant speed:
%       ell = 2 sqrt(D^2 - R^2) + R (pi - 2 acos(R/D)),  D = 0.8, R = 0.04
%       J*  = ell^2 / (4 T) = 0.0641602   (J_LB = 0.064)
%   The NLP knows none of this. Its optimum sits BELOW J* by O(h^2) (the
%   sampled constraint lets the interpolant cut the arc's corner) and
%   converges from below: 5e-4 at N = 100, 5e-5 at N = 400.
%
%   Also: the z rows are exactly zero (2D), the defects are re-evaluated,
%   and a bump of zero is refused.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));

[x0, xT, geom] = drone_layouts('headon');
T = 10;  D = 0.8;  R = 2*geom.r;
ell   = 2*sqrt(D^2 - R^2) + R*(pi - 2*acos(R/D));
Jstar = ell^2/(4*T);
ok = chk(ok, abs(Jstar - 0.0641602) < 5e-8, sprintf('closed form J* = %.7f', Jstar));

out100 = solve_drones_direct(x0, xT, geom, struct('N', 100, 'T', T, 'mode', 'hard'));
ok = chk(ok, out100.success, sprintf('N = 100 solve succeeded (%s, %d iterations, %.1f s)', out100.ipoptStatus, out100.iterations, out100.wall));
rel100 = (out100.J - Jstar)/Jstar;
ok = chk(ok, abs(rel100) < 5e-4, sprintf('N = 100: J = %.7f, (J - J*)/J* = %+.2e (bar 5e-4)', out100.J, rel100));
ok = chk(ok, rel100 <= 1e-6, 'N = 100: the discrete optimum does not exceed J* (it cuts the corner)');
ok = chk(ok, out100.J >= out100.JLB, sprintf('J >= J_LB = %.4f', out100.JLB));
ok = chk(ok, max(max(abs(out100.X(3:3:end,:)))) == 0 && max(max(abs(out100.V(3:3:end,:)))) == 0, '2D: z rows exactly zero');
ok = chk(ok, out100.maxDefect < 1e-9, sprintf('defects re-evaluated: %.1e', out100.maxDefect));
ok = chk(ok, out100.minHnodes >= -1e-9, sprintf('nodal clearance %.2e >= 0', out100.minHnodes));

out400 = solve_drones_direct(x0, xT, geom, struct('N', 400, 'T', T, 'mode', 'hard'));
rel400 = (out400.J - Jstar)/Jstar;
ok = chk(ok, out400.success && abs(rel400) < 5e-5, sprintf('N = 400: (J - J*)/J* = %+.2e (bar 5e-5)', rel400));
ok = chk(ok, abs(rel400) < abs(rel100)/4, 'second-order convergence: the N = 400 gap is under a quarter of the N = 100 gap');

try
    solve_drones_direct(x0, xT, geom, struct('bump', 0));
    ok = chk(ok, false, 'bump = 0 refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'solve_drones_direct:bump'), 'bump = 0 refused by solve_drones_direct:bump');
end

if ok, fprintf('TEST_DIRECT_ORACLE: ALL PASS\n'); else, fprintf('TEST_DIRECT_ORACLE: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 2: Run it to verify it fails**

Expected: `Unrecognized function or variable 'solve_drones_direct'`.

- [ ] **Step 3: Implement the path helper and the solver**

```matlab
function casadi_path()
%% Purpose:
%
%   Put CasADi on the path if casadi.Opti is not already there, from the
%   user's home (the cr3bp_minenergy_pmp convention). No absolute path is
%   written down.
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
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if isempty(which('casadi.Opti'))
    addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
end
if isempty(which('casadi.Opti'))
    error('casadi_path:missing', 'CasADi not found (expected ~/casadi-3.7.0)');
end
end
```

```matlab
function out = solve_drones_direct(x0, xT, geom, opts)
%% Purpose:
%
%   The multi-agent kinematic drone problem (Zhang, Liu, Darbon &
%   Karniadakis, arXiv:2610.30388, eq. 19) as a trapezoidal direct
%   transcription solved by IPOPT through CasADi:
%
%       min  integral_0^T sum_i 1/2 |v_i|^2 dt
%       s.t. xdot_i = v_i,  x_i(0) = x0_i,  x_i(T) = xT_i,
%            |x_i - x_j| >= 2 r  (mode 'hard': at every node and every
%                                 Hermite midpoint, squared form)
%            or the paper's C1 log barrier added to the running cost
%                                (mode 'barrier': the smooth twin Route 2
%                                 shoots; no inequalities)
%            |x_i - c| >= R_ball  when geom.ball is set (3D)
%
%   2D layouts (geom.d == 2) pin v_z = 0 by equality so no agent can dodge
%   in z. The initial guess is the straight line plus a tangential bump of
%   amplitude opts.bump with a common rotation sense: every circle layout
%   is degenerate at the straight line (all agents meet at the centre at
%   T/2), and the mirror image of any solution is a solution.
%
%  ASSUMPTIONS / NOTES:
%
% • The defect constraints are added FIRST so their multipliers are
%   opti.lam_g(1:3MN), read by row range (never opti.dual: canonicalised
%   orientation, see memory casadi-optidual-sign-bug).
% • The barrier expression is a second copy of barrier_shape's formulas
%   (CasADi needs if_else); out.Unodes lets a test pin it to barrier_terms.
% • A sampled path constraint holds between samples only to O(h^2); the
%   clearance gate in verify_drones_direct reports the fine-grid minimum.
%
%% Inputs:
%
%  x0, xT                   [3 x M]                 From drone_layouts
%
%  geom                     struct                  From drone_layouts
%
%  opts                     struct (optional)       .N [100] .T [10]
%                                                   .mode ['hard'|'barrier']
%                                                   .sense [+1] .bump
%                                                   [0.05] .pairs
%                                                   [pair_list(M)]
%                                                   .midpoints [true]
%                                                   .maxIter [3000]
%                                                   .print [0] .X0 .V0
%                                                   (warm start)
%
%% Outputs:
%
%  out                      struct                  .X .V [3M x N+1] .tN
%                                                   .J .JLB .success
%                                                   .ipoptStatus
%                                                   .iterations .maxDefect
%                                                   .lamDefect [3M x N]
%                                                   .minHnodes .Unodes
%                                                   .mode .N .T .geom
%                                                   .pairs .wall
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the paper's M = 4 circle, hard clearance:
     [a, b, g] = drone_layouts('paper', 4);
     o = solve_drones_direct(a, b, g);
     fprintf('M = 4: J = %.6f (J_LB %.4f, OCFormer 0.2003), %s\n', o.J, o.JLB, o.ipoptStatus);
     return
end
if nargin < 4, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
casadi_path();
import casadi.*

d = @(f, v) fieldd(opts, f, v);
N      = d('N', 100);
T      = d('T', 10);
mode   = lower(d('mode', 'hard'));
sense  = d('sense', 1);
bump   = d('bump', 0.05);
M      = geom.M;
pairs  = d('pairs', pair_list(M));
midpts = d('midpoints', true);
maxIt  = d('maxIter', 3000);
prnt   = d('print', 0);
if ~(bump > 0)
    error('solve_drones_direct:bump', 'the tangential bump must be positive: the straight line is degenerate');
end
if ~any(strcmp(mode, {'hard', 'barrier'}))
    error('solve_drones_direct:mode', 'mode must be ''hard'' or ''barrier''');
end

nX = 3*M;  nN = N + 1;  h = T/N;
tN = linspace(0, T, nN);
P  = size(pairs, 2);
t0 = tic;

%% Variables and the defects (FIRST, for the multiplier rows)
opti = Opti();
X = opti.variable(nX, nN);
V = opti.variable(nX, nN);
opti.subject_to(X(:, 2:end) - X(:, 1:end-1) - (h/2)*(V(:, 1:end-1) + V(:, 2:end)) == 0);
opti.subject_to(X(:, 1) == x0(:));
opti.subject_to(X(:, end) == xT(:));
if geom.d == 2
    opti.subject_to(V(3:3:end, :) == 0);
end

%% Cost: trapezoid of 1/2 |v|^2
L = 0.5*sum1(V.^2);                           % [1 x nN]
J = h*sum2(L) - (h/2)*(L(1) + L(end));

%% Pair geometry at the nodes and (optionally) the Hermite midpoints
ia = reshape(3*(pairs(1,:) - 1) + [1; 2; 3], [], 1);
ib = reshape(3*(pairs(2,:) - 1) + [1; 2; 3], [], 1);
Xmid = (X(:, 1:end-1) + X(:, 2:end))/2 + (h/8)*(V(:, 1:end-1) - V(:, 2:end));
if midpts, Xs = [X, Xmid]; else, Xs = X; end
Unode = [];
if P > 0
    Dp = Xs(ia, :) - Xs(ib, :);                                        % [3P x ns]
    sq = Dp(1:3:end, :).^2 + Dp(2:3:end, :).^2 + Dp(3:3:end, :).^2;   % [P x ns]
    switch mode
        case 'hard'
            opti.subject_to(sq(:) >= (2*geom.r)^2);
        case 'barrier'
            hh = sqrt(sq(:, 1:nN)) - 2*geom.r;                        % nodes only
            Ub = barrier_casadi(hh, geom.eps, geom.l);
            Unode = sum1(Ub);                                          % [1 x nN]
    end
end
if ~isempty(geom.ball)
    Db = Xs - repmat(geom.ball.c, M, size(Xs, 2));
    sqb = Db(1:3:end, :).^2 + Db(2:3:end, :).^2 + Db(3:3:end, :).^2;  % [M x ns]
    switch mode
        case 'hard'
            opti.subject_to(sqb(:) >= geom.ball.R^2);
        case 'barrier'
            hb = sqrt(sqb(:, 1:nN)) - geom.ball.R;
            Ubb = sum1(barrier_casadi(hb, geom.eps, geom.l));
            if isempty(Unode), Unode = Ubb; else, Unode = Unode + Ubb; end
    end
end
if strcmp(mode, 'barrier') && ~isempty(Unode)
    J = J + h*sum2(Unode) - (h/2)*(Unode(1) + Unode(end));
end
opti.minimize(J);

%% Initial guess: straight lines plus the tangential bump
X0 = d('X0', []);  V0 = d('V0', []);
if isempty(X0) || isempty(V0)
    [X0, V0] = bumped_lines(x0, xT, tN, bump, sense);
end
opti.set_initial(X, X0);
opti.set_initial(V, V0);

%% Solve
ip = struct('print_level', prnt, 'max_iter', maxIt, 'tol', 1e-10, ...
            'acceptable_tol', 1e-8, 'linear_solver', 'mumps');
opti.solver('ipopt', struct('print_time', false), ip);
try
    sol = opti.solve();  solved = true;
catch
    sol = opti.debug;    solved = false;
end

%% Report, with the defects re-evaluated from the returned numbers
out = struct();
out.X  = full(sol.value(X));
out.V  = full(sol.value(V));
out.tN = tN;
out.J  = full(sol.value(J));
out.JLB = lower_bound(x0, xT, T);
out.ipoptStatus = opti.stats().return_status;
out.iterations  = opti.stats().iter_count;
out.success = solved && any(strcmp(out.ipoptStatus, {'Solve_Succeeded', 'Solved_To_Acceptable_Level'}));
def = out.X(:, 2:end) - out.X(:, 1:end-1) - (h/2)*(out.V(:, 1:end-1) + out.V(:, 2:end));
out.maxDefect = max(abs(def(:)));
lamAll = full(sol.value(opti.lam_g));
out.lamDefect = reshape(lamAll(1:nX*N), nX, N);
if P > 0
    hn = clearance_terms(out.X, pairs, geom);
    out.minHnodes = min(hn(:));
else
    out.minHnodes = Inf;
end
if isempty(Unode), out.Unodes = []; else, out.Unodes = full(sol.value(Unode)); end
out.mode = mode;  out.N = N;  out.T = T;  out.geom = geom;  out.pairs = pairs;
out.sense = sense;  out.bump = bump;  out.midpoints = midpts;
out.wall = toc(t0);
end

% ------------------------------------------------------------------------
function U = barrier_casadi(hh, epsB, l)
%% Purpose:
%
%   The paper's C1 log barrier as a CasADi expression (a copy of
%   barrier_shape's two formulas; barrier_shape itself branches on
%   real(h), which has no symbolic meaning). Pinned to barrier_terms by
%   test_direct_small through out.Unodes.
%
import casadi.*
U = if_else(hh > l, -epsB*log(hh), -epsB*log(l) + (epsB/2)*(((hh - 2*l)/l).^2 - 1));
end

function [X0, V0] = bumped_lines(x0, xT, tN, bump, sense)
%% Purpose:
%
%   Straight lines at constant speed plus a tangential bump
%   A sin(pi t/T) n_perp with a common rotation sense about e3, so the
%   agents start swirling the same way instead of meeting at the centre.
%
M  = size(x0, 2);  T = tN(end);  s = tN/T;
X0 = zeros(3*M, numel(tN));  V0 = X0;
for a = 1:M
    dv = xT(:, a) - x0(:, a);
    L = sqrt(sum(dv.^2));
    if L > 0, dir = dv/L; else, dir = [1; 0; 0]; end
    perp = [-dir(2); dir(1); 0];                       % e3 x dir
    if sqrt(sum(perp.^2)) < 1e-6
        perp = [0; -dir(3); dir(2)];                   % e1 x dir, for a z-aligned line
    end
    perp = sense*perp/sqrt(sum(perp.^2));
    ia = 3*(a-1) + (1:3);
    X0(ia, :) = x0(:, a) + dv*s + bump*perp*sin(pi*s);
    V0(ia, :) = dv/T + bump*perp*(pi/T)*cos(pi*s);
end
end

function v = fieldd(s, f, v)
%% Purpose:
%
%   s.(f) if present and non-empty, else the default v.
%
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
```

- [ ] **Step 4: Run the oracle test**

Run: `matlab -batch "cd('optimal_control_examples/drone/route1_direct/tests'); test_direct_oracle"`
Expected: all `PASS`. The N = 400 solve may take a minute. If the N = 100 gap exceeds 5e-4, print `out100.minHnodes` and the fine-grid minimum (Task 2's `hermite_resample` check) before touching the bar: a gap from the sagitta is of order 2.5e-4 and sits below J*; a gap above J* means the swirl sense or the midpoint term is wrong.

- [ ] **Step 5: Commit**

```bash
cd ~/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/route1_direct/solve_drones_direct.m optimal_control_examples/drone/route1_direct/casadi_path.m optimal_control_examples/drone/route1_direct/tests/test_direct_oracle.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone route1: solve_drones_direct (CasADi/IPOPT trapezoid, hard clearance at nodes + Hermite midpoints or the paper's barrier, 2D/3D ball), with the two-agent closed-form oracle test (5e-4 at N = 100, 5e-5 at N = 400, converging from below)

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 2: `verify_drones_direct`

**Files:**
- Create: `optimal_control_examples/drone/route1_direct/verify_drones_direct.m`
- Test: `optimal_control_examples/drone/route1_direct/tests/test_direct_small.m`

**Interfaces:**
- Consumes: `solve_drones_direct` out struct, `hermite_resample`, `clearance_terms`, `barrier_terms`.
- Produces: `rep = verify_drones_direct(out, opts)` with `opts.mirror [true] .mesh [false] .nSub [10] .tolG1 [1e-2 * 2r] .tolDefect [1e-9] .tolMesh [1e-3] .tolMirror [1e-6]`; `rep` fields `.ok .why .G1 (minHfine, minHballFine) .G2 (gap = J - JLB) .G3 (status, maxDefect) .G4 (J2N, relChange, ran) .G5 (Jmirror, relDiff, ran, note) .notes`.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_direct_small()
%% Purpose:
%
%   The M = 4 paper circle through the solver and the gates:
%     1. hard mode: verify_drones_direct says ok, every gate's number is in
%        range, J is above J_LB = 0.200 and BELOW the paper's OCFormer
%        0.2003 (a hard-constraint optimum beats a barrier-biased one);
%     2. the mirror solve (opposite swirl) agrees in cost to 1e-6;
%     3. barrier mode: solves; out.Unodes equals barrier_terms at the
%        solution nodes to 1e-12 (pins the CasADi copy of the barrier);
%        its fine-grid clearance is above l; its cost is above the hard
%        optimum (the barrier can only add);
%     4. a failed solve (max_iter = 2) is NOT verified, and .why names the
%        IPOPT status;
%     5. the sphere layout (not mirror-symmetric) skips G5 with a note.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));

[x0, xT, geom] = drone_layouts('paper', 4);
oH = solve_drones_direct(x0, xT, geom, struct('mode', 'hard'));
rep = verify_drones_direct(oH, struct('mirror', true, 'mesh', false));
ok = chk(ok, rep.ok, sprintf('hard M = 4: verify ok (%s)', tern(isempty(rep.why), 'no reason recorded', rep.why)));
ok = chk(ok, oH.J > 0.200 && oH.J < 0.2003, sprintf('hard M = 4: J = %.6f in (0.200, 0.2003): above J_LB, below OCFormer', oH.J));
ok = chk(ok, rep.G1.minHfine >= -1e-2*2*geom.r, sprintf('G1 fine-grid clearance %.2e (bar %.1e)', rep.G1.minHfine, -1e-2*2*geom.r));
ok = chk(ok, rep.G3.maxDefect < 1e-9, sprintf('G3 defects %.1e', rep.G3.maxDefect));
ok = chk(ok, rep.G5.ran && rep.G5.relDiff < 1e-6, sprintf('G5 mirror: J = %.6f vs %.6f (%.1e rel)', oH.J, rep.G5.Jmirror, rep.G5.relDiff));

oB = solve_drones_direct(x0, xT, geom, struct('mode', 'barrier'));
ok = chk(ok, oB.success, sprintf('barrier M = 4: %s', oB.ipoptStatus));
Ucheck = zeros(1, oB.N + 1);
for k = 1:oB.N + 1
    Ucheck(k) = barrier_terms(oB.X(:, k), oB.pairs, geom);
end
ok = chk(ok, max(abs(Ucheck - oB.Unodes)) < 1e-12, sprintf('CasADi barrier = barrier_terms at the nodes (%.1e)', max(abs(Ucheck - oB.Unodes))));
repB = verify_drones_direct(oB, struct('mirror', false));
ok = chk(ok, repB.ok, sprintf('barrier M = 4: verify ok (%s)', tern(isempty(repB.why), '', repB.why)));
ok = chk(ok, repB.G1.minHfine > geom.l, sprintf('barrier: fine-grid clearance %.2e > l = %.0e', repB.G1.minHfine, geom.l));
ok = chk(ok, oB.J > oH.J, sprintf('barrier cost %.6f > hard cost %.6f', oB.J, oH.J));

oF = solve_drones_direct(x0, xT, geom, struct('mode', 'hard', 'maxIter', 2));
repF = verify_drones_direct(oF, struct('mirror', false));
ok = chk(ok, ~oF.success && ~repF.ok && contains(repF.why, oF.ipoptStatus), ...
         sprintf('a failed solve is not verified: %s', repF.why));

[s0, sT, gS] = drone_layouts('sphere', 6);
oS = solve_drones_direct(s0, sT, gS, struct('mode', 'hard', 'N', 60));
repS = verify_drones_direct(oS, struct('mirror', true));
ok = chk(ok, oS.success && repS.ok && ~repS.G5.ran && ~isempty(repS.G5.note), ...
         sprintf('sphere M = 6: ok, mirror skipped (%s), ball clearance %.2e', repS.G5.note, repS.G1.minHballFine));

if ok, fprintf('TEST_DIRECT_SMALL: ALL PASS\n'); else, fprintf('TEST_DIRECT_SMALL: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function s = tern(c, a, b)
%% Purpose:
%
%   a if c, else b.
%
if c, s = a; else, s = b; end
end

function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 2: Run it to verify it fails**

Expected: `Unrecognized function or variable 'verify_drones_direct'`.

- [ ] **Step 3: Implement**

```matlab
function rep = verify_drones_direct(out, opts)
%% Purpose:
%
%   The gates on one direct solve, with its own verdict:
%     G1 clearance   min over pairs (and the ball) of the gap on a nSub-times
%                    finer cubic-Hermite reconstruction from (X, V) -- a
%                    nodal constraint says nothing between nodes, and the
%                    interpolant cuts a constraint arc's corner by up to
%                    s^2/(16 r), so the bar is mesh-aware (-1e-2 * 2r) and
%                    the number is REPORTED, not hidden
%     G2 bound       J >= J_LB
%     G3 solver      IPOPT status ok and the re-evaluated defects < 1e-9
%     G4 mesh        (opts.mesh) re-solve at 2N warm-started from the
%                    interpolant; relative cost change
%     G5 mirror      (opts.mirror) re-solve with the opposite swirl; costs
%                    agree -- only meaningful on a mirror-symmetric layout
%                    (the paper's circles/rings and headon); skipped with
%                    a note otherwise
%
%% Inputs:
%
%  out                      struct                  From solve_drones_direct
%
%  opts                     struct (optional)       .mirror [true] .mesh
%                                                   [false] .nSub [10]
%                                                   .tolG1 [1e-2 * 2r]
%                                                   .tolDefect [1e-9]
%                                                   .tolMesh [1e-3]
%                                                   .tolMirror [1e-6]
%
%% Outputs:
%
%  rep                      struct                  .ok .why .G1 .G2 .G3
%                                                   .G4 .G5 .notes
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
d = @(f, v) fieldd(opts, f, v);
geom = out.geom;
doMirror = d('mirror', true);
doMesh   = d('mesh', false);
nSub     = d('nSub', 10);
tolG1    = d('tolG1', 1e-2*2*geom.r);
tolDef   = d('tolDefect', 1e-9);
tolMesh  = d('tolMesh', 1e-3);
tolMir   = d('tolMirror', 1e-6);
rep = struct('ok', false, 'why', '', 'notes', {{}});

%% G3 first: a failed solve has no gates worth reading
rep.G3 = struct('status', out.ipoptStatus, 'maxDefect', out.maxDefect, ...
                'ok', out.success && out.maxDefect < tolDef);

%% G1: the fine grid
[Xf, ~] = hermite_resample(out.X, out.V, out.tN, nSub);
if size(out.pairs, 2) > 0
    [hf, ~, ~, hbf] = clearance_terms(Xf, out.pairs, geom);
    minH = min(hf(:));
else
    [~, ~, ~, hbf] = clearance_terms(Xf, out.pairs, geom);
    minH = Inf;
end
if isempty(hbf), minHb = Inf; else, minHb = min(hbf(:)); end
rep.G1 = struct('minHfine', minH, 'minHballFine', minHb, 'nSub', nSub, 'bar', -tolG1);
switch out.mode
    case 'hard',    rep.G1.ok = minH >= -tolG1 && minHb >= -tolG1;
    case 'barrier', rep.G1.ok = true;   % reported, not gated: the barrier is the constraint
end

%% G2
rep.G2 = struct('gap', out.J - out.JLB, 'ok', out.J >= out.JLB - 1e-12);

%% G4: mesh
rep.G4 = struct('ran', false, 'J2N', NaN, 'relChange', NaN, 'ok', true);
if doMesh && out.success
    [X2, V2] = hermite_resample(out.X, out.V, out.tN, 2);
    o2 = solve_drones_direct(reshape(out.X(:,1), 3, geom.M), reshape(out.X(:,end), 3, geom.M), geom, ...
                             struct('N', 2*out.N, 'T', out.T, 'mode', out.mode, 'sense', out.sense, ...
                                    'bump', out.bump, 'pairs', out.pairs, 'midpoints', out.midpoints, ...
                                    'X0', X2, 'V0', V2));
    rep.G4.ran = true;  rep.G4.J2N = o2.J;
    rep.G4.relChange = abs(o2.J - out.J)/abs(out.J);
    rep.G4.ok = o2.success && rep.G4.relChange < tolMesh;
end

%% G5: mirror
rep.G5 = struct('ran', false, 'Jmirror', NaN, 'relDiff', NaN, 'ok', true, 'note', '');
symmetric = any(strcmp(geom.which, {'paper', 'circle', 'rings', 'headon'}));
if doMirror && out.success
    if symmetric
        om = solve_drones_direct(reshape(out.X(:,1), 3, geom.M), reshape(out.X(:,end), 3, geom.M), geom, ...
                                 struct('N', out.N, 'T', out.T, 'mode', out.mode, 'sense', -out.sense, ...
                                        'bump', out.bump, 'pairs', out.pairs, 'midpoints', out.midpoints));
        rep.G5.ran = true;  rep.G5.Jmirror = om.J;
        rep.G5.relDiff = abs(om.J - out.J)/abs(out.J);
        rep.G5.ok = om.success && rep.G5.relDiff < tolMir;
    else
        rep.G5.note = sprintf('layout ''%s'' is not mirror-symmetric: G5 skipped', geom.which);
        rep.notes{end+1} = rep.G5.note;
    end
end

%% Verdict
if ~rep.G3.ok
    rep.why = sprintf('G3: IPOPT %s, max defect %.2e', out.ipoptStatus, out.maxDefect);
elseif ~rep.G1.ok
    rep.why = sprintf('G1: fine-grid clearance %.2e below %.1e', min(minH, minHb), -tolG1);
elseif ~rep.G2.ok
    rep.why = sprintf('G2: J = %.6f below J_LB = %.6f', out.J, out.JLB);
elseif ~rep.G4.ok
    rep.why = sprintf('G4: cost changed %.2e relative at 2N', rep.G4.relChange);
elseif ~rep.G5.ok
    rep.why = sprintf('G5: mirror cost differs by %.2e relative', rep.G5.relDiff);
end
rep.ok = isempty(rep.why);
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v)
%% Purpose:
%
%   s.(f) if present and non-empty, else the default v.
%
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
```

- [ ] **Step 4: Run the test**

Expected: all `PASS`. If the hard M = 4 cost is not below 0.2003, inspect the fine-grid clearance: a cost above the paper's with a clearance well above zero means the solve stopped at a saddle (re-check the bump sense); a cost below J_LB means a constraint is missing.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/route1_direct/verify_drones_direct.m optimal_control_examples/drone/route1_direct/tests/test_direct_small.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone route1: verify_drones_direct (fine-grid clearance, bound, solver, mesh, mirror gates with its own verdict); M = 4 circle test pins the CasADi barrier to barrier_terms and beats the paper's OCFormer cost

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 3: `drone_table1`, fixtures, batch launcher

**Files:**
- Create: `optimal_control_examples/drone/route1_direct/drone_table1.m`
- Create: `optimal_control_examples/drone/route1_direct/batch/run_table1.sh`
- Create: `optimal_control_examples/drone/route1_direct/batch/job_table1.m`
- Create (generated, force-added): `optimal_control_examples/drone/route1_direct/data/drones_direct_M4.mat`, `drones_direct_M8.mat`
- Create (generated, tracked): `optimal_control_examples/drone/route1_direct/results/table1.csv`
- Test: `optimal_control_examples/drone/route1_direct/tests/test_table1_fixture.m`

**Interfaces:**
- Produces: `tbl = drone_table1(opts)` with `opts.Ms [4 8 16 32] .modes {'hard','barrier'} .mesh [true for M <= 16] .saveFixtures [true] .csv [results/table1.csv] .logFile ['']`; returns a MATLAB table with columns `M, mode, N, J_LB, J, gapPct, OCFormer, SympOCNet, minHfine, maxDefect, meshRel, mirrorRel, iters, wallSec, status`. Fixture file fields: `out` (hard), `outB` (barrier), `rep`, `repB`, `geom`, `x0`, `xT`, `paper` (struct of the paper's row).

- [ ] **Step 1: Write the failing fixture test**

```matlab
function ok = test_table1_fixture()
%% Purpose:
%
%   The committed M = 4 and M = 8 fixtures are what they claim: the stored
%   X, V re-integrate to the stored J (trapezoid, to 1e-12), the gates
%   re-derived from the stored X, V agree with the stored report, both
%   modes are present with the multipliers Route 2 needs, and the hard
%   costs beat the paper's. An edited fixture cannot sit green.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
paperJ = struct('M4', 0.2003, 'M8', 0.4036);

for M = [4 8]
    f = fullfile(here, 'data', sprintf('drones_direct_M%d.mat', M));
    ok = chk(ok, isfile(f), sprintf('fixture M = %d exists', M));
    if ~isfile(f), continue, end
    S = load(f);
    for which = {'out', 'outB'}
        o = S.(which{1});
        h = o.T/o.N;
        L = 0.5*sum(o.V.^2, 1);
        Jre = h*sum(L) - (h/2)*(L(1) + L(end));
        if strcmp(o.mode, 'barrier')
            Jre = Jre + h*sum(o.Unodes) - (h/2)*(o.Unodes(1) + o.Unodes(end));
        end
        ok = chk(ok, abs(Jre - o.J) < 1e-12*max(1, abs(o.J)), sprintf('M = %d %s: stored J re-integrates (%.1e)', M, o.mode, abs(Jre - o.J)));
        ok = chk(ok, isequal(size(o.lamDefect), [3*M o.N]) && all(isfinite(o.lamDefect(:))), sprintf('M = %d %s: defect multipliers present', M, o.mode));
        rep = verify_drones_direct(o, struct('mirror', false));
        ok = chk(ok, rep.ok, sprintf('M = %d %s: gates re-derived ok (%s)', M, o.mode, rep.why));
    end
    ok = chk(ok, S.out.J < paperJ.(sprintf('M%d', M)) && S.out.J > S.out.JLB, ...
             sprintf('M = %d: hard J = %.6f beats OCFormer %.4f, above J_LB %.4f', M, S.out.J, paperJ.(sprintf('M%d', M)), S.out.JLB));
    ok = chk(ok, S.rep.ok && S.rep.G5.ran && S.rep.G4.ran, sprintf('M = %d: stored report ran mirror and mesh and is ok', M));
end

csvFile = fullfile(here, 'results', 'table1.csv');
ok = chk(ok, isfile(csvFile), 'results/table1.csv exists');
if isfile(csvFile)
    tb = readtable(csvFile);
    ok = chk(ok, all(ismember([4 8 16 32], tb.M)), 'table1.csv has M = 4, 8, 16, 32');
    hard = tb(strcmp(tb.mode, 'hard'), :);
    ok = chk(ok, all(hard.J < hard.OCFormer) && all(hard.J >= hard.J_LB), 'every hard row beats OCFormer and respects J_LB');
end

if ok, fprintf('TEST_TABLE1_FIXTURE: ALL PASS\n'); else, fprintf('TEST_TABLE1_FIXTURE: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 2: Run it to verify it fails**

Expected: `FAIL  fixture M = 4 exists` lines and `results/table1.csv exists` FAIL.

- [ ] **Step 3: Write `drone_table1`**

```matlab
function tbl = drone_table1(opts)
%% Purpose:
%
%   Reproduce Table 1 of Zhang, Liu, Darbon & Karniadakis (arXiv:
%   2610.30388) with certified classical solves: for each M, the hard-
%   clearance optimum (the reference) and the paper's barrier form, each
%   through verify_drones_direct, written to results/table1.csv beside the
%   paper's OCFormer and SympOCNet costs, with the M = 4 and M = 8 fixtures
%   saved for Route 2. The paper's own numbers (its Table 1):
%
%       M     J_LB    OCFormer  SympOCNet
%       4     0.200   0.2003    0.2001
%       8     0.400   0.4036    0.4046
%       16    0.800   0.8307    0.8345
%       32    1.156   1.2531    1.4413
%       64    2.297   2.448     --
%       128   4.704   5.1424    13.3 (did not converge)
%
%   M = 64 and 128 are STRETCH: run them from batch/run_table1.sh with a
%   log file, never from the test suite.
%
%% Inputs:
%
%  opts                     struct (optional)       .Ms [4 8 16 32]
%                                                   .modes {'hard',
%                                                   'barrier'} .mesh
%                                                   [true for M <= 16]
%                                                   .saveFixtures [true]
%                                                   .csv [results/
%                                                   table1.csv] .logFile
%                                                   [''] .N [100]
%
%% Outputs:
%
%  tbl                      table                   One row per (M, mode)
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
d = @(f, v) fieldd(opts, f, v);
Ms     = d('Ms', [4 8 16 32]);
modes  = d('modes', {'hard', 'barrier'});
doMesh = d('mesh', []);
saveFx = d('saveFixtures', true);
csvOut = d('csv', fullfile(here, 'results', 'table1.csv'));
logF   = d('logFile', '');
N      = d('N', 100);
lg = @(varargin) logmsg(logF, sprintf(varargin{:}));

paper = containers.Map('KeyType', 'double', 'ValueType', 'any');
paper(4)   = [0.2003, 0.2001];   paper(8)   = [0.4036, 0.4046];
paper(16)  = [0.8307, 0.8345];   paper(32)  = [1.2531, 1.4413];
paper(64)  = [2.448,  NaN];      paper(128) = [5.1424, 13.3];

rows = {};
for M = Ms
    [x0, xT, geom] = drone_layouts('paper', M);
    mesh = doMesh;  if isempty(mesh), mesh = M <= 16; end
    store = struct();
    for mm = 1:numel(modes)
        mode = modes{mm};
        lg('M = %d, mode %s: solving N = %d', M, mode, N);
        o = solve_drones_direct(x0, xT, geom, struct('mode', mode, 'N', N));
        rep = verify_drones_direct(o, struct('mirror', true, 'mesh', mesh));
        pp = paper(M);
        rows(end+1, :) = {M, mode, N, o.JLB, o.J, 100*(o.J - o.JLB)/o.JLB, pp(1), pp(2), ...
                          rep.G1.minHfine, o.maxDefect, rep.G4.relChange, rep.G5.relDiff, ...
                          o.iterations, o.wall, o.ipoptStatus, rep.ok};
        lg('M = %d, mode %s: J = %.6f (J_LB %.4f, OCFormer %.4f), %s, gates %s %s', ...
           M, mode, o.J, o.JLB, pp(1), o.ipoptStatus, mat2str(rep.ok), rep.why);
        if strcmp(mode, 'hard'), store.out = o;  store.rep = rep;
        else,                    store.outB = o; store.repB = rep;
        end
    end
    if saveFx && any(M == [4 8]) && isfield(store, 'out') && isfield(store, 'outB')
        out = store.out;  outB = store.outB;  rep = store.rep;  repB = store.repB;
        paperRow = struct('OCFormer', paper(M)*[1; 0], 'SympOCNet', paper(M)*[0; 1]);
        if ~isfolder(fullfile(here, 'data')), mkdir(fullfile(here, 'data')); end
        save(fullfile(here, 'data', sprintf('drones_direct_M%d.mat', M)), ...
             'out', 'outB', 'rep', 'repB', 'geom', 'x0', 'xT', 'paperRow');
        lg('M = %d: fixture saved', M);
    end
end
tbl = cell2table(rows, 'VariableNames', {'M', 'mode', 'N', 'J_LB', 'J', 'gapPct', 'OCFormer', 'SympOCNet', ...
                 'minHfine', 'maxDefect', 'meshRel', 'mirrorRel', 'iters', 'wallSec', 'status', 'gatesOk'});
if ~isfolder(fileparts(csvOut)), mkdir(fileparts(csvOut)); end
writetable(tbl, csvOut);
lg('wrote %s', csvOut);
disp(tbl);
end

% ------------------------------------------------------------------------
function logmsg(f, s)
%% Purpose:
%
%   Append one line to the log file, or print it when no file is set.
%
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a');  fprintf(fid, '%s\n', s);  fclose(fid);
end
end

function v = fieldd(s, f, v)
%% Purpose:
%
%   s.(f) if present and non-empty, else the default v.
%
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
```

- [ ] **Step 4: Write the batch launcher (campaign discipline)**

`batch/job_table1.m`:

```matlab
%% JOB_TABLE1  Batch driver for drone_table1 (the full table, incl. stretch).
%
%  Launched by run_table1.sh in a clean MATLAB; logs to a file beside the
%  results and writes a JOB DONE marker. Edit the Ms line for the stretch
%  rows (64, 128 take a long time: thousands of pair constraints).
%
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
logFile = fullfile(fileparts(here), 'results', 'table1_batch.log');
try
    drone_table1(struct('Ms', [4 8 16 32], 'logFile', logFile));
catch err
    fid = fopen(logFile, 'a');
    fprintf(fid, 'JOB FAILED: %s\n', getReport(err, 'extended', 'hyperlinks', 'off'));
    fclose(fid);
end
fid = fopen(logFile, 'a');  fprintf(fid, 'JOB DONE\n');  fclose(fid);
```

`batch/run_table1.sh`:

```bash
#!/bin/zsh
# Launch drone_table1 in a clean MATLAB, detached, logging to results/.
# Usage: ./run_table1.sh        (then tail -f ../results/table1_batch.log)
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$HERE/../results"
nohup /Applications/MATLAB_R2026a.app/bin/matlab -batch "run('$HERE/job_table1.m')" \
    > "$HERE/../results/table1_batch.stdout" 2>&1 &
echo "launched, pid $!; log: $HERE/../results/table1_batch.log"
```

Make it executable: `chmod +x batch/run_table1.sh`. Add `results/*.log`, `results/*.stdout` to the folder's `.gitignore` (create `optimal_control_examples/drone/route1_direct/.gitignore` with those two lines); `results/table1.csv` stays tracked.

- [ ] **Step 5: Generate the table and fixtures**

Run: `cd optimal_control_examples/drone/route1_direct/batch && ./run_table1.sh`, then wait for `JOB DONE` in `results/table1_batch.log` (M = 32 with the mirror solve is the long row; expect tens of minutes). Read the log, not stdout.

Expected: `table1.csv` with 8 rows, every `gatesOk` true, every hard `J` below the paper's OCFormer column and above `J_LB`; `data/drones_direct_M4.mat` and `M8.mat` present.

- [ ] **Step 6: Run the fixture test**

Expected: all `PASS`.

- [ ] **Step 7: Commit (force-add the fixtures)**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/route1_direct/drone_table1.m optimal_control_examples/drone/route1_direct/batch optimal_control_examples/drone/route1_direct/.gitignore optimal_control_examples/drone/route1_direct/results/table1.csv optimal_control_examples/drone/route1_direct/tests/test_table1_fixture.m
/Library/Developer/CommandLineTools/usr/bin/git add -f optimal_control_examples/drone/route1_direct/data/drones_direct_M4.mat optimal_control_examples/drone/route1_direct/data/drones_direct_M8.mat
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone route1: drone_table1 reproduces the paper's Table 1 (M = 4..32, hard + barrier, all gates) into results/table1.csv; M = 4 and 8 fixtures with defect multipliers for Route 2; batch launcher under the campaign discipline

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 4: `plot_drones` and `movie_drones`

**Files:**
- Create: `optimal_control_examples/drone/route1_direct/plot_drones.m`
- Create: `optimal_control_examples/drone/route1_direct/movie_drones.m`
- Test: `optimal_control_examples/drone/route1_direct/tests/test_plot_drones.m`

**Interfaces:**
- Produces: `fig = plot_drones(out, opts)`: a 1×4 panel of the configuration at `t = 0, T/3, 2T/3, T` with the paths as faint grey curves and agents as coloured circles of radius `r` (the paper's Figure 4 style), `opts.visible [true] .file ['']` (PNG export via `exportgraphics` when set); `movie_drones(out, file, opts)`: MP4 at 1280×720 (a ÷16 frame size, memory `matlab-movie-diagonal-streaks`), `opts.fps [30] .nFrames [300]`, through `VideoWriter('MPEG-4')`.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_plot_drones()
%% Purpose:
%
%   plot_drones renders the M = 4 fixture to a PNG headlessly (no visible
%   figure), the file is non-empty, and the figure has four axes;
%   movie_drones writes a short MP4 whose frame size is 1280 x 720 (a
%   multiple of 16: anything else shears into diagonal streaks).
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
S = load(fullfile(here, 'data', 'drones_direct_M4.mat'));
png = [tempname, '.png'];  mp4 = [tempname, '.mp4'];
cleanup = onCleanup(@() cellfun(@(f) delete(f), {png, mp4}));

fig = plot_drones(S.out, struct('visible', false, 'file', png));
ax = findobj(fig, 'Type', 'axes');
ok = chk(ok, numel(ax) == 4, sprintf('four panels (%d axes)', numel(ax)));
ok = chk(ok, isfile(png) && dir(png).bytes > 1000, 'PNG written');
close(fig);

movie_drones(S.out, mp4, struct('nFrames', 12, 'fps', 6));
ok = chk(ok, isfile(mp4), 'MP4 written');
if isfile(mp4)
    vr = VideoReader(mp4);
    ok = chk(ok, vr.Width == 1280 && vr.Height == 720, sprintf('frame %dx%d (multiple of 16)', vr.Width, vr.Height));
end

if ok, fprintf('TEST_PLOT_DRONES: ALL PASS\n'); else, fprintf('TEST_PLOT_DRONES: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 2: Run it to verify it fails**

Expected: `Unrecognized function or variable 'plot_drones'`.

- [ ] **Step 3: Implement**

```matlab
function fig = plot_drones(out, opts)
%% Purpose:
%
%   The paper's Figure 4 view of one solve: four panels at t = 0, T/3,
%   2T/3, T, paths as faint grey curves, agents as circles of radius r.
%   2D layouts are drawn in the plane; 3D ones in a 3D axes with the ball
%   when present.
%
%% Inputs:
%
%  out                      struct                  From solve_drones_direct
%
%  opts                     struct (optional)       .visible [true] .file
%                                                   [''] PNG to export
%
%% Outputs:
%
%  fig                      handle                  The figure
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the M = 4 fixture:
     here = fileparts(mfilename('fullpath'));
     S = load(fullfile(here, 'data', 'drones_direct_M4.mat'));
     plot_drones(S.out);
     return
end
if nargin < 2, opts = struct(); end
vis = true;  if isfield(opts, 'visible'), vis = opts.visible; end
file = '';   if isfield(opts, 'file'), file = opts.file; end
geom = out.geom;  M = geom.M;
[Xf, ~, tf] = hermite_resample(out.X, out.V, out.tN, 5);
snaps = [0, 1/3, 2/3, 1]*out.T;
cols = lines(M);
th = linspace(0, 2*pi, 40);
fig = figure('color', [1 1 1], 'Position', [100 100 1600 420], 'Visible', tern(vis, 'on', 'off'));
for s = 1:4
    ax = subplot(1, 4, s);  hold(ax, 'on');
    [~, k] = min(abs(tf - snaps(s)));
    for a = 1:M
        ia = 3*(a-1) + (1:3);
        if geom.d == 2
            plot(ax, Xf(ia(1), :), Xf(ia(2), :), '-', 'Color', [0.8 0.8 0.8]);
            plot(ax, Xf(ia(1), k) + geom.r*cos(th), Xf(ia(2), k) + geom.r*sin(th), '-', 'Color', cols(a, :), 'LineWidth', 1.2);
        else
            plot3(ax, Xf(ia(1), :), Xf(ia(2), :), Xf(ia(3), :), '-', 'Color', [0.8 0.8 0.8]);
            plot3(ax, Xf(ia(1), k), Xf(ia(2), k), Xf(ia(3), k), 'o', 'Color', cols(a, :), 'MarkerFaceColor', cols(a, :));
        end
    end
    if geom.d == 3 && ~isempty(geom.ball)
        [sx, sy, sz] = sphere(24);
        surf(ax, geom.ball.c(1) + geom.ball.R*sx, geom.ball.c(2) + geom.ball.R*sy, geom.ball.c(3) + geom.ball.R*sz, ...
             'FaceColor', [0.6 0.6 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        view(ax, 3);
    end
    axis(ax, 'equal');  grid(ax, 'on');
    title(ax, sprintf('t = %.2f', snaps(s)));
    if geom.d == 2, xlim(ax, [-0.55 0.55]); ylim(ax, [-0.55 0.55]); end
end
sgtitle(fig, sprintf('M = %d, %s, J = %.4f (J_{LB} %.4f)', M, out.mode, out.J, out.JLB));
if ~isempty(file), exportgraphics(fig, file, 'Resolution', 120); end
end

% ------------------------------------------------------------------------
function s = tern(c, a, b)
%% Purpose:
%
%   a if c, else b.
%
if c, s = a; else, s = b; end
end
```

```matlab
function movie_drones(out, file, opts)
%% Purpose:
%
%   An MP4 of one solve: agents as circles moving along faint paths. Every
%   frame is 1280 x 720 (a multiple of 16; H.264 shears any other size
%   into diagonal streaks -- memory matlab-movie-diagonal-streaks).
%
%% Inputs:
%
%  out                      struct                  From solve_drones_direct
%
%  file                     char                    Output .mp4
%
%  opts                     struct (optional)       .fps [30] .nFrames
%                                                   [300]
%
%% Outputs:
%
%  none (writes file)
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3, opts = struct(); end
fps = 30;  if isfield(opts, 'fps'), fps = opts.fps; end
nF = 300;  if isfield(opts, 'nFrames'), nF = opts.nFrames; end
geom = out.geom;  M = geom.M;
[Xf, ~, tf] = hermite_resample(out.X, out.V, out.tN, 5);
tt = linspace(0, out.T, nF);
cols = lines(M);  th = linspace(0, 2*pi, 40);
fig = figure('color', [1 1 1], 'Visible', 'off', 'Position', [100 100 1280 720]);
ax = axes(fig);  hold(ax, 'on');  axis(ax, 'equal');  grid(ax, 'on');
if geom.d == 2, xlim(ax, [-0.6 0.6]); ylim(ax, [-0.6 0.6]); else, view(ax, 3); end
for a = 1:M
    ia = 3*(a-1) + (1:3);
    if geom.d == 2, plot(ax, Xf(ia(1), :), Xf(ia(2), :), '-', 'Color', [0.85 0.85 0.85]);
    else, plot3(ax, Xf(ia(1), :), Xf(ia(2), :), Xf(ia(3), :), '-', 'Color', [0.85 0.85 0.85]);
    end
end
hA = gobjects(1, M);
for a = 1:M
    if geom.d == 2, hA(a) = plot(ax, NaN, NaN, '-', 'Color', cols(a, :), 'LineWidth', 1.5);
    else, hA(a) = plot3(ax, NaN, NaN, NaN, 'o', 'Color', cols(a, :), 'MarkerFaceColor', cols(a, :));
    end
end
vw = VideoWriter(file, 'MPEG-4');  vw.FrameRate = fps;  open(vw);
for f = 1:nF
    Xt = interp1(tf, Xf.', tt(f)).';
    for a = 1:M
        ia = 3*(a-1) + (1:3);
        if geom.d == 2, set(hA(a), 'XData', Xt(ia(1)) + geom.r*cos(th), 'YData', Xt(ia(2)) + geom.r*sin(th));
        else, set(hA(a), 'XData', Xt(ia(1)), 'YData', Xt(ia(2)), 'ZData', Xt(ia(3)));
        end
    end
    title(ax, sprintf('t = %.2f / %.0f', tt(f), out.T));
    fr = getframe(fig);
    im = imresize(fr.cdata, [720 1280]);
    writeVideo(vw, im);
end
close(vw);  close(fig);
end
```

- [ ] **Step 4: Run the test**

Expected: all `PASS`.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/route1_direct/plot_drones.m optimal_control_examples/drone/route1_direct/movie_drones.m optimal_control_examples/drone/route1_direct/tests/test_plot_drones.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone route1: plot_drones (the paper's four-snapshot view) and movie_drones (1280x720 MP4), tested headless

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 5: Runner registration and README

**Files:**
- Modify: `optimal_control_examples/drone/run_tests.m` (append the four Route 1 tests and the `route1_direct` paths)
- Modify: `optimal_control_examples/drone/README.md` (Stage 1 row "built", the Table 1 numbers)
- Create: `optimal_control_examples/drone/route1_direct/README.md`

- [ ] **Step 1: Register the tests**

In `drone/run_tests.m`, extend the `addpath` line with `fullfile(here, 'route1_direct'), fullfile(here, 'route1_direct', 'tests')` and append to `names`:

```matlab
         'test_direct_oracle'         % the two-agent closed form (G6), N = 100 and 400
         'test_direct_small'          % M = 4 circle: all gates, barrier pinned, failure path, sphere
         'test_table1_fixture'        % the committed fixtures and the CSV
         'test_plot_drones'           % headless plot + 1280x720 movie
```

- [ ] **Step 2: Write `route1_direct/README.md`**

```markdown
# route1_direct

The paper's multi-agent kinematic drone problem (Zhang, Liu, Darbon &
Karniadakis, arXiv:2610.30388, eq. 19) as a trapezoidal NLP in CasADi/IPOPT,
with hard pairwise clearance (the reference optimum) or the paper's own C1
log barrier (the smooth twin Route 2 shoots).

| file | what |
|---|---|
| `solve_drones_direct.m` | one solve: modes `hard` / `barrier`, 2D (v_z pinned) or 3D with a ball obstacle, clearance at nodes AND Hermite midpoints, straight-line-plus-swirl initial guess, defect multipliers by row range |
| `verify_drones_direct.m` | gates: G1 fine-grid clearance, G2 bound, G3 solver + defects, G4 mesh (2N), G5 mirror swirl; `.ok/.why` |
| `drone_table1.m` | the Table 1 loop -> `results/table1.csv`, fixtures `data/drones_direct_M{4,8}.mat` |
| `plot_drones.m`, `movie_drones.m` | the paper's four-snapshot figure; MP4 |
| `batch/run_table1.sh` | detached batch run with a log file (`JOB DONE` marker); the only way to run M = 64 / 128 |
| `tests/` | `test_direct_oracle` (closed-form head-on case), `test_direct_small`, `test_table1_fixture`, `test_plot_drones` |

Results: see `results/table1.csv` (every hard row sits above J_LB and below
the paper's OCFormer cost; the fine-grid clearance minimum is printed per
row: a sampled path constraint holds between samples only to O(h^2)).
```

Fill the `drone/README.md` Stage 1 row with "built" and paste the hard-mode `M, J_LB, J, OCFormer` columns of `table1.csv` as a small table.

- [ ] **Step 3: Run the folder runner and the tree runner**

Run: `matlab -batch "cd('optimal_control_examples/drone'); exit(~run_tests())"` and `matlab -batch "cd('optimal_control_examples'); exit(~run_all_tests())"`.
Expected: `12 of 12 passed`; tree runner exit 0.

- [ ] **Step 4: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/run_tests.m optimal_control_examples/drone/README.md optimal_control_examples/drone/route1_direct/README.md
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone route1: tests registered in the folder runner (12), READMEs with the Table 1 numbers

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```
