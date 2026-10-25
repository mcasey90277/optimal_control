# Chaplygin Sleigh and the Mechanics Library Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A user-friendly MATLAB script that integrates the Chaplygin sleigh from user-set geometry and initial velocity, verifies the solution against an independent symbolic builder and an exact solution, and writes diagnostic plots and a 1280 x 720 movie; the three generic pieces (equation builder, conservation check, planar-body animator) land in a new mechanics library for the next nonholonomic systems.

**Architecture:** `lagrangianHamiltonianMechanics/mechlib/+mech/` holds system-agnostic structure (Lagrange-d'Alembert with multipliers from a symbolic Lagrangian and constraint matrix, an energy/constraint diagnostic, a planar animator). `nonholonomic/ChaplyginSleigh/` holds the sleigh's parameters, its symbolic Lagrangian for the builder, the hand-derived reduced equations the script integrates, its exact solution and derived quantities, the plot and body-shape helpers, the front-door script and a LaTeX note. The hand-derived equations and the builder check each other (the physics-oracle rule); the exact tanh solution pins the integrator.

**Tech Stack:** MATLAB R2026a (shared session through the `mcp__matlab__*` tools), Symbolic Math Toolbox (builder only), ode113, VideoWriter MPEG-4, pdflatex.

**Spec:** `docs/superpowers/specs/2026-10-25-chaplygin-sleigh-design.md` (read it first; section 3 fixes every sign convention).

## Global Constraints

- Header quartet in every MATLAB function: `%% Purpose:` / `%% Inputs:` / `%% Outputs:` / `%% Revision History:` with `%  M. Casey                                                   (c) 10/25/2026` and `%  Copyright Coorbital Inc.`, then `%% ------------------------ Begin Code Sequence ---------------------------`.
- `if nargin == 0` self-demo in every function (not scripts, not tests).
- Never `i` or `j` as variables; never `norm` (write `sqrt(sum(x.^2))`); no `%#ok` pragmas; no emojis; ASCII READMEs; READMEs carry no author signature line.
- World frame x east, y north, theta counterclockwise from x; e1 = (cos th, sin th), e2 = (-sin th, cos th); blade point r_b = r_c - a e1; a > 0 = center of mass AHEAD of the blade.
- Reduced state the script integrates: y = [x_b; y_b; th; v; w] of the BLADE point; vd = a w^2, wd = -k v w, k = m a / (J + m a^2), lambda = m J v w / (J + m a^2).
- Movies exactly 1280 x 720, MPEG-4 quality 95; `results/` git-ignored (`*.png`, `*.mp4`, `*.log`).
- Tolerances: ode113 RelTol = AbsTol = 1e-12. Bars: D1 1e-9, D2 1e-8, D3 1e-9, D4 1e-9.
- Call library functions as `mech.<fn>` after `addpath('.../mechlib')`; the sleigh never imports `oc.*`; mechlib never imports anything from a system folder.
- Symbolic Math is used by `mech.lagrange_dalembert` only; `sleigh_rhs` stays hand-written.
- Subagents never commit. The controller commits after each clean review with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`, using `/Library/Developer/CommandLineTools/usr/bin/git`; the unrelated `missiles/` deletions in the working tree are never staged.
- MATLAB tests are function files: `cd` into the `tests` folder and call the test by name through `mcp__matlab__evaluate_matlab_code`; each returns `ok` and prints `PASS`/`FAIL` lines with the measured number and its bar.

## Review Focus

1. a = 0 with w0 ~= 0 (the knife edge): v and w must stay constant and the sleigh must run a circle with lambda = m v w; k = 0 must not produce NaN anywhere (t0 = atanh(.)/0). Pinned in Task 4's oracle test.
2. v0 = w0 = 0 (at rest): C = 0; the exact solution must be the constant zero, not NaN from 0/0, and plot_sleigh's ellipse must degenerate to a point without error. Pinned in Task 4 (quantities) and Task 5 (plot).
3. An unsorted or repeated time grid handed to the animator: interp1 would throw an unhelpful error; the animator refuses with `mech:animate_planar_body:t`. Pinned in Task 3.
4. A transposed solution ([N x 2n] instead of [2n x N]) handed to the conservation check: refused with `mech:conservation_check:shape` instead of silently reading velocities as positions. Pinned in Task 2.
5. Generalized forces that depend on velocity (damping): the builder must put Q on the right-hand side so energy decreases; a builder that silently dropped Q would still pass every conservative test. Pinned in Task 1.

---

### Task 1: mechlib scaffold, `mech.lagrange_dalembert`, the folder runner

**Files:**
- Create: `lagrangianHamiltonianMechanics/mechlib/+mech/lagrange_dalembert.m`
- Create: `lagrangianHamiltonianMechanics/mechlib/tests/test_lagrange_dalembert.m`
- Create: `lagrangianHamiltonianMechanics/mechlib/README.md`
- Create: `lagrangianHamiltonianMechanics/README.md`
- Create: `lagrangianHamiltonianMechanics/run_tests.m`

**Interfaces:**
- Consumes: nothing.
- Produces: `sys = mech.lagrange_dalembert(q, qd, L, A, Q)` returning a struct with fields `f` (`@(t, y)` first-order RHS on y = [q; qd], column [2n x 1]), `lambda` (`@(q, qd)` [k x 1]), `energy` (`@(q, qd)` scalar), `constraint` (`@(q, qd)` [k x 1] = A qd), `qddSym` [n x 1] sym, `lambdaSym` [k x 1] sym, `n`, `k`. Error identifiers `mech:lagrange_dalembert:shape`, `:M`, `:A`.

- [ ] **Step 1: Write the failing test**

`lagrangianHamiltonianMechanics/mechlib/tests/test_lagrange_dalembert.m`:

```matlab
function ok = test_lagrange_dalembert()
%% Purpose:
%
%   mech.lagrange_dalembert on problems with exact answers: a free
%   particle (no constraint), a damped particle (generalized forces enter),
%   the planar pendulum written as a constrained particle (the multiplier
%   is the string tension), a knife-edge particle (v and w constant,
%   lambda = m v w), and the three refusals.
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
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                      % mechlib (the +mech package root)

% 1. free particle: qdd = 0, lambda is [0 x 1]
syms x y xd yd real
L = (xd^2 + yd^2)/2;
sys = mech.lagrange_dalembert([x; y], [xd; yd], L, sym(zeros(0, 2)));
yd0 = sys.f(0, [1; 2; 3; 4]);
ok = check(ok, isequal(size(yd0), [4 1]) && max(abs(yd0 - [3; 4; 0; 0])) < 1e-15, ...
           'free particle: f = [qd; 0], column [4 x 1]');
ok = check(ok, isequal(size(sys.lambda([1; 2], [3; 4])), [0 1]) && sys.k == 0 && sys.n == 2, ...
           'free particle: lambda is [0 x 1], k = 0, n = 2');
ok = check(ok, abs(sys.energy([1; 2], [3; 4]) - 12.5) < 1e-15, 'free particle: energy = |qd|^2/2');

% 2. damped particle: Q = -c qd must reach the right-hand side
c = 0.7;
sysD = mech.lagrange_dalembert([x; y], [xd; yd], L, sym(zeros(0, 2)), -c*[xd; yd]);
ydD = sysD.f(0, [0; 0; 3; 4]);
ok = check(ok, max(abs(ydD(3:4) - (-c*[3; 4]))) < 1e-15, 'damped particle: qdd = -c qd (Q is used)');

% 3. pendulum as a constrained particle: A = [x, y], lambda = -T/l
m = 1; g = 9.81; l = 1;
Lp = m/2*(xd^2 + yd^2) - m*g*y;
sysP = mech.lagrange_dalembert([x; y], [xd; yd], Lp, [x, y]);
v = 0.7;
lam = sysP.lambda([0; -l], [v; 0]);
T = m*g + m*v^2/l;
ok = check(ok, abs(lam + T/l) < 1e-10, sprintf('pendulum at the bottom: lambda = -T/l (%.2e)', abs(lam + T/l)));
phi = 0.6;
qg = l*[sin(phi); -cos(phi)];  qdg = v*[cos(phi); sin(phi)];
lamg = sysP.lambda(qg, qdg);
Tg = m*g*cos(phi) + m*v^2/l;
ok = check(ok, abs(lamg + Tg/l) < 1e-10, sprintf('pendulum at phi = 0.6: lambda = -T/l (%.2e)', abs(lamg + Tg/l)));
ok = check(ok, isequal(size(sysP.qddSym), [2 1]) && isequal(size(sysP.lambdaSym), [1 1]), ...
           'symbolic reduced equations have the right shapes');

% 3b. integrate two periods: energy drift and the small-amplitude period
optsO = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
phi0 = 0.3;  y0 = [l*sin(phi0); -l*cos(phi0); 0; 0];
T0 = 2*pi*sqrt(l/g);
[tt, yy] = ode113(sysP.f, linspace(0, 2*T0, 801), y0, optsO);
E = arrayfun(@(kk) sysP.energy(yy(kk, 1:2).', yy(kk, 3:4).'), 1:numel(tt));
drift = max(abs(E - E(1)))/abs(E(1));
ok = check(ok, drift < 1e-9, sprintf('pendulum: energy drift over two periods %.2e (< 1e-9)', drift));
rad = max(abs(sqrt(sum(yy(:, 1:2).^2, 2)) - l));
ok = check(ok, rad < 1e-9, sprintf('pendulum: radius drift %.2e (< 1e-9)', rad));
optsE = odeset(optsO, 'Events', @(t, y) deal(y(1), 0, 1));
[~, ~, te] = ode113(sysP.f, [0 3*T0], [l*sin(1e-3); -l*cos(1e-3); 0; 0], optsE);
per = te(2) - te(1);
ok = check(ok, abs(per - T0)/T0 < 1e-6, sprintf('pendulum: small-amplitude period %.9f vs 2 pi sqrt(l/g) %.9f (rel %.1e < 1e-6)', per, T0, abs(per - T0)/T0));

% 4. knife-edge particle on q = (x, y, th): v, w constant, lambda = m v w
syms th thd real
mk = 1.3; Jk = 0.4;
Lk = mk/2*(xd^2 + yd^2) + Jk/2*thd^2;
sysK = mech.lagrange_dalembert([x; y; th], [xd; yd; thd], Lk, [-sin(th), cos(th), 0]);
thk = 0.4; vk = 1.3; wk = 0.8;
qk = [0; 0; thk];  qdk = [vk*cos(thk); vk*sin(thk); wk];
lamk = sysK.lambda(qk, qdk);
ok = check(ok, abs(lamk - mk*vk*wk) < 1e-12, sprintf('knife edge: lambda = m v w (%.2e)', abs(lamk - mk*vk*wk)));
fk = sysK.f(0, [qk; qdk]);
vdot = fk(4)*cos(thk) + fk(5)*sin(thk) + wk*(-qdk(1)*sin(thk) + qdk(2)*cos(thk));
ok = check(ok, abs(vdot) < 1e-12 && abs(fk(6)) < 1e-12, sprintf('knife edge: vd = %.1e, wd = %.1e (both 0)', vdot, fk(6)));
ok = check(ok, abs(sysK.constraint(qk, qdk)) < 1e-15, 'knife edge: constraint handle returns A qd = 0 on an admissible velocity');

% 5. refusals
ok = check(ok, throws(@() mech.lagrange_dalembert([x; y], xd, L, sym(zeros(0, 2))), 'mech:lagrange_dalembert:shape'), ...
           'mismatched qd refused with :shape');
ok = check(ok, throws(@() mech.lagrange_dalembert([x; y], [xd; yd], xd^2/2, sym(zeros(0, 2))), 'mech:lagrange_dalembert:M'), ...
           'singular mass matrix refused with :M');
ok = check(ok, throws(@() mech.lagrange_dalembert([x; y], [xd; yd], L, [x, y; 2*x, 2*y]), 'mech:lagrange_dalembert:A'), ...
           'rank-deficient A refused with :A');

if ok, fprintf('TEST_LAGRANGE_DALEMBERT: ALL PASS\n'); else, fprintf('TEST_LAGRANGE_DALEMBERT: FAILURES\n'); end
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
%% Purpose:
%
%   Print one PASS/FAIL line and fold the condition into the verdict.
%
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && logical(cond);
end

% ------------------------------------------------------------------------
function tf = throws(fun, id)
%% Purpose:
%
%   True when fun() throws an error with identifier id.
%
tf = false;
try
    fun();
catch err
    tf = strcmp(err.identifier, id);
end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Through `mcp__matlab__evaluate_matlab_code`:
```matlab
cd('/Users/msc/Desktop/optimal_control/lagrangianHamiltonianMechanics/mechlib/tests'); test_lagrange_dalembert
```
Expected: an error that `mech.lagrange_dalembert` is not found.

- [ ] **Step 3: Write the builder**

`lagrangianHamiltonianMechanics/mechlib/+mech/lagrange_dalembert.m`:

```matlab
function sys = lagrange_dalembert(q, qd, L, A, Q)
%% Purpose:
%
%   Lagrange-d'Alembert equations for a mechanical system whose
%   constraints are linear in velocity, A(q) qd = 0 (nonholonomic, or
%   holonomic written in velocity form), from a SYMBOLIC Lagrangian. With
%   p = dL/dqd, M = dp/dqd and h = (dp/dq) qd - dL/dq, the equations
%
%       M qdd + h = Q + A' lambda,     A qdd + Ad qd = 0
%
%   are solved as the block system [M A'; A 0] [qdd; -lambda] = [Q - h;
%   -Ad qd], where Ad = sum_j (dA/dq_j) qd_j. The result is returned as
%   numeric function handles (matlabFunction) and as the symbolic reduced
%   expressions, so a document can print them. Physics never enters this
%   file: the caller supplies L, A and Q.
%
%% Inputs:
%
%  q                        [n x 1] sym             Generalized coordinates
%
%  qd                       [n x 1] sym             Generalized velocities
%
%  L                        sym scalar              Lagrangian L(q, qd)
%
%  A                        [k x n] sym             Constraint matrix,
%                                                   A(q) qd = 0; k = 0 is
%                                                   sym(zeros(0, n))
%
%  Q                        [n x 1] sym (optional)  Generalized forces
%                                                   Q(q, qd); zeros
%
%% Outputs:
%
%  sys                      struct                  .f(t, y): first-order
%                                                   RHS, y = [q; qd],
%                                                   returns [qd; qdd]
%                                                   [2n x 1]; .lambda(q,
%                                                   qd) [k x 1]; .energy(q,
%                                                   qd) = p' qd - L;
%                                                   .constraint(q, qd) = A
%                                                   qd [k x 1]; .qddSym [n
%                                                   x 1]; .lambdaSym [k x
%                                                   1]; .n; .k
%
%% References:
%   [1] Bloch, "Nonholonomic Mechanics and Control," 2nd ed., Springer,
%       2015, Sec. 1.4 (Lagrange-d'Alembert principle).
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the planar pendulum as a constrained particle; lambda = -T/l
     syms x y xd yd real
     m = 1;  g = 9.81;  l = 1;
     sys = mech.lagrange_dalembert([x; y], [xd; yd], m/2*(xd^2 + yd^2) - m*g*y, [x, y]);
     v = 0.7;
     lam = sys.lambda([0; -l], [v; 0]);
     fprintf('pendulum at the bottom, v = %.2f: lambda = %.6f, -T/l = %.6f\n', v, lam, -(m*g + m*v^2/l)/l);
     disp('reduced accelerations:');  disp(sys.qddSym);
     return
end

q = q(:);  qd = qd(:);
n = numel(q);
if nargin < 5 || isempty(Q), Q = sym(zeros(n, 1)); end
Q = Q(:);
if isempty(A), A = sym(zeros(0, n)); end
if numel(qd) ~= n || numel(Q) ~= n || ~isscalar(L) || size(A, 2) ~= n
    error('mech:lagrange_dalembert:shape', ...
          'q and qd must be [n x 1], Q [n x 1], L scalar and A [k x n] (n = %d; got qd %d, Q %d, A columns %d)', ...
          n, numel(qd), numel(Q), size(A, 2));
end
k = size(A, 1);

p = jacobian(L, qd).';                      % generalized momenta [n x 1]
M = jacobian(p, qd);                        % mass matrix [n x n]
if isAlways(simplify(det(M)) == 0, 'Unknown', 'false')
    error('mech:lagrange_dalembert:M', 'the mass matrix d2L/dqd2 is singular');
end
if k > 0 && rank(A) < k
    error('mech:lagrange_dalembert:A', 'the constraint matrix A has rank %d < k = %d', rank(A), k);
end
h  = jacobian(p, q)*qd - jacobian(L, q).';  % Euler-Lagrange terms without M qdd
Ad = sym(zeros(k, n));
for jj = 1:n
    Ad = Ad + diff(A, q(jj))*qd(jj);
end

K   = [M, A.'; A, sym(zeros(k, k))];
rhs = [Q - h; -Ad*qd];
sol = K \ rhs;
qdd = simplify(sol(1:n));
lam = simplify(-sol(n+1:end));
E   = simplify(p.'*qd - L);

fq = matlabFunction(qdd, 'Vars', {q, qd});
fE = matlabFunction(E,   'Vars', {q, qd});
if k > 0
    fl = matlabFunction(lam, 'Vars', {q, qd});
    fA = matlabFunction(A,   'Vars', {q});
    sys.lambda     = @(qq, qqd) reshape(fl(qq(:), qqd(:)), k, 1);
    sys.constraint = @(qq, qqd) reshape(fA(qq(:))*qqd(:), k, 1);
else
    sys.lambda     = @(qq, qqd) zeros(0, 1);
    sys.constraint = @(qq, qqd) zeros(0, 1);
end
sys.f         = @(t, y) [y(n+1:2*n); reshape(fq(y(1:n), y(n+1:2*n)), n, 1)];
sys.energy    = @(qq, qqd) fE(qq(:), qqd(:));
sys.qddSym    = qdd;
sys.lambdaSym = lam;
sys.n         = n;
sys.k         = k;
end
```

- [ ] **Step 4: Run the test to verify it passes**

Same call as Step 2. Expected: every line `PASS`, final line `TEST_LAGRANGE_DALEMBERT: ALL PASS`. If the pendulum period line fails by more than 1e-6, the events call is the first suspect: `te` must have at least two entries (3 periods of integration), and `direction = 1` picks the upward crossings of x only.

- [ ] **Step 5: Write the two READMEs and the runner**

`lagrangianHamiltonianMechanics/mechlib/README.md`:

```markdown
# mechlib: the mechanics library (`+mech`)

Shared machinery for the Lagrangian/Hamiltonian mechanics systems in this
tree (the Chaplygin sleigh first; rolling disk, rolling coin and knife edge
to follow). MATLAB package form: `addpath('.../mechlib')` once, then call
`mech.<function>`. It is a sibling of `oclib`, not part of it: oclib holds
optimal-control structure and stays physics-free, and when an
optimal-control problem on one of these systems appears it will consume
this library as a plant.

**Admission rule:** a function enters mechlib when it is written without
any system-specific quantity (no mass, length or named coordinate of a
particular system in its body) and carries a test on a problem with an
exact answer. System folders keep nothing generic.

## Contents

| function | what | exact-answer test |
|---|---|---|
| `mech.lagrange_dalembert` | Lagrange-d'Alembert equations with multipliers from a symbolic Lagrangian L(q, qd), constraint matrix A(q) (A qd = 0) and generalized forces Q; returns numeric handles (RHS, multipliers, energy, constraint) and the symbolic reduced equations | free particle, damped particle, pendulum as a constrained particle (lambda = -T/l, period), knife-edge particle (lambda = m v w) |
| `mech.conservation_check` | energy drift and constraint residual along a solution, with bars | exact circular motion passes; injected energy ramp and constraint violation fail with the right reason |
| `mech.animate_planar_body` | top-view MP4 of a planar rigid body: pose history, outline and markers in body coordinates, arrows and subtitle as functions of time, growing traces, fixed limits, exact 1280 x 720 frames | a 20-frame movie exists with 1280 x 720 frames; a single-point pose still gets finite limits; unsorted time refused |

## Tests

`tests/` holds one test per function; each returns `ok`, prints PASS/FAIL
lines with the measured number and its bar, and runs from a fresh clone.
Run them all with `run_tests` from `lagrangianHamiltonianMechanics/`.

## House style

Header quartet, `nargin == 0` self-demos, never `i`/`j`, never `norm`, no
Code Analyzer pragmas. Symbolic Math Toolbox is used by
`lagrange_dalembert` only.
```

`lagrangianHamiltonianMechanics/README.md`:

```markdown
# lagrangianHamiltonianMechanics

Worked mechanical systems, each in its own folder with a front-door script
that integrates the equations of motion, checks them against an
independent derivation and an exact result where one exists, and writes
static plots and a movie. Shared machinery lives in `mechlib/` (`+mech`).

| folder | system | status |
|---|---|---|
| `mechlib/` | `mech.lagrange_dalembert`, `mech.conservation_check`, `mech.animate_planar_body` | built |
| `nonholonomic/ChaplyginSleigh/` | the Chaplygin sleigh: knife-edge blade offset from the center of mass; `run_sleigh.m` | built |

Run every test with `run_tests` from this folder.
```

`lagrangianHamiltonianMechanics/run_tests.m` (Task 1 lists one test; later tasks append their names in the order given in the final list below):

```matlab
function ok = run_tests()
%% Purpose:
%
%   THE runner for lagrangianHamiltonianMechanics: every mechlib test and
%   every system test, library first, one logical verdict, and a THROW on
%   any failure so `matlab -batch "exit(~run_tests())"` returns a real
%   exit code. New systems append their tests to names below.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every test passed (also
%                                                   THROWS on failure)
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'mechlib'), fullfile(here, 'mechlib', 'tests'), ...
        fullfile(here, 'nonholonomic', 'ChaplyginSleigh'), ...
        fullfile(here, 'nonholonomic', 'ChaplyginSleigh', 'tests'));

names = {'test_lagrange_dalembert'    % the builder on exact problems
         'test_conservation_check'    % the diagnostic passes exact motion, fails injected drift
         'test_animate_planar_body'   % 1280 x 720 MP4, degenerate path, refusals
         'test_sleigh_oracle'         % hand-derived vs builder, tanh solution, energy, constraint, stability
         'test_sleigh_pipeline'};     % the script parses, params refusals, plots + movie headless

res  = false(1, numel(names));
secs = zeros(1, numel(names));
for kk = 1:numel(names)
    fprintf('\n==== %s\n', names{kk});
    t0 = tic;
    try
        v = feval(names{kk});
        if ~isscalar(v)
            error('run_tests:badVerdict', '%s did not return a scalar verdict', names{kk});
        end
        res(kk) = logical(v);
    catch err
        fprintf('  THREW  %s (%s)\n', err.message, err.identifier);
        res(kk) = false;
    end
    secs(kk) = toc(t0);
end

fprintf('\n==== SUMMARY (lagrangianHamiltonianMechanics) ====\n');
for kk = 1:numel(names)
    if res(kk), verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('  %-28s %-4s %5.0f s\n', names{kk}, verdict, secs(kk));
end
ok = all(res);
fprintf('  %d of %d passed, %.0f s total\n', nnz(res), numel(res), sum(secs));
if ~ok
    error('run_tests:failed', '%d of %d tests FAILED: %s', nnz(~res), numel(res), strjoin(names(~res), ', '));
end
end
```

In Task 1 the runner's `names` holds only `'test_lagrange_dalembert'`; Tasks 2-5 add their line as they land. The final list is the five above (spec: "lists all seven" counts the two mechlib self-demos that are not tests; the runner lists the five test functions and must print `5 of 5 passed`).

- [ ] **Step 6: Run the runner**

```matlab
cd('/Users/msc/Desktop/optimal_control/lagrangianHamiltonianMechanics'); run_tests
```
Expected: `1 of 1 passed`. Also run `mech.lagrange_dalembert` with no arguments (the demo) and confirm it prints lambda and -T/l agreeing to 6 decimals.

- [ ] **Step 7: Report**

Write the task report (status, files, test output tail, concerns). Do not commit.

---

### Task 2: `mech.conservation_check`

**Files:**
- Create: `lagrangianHamiltonianMechanics/mechlib/+mech/conservation_check.m`
- Create: `lagrangianHamiltonianMechanics/mechlib/tests/test_conservation_check.m`
- Modify: `lagrangianHamiltonianMechanics/run_tests.m` (append `'test_conservation_check'`)

**Interfaces:**
- Consumes: `sys` from Task 1 (`.n`, `.energy(q, qd)`, `.constraint(q, qd)`).
- Produces: `chk = mech.conservation_check(t, Y, sys, bars)` with `.E` [1 x N], `.driftRel`, `.resid`, `.residSeries` [1 x N], `.okEnergy`, `.okConstraint`, `.ok`, `.why`. Error identifier `mech:conservation_check:shape`.

- [ ] **Step 1: Write the failing test**

`lagrangianHamiltonianMechanics/mechlib/tests/test_conservation_check.m`:

```matlab
function ok = test_conservation_check()
%% Purpose:
%
%   mech.conservation_check on exact circular motion of a free particle
%   held on a circle (passes), an injected energy ramp (fails with the
%   energy reason), an injected constraint violation (fails with the
%   constraint reason), and the transposed-input refusal.
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
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

syms x y xd yd real
sys = mech.lagrange_dalembert([x; y], [xd; yd], (xd^2 + yd^2)/2, [x, y]);
t = linspace(0, 4*pi, 401);
Y = [cos(t); sin(t); -sin(t); cos(t)];               % exact circle, |v| = 1

chk = mech.conservation_check(t, Y, sys);
ok = check(ok, chk.ok && chk.driftRel < 1e-12 && chk.resid < 1e-12, ...
           sprintf('exact circle: drift %.1e, residual %.1e, ok = %d', chk.driftRel, chk.resid, chk.ok));
ok = check(ok, isequal(size(chk.E), [1 401]) && max(abs(chk.E - 0.5)) < 1e-12, 'E is [1 x N] and equals 1/2');
ok = check(ok, isempty(chk.why), 'why is empty on a pass');

Yr = Y;  Yr(3:4, :) = Yr(3:4, :).*(1 + 1e-6*t/t(end));   % energy ramp 2e-6 relative
chkR = mech.conservation_check(t, Yr, sys);
ok = check(ok, ~chkR.ok && ~chkR.okEnergy && chkR.okConstraint, 'injected energy ramp FAILS the energy bar only');
ok = check(ok, contains(chkR.why, 'energy'), sprintf('why names energy: "%s"', chkR.why));

Yc = Y;  Yc(3:4, :) = Yc(3:4, :) + 1e-6*Y(1:2, :);         % radial velocity 1e-6
chkC = mech.conservation_check(t, Yc, sys);
ok = check(ok, ~chkC.ok && ~chkC.okConstraint && chkC.okEnergy, 'injected radial velocity FAILS the constraint bar only');
ok = check(ok, contains(chkC.why, 'constraint'), sprintf('why names the constraint: "%s"', chkC.why));

chkB = mech.conservation_check(t, Yr, sys, struct('energy', 1e-3));
ok = check(ok, chkB.ok, 'a looser energy bar passes the ramp');

ok = check(ok, throws(@() mech.conservation_check(t, Y.', sys), 'mech:conservation_check:shape'), ...
           'transposed Y refused with :shape');
ok = check(ok, throws(@() mech.conservation_check(t(1:10), Y, sys), 'mech:conservation_check:shape'), ...
           't and Y of different lengths refused with :shape');

if ok, fprintf('TEST_CONSERVATION_CHECK: ALL PASS\n'); else, fprintf('TEST_CONSERVATION_CHECK: FAILURES\n'); end
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
%% Purpose:
%
%   Print one PASS/FAIL line and fold the condition into the verdict.
%
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && logical(cond);
end

% ------------------------------------------------------------------------
function tf = throws(fun, id)
%% Purpose:
%
%   True when fun() throws an error with identifier id.
%
tf = false;
try
    fun();
catch err
    tf = strcmp(err.identifier, id);
end
end
```

- [ ] **Step 2: Run the test to verify it fails**

```matlab
cd('/Users/msc/Desktop/optimal_control/lagrangianHamiltonianMechanics/mechlib/tests'); test_conservation_check
```
Expected: `mech.conservation_check` not found.

- [ ] **Step 3: Write the function**

`lagrangianHamiltonianMechanics/mechlib/+mech/conservation_check.m`:

```matlab
function chk = conservation_check(t, Y, sys, bars)
%% Purpose:
%
%   Energy drift and constraint residual along a solution of a system
%   built by mech.lagrange_dalembert, each against a bar, with one verdict
%   and a reason. Prints nothing; the caller prints.
%
%% Inputs:
%
%  t                        [1 x N]                 Times (used for the
%                                                   length check only)
%
%  Y                        [2n x N]                Columns [q; qd]
%
%  sys                      struct                  From
%                                                   mech.lagrange_dalembert
%                                                   (.n, .energy,
%                                                   .constraint)
%
%  bars                     struct (optional)       .energy [1e-9] relative
%                                                   drift bar; .constraint
%                                                   [1e-9] max |A qd| bar
%
%% Outputs:
%
%  chk                      struct                  .E [1 x N]; .driftRel
%                                                   = max |E - E(1)| /
%                                                   max(|E(1)|, eps);
%                                                   .resid = max over the
%                                                   grid of max |A qd|;
%                                                   .residSeries [1 x N];
%                                                   .okEnergy;
%                                                   .okConstraint; .ok;
%                                                   .why (empty on a pass)
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a particle on a circle, exact, then with an injected drift
     syms x y xd yd real
     sys = mech.lagrange_dalembert([x; y], [xd; yd], (xd^2 + yd^2)/2, [x, y]);
     t = linspace(0, 2*pi, 201);
     Y = [cos(t); sin(t); -sin(t); cos(t)];
     chk = mech.conservation_check(t, Y, sys);
     fprintf('exact circle: drift %.1e, residual %.1e, ok = %d\n', chk.driftRel, chk.resid, chk.ok);
     Y(3:4, :) = Y(3:4, :).*(1 + 1e-6*t/t(end));          % a ramp: a uniform scale would rescale E(1) too and show nothing
     chk = mech.conservation_check(t, Y, sys);
     fprintf('ramped velocity: drift %.1e, ok = %d, why = %s\n', chk.driftRel, chk.ok, chk.why);
     return
end

if nargin < 4 || isempty(bars), bars = struct(); end
bE = getf(bars, 'energy', 1e-9);
bC = getf(bars, 'constraint', 1e-9);
n  = sys.n;
N  = size(Y, 2);
if size(Y, 1) ~= 2*n || numel(t) ~= N
    error('mech:conservation_check:shape', 'Y must be [2n x N] with N = numel(t) (n = %d; got Y %d x %d, numel(t) = %d)', ...
          n, size(Y, 1), size(Y, 2), numel(t));
end

E     = zeros(1, N);
resid = zeros(1, N);
for kk = 1:N
    qk  = Y(1:n, kk);
    qdk = Y(n+1:2*n, kk);
    E(kk) = sys.energy(qk, qdk);
    c = sys.constraint(qk, qdk);
    if ~isempty(c), resid(kk) = max(abs(c)); end
end

chk.E            = E;
chk.driftRel     = max(abs(E - E(1)))/max(abs(E(1)), eps);
chk.residSeries  = resid;
chk.resid        = max(resid);
chk.okEnergy     = chk.driftRel <= bE;
chk.okConstraint = chk.resid <= bC;
chk.ok           = chk.okEnergy && chk.okConstraint;
chk.why          = '';
if ~chk.okEnergy
    chk.why = sprintf('energy drift %.2e > %.1e', chk.driftRel, bE);
end
if ~chk.okConstraint
    if isempty(chk.why), sep = ''; else, sep = '; '; end
    chk.why = sprintf('%s%sconstraint residual %.2e > %.1e', chk.why, sep, chk.resid, bC);
end
end

% ------------------------------------------------------------------------
function v = getf(s, name, default)
%% Purpose:
%
%   Field with a default.
%
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Same call as Step 2. Expected: all `PASS`, `TEST_CONSERVATION_CHECK: ALL PASS`.

- [ ] **Step 5: Register and rerun the runner**

Append `'test_conservation_check'` to `names` in `run_tests.m`; run it: `2 of 2 passed`. Add the README row if Task 1's README did not already carry it (it does). Report; do not commit.

---

### Task 3: `mech.animate_planar_body`

**Files:**
- Create: `lagrangianHamiltonianMechanics/mechlib/+mech/animate_planar_body.m`
- Create: `lagrangianHamiltonianMechanics/mechlib/tests/test_animate_planar_body.m`
- Modify: `lagrangianHamiltonianMechanics/run_tests.m` (append `'test_animate_planar_body'`)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `file = mech.animate_planar_body(t, pose, body, stem, opts)`. `pose` [3 x N] (x, y, th) of the reference point; `body.outline` [2 x P] body coordinates; `body.markers` struct array (`.name`, `.pos` [2 x 1], `.color` [1 x 3]); `body.segments` struct array (`.pts` [2 x 2] body-coordinate endpoints as columns, `.color`, `.width`); `body.arrows` struct array (`.name`, `.fun` = `@(tk)` returning [base; vec] [4 x 1] in WORLD coordinates at time tk, `.color`); `body.subtitle` = `@(tk)` returning a char line (optional). `opts.fps` [30], `.seconds` [t(end) - t(1)], `.traces` [2 x N x T] world paths on the grid t (default: the reference point), `.title` [''], `.traceColors` [T x 3]. Returns the MP4 path. Error identifier `mech:animate_planar_body:t` (unsorted, repeated, or fewer than two times) and `:shape`.

Note on the spec: section 4.3 wrote arrows and subtitle as functions of the frame INDEX. The plan uses the movie TIME, because the caller does not know the movie's own grid; the spec is amended to match.

- [ ] **Step 1: Write the failing test**

`lagrangianHamiltonianMechanics/mechlib/tests/test_animate_planar_body.m`:

```matlab
function ok = test_animate_planar_body()
%% Purpose:
%
%   mech.animate_planar_body writes a 20-frame MP4 with exact 1280 x 720
%   frames from a rectangle sliding on a circle (with a marker, a segment,
%   an arrow and a subtitle), survives a single-point pose, and refuses an
%   unsorted time grid.
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
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

t = linspace(0, 2, 41);
pose = [cos(pi*t); sin(pi*t); pi*t + pi/2];
body.outline  = [-0.15 0.15 0.15 -0.15; -0.08 -0.08 0.08 0.08];
body.markers  = struct('name', 'nose', 'pos', [0.15; 0], 'color', [1 0.3 0.3]);
body.segments = struct('pts', [-0.1 0.1; 0 0], 'color', [0.5 0.9 1], 'width', 3);
body.arrows   = struct('name', 'v', 'fun', @(tk) [cos(pi*tk); sin(pi*tk); -0.3*sin(pi*tk); 0.3*cos(pi*tk)], 'color', [1 0.75 0.3]);
body.subtitle = @(tk) sprintf('t = %.2f s', tk);
stem = fullfile(tempdir, 'test_animate_planar_body');
file = mech.animate_planar_body(t, pose, body, stem, struct('fps', 10, 'seconds', 2, 'title', 'test'));
ok = check(ok, exist(file, 'file') == 2, sprintf('MP4 written: %s', file));
d = dir(file);
ok = check(ok, d.bytes > 10000, sprintf('MP4 is %d bytes (> 10000)', d.bytes));
vr = VideoReader(file);
fr = readFrame(vr);
ok = check(ok, isequal(size(fr, 1:2), [720 1280]), sprintf('frame is %d x %d (720 x 1280)', size(fr, 1), size(fr, 2)));
ok = check(ok, abs(vr.NumFrames - 20) <= 1, sprintf('%d frames (20 asked)', vr.NumFrames));

file1 = mech.animate_planar_body([0 1], [0 0; 0 0; 0 0], struct('outline', body.outline), fullfile(tempdir, 'test_animate_point'), struct('fps', 5, 'seconds', 0.4));
ok = check(ok, exist(file1, 'file') == 2, 'a single-point pose (zero path extent) still writes a movie with finite limits');

ok = check(ok, throws(@() mech.animate_planar_body([0 2 1], pose(:, 1:3), body, stem), 'mech:animate_planar_body:t'), ...
           'unsorted t refused with :t');
ok = check(ok, throws(@() mech.animate_planar_body(t, pose(1:2, :), body, stem), 'mech:animate_planar_body:shape'), ...
           'pose with two rows refused with :shape');

if ok, fprintf('TEST_ANIMATE_PLANAR_BODY: ALL PASS\n'); else, fprintf('TEST_ANIMATE_PLANAR_BODY: FAILURES\n'); end
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
%% Purpose:
%
%   Print one PASS/FAIL line and fold the condition into the verdict.
%
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && logical(cond);
end

% ------------------------------------------------------------------------
function tf = throws(fun, id)
%% Purpose:
%
%   True when fun() throws an error with identifier id.
%
tf = false;
try
    fun();
catch err
    tf = strcmp(err.identifier, id);
end
end
```

- [ ] **Step 2: Run the test to verify it fails**

```matlab
cd('/Users/msc/Desktop/optimal_control/lagrangianHamiltonianMechanics/mechlib/tests'); test_animate_planar_body
```
Expected: `mech.animate_planar_body` not found.

- [ ] **Step 3: Write the animator**

`lagrangianHamiltonianMechanics/mechlib/+mech/animate_planar_body.m` (the frame discipline is `drone/quad_flight/movie_quad.m`'s: fixed limits, dark theme, every frame forced to 1280 x 720 because H.264 shears other sizes into diagonal streaks):

```matlab
function file = animate_planar_body(t, pose, body, stem, opts)
%% Purpose:
%
%   A 1280 x 720 top-view MP4 of a planar rigid body: its outline, markers
%   and decorative segments given once in body coordinates and carried by
%   the pose (x, y, theta) of a reference point; arrows and a subtitle
%   supplied by the caller as functions of time; traces that grow with
%   time; a fixed camera and fixed, equal-aspect limits covering the whole
%   motion (a frame that rescales reads as motion that is not there).
%   Frames are forced to an exact 1280 x 720: H.264 shears a frame whose
%   size is not a multiple of 16 into diagonal streaks. Knows no system.
%
%% Inputs:
%
%  t                        [1 x N]                 Times (s), strictly
%                                                   increasing, N >= 2
%
%  pose                     [3 x N]                 (x, y, theta) of the
%                                                   reference point in the
%                                                   world frame; theta
%                                                   counterclockwise from x
%
%  body                     struct                  .outline [2 x P] body
%                                                   coordinates (closed by
%                                                   this function);
%                                                   .markers struct array
%                                                   (.name, .pos [2 x 1],
%                                                   .color); .segments
%                                                   struct array (.pts [2
%                                                   x 2] endpoints as
%                                                   columns, .color,
%                                                   .width); .arrows struct
%                                                   array (.name, .fun =
%                                                   @(tk) -> [base; vec]
%                                                   [4 x 1] WORLD, .color);
%                                                   .subtitle = @(tk) ->
%                                                   char. All but .outline
%                                                   optional
%
%  stem                     char                    Output path without
%                                                   extension
%
%  opts                     struct (optional)       .fps [30]; .seconds
%                                                   [t(end) - t(1)];
%                                                   .traces [2 x N x T]
%                                                   world paths on t
%                                                   [default: the
%                                                   reference point];
%                                                   .traceColors [T x 3];
%                                                   .title ['']
%
%% Outputs:
%
%  file                     char                    The MP4 written
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a rectangle sliding on a circle, two seconds, to tempdir
     t = linspace(0, 2, 61);
     pose = [cos(pi*t); sin(pi*t); pi*t + pi/2];
     body.outline = [-0.15 0.15 0.15 -0.15; -0.08 -0.08 0.08 0.08];
     body.markers = struct('name', 'nose', 'pos', [0.15; 0], 'color', [1 0.3 0.3]);
     body.subtitle = @(tk) sprintf('t = %.2f s', tk);
     file = mech.animate_planar_body(t, pose, body, fullfile(tempdir, 'demo_planar_body'), struct('title', 'demo'));
     fprintf('wrote %s\n', file);
     return
end

if nargin < 5 || isempty(opts), opts = struct(); end
t = t(:).';
N = numel(t);
if N < 2 || any(diff(t) <= 0)
    error('mech:animate_planar_body:t', 't must be strictly increasing with at least two entries (got %d)', N);
end
if ~isequal(size(pose), [3 N])
    error('mech:animate_planar_body:shape', 'pose must be [3 x N] with N = numel(t) = %d (got %d x %d)', N, size(pose, 1), size(pose, 2));
end
if ~isfield(body, 'outline') || size(body.outline, 1) ~= 2
    error('mech:animate_planar_body:shape', 'body.outline must be [2 x P]');
end
markers  = getf(body, 'markers',  struct('name', {}, 'pos', {}, 'color', {}));
segments = getf(body, 'segments', struct('pts', {}, 'color', {}, 'width', {}));
arrows   = getf(body, 'arrows',   struct('name', {}, 'fun', {}, 'color', {}));
subtitle = getf(body, 'subtitle', @(tk) sprintf('t = %.2f s', tk));
fps      = getf(opts, 'fps', 30);
seconds  = getf(opts, 'seconds', t(end) - t(1));
ttl      = getf(opts, 'title', '');
traces   = getf(opts, 'traces', reshape(pose(1:2, :), 2, N, 1));
nTr      = size(traces, 3);
trColors = getf(opts, 'traceColors', repmat([0.45 0.70 1.00], nTr, 1));

nFrame = max(2, round(fps*seconds));
tMov   = linspace(t(1), t(end), nFrame);
P      = interp1(t, pose.', tMov, 'pchip').';                 % [3 x nFrame]
Tr     = zeros(2, nFrame, nTr);
for tt = 1:nTr
    Tr(:, :, tt) = interp1(t, traces(:, :, tt).', tMov, 'pchip').';
end

% body size and fixed limits: every trace, every reference point, plus the
% body's own reach, with a margin; equal aspect at 16:9 so the frame fills
pts   = [body.outline, [markers.pos], reshape([segments.pts], 2, [])];
if isempty(pts), pts = zeros(2, 1); end
reach = max([sqrt(sum(pts.^2, 1)), eps]);
allXY = [reshape(Tr, 2, []), P(1:2, :)];
lo    = min(allXY, [], 2) - 1.2*reach;
hi    = max(allXY, [], 2) + 1.2*reach;
span  = max(hi - lo);
lo    = lo - 0.05*span;  hi = hi + 0.05*span;
mid   = (lo + hi)/2;
half  = max((hi - lo)/2, 1.5*reach);
if half(1)/half(2) < 16/9, half(1) = half(2)*16/9; else, half(2) = half(1)*9/16; end
lim   = [mid - half, mid + half];                             % [2 x 2] rows x, y

fig = figure('Color', 'k', 'Position', [100 100 1280 720], 'InvertHardcopy', 'off', 'Visible', 'off');
ax  = axes(fig, 'Position', [0.06 0.08 0.88 0.80]);
frames(nFrame) = struct('cdata', [], 'colormap', []);
closedOutline = [body.outline, body.outline(:, 1)];
for kk = 1:nFrame
    cla(ax);  hold(ax, 'on');
    for tt = 1:nTr
        plot(ax, Tr(1, 1:kk, tt), Tr(2, 1:kk, tt), '-', 'Color', [trColors(tt, :) 0.8], 'LineWidth', 1.6);
    end
    c  = P(1:2, kk);  th = P(3, kk);
    R  = [cos(th), -sin(th); sin(th), cos(th)];
    O  = c + R*closedOutline;
    patch(ax, O(1, :), O(2, :), [1 0.75 0.3], 'FaceAlpha', 0.35, 'EdgeColor', [1 0.75 0.3], 'LineWidth', 2);
    for ss = 1:numel(segments)
        S = c + R*segments(ss).pts;
        plot(ax, S(1, :), S(2, :), '-', 'Color', segments(ss).color, 'LineWidth', segments(ss).width);
    end
    for mm = 1:numel(markers)
        Mk = c + R*markers(mm).pos;
        plot(ax, Mk(1), Mk(2), 'o', 'MarkerFaceColor', markers(mm).color, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
    end
    for aa = 1:numel(arrows)
        bv = arrows(aa).fun(tMov(kk));
        quiver(ax, bv(1), bv(2), bv(3), bv(4), 0, 'Color', arrows(aa).color, 'LineWidth', 2, 'MaxHeadSize', 0.6);
    end
    set(ax, 'Color', 'k', 'XColor', [0.75 0.75 0.75], 'YColor', [0.75 0.75 0.75], ...
            'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.35, 'XLim', lim(1, :), 'YLim', lim(2, :), ...
            'DataAspectRatio', [1 1 1]);
    grid(ax, 'on');
    xlabel(ax, 'x  [m]', 'Color', [0.8 0.8 0.8]);
    ylabel(ax, 'y  [m]', 'Color', [0.8 0.8 0.8]);
    annotation(fig, 'textbox', [0 0.935 1 0.05], 'String', ttl, 'Color', 'w', 'FontSize', 18, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'Interpreter', 'none');
    annotation(fig, 'textbox', [0 0.895 1 0.04], 'String', subtitle(tMov(kk)), 'Color', [0.72 0.72 0.72], ...
        'FontSize', 13, 'FontName', 'FixedWidth', 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'Interpreter', 'none');
    drawnow;
    frames(kk) = getframe(fig);
    delete(findall(fig, 'Type', 'textboxshape'));
end
close(fig);

for kk = 1:nFrame
    if ~isequal(size(frames(kk).cdata, 1:2), [720 1280])
        frames(kk).cdata = imresize(frames(kk).cdata, [720 1280]);
    end
end
file = [stem '.mp4'];
v = VideoWriter(file, 'MPEG-4');
v.FrameRate = fps;  v.Quality = 95;
open(v);
for kk = 1:nFrame, writeVideo(v, frames(kk)); end
close(v);
end

% ------------------------------------------------------------------------
function v = getf(s, name, default)
%% Purpose:
%
%   Field with a default.
%
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Same call as Step 2. Expected: all `PASS`. If `getframe` on an invisible figure returns a frame of the wrong size, the `imresize` pass covers it; if `getframe` returns empty, set the figure `'Visible', 'on'` (movie_quad does) and re-run.

- [ ] **Step 5: Register and rerun the runner**

Append `'test_animate_planar_body'` to `names`; `run_tests` prints `3 of 3 passed`. Report; do not commit.

---

### Task 4: the sleigh physics and its oracle test

**Files:**
- Create: `lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/sleigh_params.m`
- Create: `lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/sleigh_lagrangian.m`
- Create: `lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/sleigh_rhs.m`
- Create: `lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/sleigh_quantities.m`
- Create: `lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/tests/test_sleigh_oracle.m`
- Modify: `lagrangianHamiltonianMechanics/run_tests.m` (append `'test_sleigh_oracle'`)

**Interfaces:**
- Consumes: `mech.lagrange_dalembert` (Task 1), `mech.conservation_check` (Task 2).
- Produces: `p = sleigh_params(m, J, a)` with `.m .J .a .Jb .k`; `[q, qd, L, A] = sleigh_lagrangian(p)`; `yd = sleigh_rhs(t, y, p)` vectorized over columns of y = [x_b; y_b; th; v; w]; `out = sleigh_quantities(t, Y, p)` with `.rb .rc` [2 x N], `.E .lambda .vExact .wExact` [1 x N], `.C`, `.asymptote`, `.Yfull` [6 x N] = [q; qd], `.Yfulld` [6 x N] = [qd; qdd] from the hand-derived equations. Error identifier `sleigh_params:value`.

- [ ] **Step 1: Write the failing test**

`lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/tests/test_sleigh_oracle.m`:

```matlab
function ok = test_sleigh_oracle()
%% Purpose:
%
%   The sleigh's hand-derived reduced equations against three independent
%   references: the symbolic Lagrange-d'Alembert builder evaluated on the
%   same solution (D1), the exact tanh solution (D2), energy conservation
%   and the constraint residual (D3, D4) -- for both signs of a, both signs
%   of v0, the knife edge a = 0, the straight ride w0 = 0 and rest. Also
%   the stability statement (D5) and the center-of-mass mapping.
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
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(here)));          % lagrangianHamiltonianMechanics
addpath(fileparts(here), fullfile(root, 'mechlib'));

odeOpts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
t = 0:0.01:8;
w0 = 2;

for a = [0.2, -0.2, 0]
    p = sleigh_params(1, 0.1, a);
    [q, qd, L, A] = sleigh_lagrangian(p);
    sys = mech.lagrange_dalembert(q, qd, L, A);
    for v0 = [1, -1]
        tag = sprintf('a = %+.1f, v0 = %+.0f', a, v0);
        [~, Yt] = ode113(@(tt, yy) sleigh_rhs(tt, yy, p), t, [0; 0; 0; v0; w0], odeOpts);
        Y = Yt.';
        out = sleigh_quantities(t, Y, p);

        % D1: hand-derived derivative vs the builder, on the full state
        d1 = 0;
        for kk = 1:20:numel(t)
            d1 = max(d1, max(abs(sys.f(t(kk), out.Yfull(:, kk)) - out.Yfulld(:, kk))));
        end
        ok = check(ok, d1 <= 1e-9, sprintf('%s: D1 hand-derived vs builder %.2e (<= 1e-9)', tag, d1));

        % D2: exact solution
        d2v = max(abs(Y(4, :) - out.vExact));  d2w = max(abs(Y(5, :) - out.wExact));
        ok = check(ok, d2v <= 1e-8 && d2w <= 1e-8, sprintf('%s: D2 exact solution, v %.2e, w %.2e (<= 1e-8)', tag, d2v, d2w));

        % D3, D4
        chk = mech.conservation_check(t, out.Yfull, sys);
        ok = check(ok, chk.ok, sprintf('%s: D3 energy drift %.2e, D4 constraint %.2e (<= 1e-9) %s', tag, chk.driftRel, chk.resid, chk.why));

        % mapping: the center-of-mass velocity satisfies the constraint
        ok = check(ok, chk.resid <= 1e-12, sprintf('%s: center-of-mass velocity satisfies A qd = 0 (%.1e)', tag, chk.resid));

        % D5 and the sign structure
        ok = check(ok, all(sign(Y(5, :)) == sign(w0)), sprintf('%s: w never changes sign', tag));
        if a ~= 0
            ok = check(ok, abs(Y(4, end) - out.asymptote) < 1e-3, ...
                       sprintf('%s: v(end) = %.6f -> sign(a) sqrt(C) = %.6f (|diff| < 1e-3)', tag, Y(4, end), out.asymptote));
            vs = Y(4, Y(4, :) ~= 0);                  % ignore exact-zero samples: (-, 0, +) is ONE crossing
            nFlip = nnz(diff(sign(vs)) ~= 0);
            if a*v0 > 0
                ok = check(ok, nFlip == 0, sprintf('%s: center of mass leads, v keeps its sign (flips %d)', tag, nFlip));
            else
                ok = check(ok, nFlip == 1, sprintf('%s: center of mass trails, the signed speed crosses zero exactly once (flips %d)', tag, nFlip));
            end
            dth = abs(Y(3, end) - out.thetaInf);
            ok = check(ok, dth < 1e-4, sprintf('%s: total heading change %.6f vs exact sign(w0) (b/|a|) acos(sign(a) v0/S) = %.6f (|diff| %.1e < 1e-4)', ...
                       tag, Y(3, end) - Y(3, 1), out.thetaInf - Y(3, 1), dth));
        else
            cv = max(abs(Y(4, :) - v0));  cw = max(abs(Y(5, :) - w0));
            ok = check(ok, cv < 1e-10 && cw < 1e-10, sprintf('%s: knife edge, v and w constant (%.1e, %.1e)', tag, cv, cw));
            ok = check(ok, max(abs(out.lambda - p.m*Y(4, :).*Y(5, :))) < 1e-12, sprintf('%s: knife edge, lambda = m v w', tag));
            rad = sqrt(sum((out.rb - [v0/w0*(-sin(0)); v0/w0]).^2, 1));   % circle centre (0, v0/w0) for the start at the origin heading east
            ok = check(ok, max(abs(rad - abs(v0/w0))) < 1e-9, sprintf('%s: knife edge runs a circle of radius |v0/w0| (%.1e)', tag, max(abs(rad - abs(v0/w0)))));
        end
    end
end

% w0 = 0: straight line at speed v0 for a ~= 0; rest when v0 = w0 = 0
p = sleigh_params(1, 0.1, 0.2);
[~, Ys] = ode113(@(tt, yy) sleigh_rhs(tt, yy, p), t, [0; 0; 0; 1; 0], odeOpts);
outS = sleigh_quantities(t, Ys.', p);
ok = check(ok, max(abs(Ys(:, 2))) < 1e-12 && max(abs(Ys(:, 3))) < 1e-12 && max(abs(Ys(:, 4) - 1)) < 1e-12, ...
           'w0 = 0: straight line, th constant, v constant');
ok = check(ok, max(abs(outS.vExact - 1)) < 1e-15 && all(outS.wExact == 0) && all(isfinite(outS.vExact)), ...
           'w0 = 0: exact solution is the constant, not NaN');
outR = sleigh_quantities(t, zeros(5, numel(t)), p);
ok = check(ok, all(isfinite(outR.vExact)) && all(isfinite(outR.wExact)) && outR.C == 0 && all(outR.E == 0), ...
           'v0 = w0 = 0: C = 0, exact solution finite (rest), E = 0');

% trailing launch with w0 = 0 is an (unstable) equilibrium: it never reverses
[~, Yu] = ode113(@(tt, yy) sleigh_rhs(tt, yy, p), t, [0; 0; 0; -1; 0], odeOpts);
ok = check(ok, max(abs(Yu(:, 4) + 1)) < 1e-12 && max(abs(Yu(:, 5))) < 1e-12, ...
           'a = +0.2, v0 = -1, w0 = 0: the trailing straight ride is an exact (unstable) equilibrium, v stays -1');

% tiny w0: the exact solution must stay finite (the atanh form would give NaN because C rounds to v0^2)
tL = 0:0.05:25;
[~, Yt] = ode113(@(tt, yy) sleigh_rhs(tt, yy, p), tL, [0; 0; 0; -1; 1e-10], odeOpts);
outT = sleigh_quantities(tL, Yt.', p);
ok = check(ok, all(isfinite(outT.vExact)) && all(isfinite(outT.wExact)) && abs(outT.wExact(1) - 1e-10) < 1e-25, ...
           'w0 = 1e-10: exact solution finite, wExact(1) = w0');
kz = find(Yt(:, 4) > 0, 1);
ok = check(ok, ~isempty(kz) && abs(tL(kz) - 17.29) < 0.1 && max(abs(Yt(:, 4).' - outT.vExact)) < 1e-6, ...
           sprintf('w0 = 1e-10: the trailing sleigh reverses near 17.29 s (first v > 0 at %.2f s), exact solution tracks it', tL(max(kz, 1))));

ok = check(ok, throws(@() sleigh_params(0, 0.1, 0.2), 'sleigh_params:value') && ...
               throws(@() sleigh_params(1, -1, 0.2), 'sleigh_params:value') && ...
               throws(@() sleigh_params(1, 0.1, NaN), 'sleigh_params:value'), ...
           'sleigh_params refuses m <= 0, J <= 0, non-finite a with :value');
p0 = sleigh_params();
ok = check(ok, p0.m == 1 && p0.J == 0.1 && p0.a == 0.2 && abs(p0.Jb - 0.14) < 1e-15 && abs(p0.k - 0.2/0.14) < 1e-15, ...
           'sleigh_params defaults and derived constants');

if ok, fprintf('TEST_SLEIGH_ORACLE: ALL PASS\n'); else, fprintf('TEST_SLEIGH_ORACLE: FAILURES\n'); end
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
%% Purpose:
%
%   Print one PASS/FAIL line and fold the condition into the verdict.
%
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && logical(cond);
end

% ------------------------------------------------------------------------
function tf = throws(fun, id)
%% Purpose:
%
%   True when fun() throws an error with identifier id.
%
tf = false;
try
    fun();
catch err
    tf = strcmp(err.identifier, id);
end
end
```

(The knife-edge circle: starting at the origin heading east with speed v0 and yaw rate w0, the blade runs the circle of radius v0/w0 centred at (0, v0/w0); for v0 = -1, w0 = 2 the centre is (0, -1/2) and the radius 1/2; `abs(v0/w0)` covers both.)

- [ ] **Step 2: Run the test to verify it fails**

```matlab
cd('/Users/msc/Desktop/optimal_control/lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/tests'); test_sleigh_oracle
```
Expected: `sleigh_params` not found.

- [ ] **Step 3: Write the four physics files**

`sleigh_params.m`:

```matlab
function p = sleigh_params(m, J, a)
%% Purpose:
%
%   The Chaplygin sleigh's constants and the two derived ones every other
%   file uses: Jb = J + m a^2 (inertia about the blade point) and k = m a
%   / Jb (the decay rate in wd = -k v w). a is SIGNED: a > 0 puts the
%   center of mass ahead of the blade in the forward (+e1) direction.
%
%% Inputs:
%
%  m                        [scalar]                Mass (kg), > 0 [1]
%
%  J                        [scalar]                Inertia about the
%                                                   center of mass (kg
%                                                   m^2), > 0 [0.1]
%
%  a                        [scalar]                Blade to center of
%                                                   mass (m), signed, finite
%                                                   [0.2]
%
%% Outputs:
%
%  p                        struct                  .m .J .a .Jb .k
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the defaults
     p = sleigh_params(1, 0.1, 0.2);
     disp(p);
     return
end
if ~(isscalar(m) && isfinite(m) && m > 0 && isscalar(J) && isfinite(J) && J > 0 && isscalar(a) && isfinite(a))
    error('sleigh_params:value', 'need scalar m > 0, J > 0 and finite a (got m = %g, J = %g, a = %g)', m, J, a);
end
p.m  = m;
p.J  = J;
p.a  = a;
p.Jb = J + m*a^2;
p.k  = m*a/p.Jb;
end
```

`sleigh_lagrangian.m`:

```matlab
function [q, qd, L, A] = sleigh_lagrangian(p)
%% Purpose:
%
%   The sleigh as mech.lagrange_dalembert wants it: generalized
%   coordinates q = (x_c, y_c, th) of the CENTER OF MASS, the free-body
%   Lagrangian L = m/2 |rd_c|^2 + J/2 thd^2 (no potential), and the
%   knife-edge constraint at the blade point r_b = r_c - a e1: zero
%   sideways velocity, rd_b . e2 = -xd_c sin th + yd_c cos th - a thd = 0,
%   so A = [-sin th, cos th, -a]. Nothing is solved here.
%
%% Inputs:
%
%  p                        struct                  From sleigh_params
%
%% Outputs:
%
%  q                        [3 x 1] sym             [xc; yc; th]
%
%  qd                       [3 x 1] sym             [xcd; ycd; thd]
%
%  L                        sym                     Lagrangian
%
%  A                        [1 x 3] sym             Constraint row
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: print the objects for the defaults
     [q, qd, L, A] = sleigh_lagrangian(sleigh_params());
     disp(q.');  disp(qd.');  disp(L);  disp(A);
     return
end
syms xc yc th xcd ycd thdot real            % not `thd`: a Signal Processing Toolbox function shadows that name
q  = [xc; yc; th];
qd = [xcd; ycd; thdot];
L  = p.m/2*(xcd^2 + ycd^2) + p.J/2*thdot^2;
A  = [-sin(th), cos(th), -p.a];
end
```

`sleigh_rhs.m`:

```matlab
function yd = sleigh_rhs(t, y, p)
%% Purpose:
%
%   The hand-derived reduced equations of the Chaplygin sleigh on the
%   BLADE-point state y = [x_b; y_b; th; v; w]: v is the forward speed of
%   the blade (equal to the forward speed of the center of mass), w the
%   yaw rate. From m rdd_c = lambda e2 and J wd = -a lambda with rd_c = v
%   e1 + a w e2:
%
%       xd_b = v cos th,  yd_b = v sin th,  thd = w,
%       vd = a w^2,       wd = -k v w,      k = m a / (J + m a^2).
%
%   No symbolic call: mech.lagrange_dalembert is this file's independent
%   check (test_sleigh_oracle D1). Vectorized over columns of y.
%
%% Inputs:
%
%  t                        [scalar]                Time (unused;
%                                                   autonomous)
%
%  y                        [5 x N]                 [x_b; y_b; th; v; w]
%
%  p                        struct                  From sleigh_params
%
%% Outputs:
%
%  yd                       [5 x N]                 Time derivative
%
%% References:
%   [1] Bloch, "Nonholonomic Mechanics and Control," 2nd ed., Springer,
%       2015, Sec. 1.7 (the Chaplygin sleigh).
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: five seconds from the defaults
     p = sleigh_params();
     [tt, yy] = ode113(@(a, b) sleigh_rhs(a, b, p), [0 5], [0; 0; 0; 1; 2], odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
     fprintf('t = %.1f s: v = %.6f, w = %.6f\n', tt(end), yy(end, 4), yy(end, 5));
     yd = sleigh_rhs(0, yy(end, :).', p);
     return
end
th = y(3, :);  v = y(4, :);  w = y(5, :);
yd = [v.*cos(th); v.*sin(th); w; p.a*w.^2; -p.k*v.*w];
end
```

`sleigh_quantities.m`:

```matlab
function out = sleigh_quantities(t, Y, p)
%% Purpose:
%
%   Everything derived from a reduced solution: the blade and center-of-
%   mass paths, the energy E = m/2 v^2 + Jb/2 w^2, the blade's constraint
%   force lambda = m J v w / Jb, the invariant C = v^2 + (Jb/m) w^2, the
%   exact solution in its well-conditioned form (b = sqrt(Jb/m), S =
%   sqrt(C), alpha = asinh(v0/(b |w0|)), z = alpha + k S t: v = S tanh z,
%   w = sign(w0) (S/b) sech z; constants when a = 0 or w0 = 0), the
%   asymptote sign(a) S (v0 when a = 0 or w0 = 0), the exact total heading
%   change theta_inf = theta_0 + sign(w0) (b/|a|) acos(sign(a) v0/S), and
%   the full-state arrays [q; qd] and [qd; qdd] of the center-of-mass
%   description that the symbolic builder is checked against.
%
%% Inputs:
%
%  t                        [1 x N]                 Times (s)
%
%  Y                        [5 x N]                 Reduced solution
%                                                   columns [x_b; y_b;
%                                                   th; v; w]
%
%  p                        struct                  From sleigh_params
%
%% Outputs:
%
%  out                      struct                  .rb .rc [2 x N]; .E
%                                                   .lambda .vExact
%                                                   .wExact [1 x N]; .C;
%                                                   .asymptote; .thetaInf
%                                                   (NaN when a = 0 or
%                                                   w0 = 0); .Yfull [6
%                                                   x N] = [xc; yc; th;
%                                                   xcd; ycd; thd];
%                                                   .Yfulld [6 x N] =
%                                                   [qd; qdd] from
%                                                   sleigh_rhs
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the defaults, five seconds
     p = sleigh_params();
     t = 0:0.01:5;
     [~, yy] = ode113(@(a, b) sleigh_rhs(a, b, p), t, [0; 0; 0; 1; 2], odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
     out = sleigh_quantities(t, yy.', p);
     fprintf('C = %.6f, asymptote %.6f, v(end) %.6f, max |v - vExact| %.1e, energy spread %.1e\n', ...
             out.C, out.asymptote, yy(end, 4), max(abs(yy(:, 4).' - out.vExact)), max(out.E) - min(out.E));
     return
end
t  = t(:).';
th = Y(3, :);  v = Y(4, :);  w = Y(5, :);
e1 = [cos(th); sin(th)];
e2 = [-sin(th); cos(th)];

out.rb     = Y(1:2, :);
out.rc     = out.rb + p.a*e1;
vc         = v.*e1 + p.a*w.*e2;                      % center-of-mass velocity
out.E      = p.m/2*v.^2 + p.Jb/2*w.^2;
out.lambda = p.m*p.J*v.*w/p.Jb;
out.C      = v(1)^2 + (p.Jb/p.m)*w(1)^2;

yd = sleigh_rhs(t, Y, p);
vd = yd(4, :);  wd = yd(5, :);
ac = (vd - p.a*w.^2).*e1 + (v.*w + p.a*wd).*e2;     % d/dt (v e1 + a w e2)
out.Yfull  = [out.rc; th; vc; w];
out.Yfulld = [vc; w; ac; wd];

s = t - t(1);
if p.a ~= 0 && w(1) ~= 0
    % the well-conditioned form: alpha = asinh(v0/(b |w0|)) never sees the
    % rounded ratio v0/sqrt(C) = +-1 that the atanh form produces for tiny w0
    b  = sqrt(p.Jb/p.m);
    S  = sqrt(out.C);
    al = asinh(v(1)/(b*abs(w(1))));
    z  = al + p.k*S*s;
    out.vExact    = S*tanh(z);
    out.wExact    = sign(w(1))*(S/b)*sech(z);
    out.asymptote = sign(p.a)*S;
    out.thetaInf  = th(1) + sign(w(1))*(b/abs(p.a))*acos(sign(p.a)*v(1)/S);   % total heading change, exact
else
    out.vExact    = v(1)*ones(size(s));              % knife edge or straight ride: constants
    out.wExact    = w(1)*ones(size(s));
    out.asymptote = v(1);
    out.thetaInf  = NaN;                             % the heading keeps turning (a = 0) or never changes (w0 = 0)
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Same call as Step 2. Expected: all `PASS` for the six (a, v0) cases plus the special cases. If D1 fails at ~1e-8 rather than ~1e-15, the builder's `qdd` and the hand-derived `ac` disagree on a sign: re-derive `ac` from rd_c = v e1 + a w e2 (de1/dt = w e2, de2/dt = -w e1) before touching either file. If D2 fails for a = -0.2, check that `p.k` carries the sign of a (it must) and that `t0` uses `atanh(v0/sC)` with v0 signed.

- [ ] **Step 5: Register and rerun the runner**

Append `'test_sleigh_oracle'`; `run_tests` prints `4 of 4 passed`. Report; do not commit.

---

### Task 5: body shape, plots, the front-door script, pipeline test, folder README

**Files:**
- Create: `lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/sleigh_body.m`
- Create: `lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/plot_sleigh.m`
- Create: `lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/run_sleigh.m`
- Create: `lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/tests/test_sleigh_pipeline.m`
- Create: `lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/README.md`
- Create: `lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/.gitignore` (`results/`)
- Modify: `lagrangianHamiltonianMechanics/run_tests.m` (append `'test_sleigh_pipeline'`)

**Interfaces:**
- Consumes: Task 4's five functions, `mech.lagrange_dalembert`, `mech.conservation_check`, `mech.animate_planar_body`.
- Produces: `body = sleigh_body(p)` (outline, markers com/blade, the blade segment; reference point = the blade); `files = plot_sleigh(t, Y, out, p, info, outDir)` writing `<info.label>_path.png` and `<info.label>_diag.png`, `info.visible` default false; the script `run_sleigh.m`.

- [ ] **Step 1: Write the failing test**

`tests/test_sleigh_pipeline.m`:

```matlab
function ok = test_sleigh_pipeline()
%% Purpose:
%
%   The pieces run_sleigh.m calls, headless and into a temp folder: the
%   body shape, plot_sleigh's two PNGs (also for the degenerate rest
%   case), a 20-frame movie through mech.animate_planar_body with the
%   sleigh's arrows and subtitle, and a parse check of the script itself
%   (scripts take no arguments, so the script is run by hand once).
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
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(mfilename('fullpath'));
sdir = fileparts(here);
root = fileparts(fileparts(sdir));
addpath(sdir, fullfile(root, 'mechlib'));

p = sleigh_params(1, 0.1, 0.2);
t = 0:0.02:4;
[~, yy] = ode113(@(a, b) sleigh_rhs(a, b, p), t, [0; 0; 0; 1; 2], odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
Y = yy.';
out = sleigh_quantities(t, Y, p);

body = sleigh_body(p);
ok = check(ok, size(body.outline, 1) == 2 && size(body.outline, 2) >= 4, 'sleigh_body: outline is [2 x P], P >= 4');
names = {body.markers.name};
ok = check(ok, any(strcmp(names, 'com')) && any(strcmp(names, 'blade')), 'sleigh_body: markers com and blade');
cm = body.markers(strcmp(names, 'com')).pos;
ok = check(ok, max(abs(cm - [p.a; 0])) < 1e-15, 'sleigh_body: center of mass at (a, 0) in blade coordinates');
ok = check(ok, numel(body.segments) >= 1, 'sleigh_body: the blade segment is present');

outDir = fullfile(tempdir, 'test_sleigh_pipeline');
if ~exist(outDir, 'dir'), mkdir(outDir); end
files = plot_sleigh(t, Y, out, p, struct('label', 'tsp', 'visible', false), outDir);
ok = check(ok, numel(files) == 2 && all(cellfun(@(f) exist(f, 'file') == 2, files)), 'plot_sleigh writes its two PNGs');
outR = sleigh_quantities(t, zeros(5, numel(t)), p);
filesR = plot_sleigh(t, zeros(5, numel(t)), outR, p, struct('label', 'tsp_rest', 'visible', false), outDir);
ok = check(ok, numel(filesR) == 2, 'plot_sleigh survives the rest case (C = 0, degenerate ellipse)');

arrows = sleigh_arrows(t, Y, out, p);
body.arrows   = arrows;
body.subtitle = @(tk) sprintf('t = %5.2f s', tk);
file = mech.animate_planar_body(t, [out.rb; Y(3, :)], body, fullfile(outDir, 'tsp'), ...
                                struct('fps', 10, 'seconds', 2, 'traces', cat(3, out.rc, out.rb), 'title', 'test'));
ok = check(ok, exist(file, 'file') == 2, sprintf('movie written through mech.animate_planar_body: %s', file));

msgs = checkcode(fullfile(sdir, 'run_sleigh.m'), '-id');
bad = msgs(arrayfun(@(mm) startsWith(mm.id, 'E') || contains(mm.message, 'Parse error'), msgs));
ok = check(ok, isempty(bad), sprintf('run_sleigh.m parses (%d checkcode errors)', numel(bad)));

if ok, fprintf('TEST_SLEIGH_PIPELINE: ALL PASS\n'); else, fprintf('TEST_SLEIGH_PIPELINE: FAILURES\n'); end
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
%% Purpose:
%
%   Print one PASS/FAIL line and fold the condition into the verdict.
%
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && logical(cond);
end
```

The test calls `sleigh_arrows(t, Y, out, p)`: a small function this task also creates (see Step 3), so the script and the test build the velocity and constraint-force arrows from one place.

- [ ] **Step 2: Run the test to verify it fails**

```matlab
cd('/Users/msc/Desktop/optimal_control/lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/tests'); test_sleigh_pipeline
```
Expected: `sleigh_body` not found.

- [ ] **Step 3: Write sleigh_body, sleigh_arrows, plot_sleigh**

`sleigh_body.m`:

```matlab
function body = sleigh_body(p)
%% Purpose:
%
%   The sleigh's drawing in BLADE-point body coordinates for
%   mech.animate_planar_body: a pointed hull whose length scales with
%   max(|a|, 0.1) and reaches from behind the blade to ahead of the center
%   of mass, the blade as a thick segment along e1 through the blade
%   point, and two markers (center of mass at (a, 0), blade at the
%   origin). The reference point the animator carries is the blade.
%
%% Inputs:
%
%  p                        struct                  From sleigh_params
%
%% Outputs:
%
%  body                     struct                  .outline [2 x 5];
%                                                   .markers (com, blade);
%                                                   .segments (the blade)
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: draw the defaults
     body = sleigh_body(sleigh_params());
     figure;  O = [body.outline, body.outline(:, 1)];
     plot(O(1, :), O(2, :), 'k-', body.segments(1).pts(1, :), body.segments(1).pts(2, :), 'c-', 'LineWidth', 2);
     hold on;  plot(body.markers(1).pos(1), body.markers(1).pos(2), 'ro', body.markers(2).pos(1), body.markers(2).pos(2), 'bo');
     axis equal;  grid on;  title('sleigh body, blade at the origin, +x forward');
     return
end
Lb = max(abs(p.a), 0.1);
xf = max(p.a, 0) + 0.5*Lb;             % nose
xr = min(p.a, 0) - 0.5*Lb;             % tail
hw = 0.3*Lb;                           % half width
body.outline  = [xf, xf - 0.3*Lb, xr, xr, xf - 0.3*Lb; 0, hw, 0.7*hw, -0.7*hw, -hw];
body.markers  = struct('name', {'com', 'blade'}, 'pos', {[p.a; 0], [0; 0]}, ...
                       'color', {[1 0.45 0.2], [0.5 0.9 1]});
body.segments = struct('pts', [-0.25*Lb, 0.25*Lb; 0, 0], 'color', [0.5 0.9 1], 'width', 4);
end
```

`sleigh_arrows.m`:

```matlab
function arrows = sleigh_arrows(t, Y, out, p)
%% Purpose:
%
%   The two arrows the movie draws, as functions of time for
%   mech.animate_planar_body: the center-of-mass velocity at the center of
%   mass and the blade's constraint force lambda e2 at the blade point.
%   Each is scaled so the largest arrow over the run is one body length,
%   max(|a|, 0.1), and interpolated onto the movie's time with pchip.
%
%% Inputs:
%
%  t                        [1 x N]                 Times (s)
%
%  Y                        [5 x N]                 Reduced solution
%
%  out                      struct                  From sleigh_quantities
%
%  p                        struct                  From sleigh_params
%
%% Outputs:
%
%  arrows                   struct array [1 x 2]    .name, .fun = @(tk) ->
%                                                   [base; vec] world, .color
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the defaults, print both arrows at t = 1
     p = sleigh_params();  t = 0:0.01:3;
     [~, yy] = ode113(@(a, b) sleigh_rhs(a, b, p), t, [0; 0; 0; 1; 2], odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
     out = sleigh_quantities(t, yy.', p);
     arrows = sleigh_arrows(t, yy.', out, p);
     fprintf('%s at t = 1: %s\n', arrows(1).name, mat2str(arrows(1).fun(1).', 4));
     fprintf('%s at t = 1: %s\n', arrows(2).name, mat2str(arrows(2).fun(1).', 4));
     return
end
Lb = max(abs(p.a), 0.1);
th = Y(3, :);
e2 = [-sin(th); cos(th)];
vc = out.Yfull(4:5, :);
F  = out.lambda.*e2;
sv = Lb/max(max(sqrt(sum(vc.^2, 1))), eps);
sf = Lb/max(max(abs(out.lambda)), eps);
vcI = @(tk) interp1(t, vc.', tk, 'pchip').';
FI  = @(tk) interp1(t, F.', tk, 'pchip').';
rcI = @(tk) interp1(t, out.rc.', tk, 'pchip').';
rbI = @(tk) interp1(t, out.rb.', tk, 'pchip').';
arrows = struct('name', {'velocity', 'constraint force'}, ...
                'fun',  {@(tk) [rcI(tk); sv*vcI(tk)], @(tk) [rbI(tk); sf*FI(tk)]}, ...
                'color', {[1 0.75 0.3], [0.9 0.4 0.9]});
end
```

`plot_sleigh.m`:

```matlab
function files = plot_sleigh(t, Y, out, p, info, outDir)
%% Purpose:
%
%   The static diagnostics: a path figure (center-of-mass and blade
%   traces, heading ticks every 0.5 s, start and end markers, equal
%   aspect) and a 2 x 3 panel (v(t) with the exact solution, the
%   asymptote and lambda(t) on a right axis; w(t) with the exact
%   solution; th(t); energy drift; constraint residual |A qd| from the
%   center-of-mass description; the (v, w) phase portrait on its energy
%   ellipse). Two PNGs; headless by default.
%
%% Inputs:
%
%  t                        [1 x N]                 Times (s)
%
%  Y                        [5 x N]                 Reduced solution
%
%  out                      struct                  From sleigh_quantities
%
%  p                        struct                  From sleigh_params
%
%  info                     struct                  .label; .visible
%                                                   [false]; .resid [1 x
%                                                   N] constraint residual
%                                                   (optional; computed
%                                                   from out.Yfull when
%                                                   absent)
%
%  outDir                   char                    Folder for the PNGs
%
%% Outputs:
%
%  files                    cell [1 x 2]            The PNG paths
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the defaults, three seconds, into tempdir
     p = sleigh_params();  t = 0:0.01:3;
     [~, yy] = ode113(@(a, b) sleigh_rhs(a, b, p), t, [0; 0; 0; 1; 2], odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
     out = sleigh_quantities(t, yy.', p);
     files = plot_sleigh(t, yy.', out, p, struct('label', 'demo', 'visible', true), tempdir);
     fprintf('wrote %s and %s\n', files{1}, files{2});
     return
end
if ~isfield(info, 'visible') || isempty(info.visible), info.visible = false; end
if info.visible, vis = 'on'; else, vis = 'off'; end
th = Y(3, :);  v = Y(4, :);  w = Y(5, :);
if isfield(info, 'resid') && ~isempty(info.resid)
    resid = info.resid;
else
    vc = out.Yfull(4:5, :);
    resid = abs(-vc(1, :).*sin(th) + vc(2, :).*cos(th) - p.a*w);
end
E0 = out.E(1);
drift = (out.E - E0)/max(abs(E0), eps);

% 1. the path
f1 = figure('Visible', vis, 'Color', 'w', 'Position', [100 100 900 700]);
ax = axes(f1);  hold(ax, 'on');
plot(ax, out.rc(1, :), out.rc(2, :), '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.8);
plot(ax, out.rb(1, :), out.rb(2, :), '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
tick = 0:0.5:t(end);
kt = arrayfun(@(tt) find(t >= tt, 1), tick);
Lb = max(abs(p.a), 0.1);
quiver(ax, out.rb(1, kt), out.rb(2, kt), Lb*cos(th(kt)), Lb*sin(th(kt)), 0, 'Color', [0.4 0.4 0.4], 'MaxHeadSize', 0.8);
plot(ax, out.rb(1, 1), out.rb(2, 1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(ax, out.rb(1, end), out.rb(2, end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
axis(ax, 'equal');  grid(ax, 'on');
xlabel(ax, 'x  [m]');  ylabel(ax, 'y  [m]');
title(ax, sprintf('%s: a = %+.3g m, v0 = %.3g m/s, w0 = %.3g rad/s (ticks every 0.5 s)', info.label, p.a, v(1), w(1)), 'Interpreter', 'none');
legend(ax, {'center of mass', 'blade', 'heading', 'start', 'end'}, 'Location', 'best');
files{1} = fullfile(outDir, [info.label '_path.png']);
exportgraphics(f1, files{1}, 'Resolution', 150);
close(f1);

% 2. the diagnostics
f2 = figure('Visible', vis, 'Color', 'w', 'Position', [100 100 1400 800]);
tl = tiledlayout(f2, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile(tl);  hold(ax, 'on');
yyaxis(ax, 'left');
plot(ax, t, v, 'b-', 'LineWidth', 1.6);  plot(ax, t, out.vExact, 'k--');
yline(ax, out.asymptote, ':', 'Color', [0.3 0.3 0.3]);
ylabel(ax, 'v  [m/s]');
yyaxis(ax, 'right');  plot(ax, t, out.lambda, '-', 'Color', [0.9 0.4 0.9]);  ylabel(ax, '\lambda  [N]');
xlabel(ax, 't  [s]');  grid(ax, 'on');  title(ax, 'forward speed (exact dashed, asymptote dotted), blade force');
ax = nexttile(tl);  hold(ax, 'on');
plot(ax, t, w, 'b-', 'LineWidth', 1.6);  plot(ax, t, out.wExact, 'k--');
xlabel(ax, 't  [s]');  ylabel(ax, '\omega  [rad/s]');  grid(ax, 'on');  title(ax, 'yaw rate (exact dashed)');
ax = nexttile(tl);
plot(ax, t, th, 'b-', 'LineWidth', 1.6);
xlabel(ax, 't  [s]');  ylabel(ax, '\theta  [rad]');  grid(ax, 'on');  title(ax, 'heading');
ax = nexttile(tl);
plot(ax, t, drift, 'b-');
xlabel(ax, 't  [s]');  ylabel(ax, '(E - E_0)/E_0');  grid(ax, 'on');  title(ax, sprintf('energy drift, max %.1e', max(abs(drift))));
ax = nexttile(tl);
semilogy(ax, t, max(resid, 1e-20), 'b-');
xlabel(ax, 't  [s]');  ylabel(ax, '|A \dot{q}|');  grid(ax, 'on');  title(ax, sprintf('constraint residual, max %.1e', max(resid)));
ax = nexttile(tl);  hold(ax, 'on');
sC = sqrt(max(out.C, 0));  s = linspace(0, 2*pi, 361);
plot(ax, sC*cos(s), sC*sqrt(p.m/p.Jb)*sin(s), '-', 'Color', [0.6 0.6 0.6]);
plot(ax, v, w, 'b-', 'LineWidth', 1.6);
plot(ax, v(1), w(1), 'go', 'MarkerFaceColor', 'g');  plot(ax, v(end), w(end), 'rs', 'MarkerFaceColor', 'r');
axis(ax, 'equal');  grid(ax, 'on');
xlabel(ax, 'v  [m/s]');  ylabel(ax, '\omega  [rad/s]');  title(ax, 'phase portrait on the energy ellipse');
title(tl, sprintf('%s: m = %.3g kg, J = %.3g kg m^2, a = %+.3g m', info.label, p.m, p.J, p.a), 'Interpreter', 'none');
files{2} = fullfile(outDir, [info.label '_diag.png']);
exportgraphics(f2, files{2}, 'Resolution', 150);
close(f2);
end
```

- [ ] **Step 4: Write the script**

`run_sleigh.m`:

```matlab
%% Purpose:
%
%   THE front door for the Chaplygin sleigh. Set the geometry and the
%   initial velocity in section 1, run, and get: the integrated motion
%   (ode113 at 1e-12), five printed diagnostics with PASS/FAIL (hand-
%   derived equations against the symbolic Lagrange-d'Alembert builder,
%   the exact tanh solution, energy, the constraint, the stability
%   statement), two PNG figures and a 1280 x 720 MP4 in results/ beside
%   this script. The derivation and the conventions are in
%   doc/chaplygin_sleigh_note.pdf; the sign of a is the whole story:
%   a > 0 puts the center of mass AHEAD of the blade and straight motion
%   is then asymptotically stable, a < 0 (or v0 < 0) makes it unstable and
%   the sleigh turns around.
%
%% What the user sets (section 1):
%
%   label, m, J, a, v0, w0, tEnd, dtOut, movieOpts.
%
%% What it writes:
%
%   results/<label>_path.png, results/<label>_diag.png, results/<label>.mp4
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% 1. USER SETTINGS
label = 'sleigh_forward';
m     = 1.0;     % mass (kg)
J     = 0.1;     % inertia about the center of mass (kg m^2)
a     = 0.2;     % blade to center of mass (m); a > 0 = center of mass ahead of the blade
v0    = 1.0;     % initial forward speed along the blade (m/s)
w0    = 2.0;     % initial yaw rate (rad/s)
tEnd  = 10;      % final time (s)
dtOut = 0.01;    % output step (s)
movieOpts = struct('fps', 30, 'seconds', tEnd);

%% 2. SETUP
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'mechlib'));
outDir = fullfile(here, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
assert(tEnd > 0 && dtOut > 0 && dtOut <= tEnd, 'run_sleigh:time', 'need 0 < dtOut <= tEnd');
p = sleigh_params(m, J, a);
C = v0^2 + (p.Jb/p.m)*w0^2;
fprintf('sleigh: m %.3g kg, J %.3g kg m^2, a %+.3g m -> Jb %.4g, k %.4g 1/m; C %.6g, v -> %+.6g m/s\n', ...
        m, J, a, p.Jb, p.k, C, sign(a)*sqrt(C));

%% 3. INTEGRATE (blade-point state [x_b; y_b; th; v; w], blade at the origin heading east)
t  = 0:dtOut:tEnd;
y0 = [0; 0; 0; v0; w0];
[~, yy] = ode113(@(tt, y) sleigh_rhs(tt, y, p), t, y0, odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
Y   = yy.';
out = sleigh_quantities(t, Y, p);
fprintf('integrated %d steps to %.2f s: v(end) %+.6f m/s, w(end) %+.3e rad/s, path length %.3f m\n', ...
        numel(t) - 1, t(end), Y(4, end), Y(5, end), sum(sqrt(sum(diff(out.rb, 1, 2).^2, 1))));

%% 4. THE SYMBOLIC BUILDER (the independent derivation)
[q, qd, L, A] = sleigh_lagrangian(p);
sys = mech.lagrange_dalembert(q, qd, L, A);
fprintf('builder: reduced accelerations and lambda formed symbolically (n = %d, k = %d)\n', sys.n, sys.k);

%% 5. DIAGNOSTICS
d1 = 0;
for kk = 1:numel(t)
    d1 = max(d1, max(abs(sys.f(t(kk), out.Yfull(:, kk)) - out.Yfulld(:, kk))));
end
verdict('D1 hand-derived vs builder', d1, 1e-9);
d2 = max([abs(Y(4, :) - out.vExact), abs(Y(5, :) - out.wExact)]);
verdict('D2 exact solution (v, w)', d2, 1e-8);
chk = mech.conservation_check(t, out.Yfull, sys);
verdict('D3 energy drift (rel)', chk.driftRel, 1e-9);
verdict('D4 constraint residual |A qd|', chk.resid, 1e-9);
leads = a*Y(4, end) > 0;
settled = abs(Y(5, end)) <= 1e-3*max(abs(w0), eps);
if a == 0
    fprintf('  NOTE  D5 a = 0: the knife edge; v and w stay constant, the blade runs a circle of radius %.4g m\n', abs(v0/max(abs(w0), eps)));
elseif w0 == 0 && a*v0 > 0
    fprintf('  NOTE  D5 w0 = 0: a straight ride at v = %.3g m/s; with the center of mass leading (a v > 0) it is stable to yaw perturbations\n', v0);
elseif w0 == 0
    fprintf('  NOTE  D5 w0 = 0: a straight ride at v = %.3g m/s with the center of mass TRAILING (a v < 0): an exact but unstable equilibrium, any yaw would grow and reverse it\n', v0);
elseif leads && settled
    fprintf('  PASS  D5 asymptotics: v -> sign(a) sqrt(C) = %+.6f (v(end) %+.6f), |w(end)| = %.1e, heading change %.4f rad vs exact %.4f; straight motion with the center of mass leading is stable to yaw (one of a line of equilibria: it attracts transversely, on its energy ellipse)\n', ...
            out.asymptote, Y(4, end), abs(Y(5, end)), Y(3, end) - Y(3, 1), out.thetaInf - Y(3, 1));
elseif leads
    fprintf('  NOTE  D5 asymptotics: center of mass leads at the end but |w(end)| = %.1e has not settled below 1e-3 |w0| by t = %.3g s (lengthen tEnd)\n', abs(Y(5, end)), tEnd);
else
    fprintf('  NOTE  D5 asymptotics: the center of mass still trails (a v < 0) at t = %.3g s; this is the unstable branch, |w| grows until the signed speed crosses zero (at t = %.3g s) and the sleigh then departs along sign(a) (lengthen tEnd to see it)\n', ...
            tEnd, -atanh(v0/sqrt(C))/(p.k*sqrt(C)));
end

%% 6. PLOTS AND MOVIE
files = plot_sleigh(t, Y, out, p, struct('label', label, 'visible', false, 'resid', chk.residSeries), outDir);
fprintf('wrote %s\n       %s\n', files{1}, files{2});
body = sleigh_body(p);
body.arrows   = sleigh_arrows(t, Y, out, p);
lamI = @(tk) interp1(t, out.lambda, tk, 'pchip');
vI   = @(tk) interp1(t, Y(4, :), tk, 'pchip');
wI   = @(tk) interp1(t, Y(5, :), tk, 'pchip');
EI   = @(tk) interp1(t, (out.E - out.E(1))/max(abs(out.E(1)), eps), tk, 'pchip');
body.subtitle = @(tk) sprintf('t =%6.2f s   v =%7.3f m/s   w =%7.3f rad/s   lambda =%7.3f N   dE/E0 = %+.1e', tk, vI(tk), wI(tk), lamI(tk), EI(tk));
mfile = mech.animate_planar_body(t, [out.rb; Y(3, :)], body, fullfile(outDir, label), ...
                                 struct('fps', movieOpts.fps, 'seconds', movieOpts.seconds, ...
                                        'traces', cat(3, out.rc, out.rb), 'traceColors', [0.85 0.33 0.1; 0.45 0.70 1.00], ...
                                        'title', sprintf('Chaplygin sleigh: a = %+.3g m, v0 = %.3g m/s, w0 = %.3g rad/s', a, v0, w0)));
fprintf('wrote %s\n', mfile);

% ------------------------------------------------------------------------
function verdict(name, value, bar)
%% Purpose:
%
%   One PASS/FAIL line with the measured value and its bar.
%
if value <= bar, s = 'PASS'; else, s = 'FAIL'; end
fprintf('  %s  %s: %.2e (bar %.0e)\n', s, name, value, bar);
end
```

(Scripts may define local functions at the end of the file since R2016b; `verdict` lives there.)

- [ ] **Step 5: README and .gitignore**

`.gitignore`:
```
results/
```

`README.md`:

```markdown
# Chaplygin sleigh

A planar body with a knife-edge blade a signed distance `a` from its center
of mass: the blade forbids sideways motion, a nonholonomic constraint. The
first system of the `lagrangianHamiltonianMechanics` series; the
derivation, conventions and the exact solution are in
`doc/chaplygin_sleigh_note.pdf`.

## Run it

Open `run_sleigh.m`, edit section 1 (`m`, `J`, `a`, `v0`, `w0`, `tEnd`,
`dtOut`, movie options), run. It prints five diagnostics with PASS/FAIL,
writes `results/<label>_path.png`, `results/<label>_diag.png` and
`results/<label>.mp4` (1280 x 720). `results/` is git-ignored.

Try `a = -0.2` (or `v0 = -1`): the center of mass trails the blade, the
straight ride is unstable to yaw, |w| grows until the signed speed crosses
zero, and the sleigh departs the other way. With `w0 = 0` nothing happens
in either case: a straight ride is an exact equilibrium, stable or not.

## What is checked, every run

| line | what | reference |
|---|---|---|
| D1 | the hand-derived reduced equations (`sleigh_rhs`) against the symbolic Lagrange-d'Alembert builder (`mech.lagrange_dalembert` on `sleigh_lagrangian`) | independent derivation |
| D2 | v(t) and w(t) against the exact solution v = S tanh z, w = sign(w0) (S/b) sech z, z = asinh(v0/(b abs(w0))) + k S t | closed form |
| D3 | energy E = m/2 v^2 + (J + m a^2)/2 w^2 | conserved |
| D4 | constraint residual A qd on the reconstructed center-of-mass state | a consistency check of the mapping (it is zero by construction); D1-D3 are the independent evidence |
| D5 | the asymptote sign(a) sqrt(C), the decay of w, the total heading change sign(w0) (b/abs(a)) acos(sign(a) v0/S) | stability statement |

Straight motion is a LINE of equilibria (v, 0): the one with the center of
mass leading (a v > 0) is Lyapunov stable and attracts transversely (on its
energy ellipse it is asymptotically stable); the trailing one (a v < 0) is
unstable to yaw. No energy is dissipated: the constraint force does no
work, yaw energy is traded into translation (and the other way first, on a
trailing launch).

## Conventions (binding)

World x east, y north, theta counterclockwise from x; e1 along the blade
(forward), e2 to the left. q = (x_c, y_c, theta) of the center of mass;
blade point r_b = r_c - a e1; `a > 0` = center of mass ahead of the blade.
Reduced state integrated: `[x_b; y_b; theta; v; w]` of the blade point.

## Files

| file | what |
|---|---|
| `run_sleigh.m` | the front door |
| `sleigh_params.m` | m, J, a and the derived Jb, k |
| `sleigh_lagrangian.m` | symbolic L and A for the builder |
| `sleigh_rhs.m` | the hand-derived reduced equations |
| `sleigh_quantities.m` | paths, energy, blade force, exact solution, full-state arrays |
| `sleigh_body.m`, `sleigh_arrows.m` | the drawing and the movie's arrows |
| `plot_sleigh.m` | the two static figures |
| `tests/test_sleigh_oracle.m` | D1-D5 for both signs of a and v0, the knife edge, the straight ride, rest |
| `tests/test_sleigh_pipeline.m` | plots, movie and the script's parse, headless |
| `doc/chaplygin_sleigh_note.tex` | the derivation |

Tests run through `run_tests` in `lagrangianHamiltonianMechanics/`.
```

- [ ] **Step 6: Run the test, the runner, and the script once**

```matlab
cd('/Users/msc/Desktop/optimal_control/lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/tests'); test_sleigh_pipeline
cd('/Users/msc/Desktop/optimal_control/lagrangianHamiltonianMechanics'); run_tests
```
Expected: all `PASS`; `5 of 5 passed`. Then run the script once with the defaults through `mcp__matlab__run_matlab_file` on `run_sleigh.m`; expected D1-D4 PASS and a D5 PASS line (with a = 0.2, v0 = 1, w0 = 2 and tEnd = 10, w decays below 2e-3 by t = 10: k sqrt(C) = (0.2/0.14) sqrt(1 + 0.14*4) ~ 1.78 /s, so sech falls to 1e-3 in about 4.3 s). Paste the printed block and the three output paths into the report. Look at `results/sleigh_forward_path.png` with the Read tool and confirm the path curls then straightens; look at one movie frame through VideoReader if in doubt. Do not commit.

---

### Task 6: the LaTeX note

**Files:**
- Create: `lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/doc/chaplygin_sleigh_note.tex`
- Create: `lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/doc/figures/sleigh_forward_path.png`, `sleigh_forward_diag.png`, `sleigh_backward_path.png` (copied from `results/` after running the script twice: defaults, then `label = 'sleigh_backward'`, `v0 = -1`)
- Create: `.../doc/chaplygin_sleigh_note.pdf` (compiled, committed)

**Interfaces:**
- Consumes: `sys.qddSym`, `sys.lambdaSym` printed through `latex()` in the shared session for section 3's cross-check; the D1-D5 lines from Task 5's run.

- [ ] **Step 1: Generate the figures and the symbolic strings**

Run `run_sleigh.m` with the defaults, then with `label = 'sleigh_backward'; v0 = -1;` (restore the file to the defaults afterwards; `git diff` on it must be empty). Copy the three PNGs into `doc/figures/`. In the session:
```matlab
p = sleigh_params(); [q, qd, L, A] = sleigh_lagrangian(p); sys = mech.lagrange_dalembert(q, qd, L, A);
disp(latex(sys.qddSym)); disp(latex(sys.lambdaSym));
```
and paste the two strings into section 3 (edit them to the note's symbols: `\dot x_c` etc., keeping the structure the builder produced).

- [ ] **Step 2: Write the note**

`doc/chaplygin_sleigh_note.tex`:

```latex
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb,bm}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{hyperref}
\usepackage{tikz}
\newcommand{\e}{\mathbf{e}}
\newcommand{\rv}{\mathbf{r}}
\title{The Chaplygin Sleigh:\\ Lagrange--d'Alembert Equations, Reduction, and the Exact Solution}
\author{Michael Casey}
\date{October 25, 2026}
\begin{document}
\maketitle

\begin{abstract}
A planar rigid body rides on a knife-edge blade offset from its center of mass. The blade cannot slip sideways: a constraint linear in velocity that fails the Frobenius integrability condition, so the system is nonholonomic and the usual Lagrange equations do not apply. This note derives the equations of motion by the Lagrange--d'Alembert principle with a multiplier, reduces them to two first-order equations in the blade's forward speed and the yaw rate, solves those in closed form (speed, yaw rate and the total heading change), and shows that among the straight motions, which form a line of equilibria, those with the center of mass leading the blade are stable to yaw and attract their neighbours while those with it trailing are unstable. The MATLAB script \texttt{run\_sleigh.m} checks every step numerically.
\end{abstract}

\section{Setup and conventions}
World frame: $x$ east, $y$ north, heading $\theta$ counterclockwise from $x$. Body axes
\[
\e_1 = (\cos\theta,\ \sin\theta),\qquad \e_2 = (-\sin\theta,\ \cos\theta),\qquad \dot\e_1 = \omega\,\e_2,\quad \dot\e_2 = -\omega\,\e_1,\quad \omega = \dot\theta .
\]
The body has mass $m$ and inertia $J$ about its center of mass $\rv_c = (x_c, y_c)$. The blade touches the plane at
\[
\rv_b = \rv_c - a\,\e_1 ,
\]
so $a > 0$ puts the center of mass \emph{ahead} of the blade in the forward direction (Figure~\ref{fig:body}). The blade is a knife edge aligned with $\e_1$: it may slide along $\e_1$ and turn, but it cannot move along $\e_2$.

\begin{figure}[h]
\centering
\begin{tikzpicture}[scale=2.2]
  \draw[->, gray] (-0.6,0) -- (1.4,0) node[right] {$x$};
  \draw[->, gray] (0,-0.5) -- (0,0.8) node[above] {$y$};
  \begin{scope}[rotate=25]
    \draw[thick, fill=orange!20] (0.7,0) -- (0.55,0.18) -- (-0.4,0.12) -- (-0.4,-0.12) -- (0.55,-0.18) -- cycle;
    \draw[very thick, cyan!70!black] (-0.2,0) -- (0.2,0);
    \fill (0,0) circle (0.025) node[below left] {$\rv_b$};
    \fill[red] (0.35,0) circle (0.03) node[above] {$\rv_c$};
    \draw[->, thick] (0.35,0) -- (0.85,0) node[right] {$\e_1$};
    \draw[->, thick] (0.35,0) -- (0.35,0.45) node[above] {$\e_2$};
    \draw[<->] (0,-0.28) -- (0.35,-0.28) node[midway, below] {$a$};
  \end{scope}
\end{tikzpicture}
\caption{The sleigh. The blade (cyan) sits at $\rv_b$ and is aligned with $\e_1$; the center of mass (red) is $a\,\e_1$ ahead of it.}
\label{fig:body}
\end{figure}

\section{Lagrange--d'Alembert with the multiplier}
Generalized coordinates $q = (x_c, y_c, \theta)$. There is no potential, so
\[
L = \tfrac{m}{2}\left(\dot x_c^2 + \dot y_c^2\right) + \tfrac{J}{2}\dot\theta^2,\qquad M = \mathrm{diag}(m, m, J).
\]
The blade's velocity is $\dot\rv_b = \dot\rv_c - a\omega\,\e_2$, and the knife edge demands $\dot\rv_b\cdot\e_2 = 0$:
\begin{equation}
A(q)\,\dot q = 0,\qquad A = \begin{bmatrix} -\sin\theta & \cos\theta & -a \end{bmatrix}.
\label{eq:constraint}
\end{equation}
Is this constraint holonomic in disguise? A one-form $\alpha$ defines a holonomic constraint when some integrating factor makes it exact, and Frobenius' theorem says that happens if and only if $\alpha\wedge d\alpha = 0$. Here $\alpha = -\sin\theta\,dx_c + \cos\theta\,dy_c - a\,d\theta$, $d\alpha = -\cos\theta\,d\theta\wedge dx_c - \sin\theta\,d\theta\wedge dy_c$, and
\[
\alpha\wedge d\alpha = -\,dx_c\wedge dy_c\wedge d\theta \neq 0
\]
for every $a$, including $a = 0$. So no coordinate change turns \eqref{eq:constraint} into a condition on positions: the sleigh has three configuration degrees of freedom but only two velocity degrees of freedom. (Non-closedness alone would not prove this: $e^x\,dy$ is not closed, yet its zero set is $\dot y = 0$.) The Lagrange--d'Alembert principle states that the constraint force does no work on admissible velocities, which puts it in the row space of $A$:
\begin{equation}
\frac{d}{dt}\frac{\partial L}{\partial \dot q} - \frac{\partial L}{\partial q} = A^{\!\top}\lambda
\quad\Longrightarrow\quad
m\ddot x_c = -\lambda\sin\theta,\qquad m\ddot y_c = \lambda\cos\theta,\qquad J\ddot\theta = -a\lambda .
\label{eq:LdA}
\end{equation}
$\lambda$ is the sideways force the ice exerts on the blade, along $\e_2$, and $-a\lambda$ is its moment about the center of mass. In general, with $p = \partial L/\partial\dot q$, $M = \partial p/\partial\dot q$ and $h = (\partial p/\partial q)\,\dot q - \partial L/\partial q$ (the Euler--Lagrange terms left after removing $M\ddot q$, for an autonomous $L$), and with generalized forces $Q$, the equations are $M\ddot q + h = Q + A^{\!\top}\lambda$. Together with the time derivative of the constraint, $A\ddot q + \dot A\dot q = 0$, the unknowns $(\ddot q, \lambda)$ satisfy the block system
\[
\begin{bmatrix} M & A^{\!\top} \\ A & 0 \end{bmatrix}
\begin{bmatrix} \ddot q \\ -\lambda \end{bmatrix} =
\begin{bmatrix} Q - h \\ -\dot A\dot q \end{bmatrix},
\]
which is what \texttt{mech.lagrange\_dalembert} solves symbolically for any $L$, $A$ and $Q$ with $M$ nonsingular and $A$ of full row rank. For the sleigh $Q = 0$ and $h = 0$ (Cartesian coordinates, no potential), and the top block reduces to \eqref{eq:LdA}.

\section{Reduction to the blade's speed and the yaw rate}
Define the quasi-velocities $v = \dot\rv_b\cdot\e_1 = \dot x_c\cos\theta + \dot y_c\sin\theta$ (the forward speed, the same for blade and center of mass) and $\omega = \dot\theta$. The constraint gives $\dot\rv_c = v\,\e_1 + a\omega\,\e_2$, and differentiating,
\[
\ddot\rv_c = (\dot v - a\omega^2)\,\e_1 + (v\omega + a\dot\omega)\,\e_2 .
\]
Comparing with $m\ddot\rv_c = \lambda\,\e_2$ from \eqref{eq:LdA}: the $\e_1$ component gives $\dot v = a\omega^2$ and the $\e_2$ component $\lambda = m(v\omega + a\dot\omega)$. Substituting into $J\dot\omega = -a\lambda$ and collecting $\dot\omega$,
\begin{equation}
\boxed{\;\dot v = a\,\omega^2,\qquad \dot\omega = -k\,v\,\omega,\qquad k = \frac{m a}{J + m a^2},\qquad \lambda = \frac{m J}{J + m a^2}\,v\,\omega\;}
\label{eq:reduced}
\end{equation}
with the kinematics $\dot x_b = v\cos\theta$, $\dot y_b = v\sin\theta$, $\dot\theta = \omega$. The script integrates this five-state system (\texttt{sleigh\_rhs}). The symbolic builder, given only $L$ and $A$, returns for the defaults $m = 1$, $J = 0.1$, $a = 0.2$
\[
\ddot q = \text{(paste the \texttt{latex(sys.qddSym)} string here)},\qquad
\lambda = \text{(paste the \texttt{latex(sys.lambdaSym)} string here)},
\]
which the script compares pointwise with \eqref{eq:reduced} along the solution (line D1).

\section{Energy, the invariant ellipse, and the exact solution}
The kinetic energy about the center of mass, with $|\dot\rv_c|^2 = v^2 + a^2\omega^2$, is
\begin{equation}
E = \tfrac{m}{2}v^2 + \tfrac{1}{2}(J + m a^2)\,\omega^2 ,
\label{eq:energy}
\end{equation}
and \eqref{eq:reduced} conserves it: $\dot E = m v\,a\omega^2 - (J + ma^2)\,k\,v\,\omega^2 = 0$ since $(J + ma^2)k = ma$. So the motion in the $(v, \omega)$ plane stays on the ellipse
\[
v^2 + \frac{J + m a^2}{m}\,\omega^2 = C,\qquad C = v_0^2 + \frac{J + m a^2}{m}\,\omega_0^2 .
\]
Eliminating $\omega^2$ from $\dot v = a\omega^2$ gives $\dot v = k\,(C - v^2)$, a Riccati equation. For $a \neq 0$ and $\omega_0 \neq 0$, write $b = \sqrt{(J + ma^2)/m}$, $S = \sqrt{C}$ and
\begin{equation}
z(t) = \alpha + kS\,t,\qquad \alpha = \operatorname{arsinh}\!\frac{v_0}{b\,|\omega_0|},\qquad
v(t) = S\tanh z,\qquad \omega(t) = \operatorname{sign}(\omega_0)\,\frac{S}{b}\,\operatorname{sech} z .
\label{eq:exact}
\end{equation}
At $t = 0$ this returns $v_0$ and $\omega_0$ (since $\tanh\alpha = v_0/S$ and $\cosh\alpha = S/(b|\omega_0|)$), and it is the form to evaluate numerically: the equivalent $v = S\tanh(kS(t - t_0))$ with $t_0 = -\operatorname{artanh}(v_0/S)/(kS)$ fails for tiny $\omega_0$, where $C$ rounds to $v_0^2$ and $\operatorname{artanh}(\pm 1)$ is infinite. $\omega$ never changes sign (its equation is linear in $\omega$; equivalently $\omega = \omega_0\exp(-k\int_0^t v)$), and $v \to \operatorname{sign}(a)\,S$. The heading also integrates exactly, since $\int\operatorname{sech} z\,dz$ is the Gudermannian:
\begin{equation}
\theta_\infty - \theta_0 = \operatorname{sign}(\omega_0)\,\frac{b}{|a|}\,\arccos\!\frac{\operatorname{sign}(a)\,v_0}{S},
\label{eq:heading}
\end{equation}
which can exceed several revolutions when $|a|$ is small. Two special cases fall outside \eqref{eq:exact}: $a = 0$ is the knife edge, $k = 0$, $v$ and $\omega$ constant and the blade runs a circle of radius $|v_0/\omega_0|$ with $\lambda = m v\omega$ the centripetal force; $\omega_0 = 0$ is the straight ride at constant $v_0$ for every $a$, whether or not it is stable. Line D2 of the script measures the integrator against \eqref{eq:exact}, D3 against \eqref{eq:energy}, D5 against \eqref{eq:heading}; D4 evaluates \eqref{eq:constraint} on the reconstructed center-of-mass velocity $\dot\rv_c = v\e_1 + a\omega\e_2$, where it vanishes identically, so D4 is a consistency check of that mapping rather than independent evidence.

\section{Stability of straight motion}
Straight motion ($\omega = 0$, $v = V$) is a whole line of equilibria of \eqref{eq:reduced}, one for each speed. Linearizing at $(V, 0)$,
\[
\begin{bmatrix}\delta\dot v\\ \delta\dot\omega\end{bmatrix} = \begin{bmatrix} 0 & 0\\ 0 & -kV\end{bmatrix}\begin{bmatrix}\delta v\\ \delta\omega\end{bmatrix},
\]
and the zero eigenvalue is real: a perturbation $\delta v$ along the line is itself another equilibrium and never decays, so no single straight motion is asymptotically stable in the $(v, \omega)$ plane. The transverse eigenvalue $-kV$ is what the sign of $a$ decides. When $aV > 0$, the center of mass ahead of the blade in the direction of travel, a yaw perturbation decays and the motion converges to a nearby straight ride at the speed its energy allows: the equilibrium is Lyapunov stable and transversely attracting, and on its own energy ellipse it is asymptotically stable. When $aV < 0$ it is unstable to yaw. The exact solution says the same globally: for $a \neq 0$ and any $\omega_0 \neq 0$, however small, the sleigh ends up moving at $\operatorname{sign}(a)\,S$. A sleigh launched with the center of mass trailing therefore sees $|\omega|$ \emph{grow} (to $S/b$ at the moment $v = 0$), its signed speed cross zero exactly once, at $t_{\rm rev} = -\operatorname{artanh}(v_0/S)/(kS)$, and then departs along the stable direction, turning in total by \eqref{eq:heading} (Figure~\ref{fig:paths}). With $\omega_0 = 0$ exactly, nothing happens: the trailing straight ride is an exact, unstable equilibrium.

No energy is dissipated in this model. The constraint force does no work in total, $\lambda\e_2\cdot\dot\rv_c + (-a\lambda)\omega = a\lambda\omega - a\lambda\omega = 0$, and the kinetic energy \eqref{eq:energy} is conserved. What changes is its split: the rotational part $\tfrac{J}{2}\omega^2$ obeys $\tfrac{d}{dt}\tfrac{J}{2}\omega^2 = -Jk\,v\,\omega^2$, so it decreases when $av > 0$ and increases when $av < 0$, the translational part changing by the opposite amount. Yaw is traded into forward motion through the constraint, reversibly, not dissipated.

\begin{figure}[h]
\centering
\includegraphics[width=0.48\textwidth]{figures/sleigh_forward_path.png}\hfill
\includegraphics[width=0.48\textwidth]{figures/sleigh_backward_path.png}
\caption{Paths of the center of mass (orange) and the blade (blue) for $a = 0.2$~m, $\omega_0 = 2$~rad/s and $v_0 = +1$~m/s (left, the center of mass leads: the heading turns by $1.20$~rad in total and the ride straightens) and $v_0 = -1$~m/s (right, it trails: $|\omega|$ grows, the signed speed crosses zero at $0.617$~s, and the sleigh departs the other way after turning $4.68$~rad in total). Grey ticks mark the heading every 0.5~s.}
\label{fig:paths}
\end{figure}

\section{Numerical check}
Figure~\ref{fig:diag} is the diagnostics panel of \texttt{run\_sleigh.m} for the defaults. The printed verdicts of that run were (paste the five D lines of the run here, verbatim):
\begin{verbatim}
  PASS  D1 hand-derived vs builder: ...
  PASS  D2 exact solution (v, w): ...
  PASS  D3 energy drift (rel): ...
  PASS  D4 constraint residual |A qd|: ...
  PASS  D5 asymptotics: ...
\end{verbatim}

\begin{figure}[h]
\centering
\includegraphics[width=\textwidth]{figures/sleigh_forward_diag.png}
\caption{Forward speed with the exact solution and the asymptote, yaw rate, heading, energy drift, constraint residual, and the $(v,\omega)$ phase portrait on its energy ellipse.}
\label{fig:diag}
\end{figure}

\section*{References}
\begin{enumerate}
\item A. M. Bloch, \emph{Nonholonomic Mechanics and Control}, 2nd ed., Springer, 2015, Sections 1.4 and 1.7.
\item Ju. I. Neimark and N. A. Fufaev, \emph{Dynamics of Nonholonomic Systems}, AMS Translations of Mathematical Monographs 33, 1972.
\item S. A. Chaplygin, ``On the theory of motion of nonholonomic systems. The reducing-multiplier theorem,'' \emph{Mat. Sbornik} 28, 1911 (English translation in \emph{Regular and Chaotic Dynamics} 13, 2008).
\end{enumerate}
\end{document}
```

Replace the two "paste" placeholders in section 3 and the verbatim block in section 6 with the real strings and lines from Step 1 and Task 5's run before compiling; the task is not done while either placeholder text remains.

- [ ] **Step 3: Compile and clean**

```bash
cd lagrangianHamiltonianMechanics/nonholonomic/ChaplyginSleigh/doc && /Library/TeX/texbin/pdflatex -interaction=nonstopmode chaplygin_sleigh_note.tex > /dev/null && /Library/TeX/texbin/pdflatex -interaction=nonstopmode chaplygin_sleigh_note.tex | grep -i "warning\|error" ; rm -f *.aux *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz; ls -la chaplygin_sleigh_note.pdf
```
Expected: a PDF of about six pages, no undefined references, no overfull-box warnings above 10 pt, no leftover aux files. Read the PDF (Read tool, pages 1-6) and confirm the two paste placeholders are gone and both figures render.

- [ ] **Step 4: Report**

Status, the page count, the figure list, and any TeX warnings. Do not commit.

---

## Self-review

**Spec coverage.** Section 2 layout: Tasks 1 (mechlib scaffold, runner, both top READMEs), 4-5 (sleigh files, README, .gitignore), 6 (doc). Section 3 physics: sleigh_rhs, sleigh_quantities, the test's checks on k's sign, the circle, the straight ride, rest. Section 4.1-4.3: Tasks 1-3 with every field and identifier named. Section 5 modules: all seven sleigh files (sleigh_arrows added as the home of the movie arrows so the script and the pipeline test share them; recorded as an addition to the spec's file table). 5.1 script: section 1 verbatim, pipeline steps 1-5, D1-D5 bars, D5 as a verdict not a gate. 5.2 tests: five test functions; the spec's "lists all seven" was a miscount (it counted the two mechlib self-demos), the runner prints `5 of 5`. 5.3 note: all seven sections, three references, compiled and committed, no signature line. Section 6 house style: in Global Constraints.

**Placeholder scan.** The note carries two deliberate "paste here" spots and a verbatim block that Task 6 must fill from real output; the step says the task is not done while they remain. No other placeholders.

**Type consistency.** `sys` fields (`f`, `lambda`, `energy`, `constraint`, `qddSym`, `lambdaSym`, `n`, `k`) are used identically in Tasks 2, 4, 5. `out` fields (`rb`, `rc`, `E`, `lambda`, `C`, `asymptote`, `vExact`, `wExact`, `Yfull`, `Yfulld`) match between sleigh_quantities, the oracle test, plot_sleigh, sleigh_arrows and the script. `body` fields (`outline`, `markers`, `segments`, `arrows`, `subtitle`) match between the animator, its test, sleigh_body and the script. `chk` fields (`ok`, `driftRel`, `resid`, `residSeries`, `why`, `okEnergy`, `okConstraint`) match between Task 2 and its consumers.

**Review Focus.** 1 (a = 0) and 2 (rest) are in test_sleigh_oracle and test_sleigh_pipeline; 3 (unsorted t) in test_animate_planar_body; 4 (transposed Y) in test_conservation_check; 5 (velocity-dependent Q) in test_lagrange_dalembert.
