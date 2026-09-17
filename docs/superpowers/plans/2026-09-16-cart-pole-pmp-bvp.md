# Cart-pole PMP-BVP Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Solve the cart-pole swing-up indirectly — form its Pontryagin boundary-value problem and shoot it with `ms_bvp` — so the engine has a consumer with no orbit in it, and can move into `oclib/+oc`.

**Architecture:** Fixed final time, minimum effort, unconstrained control, so the PMP control is the smooth `u* = -lam'G/2` and the BVP is square: four unknown initial costates against four terminal state conditions. The 8-state PMP flow is propagated with `ode113`; its 8x8 Jacobian (for the multiple-shooting STM) comes from a complex step, which is only safe because the 4x4 state Jacobian inside the costate equation is generated symbolically rather than by a nested complex step. Seeded from a stored direct solution's multipliers through `oc.duals_to_costates`.

**Tech Stack:** MATLAB R2026a, Symbolic Math Toolbox (code generation, one-off), `ode113`, `costate_common/ms_bvp` then `oclib/+oc/ms_bvp`, `oclib/+oc/duals_to_costates`, `oclib/+oc/fly_control`. No CasADi, no pumpkyn.

**Spec:** `docs/superpowers/specs/2026-09-16-cart-pole-pmp-bvp-design.md`

## Global Constraints

- Problem constants, exactly the existing example's: `m1 = 5` kg, `m2 = 1` kg, `L = 2` m, `g = 9.8` m/s^2, `t_f = 5` s, `x(0) = [0;0;0;0]`, `x(t_f) = [0;pi;0;0]`.
- House MATLAB style: `%% Purpose / Inputs / Outputs / Revision History` header quartet closed by the `Begin Code Sequence` divider; **no** `%#ok` pragmas; **never** `i` or `j` as loop variables; `nargin == 0` self-demo where a demo means anything; norms as `sqrt(sum(...))`, never `norm`, in any file that must stay complex-step safe.
- Every new function file gets a test under `optimal_control_examples/ex3_cart_pole_pmp/tests/`. Tests are functions returning `ok` (logical) and printing one `PASS`/`FAIL` line per check, following `orbit_transfer/costate_common/tests/test_ms_bvp_fixedtf.m`.
- Nothing under `optimal_control_examples/ex2_cart_pole_swing_up/` is modified.
- Run MATLAB headlessly as `/Applications/MATLAB_R2026a.app/bin/matlab -batch "run('<script>')"`; the shared MCP session is not required by any step.
- A deviation from the spec is recorded in the plan and in the code's header, never applied silently.

**Deviation already recorded (Task 3):** the spec says the PMP field's Jacobian comes from a complex step. Applied to the whole 8-state field that nests a complex step inside a complex step, which corrupts the inner derivative without any warning. The state Jacobian is therefore generated symbolically once and committed; the outer 8x8 complex step in Task 5 is then legitimate.

---

### Task 1: The direct-solution fixture

**Files:**
- Create: `optimal_control_examples/ex3_cart_pole_pmp/gen_direct_ref.m`
- Create: `optimal_control_examples/ex3_cart_pole_pmp/data/cartpole_direct_ref.mat` (produced by the above)
- Test: `optimal_control_examples/ex3_cart_pole_pmp/tests/test_direct_ref.m`

**Interfaces:**
- Consumes: nothing.
- Produces: the fixture `.mat` with variables `tN` `[1 x 201]`, `X` `[4 x 201]`, `U` `[1 x 201]`, `muDefect` `[4 x 200]`, `J` (scalar), `p` (struct with `m1 m2 L g`), `tf` (scalar), `uMax` (scalar). Tasks 6 and 7 read it.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_direct_ref()
%% Purpose:
%
%   The committed direct-solution fixture is a real solution of the problem
%   the indirect demo solves: right shape, boundary conditions met, its own
%   trapezoidal defects small, and a control that stayed inside the relaxed
%   bound (so it is the UNCONSTRAINED optimum the PMP solve is compared to).
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
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
R = load(fullfile(here, 'data', 'cartpole_direct_ref.mat'));

ok = chk(ok, isequal(size(R.X), [4 201]) && isequal(size(R.U), [1 201]) ...
             && isequal(size(R.muDefect), [4 200]) && numel(R.tN) == 201, ...
         'shapes: X [4 x 201], U [1 x 201], muDefect [4 x 200]');
ok = chk(ok, R.tN(1) == 0 && abs(R.tN(end) - 5) < 1e-12 && abs(R.tf - 5) < 1e-12, ...
         'the grid spans the fixed 5 s horizon');
ok = chk(ok, max(abs(R.X(:,1) - [0;0;0;0])) < 1e-8, 'starts at rest, pendulum down');
ok = chk(ok, max(abs(R.X(:,end) - [0;pi;0;0])) < 1e-6, ...
         sprintf('ends at rest, pendulum up (miss %.1e)', max(abs(R.X(:,end) - [0;pi;0;0]))));
ok = chk(ok, max(abs(R.U)) < 0.9*R.uMax, ...
         sprintf('the control stayed off the relaxed bound: max |u| = %.1f N of %.0f N', ...
                 max(abs(R.U)), R.uMax));

% its own trapezoidal defects, recomputed here from the stored nodes
h = diff(R.tN);
F = zeros(4, 201);
for k = 1:201
    F(:,k) = ref_field(R.X(:,k), R.U(k), R.p);
end
d = R.X(:,2:end) - R.X(:,1:end-1) - (h/2).*(F(:,2:end) + F(:,1:end-1));
ok = chk(ok, max(abs(d(:))) < 1e-6, sprintf('trapezoidal defects small: %.1e', max(abs(d(:)))));
ok = chk(ok, R.J > 0 && isfinite(R.J), sprintf('a finite positive cost: J = %.6f', R.J));

if ok, fprintf('TEST_DIRECT_REF: ALL PASS\n');
else,  fprintf('TEST_DIRECT_REF: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function dx = ref_field(x, u, p)
%% Purpose:
%
%   The example's own dynamics, written out here so the fixture is checked
%   against the problem statement rather than against the code under test.
%
q2 = x(2);  q1d = x(3);  q2d = x(4);
s = sin(q2);  c = cos(q2);
D1 = p.m1 + p.m2*(1 - c^2);
D2 = p.L*(p.m1 + p.m2)*(1 - (p.m2/(p.m1 + p.m2))*c^2);
dx = [q1d; q2d;
      (p.L*p.m2*s*q2d^2 + u + p.m2*p.g*c*s)/D1;
      (p.L*p.m2*c*s*q2d^2 + u*c + (p.m1 + p.m2)*p.g*s)/D2];
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

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = test_direct_ref()"
```

Expected: FAIL — `data/cartpole_direct_ref.mat` does not exist (`load` errors).

- [ ] **Step 3: Write the generator**

`gen_direct_ref.m`: N = 200 trapezoid segments, decision vector `[q1; q2; q1dot; q2dot; u]` stacked by variable, `fmincon` SQP with the relaxed bound `uMax = 2000`, objective the trapezoidal quadrature of `u^2`, constraints the boundary conditions plus the defects. Harvest `lambda.eqnonlin` for the defect rows.

```matlab
function gen_direct_ref()
%% Purpose:
%
%   Produce the committed fixture the indirect demo seeds from and compares
%   against: one direct trapezoidal solution of the cart-pole swing-up at a
%   RELAXED force bound, so it approximates the unconstrained optimum the
%   PMP-BVP solves. Run once; the .mat is committed.
%
%   The formulation is the existing example's
%   (ex2_cart_pole_swing_up/try2), re-expressed as a function so the
%   multipliers can be harvested.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none (writes data/cartpole_direct_ref.mat)
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
p  = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
tf = 5;  N = 200;  nN = N + 1;  uMax = 2000;
tN = linspace(0, tf, nN);
h  = tN(2) - tN(1);

% seed: the example's "good initial guess"
q1  = 0.3*sin(2*pi*tN/tf);
q2  = pi*tN/tf;
q1d = (pi/tf)*cos(pi*tN/tf);
q2d = (pi/tf)*ones(1, nN);
u0  = 40*(sin(2*pi*tN/tf) + 0.5*sin(4*pi*tN/tf));
Z0  = [q1, q2, q1d, q2d, u0].';

obj = @(Z) trap_cost(Z, nN, h);
nlc = @(Z) deal([], defects(Z, nN, h, p));
lb  = -inf(5*nN, 1);  ub = inf(5*nN, 1);
lb(4*nN+1:end) = -uMax;  ub(4*nN+1:end) = uMax;

opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
    'MaxIterations', 5000, 'MaxFunctionEvaluations', 1e7, ...
    'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-8, ...
    'SpecifyObjectiveGradient', false);
[Z, J, exitflag, ~, lambda] = fmincon(obj, Z0, [], [], [], [], lb, ub, nlc, opts);
assert(exitflag > 0, 'gen_direct_ref:fmincon', 'fmincon exit flag %d', exitflag);

X = [Z(1:nN).'; Z(nN+1:2*nN).'; Z(2*nN+1:3*nN).'; Z(3*nN+1:4*nN).'];
U = Z(4*nN+1:5*nN).';
% the equality multipliers: 4 boundary rows at each end, then 4 per interval
mu = lambda.eqnonlin(:);
muDefect = reshape(mu(9:end), 4, N);

if ~isfolder(fullfile(here, 'data')), mkdir(fullfile(here, 'data')); end
save(fullfile(here, 'data', 'cartpole_direct_ref.mat'), ...
     'tN', 'X', 'U', 'muDefect', 'J', 'p', 'tf', 'uMax', 'N');
fprintf('gen_direct_ref: J = %.6f, max|u| = %.2f N, terminal miss %.2e\n', ...
        J, max(abs(U)), max(abs(X(:,end) - [0;pi;0;0])));
end

% ------------------------------------------------------------------------
function J = trap_cost(Z, nN, h)
%% Purpose:
%
%   Trapezoidal quadrature of the running cost u^2.
%
u = Z(4*nN+1:5*nN);
J = h*(sum(u.^2) - 0.5*(u(1)^2 + u(end)^2));
end

function ceq = defects(Z, nN, h, p)
%% Purpose:
%
%   Boundary conditions (8 rows) then the trapezoidal defects (4 per
%   interval), in that order -- the order gen_direct_ref harvests by.
%
X = [Z(1:nN).'; Z(nN+1:2*nN).'; Z(2*nN+1:3*nN).'; Z(3*nN+1:4*nN).'];
u = Z(4*nN+1:5*nN).';
F = zeros(4, nN);
for k = 1:nN
    F(:,k) = ref_field(X(:,k), u(k), p);
end
d = X(:,2:end) - X(:,1:end-1) - (h/2)*(F(:,2:end) + F(:,1:end-1));
ceq = [X(:,1) - [0;0;0;0]; X(:,end) - [0;pi;0;0]; d(:)];
end

function dx = ref_field(x, u, p)
%% Purpose:
%
%   The example's dynamics: cart and pendulum accelerations with the force.
%
q2 = x(2);  q1d = x(3);  q2d = x(4);
s = sin(q2);  c = cos(q2);
D1 = p.m1 + p.m2*(1 - c^2);
D2 = p.L*(p.m1 + p.m2)*(1 - (p.m2/(p.m1 + p.m2))*c^2);
dx = [q1d; q2d;
      (p.L*p.m2*s*q2d^2 + u + p.m2*p.g*c*s)/D1;
      (p.L*p.m2*c*s*q2d^2 + u*c + (p.m1 + p.m2)*p.g*s)/D2];
end
```

- [ ] **Step 4: Generate the fixture and run the test**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp
/Applications/MATLAB_R2026a.app/bin/matlab -batch "gen_direct_ref"
cd tests && /Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = test_direct_ref(); exit(~ok)"
```

Expected: PASS on all seven checks. **If `max|u|` is at the bound**, the relaxed solve is not the unconstrained optimum — raise `uMax` and regenerate; do not relax the test.

- [ ] **Step 5: Commit**

```bash
git add optimal_control_examples/ex3_cart_pole_pmp/gen_direct_ref.m \
        optimal_control_examples/ex3_cart_pole_pmp/data/cartpole_direct_ref.mat \
        optimal_control_examples/ex3_cart_pole_pmp/tests/test_direct_ref.m
git commit -m "cart-pole PMP: the direct-solution fixture and its validation"
```

---

### Task 2: The dynamics, split as drift plus control column

**Files:**
- Create: `optimal_control_examples/ex3_cart_pole_pmp/cartpole_field.m`
- Test: `optimal_control_examples/ex3_cart_pole_pmp/tests/test_cartpole_field.m`

**Interfaces:**
- Consumes: nothing.
- Produces: `[F, G] = cartpole_field(x, p)` with `x` `[4 x 1]`, `p` a struct with `m1 m2 L g`, `F` `[4 x 1]`, `G` `[4 x 1]`. `F + G*u` is the state derivative. Used by Tasks 3, 4, 6.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_cartpole_field()
%% Purpose:
%
%   The indirect demo must solve the SAME problem as the direct example, so
%   this pins cartpole_field against that example's own helpers:
%     1. F + G*u reproduces cart_accel / pendulum_accel at a scatter of
%        states and controls;
%     2. the split is affine: G is the exact d(xdot)/du;
%     3. the first two rows are the kinematics;
%     4. it is complex-step safe (an imaginary perturbation propagates).
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
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
addpath(fullfile(fileparts(here), 'ex2_cart_pole_swing_up', 'try2'));
p = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);

rngWas = rng(11);  restore = onCleanup(@() rng(rngWas));
worst = 0;  worstAffine = 0;
for k = 1:25
    x = [2*randn; pi*randn; randn; 2*randn];
    u = 50*randn;
    [F, G] = cartpole_field(x, p);
    dx = F + G*u;
    a1 = cart_accel(x(2), x(4), u, p.L, p.m1, p.m2, p.g);
    a2 = pendulum_accel(x(2), x(4), u, p.L, p.m1, p.m2, p.g);
    worst = max(worst, max(abs(dx(3:4) - [a1; a2])));
    % affine in u: the difference quotient in u is exactly G
    [F2, ~] = cartpole_field(x, p);
    worstAffine = max(worstAffine, max(abs((F2 + G*(u+1)) - (dx + G))));
end
ok = chk(ok, worst < 1e-12, sprintf('F + G*u equals the example helpers (worst %.1e)', worst));
ok = chk(ok, worstAffine < 1e-12, sprintf('the split is exactly affine in u (worst %.1e)', worstAffine));

x = [0.3; 1.1; -0.7; 0.4];
[F, G] = cartpole_field(x, p);
ok = chk(ok, isequal(F(1:2), x(3:4)) && isequal(G(1:2), [0;0]), ...
         'rows 1-2 are the kinematics, with no control');

h = 1e-20;
[Fc, Gc] = cartpole_field(x + [0; 1i*h; 0; 0], p);
ok = chk(ok, any(imag(Fc) ~= 0) && any(imag(Gc) ~= 0) && all(isfinite([Fc; Gc])), ...
         'complex-step safe: an imaginary perturbation propagates through both outputs');

if ok, fprintf('TEST_CARTPOLE_FIELD: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_FIELD: FAILURE (see lines above)\n');
end
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

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = test_cartpole_field()"
```

Expected: FAIL — `Unrecognized function or variable 'cartpole_field'`.

- [ ] **Step 3: Write the implementation**

```matlab
function [F, G] = cartpole_field(x, p)
%% Purpose:
%
%   The cart-pole dynamics, split so the control enters affinely:
%   xdot = F(x) + G(x) u. The split is what makes the Pontryagin control
%   explicit -- with a running cost u^2, stationarity gives u* = -lam'G/2 --
%   so this file, not the PMP field, is the one home for the physics.
%
%   Same equations and constants as the direct example
%   (ex2_cart_pole_swing_up/try2/cart_accel, pendulum_accel), which
%   tests/test_cartpole_field pins.
%
%  ASSUMPTIONS / NOTES:
%
% • COMPLEX-STEP SAFE, and must stay so: no abs, no norm, no max/min, no
%   branch on a state value. The STM's Jacobian is taken by complex step
%   through this function.
%
%% Inputs:
%
%  x                        [4 x 1]                 [q1; q2; q1dot; q2dot]:
%                                                   cart position (m),
%                                                   pendulum angle (rad, 0 =
%                                                   hanging), and their rates
%
%  p                        struct                  .m1 cart mass (kg), .m2
%                                                   bob mass (kg), .L
%                                                   pendulum length (m), .g
%                                                   gravity (m/s^2)
%
%% Outputs:
%
%  F                        [4 x 1]                 Drift: xdot at u = 0
%
%  G                        [4 x 1]                 Control column: d(xdot)/du
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the drift and control column at the hanging equilibrium:
     pd = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
     [F, G] = cartpole_field([0; 0; 0; 0], pd);
     fprintf('at rest, hanging:  F = [%s],  G = [%s]\n', ...
             sprintf('%.4f ', F), sprintf('%.4f ', G));
     return
end

q1d = x(3);  q2d = x(4);
s = sin(x(2));  c = cos(x(2));
D1 = p.m1 + p.m2*(1 - c^2);
D2 = p.L*(p.m1 + p.m2)*(1 - (p.m2/(p.m1 + p.m2))*c^2);

F = [q1d;
     q2d;
     (p.L*p.m2*s*q2d^2 + p.m2*p.g*c*s)/D1;
     (p.L*p.m2*c*s*q2d^2 + (p.m1 + p.m2)*p.g*s)/D2];
G = [0; 0; 1/D1; c/D2];
end
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = test_cartpole_field(); exit(~ok)"
```

Expected: PASS on all four checks.

- [ ] **Step 5: Commit**

```bash
git add optimal_control_examples/ex3_cart_pole_pmp/cartpole_field.m \
        optimal_control_examples/ex3_cart_pole_pmp/tests/test_cartpole_field.m
git commit -m "cart-pole PMP: the dynamics split as drift plus control column"
```

---

### Task 3: The state Jacobian, generated symbolically

**Files:**
- Create: `optimal_control_examples/ex3_cart_pole_pmp/gen_state_jac.m`
- Create: `optimal_control_examples/ex3_cart_pole_pmp/cartpole_state_jac.m` (generated by the above, then given a house header by hand)
- Test: `optimal_control_examples/ex3_cart_pole_pmp/tests/test_cartpole_state_jac.m`

**Interfaces:**
- Consumes: `cartpole_field` (Task 2) for the test's cross-check only.
- Produces: `A = cartpole_state_jac(x, u, p)`, `A` `[4 x 4]`, the exact `d(F + G u)/dx` at fixed `u`. Used by Task 4.

**Why generated, not complex-stepped:** the costate equation needs this Jacobian, so it sits INSIDE the PMP field. Task 5 takes a complex step through that field to build the STM's 8x8 Jacobian. A complex step inside a complex step corrupts the inner derivative silently — both use the same imaginary channel. Generating this one symbolically breaks the nesting.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_cartpole_state_jac()
%% Purpose:
%
%   The generated state Jacobian must be the derivative it claims to be, and
%   must be usable INSIDE a complex step (that is the whole reason it is
%   generated rather than complex-stepped):
%     1. it matches a complex-step derivative of cartpole_field to 1e-12;
%     2. it matches central finite differences to 1e-7 (an independent
%        method, in case the complex step and the symbolic derivation share
%        a mistake in the field);
%     3. it is real and finite for real inputs, and complex-step clean: a
%        complex-perturbed x gives a finite A with no NaN.
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
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
p = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);

rngWas = rng(3);  restore = onCleanup(@() rng(rngWas));
worstCS = 0;  worstFD = 0;
for k = 1:15
    x = [randn; 2*randn; randn; randn];
    u = 30*randn;
    A = cartpole_state_jac(x, u, p);

    Acs = zeros(4);
    hcs = 1e-20;
    for col = 1:4
        xp = x;  xp(col) = xp(col) + 1i*hcs;
        [Fp, Gp] = cartpole_field(xp, p);
        Acs(:,col) = imag(Fp + Gp*u)/hcs;
    end
    worstCS = max(worstCS, max(abs(A(:) - Acs(:))));

    Afd = zeros(4);
    hfd = 1e-6;
    for col = 1:4
        xp = x;  xp(col) = xp(col) + hfd;
        xm = x;  xm(col) = xm(col) - hfd;
        [Fp, Gp] = cartpole_field(xp, p);
        [Fm, Gm] = cartpole_field(xm, p);
        Afd(:,col) = ((Fp + Gp*u) - (Fm + Gm*u))/(2*hfd);
    end
    worstFD = max(worstFD, max(abs(A(:) - Afd(:))));
end
ok = chk(ok, worstCS < 1e-12, sprintf('matches a complex step of the field (worst %.1e)', worstCS));
ok = chk(ok, worstFD < 1e-7,  sprintf('matches central differences (worst %.1e)', worstFD));

x = [0.2; 1.0; -0.3; 0.6];
A = cartpole_state_jac(x, 12, p);
ok = chk(ok, isreal(A) && all(isfinite(A(:))), 'real and finite for real inputs');
Ac = cartpole_state_jac(x + [0; 1i*1e-20; 0; 0], 12, p);
ok = chk(ok, all(isfinite(Ac(:))), 'survives a complex-perturbed state (no nested-step damage)');

if ok, fprintf('TEST_CARTPOLE_STATE_JAC: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_STATE_JAC: FAILURE (see lines above)\n');
end
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

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = test_cartpole_state_jac()"
```

Expected: FAIL — `Unrecognized function or variable 'cartpole_state_jac'`.

- [ ] **Step 3: Write the generator and run it**

```matlab
function gen_state_jac()
%% Purpose:
%
%   Generate cartpole_state_jac.m: the exact 4x4 d(F + G u)/dx, derived with
%   the Symbolic Math Toolbox and written out with matlabFunction. Run once;
%   the generated file is committed.
%
%   It is generated rather than taken by complex step because the costate
%   equation that uses it sits INSIDE the PMP field, and the STM's Jacobian
%   complex-steps through that field. A complex step within a complex step
%   shares the imaginary channel and silently corrupts the inner derivative.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none (writes cartpole_state_jac.m beside this file)
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
syms q1 q2 q1d q2d u m1 m2 L g real
s = sin(q2);  c = cos(q2);
D1 = m1 + m2*(1 - c^2);
D2 = L*(m1 + m2)*(1 - (m2/(m1 + m2))*c^2);
f = [q1d;
     q2d;
     (L*m2*s*q2d^2 + u + m2*g*c*s)/D1;
     (L*m2*c*s*q2d^2 + u*c + (m1 + m2)*g*s)/D2];
A = simplify(jacobian(f, [q1; q2; q1d; q2d]));

out = fullfile(here, 'cartpole_state_jac_raw.m');
matlabFunction(A, 'File', out, 'Vars', {[q1; q2; q1d; q2d], u, [m1; m2; L; g]}, ...
               'Outputs', {'A'}, 'Optimize', true);
fprintf('gen_state_jac: wrote %s -- wrap it as cartpole_state_jac(x, u, p)\n', out);
end
```

Then hand-write `cartpole_state_jac.m` as the house-styled wrapper around the generated body: the header quartet, the `nargin == 0` demo, `p` unpacked into the `[m1; m2; L; g]` vector the generated code expects, and the generated expression pasted in as the body. Delete `cartpole_state_jac_raw.m` once its body has been moved in; the generator stays so the derivation is reproducible.

- [ ] **Step 4: Run the test to verify it passes**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp
/Applications/MATLAB_R2026a.app/bin/matlab -batch "gen_state_jac"
cd tests && /Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = test_cartpole_state_jac(); exit(~ok)"
```

Expected: PASS on all four checks.

- [ ] **Step 5: Commit**

```bash
git add optimal_control_examples/ex3_cart_pole_pmp/gen_state_jac.m \
        optimal_control_examples/ex3_cart_pole_pmp/cartpole_state_jac.m \
        optimal_control_examples/ex3_cart_pole_pmp/tests/test_cartpole_state_jac.m
git commit -m "cart-pole PMP: exact state Jacobian, generated symbolically"
```

---

### Task 4: The PMP field

**Files:**
- Create: `optimal_control_examples/ex3_cart_pole_pmp/cartpole_pmp_rhs.m`
- Test: `optimal_control_examples/ex3_cart_pole_pmp/tests/test_cartpole_pmp_rhs.m`

**Interfaces:**
- Consumes: `cartpole_field` (Task 2), `cartpole_state_jac` (Task 3).
- Produces: `[dy, u] = cartpole_pmp_rhs(y, p)` with `y = [x; lam]` `[8 x 1]`, `dy` `[8 x 1]`, `u` the scalar PMP control. Used by Tasks 5 and 6.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_cartpole_pmp_rhs()
%% Purpose:
%
%   The PMP field must BE the Pontryagin conditions, not merely resemble
%   them:
%     1. the control is the stationary point: 2u + lam'G = 0 exactly;
%     2. it is the MINIMISER, not just a stationary point: H at u* is below
%        H at u* +/- d for a scatter of d (Legendre, d2H/du2 = 2 > 0);
%     3. the state rows are F + G u*;
%     4. the costate rows are -dH/dx, checked against a complex step of H
%        taken with u HELD FIXED (which is legitimate: dH/du = 0 at u*);
%     5. zero costates give zero costate rates and a drifting state (the
%        uncontrolled flow), so lam = 0 is the uncontrolled solution.
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
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
p = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);

rngWas = rng(5);  restore = onCleanup(@() rng(rngWas));
worstStat = 0;  worstState = 0;  worstCostate = 0;  minGap = inf;
for k = 1:15
    y = [randn; 2*randn; randn; randn; 3*randn; 3*randn; 3*randn; 3*randn];
    x = y(1:4);  lam = y(5:8);
    [dy, u] = cartpole_pmp_rhs(y, p);
    [F, G] = cartpole_field(x, p);

    worstStat  = max(worstStat, abs(2*u + lam.'*G));
    worstState = max(worstState, max(abs(dy(1:4) - (F + G*u))));

    Hstar = u^2 + lam.'*(F + G*u);
    for d = [-1 -0.1 0.1 1]
        Hd = (u+d)^2 + lam.'*(F + G*(u+d));
        minGap = min(minGap, Hd - Hstar);
    end

    % -dH/dx by complex step, u held fixed
    dHdx = zeros(4,1);
    hcs = 1e-20;
    for col = 1:4
        xp = x;  xp(col) = xp(col) + 1i*hcs;
        [Fp, Gp] = cartpole_field(xp, p);
        Hp = u^2 + lam.'*(Fp + Gp*u);
        dHdx(col) = imag(Hp)/hcs;
    end
    worstCostate = max(worstCostate, max(abs(dy(5:8) + dHdx)));
end
ok = chk(ok, worstStat < 1e-12, sprintf('u is stationary: max |2u + lam''G| = %.1e', worstStat));
ok = chk(ok, minGap > 0, sprintf('and a MINIMISER: worst H(u*+d) - H(u*) = %.3e > 0', minGap));
ok = chk(ok, worstState < 1e-12, sprintf('state rows are F + G u* (worst %.1e)', worstState));
ok = chk(ok, worstCostate < 1e-9, sprintf('costate rows are -dH/dx (worst %.1e)', worstCostate));

y0 = [0.1; 0.2; 0; 0; 0; 0; 0; 0];
[dy0, u0] = cartpole_pmp_rhs(y0, p);
ok = chk(ok, u0 == 0 && all(dy0(5:8) == 0) && any(dy0(1:4) ~= 0), ...
         'lam = 0 gives u = 0, still costates, and the uncontrolled drift');

if ok, fprintf('TEST_CARTPOLE_PMP_RHS: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_PMP_RHS: FAILURE (see lines above)\n');
end
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

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = test_cartpole_pmp_rhs()"
```

Expected: FAIL — `Unrecognized function or variable 'cartpole_pmp_rhs'`.

- [ ] **Step 3: Write the implementation**

```matlab
function [dy, u] = cartpole_pmp_rhs(y, p)
%% Purpose:
%
%   The Pontryagin field of the cart-pole minimum-effort swing-up at fixed
%   final time: the 8-state flow y = [x; lam] whose boundary-value problem
%   IS the optimality condition.
%
%   With J = integral u^2 dt and xdot = F(x) + G(x) u, the Hamiltonian is
%       H = u^2 + lam'(F + G u),
%   so dH/du = 2u + lam'G = 0 gives the control in closed form,
%       u* = -lam'G / 2,
%   and d2H/du2 = 2 > 0 makes it the minimiser everywhere -- there is no
%   switching structure, which is exactly why this problem is the clean
%   demonstration that the shooting engine is generic. The costate obeys
%   lamdot = -dH/dx at u*; the running cost has no state dependence, so that
%   is -A(x,u*)' lam with A the state Jacobian.
%
%  ASSUMPTIONS / NOTES:
%
% • Complex-step safe: A comes from the GENERATED cartpole_state_jac, never
%   from a complex step, so a caller may complex-step through this function
%   (that is how the STM's Jacobian is built).
% • u* is substituted before the costate equation is evaluated. Doing so is
%   legitimate precisely because dH/du = 0 there: the total and partial
%   x-derivatives of H agree at the stationary control.
%
%% Inputs:
%
%  y                        [8 x 1]                 [x; lam]: the four
%                                                   states then the four
%                                                   costates
%
%  p                        struct                  .m1 .m2 .L .g
%
%% Outputs:
%
%  dy                       [8 x 1]                 [xdot; lamdot]
%
%  u                        double                  The PMP control at y
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the field and its control at a representative point:
     pd = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
     [dy, u] = cartpole_pmp_rhs([0; 0.5; 0; 0; 1; -2; 0.5; 0.3], pd);
     fprintf('u* = %.6f,  dy = [%s]\n', u, sprintf('%.4f ', dy));
     return
end

x = y(1:4);  lam = y(5:8);
[F, G] = cartpole_field(x, p);
u = -(lam.'*G)/2;
A = cartpole_state_jac(x, u, p);
dy = [F + G*u; -A.'*lam];
end
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = test_cartpole_pmp_rhs(); exit(~ok)"
```

Expected: PASS on all five checks.

- [ ] **Step 5: Commit**

```bash
git add optimal_control_examples/ex3_cart_pole_pmp/cartpole_pmp_rhs.m \
        optimal_control_examples/ex3_cart_pole_pmp/tests/test_cartpole_pmp_rhs.m
git commit -m "cart-pole PMP: the 8-state Pontryagin field"
```

---

### Task 5: The propagator with its state-transition matrix

**Files:**
- Create: `optimal_control_examples/ex3_cart_pole_pmp/cartpole_pmp_prop.m`
- Test: `optimal_control_examples/ex3_cart_pole_pmp/tests/test_cartpole_pmp_prop.m`

**Interfaces:**
- Consumes: `cartpole_pmp_rhs` (Task 4).
- Produces: `[yEnd, PHI] = cartpole_pmp_prop(dt, y0, needSTM, p)` — `yEnd` `[8 x 1]`, `PHI` `[8 x 8]` when `needSTM` else `[]`. This is the closure `ms_bvp` calls as `prob.prop`; Task 6 binds `p` into it.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_cartpole_pmp_prop()
%% Purpose:
%
%   The propagator and the STM it hands the shooting engine:
%     1. the STM matches central finite differences of the propagated map,
%        column by column (the STM is what the Newton step is built from, so
%        a wrong one shows up as slow convergence, not as a wrong answer);
%     2. PHI(0) = I;
%     3. the group property: propagating dt in one call equals two calls of
%        dt/2, and the STMs multiply;
%     4. it THROWS rather than returning junk when the flow blows up (the
%        engine's contract);
%     5. with no STM requested the state still matches the STM run.
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
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
p = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
y0 = [0; 0; 0; 0; 0.8; -1.5; 0.4; 0.2];
dt = 0.6;

[yE, PHI] = cartpole_pmp_prop(dt, y0, true, p);
ok = chk(ok, isequal(size(PHI), [8 8]) && numel(yE) == 8, 'shapes: yEnd [8 x 1], PHI [8 x 8]');

PHIfd = zeros(8);
h = 1e-6;
for col = 1:8
    yp = y0;  yp(col) = yp(col) + h;
    ym = y0;  ym(col) = ym(col) - h;
    PHIfd(:,col) = (cartpole_pmp_prop(dt, yp, false, p) - cartpole_pmp_prop(dt, ym, false, p))/(2*h);
end
rel = max(abs(PHI(:) - PHIfd(:)))/max(abs(PHIfd(:)));
ok = chk(ok, rel < 1e-6, sprintf('STM matches finite differences: relative %.1e', rel));

[~, PHI0] = cartpole_pmp_prop(0, y0, true, p);
ok = chk(ok, max(abs(PHI0(:) - reshape(eye(8), 64, 1))) < 1e-12, 'PHI(0) = I');

[yH, PHIa] = cartpole_pmp_prop(dt/2, y0, true, p);
[yF, PHIb] = cartpole_pmp_prop(dt/2, yH, true, p);
ok = chk(ok, max(abs(yF - yE)) < 1e-9, sprintf('one step equals two half steps (%.1e)', max(abs(yF - yE))));
ok = chk(ok, max(abs(reshape(PHIb*PHIa - PHI, 64, 1))) < 1e-7, 'and their STMs multiply');

yNo = cartpole_pmp_prop(dt, y0, false, p);
ok = chk(ok, max(abs(yNo - yE)) < 1e-9, 'the state is the same with and without the STM');

threw = false;
try
    cartpole_pmp_prop(50, [0; 0; 0; 0; 1e6; 1e6; 1e6; 1e6], false, p);
catch
    threw = true;
end
ok = chk(ok, threw, 'a blown-up flow THROWS rather than returning junk');

if ok, fprintf('TEST_CARTPOLE_PMP_PROP: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_PMP_PROP: FAILURE (see lines above)\n');
end
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

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = test_cartpole_pmp_prop()"
```

Expected: FAIL — `Unrecognized function or variable 'cartpole_pmp_prop'`.

- [ ] **Step 3: Write the implementation**

```matlab
function [yEnd, PHI] = cartpole_pmp_prop(dt, y0, needSTM, p)
%% Purpose:
%
%   Propagate the cart-pole PMP flow, with the 8x8 state-transition matrix
%   on request -- the propagator contract costate_common/ms_bvp expects as
%   prob.prop.
%
%   The variational equations are integrated ALONGSIDE the state (72
%   equations in all), with the Jacobian taken by COMPLEX STEP through
%   cartpole_pmp_rhs. That is safe here only because the state Jacobian
%   inside the field is generated symbolically: a complex step inside a
%   complex step shares the imaginary channel and silently corrupts the
%   inner derivative.
%
%  ASSUMPTIONS / NOTES:
%
% • Tolerances 1e-12 / 1e-14: the STM feeds a Newton step, and an integrator
%   sloppier than the step it informs turns quadratic convergence into a
%   crawl.
% • THROWS on a non-finite result, per the engine's contract -- the solver
%   converts a throw into a rejected iterate and backtracks.
%
%% Inputs:
%
%  dt                       double                  Propagation time (s);
%                                                   0 returns y0 and I
%
%  y0                       [8 x 1]                 [x; lam] at the start
%
%  needSTM                  logical                 Return PHI as well
%
%  p                        struct                  .m1 .m2 .L .g
%
%% Outputs:
%
%  yEnd                     [8 x 1]                 The propagated state
%
%  PHI                      [8 x 8] or []           d(yEnd)/d(y0)
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: half a second of the flow, with its STM conditioning:
     pd = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
     [yE, PH] = cartpole_pmp_prop(0.5, [0;0;0;0; 0.8; -1.5; 0.4; 0.2], true, pd);
     fprintf('y(0.5) = [%s]\n  cond(PHI) = %.3e\n', sprintf('%.4f ', yE), cond(PH));
     return
end

if dt == 0
    yEnd = y0(:);
    PHI = [];  if needSTM, PHI = eye(8); end
    return
end

oo = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
if needSTM
    z0 = [y0(:); reshape(eye(8), 64, 1)];
    [~, Z] = ode113(@(t, z) rhs_with_stm(z, p), [0 dt], z0, oo);
    zEnd = Z(end,:).';
    yEnd = zEnd(1:8);
    PHI  = reshape(zEnd(9:72), 8, 8);
else
    [~, Y] = ode113(@(t, y) cartpole_pmp_rhs(y, p), [0 dt], y0(:), oo);
    yEnd = Y(end,:).';
    PHI  = [];
end

if ~all(isfinite(yEnd)) || (needSTM && ~all(isfinite(PHI(:))))
    error('cartpole_pmp_prop:collapse', ...
          'the PMP flow left the finite range over dt = %g', dt);
end
end

% ------------------------------------------------------------------------
function dz = rhs_with_stm(z, p)
%% Purpose:
%
%   The state and its variational equations: zdot = [f(y); A(y) PHI], with
%   A taken by complex step through the PMP field.
%
y = z(1:8);
PHI = reshape(z(9:72), 8, 8);
dy = cartpole_pmp_rhs(y, p);
A = zeros(8);
h = 1e-20;
for col = 1:8
    yp = y;  yp(col) = yp(col) + 1i*h;
    A(:,col) = imag(cartpole_pmp_rhs(yp, p))/h;
end
dz = [dy; reshape(A*PHI, 64, 1)];
end
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = test_cartpole_pmp_prop(); exit(~ok)"
```

Expected: PASS on all seven checks. **If the throw check fails** because the integrator merely stalls instead of producing non-finite values, add `'MaxStep'`-independent protection by wrapping the `ode113` call in a try/catch that rethrows as `cartpole_pmp_prop:collapse`; do not delete the check.

- [ ] **Step 5: Commit**

```bash
git add optimal_control_examples/ex3_cart_pole_pmp/cartpole_pmp_prop.m \
        optimal_control_examples/ex3_cart_pole_pmp/tests/test_cartpole_pmp_prop.m
git commit -m "cart-pole PMP: propagator with complex-step variational STM"
```

---

### Task 6: The solve

**Files:**
- Create: `optimal_control_examples/ex3_cart_pole_pmp/run_cartpole_pmp.m`
- Create: `optimal_control_examples/ex3_cart_pole_pmp/data/cartpole_pmp_ref.mat` (written once, in Step 4)
- Test: `optimal_control_examples/ex3_cart_pole_pmp/tests/test_cartpole_pmp.m`
- Read: `orbit_transfer/costate_common/ms_bvp.m` (the engine; interface in its header)

**Interfaces:**
- Consumes: Tasks 1-5, plus `ms_bvp` and `oc.duals_to_costates`.
- Produces: `out = run_cartpole_pmp(opts)` with fields `.lam0` `[4 x 1]`, `.J` (scalar), `.missTerminal` (scalar), `.missFlown` (scalar), `.statMax` (scalar), `.t` `[1 x M]`, `.X` `[4 x M]`, `.U` `[1 x M]`, `.info` (the `ms_bvp` info struct). `opts` fields: `.K` (segments, default 8), `.plot` (default true when called with no output), `.engine` (function handle, default `@ms_bvp` — Task 7 passes `@oc.ms_bvp`).

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_cartpole_pmp()
%% Purpose:
%
%   The whole indirect solve, and the claims the spec makes for it:
%     1. it converges, and the terminal state is hit to 1e-9;
%     2. stationarity holds ALONG the trajectory, not just at the ends;
%     3. flying the recovered control through the true dynamics
%        (oc.fly_control) arrives where the solve says -- the G1b idea;
%     4. the cost agrees with the direct fixture to 1% (the gap is the
%        direct method's discretization error; the indirect solve is the
%        more accurate of the two);
%     5. the control profiles agree in shape to 2% RMS;
%     6. regression: lam0 reproduces the stored reference to 1e-8;
%     7. the same answer comes back with a different segment count -- the
%        solution is the problem's, not the discretization's.
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
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
root = fileparts(fileparts(here));                 % optimal_control
addpath(here, fullfile(root, 'oclib'), fullfile(root, 'orbit_transfer', 'costate_common'));

out = run_cartpole_pmp(struct('K', 8, 'plot', false));
R = load(fullfile(here, 'data', 'cartpole_direct_ref.mat'));

ok = chk(ok, out.info.converged, 'ms_bvp reports converged');
ok = chk(ok, out.missTerminal < 1e-9, sprintf('terminal miss %.2e < 1e-9', out.missTerminal));
ok = chk(ok, out.statMax < 1e-8, sprintf('stationarity along the arc: max |2u + lam''G| = %.2e', out.statMax));
ok = chk(ok, out.missFlown < 1e-6, sprintf('flown control arrives: %.2e', out.missFlown));

relJ = abs(out.J - R.J)/R.J;
ok = chk(ok, relJ < 0.01, sprintf('cost agrees with the direct solution: J = %.6f vs %.6f (%.3f%%)', ...
                                  out.J, R.J, 100*relJ));

Uref = interp1(R.tN, R.U, out.t, 'linear');
rmsU = sqrt(mean((out.U - Uref).^2))/sqrt(mean(Uref.^2));
ok = chk(ok, rmsU < 0.02, sprintf('control profiles agree: %.3f%% RMS', 100*rmsU));

ref = load(fullfile(here, 'data', 'cartpole_pmp_ref.mat'));
ok = chk(ok, max(abs(out.lam0 - ref.lam0)) < 1e-8, ...
         sprintf('regression: lam0 reproduces (%.1e)', max(abs(out.lam0 - ref.lam0))));

out16 = run_cartpole_pmp(struct('K', 16, 'plot', false));
ok = chk(ok, max(abs(out16.lam0 - out.lam0)) < 1e-7, ...
         sprintf('K = 16 finds the same extremal (%.1e)', max(abs(out16.lam0 - out.lam0))));

if ok, fprintf('TEST_CARTPOLE_PMP: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_PMP: FAILURE (see lines above)\n');
end
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

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = test_cartpole_pmp()"
```

Expected: FAIL — `Unrecognized function or variable 'run_cartpole_pmp'`.

- [ ] **Step 3: Write the implementation**

```matlab
function out = run_cartpole_pmp(opts)
%% Purpose:
%
%   Solve the cart-pole minimum-effort swing-up INDIRECTLY: build the
%   Pontryagin boundary-value problem and shoot it with the shared multiple-
%   shooting engine (costate_common/ms_bvp, or oclib's oc.ms_bvp).
%
%   Four unknowns -- lam(0) -- against four terminal conditions, at fixed
%   final time. The seed comes from the committed direct solution's defect
%   multipliers, mapped to costates by oc.duals_to_costates with the
%   TRAPEZOID station rule (the CR3BP campaigns use Hermite-Simpson, so this
%   demo exercises a rule the orbit work never does). Its primer sign vote
%   does not apply to a scalar force, so the sign is resolved here against
%   the direct solve's own control.
%
%   This file is also the second TOP-LEVEL consumer of the oclib package,
%   which is what admits ms_bvp to it: nothing here is an orbit, a CR3BP
%   quantity, or a pumpkyn call.
%
%% Inputs:
%
%  opts                     struct (optional)       .K segments [8], .plot
%                                                   [true when nargout = 0],
%                                                   .engine solver handle
%                                                   [@ms_bvp]
%
%% Outputs:
%
%  out                      struct                  .lam0, .J, .missTerminal,
%                                                   .missFlown, .statMax,
%                                                   .t, .X, .U, .info
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(here, fullfile(root, 'oclib'), fullfile(root, 'orbit_transfer', 'costate_common'));
d = @(f, v) fieldd(opts, f, v);
K      = d('K', 8);
engine = d('engine', @ms_bvp);
doPlot = d('plot', nargout == 0);

R = load(fullfile(here, 'data', 'cartpole_direct_ref.mat'));
p = R.p;  tf = R.tf;
xf = [0; pi; 0; 0];

%% Seed: the direct solve's multipliers become costates
%  duals_to_costates resolves the GLOBAL SIGN by a primer vote over 3-vector
%  thrust directions. This problem's control is a scalar force, so there is
%  no such vote: the mapping is asked for the costates unsigned (no uDir),
%  and the sign is fixed here by the same principle in scalar form -- the
%  implied control u = -lam'G/2 must agree in sign with the control the
%  direct solve actually used.
[lamS, tS, dg] = oc.duals_to_costates(struct('scheme', 'trapezoid', 'mu', R.muDefect, ...
                                             'tNodes', R.tN, 'velRows', 3:4));
uImplied = zeros(1, numel(tS));
Xs = interp1(R.tN, R.X.', tS, 'pchip').';
for k = 1:numel(tS)
    [~, Gk] = cartpole_field(Xs(:,k), p);
    uImplied(k) = -(lamS(:,k).'*Gk)/2;
end
uDirect = interp1(R.tN, R.U, tS, 'pchip');
if sum(uImplied.*uDirect) < 0, lamS = -lamS; end
assert(sum(uImplied.*uDirect) ~= 0, 'run_cartpole_pmp:sign', ...
       'the implied and solved controls are orthogonal: the seed carries no sign information');
tGrid = linspace(0, tf, K+1);
Xg   = interp1(R.tN, R.X.', tGrid, 'pchip').';
Lg   = interp1(tS, lamS.', tGrid, 'pchip', 'extrap').';
seed = struct('tf', tf, 'tGrid', tGrid, 'Y', [Xg; Lg]);
seed.Y(1:4,1) = [0; 0; 0; 0];                  % the fixed departure state

%% The problem, as three closures
prob = struct('ny', 8, 'freeIdx0', 5:8, ...
    'prop', @(dt, y0, needSTM) cartpole_pmp_prop(dt, y0, needSTM, p), ...
    'rhs',  @(y) cartpole_pmp_rhs(y, p), ...
    'terminal', @(y, needJ) terminalFcn(y, xf, needJ));

[~, info] = engine(prob, seed, struct('fixedTf', true, 'tolR', 1e-12, 'maxIter', 60));

%% Report: fly the answer, measure what the gates measure
lam0 = info.Y(5:8, 1);
[t, Y] = ode113(@(tt, y) cartpole_pmp_rhs(y, p), linspace(0, tf, 501), ...
                [0; 0; 0; 0; lam0], odeset('RelTol', 1e-12, 'AbsTol', 1e-14));
X = Y(:,1:4).';  Lam = Y(:,5:8).';
U = zeros(1, numel(t));  stat = zeros(1, numel(t));
for k = 1:numel(t)
    [~, G] = cartpole_field(X(:,k), p);
    U(k)    = -(Lam(:,k).'*G)/2;
    stat(k) = abs(2*U(k) + Lam(:,k).'*G);
end
J = trapz(t, U.^2);

uOf = @(tt) interp1(t, U, min(max(tt, t(1)), t(end)), 'pchip');
zEnd = oc.fly_control([0; 0; 0; 0], [0 tf], ...
    @(tt, x) flyRhs(x, uOf(tt), p), struct('mode', 'span', 'solver', @ode113));

out = struct('lam0', lam0, 'J', J, ...
    'missTerminal', max(abs(X(:,end) - xf)), ...
    'missFlown', max(abs(zEnd - xf)), 'statMax', max(stat), ...
    't', t.', 'X', X, 'U', U, 'info', info);

if doPlot, plotAgainstDirect(out, R); end
if nargout == 0
    fprintf(['cart-pole PMP-BVP: J = %.6f (direct %.6f), terminal miss %.2e, ', ...
             'flown miss %.2e, stationarity %.2e\n'], ...
            out.J, R.J, out.missTerminal, out.missFlown, out.statMax);
    clear out
end
end

% ------------------------------------------------------------------------
function [g, dgdy] = terminalFcn(y, xf, needJ)
%% Purpose:
%
%   The four terminal conditions x(tf) = xf, and their Jacobian.
%
g = y(1:4) - xf;
dgdy = [];
if needJ, dgdy = [eye(4), zeros(4)]; end
end

function dx = flyRhs(x, u, p)
%% Purpose:
%
%   The TRUE dynamics with a supplied control -- what the flown check flies.
%
[F, G] = cartpole_field(x, p);
dx = F + G*u;
end

function plotAgainstDirect(out, R)
%% Purpose:
%
%   The indirect solution against the direct one: states and control.
%
figure('color', [1 1 1]);
subplot(2,1,1); plot(out.t, out.X(2,:), 'k', R.tN, R.X(2,:), 'r--'); grid on
ylabel('q_2 (rad)'); legend({'indirect (PMP-BVP)', 'direct (collocation)'}, 'Location', 'best');
title('cart-pole swing-up: indirect vs direct');
subplot(2,1,2); plot(out.t, out.U, 'k', R.tN, R.U, 'r--'); grid on
xlabel('t (s)'); ylabel('u (N)');
end

function v = fieldd(s, f, v0)
%% Purpose:
%
%   s.(f) if present, else the default.
%
if isfield(s, f), v = s.(f); else, v = v0; end
end
```

- [ ] **Step 4: Solve once, write the regression reference, then run the test**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp
/Applications/MATLAB_R2026a.app/bin/matlab -batch "out = run_cartpole_pmp(struct('plot',false)); lam0 = out.lam0; J = out.J; save('data/cartpole_pmp_ref.mat','lam0','J'); fprintf('lam0 = [%s], J = %.8f\n', sprintf('%.12g ', lam0), J)"
cd tests && /Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = test_cartpole_pmp(); exit(~ok)"
```

Expected: PASS on all eight checks.

**If the seed does not converge** (the spec's named risk), switch to the fallback before touching any tolerance: solve with `xf = [0; theta; 0; 0]` for `theta = [0.5 1.0 1.5 2.0 2.5 3.0 pi]`, each solve's `lam0` seeding the next, and keep the homotopy in `run_cartpole_pmp` behind `opts.homotopy` (default true only if the direct seed fails). Record which route was used in the folder README and amend the spec's Seeding section — the spec says this out loud so the fallback is a decision, not a scramble.

- [ ] **Step 5: Commit**

```bash
git add optimal_control_examples/ex3_cart_pole_pmp/run_cartpole_pmp.m \
        optimal_control_examples/ex3_cart_pole_pmp/data/cartpole_pmp_ref.mat \
        optimal_control_examples/ex3_cart_pole_pmp/tests/test_cartpole_pmp.m
git commit -m "cart-pole PMP: the indirect solve, seeded from the direct solution's duals"
```

---

### Task 7: Promote ms_bvp into oclib

**Files:**
- Move: `orbit_transfer/costate_common/ms_bvp.m` to `oclib/+oc/ms_bvp.m`
- Create: `orbit_transfer/costate_common/ms_bvp.m` (delegate)
- Modify: `oclib/README.md` (contents table, roadmap item 3)
- Modify: `orbit_transfer/costate_common/README.md` (the `ms_bvp` row), `orbit_transfer/costate_common/TODO.md` (step 13)
- Modify: `optimal_control_examples/ex3_cart_pole_pmp/run_cartpole_pmp.m:engine default` → `@oc.ms_bvp`
- Test: existing suites, listed in the gates below

**Interfaces:**
- Consumes: everything above.
- Produces: `oc.ms_bvp(prob, seed, opts)` with the signature unchanged, and `ms_bvp` in costate_common forwarding to it.

- [ ] **Step 1: Record the baseline the move must reproduce**

```bash
cd /private/tmp/claude-501/-Users-msc-Desktop-optimal-control/*/scratchpad
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd /Users/msc/Desktop/optimal_control/orbit_transfer/costate_common; ok = golden_cells(); fprintf('GOLDEN %d\n', ok)" | tail -25
```

Expected: 20/20, and the printed residuals/iterations kept for the after-comparison.

- [ ] **Step 2: Move the file and write the delegate**

```bash
cd /Users/msc/Desktop/optimal_control
git mv orbit_transfer/costate_common/ms_bvp.m oclib/+oc/ms_bvp.m
```

Then create `orbit_transfer/costate_common/ms_bvp.m`:

```matlab
function [p, info] = ms_bvp(prob, seed, opts)
%% Purpose:
%
%   DELEGATE since 2026-09-16: the generic multiple-shooting engine was
%   promoted to the cross-folder library, oclib/+oc/ms_bvp, when the
%   cart-pole PMP-BVP demo (optimal_control_examples/ex3_cart_pole_pmp) became
%   its second TOP-LEVEL consumer. This file keeps every costate_common
%   caller working; the contract, the maths and the tests live with the
%   implementation.
%
%% Inputs:
%
%  prob, seed, opts                                 see oc.ms_bvp
%
%% Outputs:
%
%  p, info                                          see oc.ms_bvp
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if isempty(which('oc.ms_bvp'))
    addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'oclib'));
end
if nargin < 3, opts = struct(); end
[p, info] = oc.ms_bvp(prob, seed, opts);
end
```

- [ ] **Step 3: Point the demo at the package**

In `run_cartpole_pmp.m`, change the default engine to `@oc.ms_bvp` so the top-level consumer reaches the package directly rather than through costate_common:

```matlab
engine = d('engine', @oc.ms_bvp);
```

- [ ] **Step 4: Run every gate**

```bash
cd /Users/msc/Desktop/optimal_control/orbit_transfer/costate_common/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "for t = {'test_ms_bvp_extra','test_ms_bvp_fixedtf','test_ms_tfmin_hom','test_arclength_ms','test_folder_rules'}, ok = feval(t{1}); fprintf('%s %d\n', t{1}, ok); end"
cd ../ && /Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = golden_cells(); fprintf('GOLDEN %d\n', ok)" | tail -25
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "ok = test_cartpole_pmp(); exit(~ok)"
```

Then the four-campaign ladder re-solve, the same gate step 12 used:

```bash
cd /private/tmp/claude-501/-Users-msc-Desktop-optimal-control/*/scratchpad
GATE_PHASE=ms_bvp_promotion /Applications/MATLAB_R2026a.app/bin/matlab -batch "run('gate_job.m')" | grep -E "^== |max\|dTF\|"
```

Expected: every test green; `golden_cells` 20/20 with the SAME residuals and iteration counts as Step 1; the four campaigns at `max|dTF|` 4.55e-14 (HALO), 1.13e-14 (DPO), 3.48e-14 (HALO_HALO), 2.72e-14 (GTO), as in step 12.

**If any number moves**, stop and report it. A delegate that changes a result is not a delegate.

- [ ] **Step 5: Update the documentation and commit**

`oclib/README.md`: add the `oc.ms_bvp` row (consumers: `orbit_transfer` via the costate_common delegate, and `optimal_control_examples/ex3_cart_pole_pmp` directly; equivalence gate: golden_cells 20/20 plus the four-campaign re-solve plus the cart-pole test), and mark roadmap item 3's `ms_bvp` half done, leaving `ms_conjugate_test` open with its reason.

`costate_common/README.md`: change the `ms_bvp.m` row to say DELEGATE, pointing at `../../oclib/+oc/ms_bvp`, exactly as the `duals_to_costates` row reads.

`costate_common/TODO.md`: mark step 13's `ms_bvp` done with the date and gates; leave the rest of the list (`ms_conjugate_test`, `arclength_ms`, `newton_fixed_q`, `conj_resolve`, `lift_space_dim`) open, still blocked on a second top-level consumer.

```bash
cd /Users/msc/Desktop/optimal_control
git add oclib/+oc/ms_bvp.m oclib/README.md \
        orbit_transfer/costate_common/ms_bvp.m \
        orbit_transfer/costate_common/README.md \
        orbit_transfer/costate_common/TODO.md \
        optimal_control_examples/ex3_cart_pole_pmp/run_cartpole_pmp.m
git commit -m "oclib: promote ms_bvp, with the cart-pole PMP-BVP as its second consumer"
```

---

## Notes for the executor

- **Run order matters**: Tasks 1-5 are independent of the orbit library; Task 6 is the first to need `costate_common` and `oclib` on the path; Task 7 touches the library that five campaigns run, so its gates are not optional.
- **Do not widen a tolerance to make a check pass.** Every numeric bound in these tests came from the spec, and each has a named reason. If one fails, the finding is the result — report it.
- **The scalar-control sign** (Task 6): if the assert about orthogonal
  controls ever fires, the seed is uninformative and the homotopy fallback
  is the route, not a guessed sign.
- **The demo must stay campaign-free**: no file under `ex3_cart_pole_pmp/` may reference an orbit, a CR3BP quantity, pumpkyn, or any campaign folder. That is the property being demonstrated.
