# Drone Quadrotor Upgrade Implementation Plan (Stage 3)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Put the library's 6-DOF rigid-body quadrotor (and the point-mass middle rung) under the same path-planning task: a direct CasADi/IPOPT solve for M vehicles with clearance, rotor and attitude boxes, gated by the local residual and the fine-grid clearance; an indirect PMP-BVP for one vehicle seeded from the direct duals and conjugate-tested; a study script and a comparison table across the three models.

**Architecture:** `optimal_control_examples/drone/quadrotor/`: `quad_layouts` (hover-to-hover boundary states), `solve_quads_direct` (one NLP for `pointmass` or `rigid6dof`, M vehicles, the dynamics as a CasADi Function built by calling `drone_field` on symbolic inputs so there is one copy of the physics), `verify_quads_direct` (gates), `quad_pmp_rhs`/`quad_pmp_prop` (24-state PMP flow, `f* = −B'λ/2`, complex-step STM as the cart-pole does), `run_quad_pmp` (front door), `quad_minenergy_study.m`, `movie_quads`, `gen_quad_fixtures`, `quad_compare`.

**Tech Stack:** MATLAB R2026a, CasADi 3.7.0 (via `route1_direct/casadi_path`), `oclib`, Symbolic Math Toolbox (already used by Stage 0's generator).

**Spec:** `docs/superpowers/specs/2026-10-23-drone-multiagent-design.md` (Stage 3 section). Stages 0-2 are DONE before this plan starts. RE-CHECK before execution: (a) whether `drone_field` runs on CasADi `MX` inputs unchanged (Task 2 pins it numerically; if it does not, the fix goes INTO `drone_field`, keeping it complex-step safe, never a second copy); (b) the rotor and angle boxes are inactive on the M = 1 solution (Task 2 reports it; if active, the indirect comparison X2 is skipped with a note, per the spec); (c) Task 5's study script and its test are specified by reference to the Route 2 study (deliberately: their exact tolerances come out of Stage 2's measured numbers) and MUST be expanded into full code at this re-check before Task 5 is dispatched.

## Global Constraints

- House MATLAB style (header quartet, demos, no `i`/`j`, no `norm`, no `%#ok`, no emojis); relative paths; no `orbit_transfer`.
- Verbatim numbers: `T = 4`, `N = 100`, `r = 0.06`, M = 1 line `(−0.4, 0, 0.5) → (0.4, 0, 0.5)`, M = 4 circle radius 0.4 at `z = 0.5`, hover boundary (`v = 0`, `Θ = 0`, `ω = 0`), rotor box `[0, fMax]`, angle box `±thetaMax` on rows 7-8, point-mass `|a| ≤ aMax` (quadratic form); IPOPT `tol 1e-10`, `acceptable_tol 1e-8`, `max_iter 3000`.
- Gates: G3 defects `< 1e-9`; G1-local `oc.local_residual` max per-interval miss `< 1e-4` (reported with the mesh gate; trapezoid at `h = 0.04 s`); G1 clearance fine-grid `≥ −1e-2·2r`; bounds margin reported; G4 mesh `< 1e-3`. Indirect: N1 `|R| < 1e-10`, N2 `< 1e-9`, N3 `< 1e-6`, N4 `< 1e-9`, N5 `< 1e-9`, X2 `|J − J_direct|/J_direct < 1e-2` only when the direct bounds are inactive, X3 `< 1e-8`, S2 not FAIL.
- Defects first in the `Opti`; multipliers by `opti.lam_g` row range; never `opti.dual`.
- Subagents never commit; commit messages end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`; do not stage the unrelated `missiles/` deletions. Fixtures force-added with generators.

## Review Focus

1. Euler-angle wrap: a direct solution that yaws by 2π is a different point in state space with the same physics; the hover boundary pins `ψ(T) = 0`, and the angle box keeps pitch/roll away from the singularity; the test asserts `|Θ| ≤ thetaMax` on the fine grid, not only at nodes. (Task 2.)
2. The CasADi dynamics must be the library's: `fdyn` evaluated at 20 random `(x, u)` must equal `drone_field` to 1e-14, for both models. (Task 2 test.)
3. The PMP sign: the covector seed's global sign is resolved by correlating the implied `f* = −B'λ/2` with the direct `U`; a correlation under 0.5 is refused by identifier. (Task 4 test with a negated seed.)
4. Inactive-bounds fairness: `run_quad_pmp` compares costs only when `out.boundsActive` of the fixture is false; otherwise `.notes` says the comparison is between different problems. (Task 4.)
5. The complex-step STM through `quad_pmp_rhs` requires that no complex step happens INSIDE it: `drone_state_jac` is generated, so this holds; the propagator test checks the STM against a complex step through the whole flow (an independent differentiation). (Task 3.)

---

### Task 1: `quad_layouts`

**Files:**
- Create: `optimal_control_examples/drone/quadrotor/quad_layouts.m`
- Test: `optimal_control_examples/drone/quadrotor/tests/test_quad_layouts.m`

**Interfaces:**
- Produces: `[x0, xT, geom] = quad_layouts(which, model)` with `which ∈ {'line','circle'}`, `model ∈ {'pointmass','rigid6dof'}`; `x0, xT` [nx x M] full hover/rest states; `geom = struct('M','r',0.06,'d',3,'which','T',4,'eps',1e-5,'l',1e-4,'ball',[])`.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_quad_layouts()
%% Purpose:
%
%   quad_layouts returns hover states of the right size for each model,
%   the line and circle positions the spec names, r = 0.06, T = 4, and
%   refuses an unknown layout by identifier.
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

[x0, xT, g] = quad_layouts('line', 'rigid6dof');
ok = chk(ok, isequal(size(x0), [12 1]) && isequal(x0, [-0.4; 0; 0.5; zeros(9, 1)]) && isequal(xT, [0.4; 0; 0.5; zeros(9, 1)]), 'rigid line: hover-to-hover states');
ok = chk(ok, g.M == 1 && g.r == 0.06 && g.T == 4 && g.d == 3 && isempty(g.ball), 'geom: M = 1, r = 0.06, T = 4');
[x0, xT, g] = quad_layouts('circle', 'rigid6dof');
ok = chk(ok, isequal(size(x0), [12 4]) && g.M == 4 && all(abs(sqrt(sum(x0(1:2, :).^2, 1)) - 0.4) < 1e-15) && all(x0(3, :) == 0.5), 'rigid circle: 4 vehicles, radius 0.4, z = 0.5');
ok = chk(ok, max(max(abs(xT(1:2, :) + x0(1:2, :)))) < 1e-15 && all(xT(3, :) == 0.5) && ~any(any(x0(4:12, :))) && ~any(any(xT(4:12, :))), 'rigid circle: antipodal in the plane, at rest, level');
[x0, xT, g] = quad_layouts('circle', 'pointmass');
ok = chk(ok, isequal(size(x0), [6 4]) && ~any(any(x0(4:6, :))) && ~any(any(xT(4:6, :))) && g.M == 4, 'pointmass circle: 6-state rest-to-rest');
try
    quad_layouts('spiral', 'rigid6dof');
    ok = chk(ok, false, 'unknown layout refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'quad_layouts:which'), 'unknown layout refused by quad_layouts:which');
end

if ok, fprintf('TEST_QUAD_LAYOUTS: ALL PASS\n'); else, fprintf('TEST_QUAD_LAYOUTS: FAILURE\n'); end
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

Expected: `Unrecognized function or variable 'quad_layouts'`.

- [ ] **Step 3: Implement**

```matlab
function [x0, xT, geom] = quad_layouts(which, model)
%% Purpose:
%
%   Boundary states for the quadrotor problems: hover-to-hover (rigid6dof:
%   v = 0, Theta = 0, omega = 0; pointmass: v = 0) at the spec's positions.
%
%   'line'    M = 1, (-0.4, 0, 0.5) -> (0.4, 0, 0.5)
%   'circle'  M = 4 on a circle of radius 0.4 at z = 0.5, antipodal targets
%
%   geom carries r = 0.06 m (Crazyflie span plus margin), T = 4 s (at the
%   paper's 10 s a quadrotor barely leaves hover), d = 3, and the barrier
%   constants for a future barrier run.
%
%% Inputs:
%
%  which                    char                    'line' | 'circle'
%
%  model                    char                    'pointmass' | 'rigid6dof'
%
%% Outputs:
%
%  x0, xT                   [nx x M]                Full boundary states
%
%  geom                     struct                  .M .r .d .which .T
%                                                   .eps .l .ball
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo:
     [a, b, g] = quad_layouts('circle', 'rigid6dof');
     fprintf('M = %d, r = %.2f, T = %g; p0 of vehicle 1 = [%s]\n', g.M, g.r, g.T, sprintf('%.2f ', a(1:3, 1)));
     disp(b(1:3, :));
     return
end

p = drone_params(model);
switch lower(which)
    case 'line'
        P0 = [-0.4; 0; 0.5];  PT = [0.4; 0; 0.5];
    case 'circle'
        ang = 2*pi*(0:3)/4;
        P0 = [0.4*cos(ang); 0.4*sin(ang); 0.5*ones(1, 4)];
        PT = [-P0(1:2, :); P0(3, :)];
    otherwise
        error('quad_layouts:which', 'unknown layout ''%s''', which);
end
M = size(P0, 2);
x0 = zeros(p.nx, M);  xT = zeros(p.nx, M);
x0(1:3, :) = P0;  xT(1:3, :) = PT;
geom = struct('M', M, 'r', 0.06, 'd', 3, 'which', lower(which), 'T', 4, 'eps', 1e-5, 'l', 1e-4, 'ball', []);
end
```

- [ ] **Step 4: Run the test** — all `PASS`.

- [ ] **Step 5: Commit**

```bash
cd ~/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/quadrotor/quad_layouts.m optimal_control_examples/drone/quadrotor/tests/test_quad_layouts.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone quadrotor: quad_layouts (hover-to-hover line and circle boundary states), with test

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 2: `solve_quads_direct`, `verify_quads_direct`, fixtures

**Files:**
- Create: `optimal_control_examples/drone/quadrotor/solve_quads_direct.m`
- Create: `optimal_control_examples/drone/quadrotor/verify_quads_direct.m`
- Create: `optimal_control_examples/drone/quadrotor/gen_quad_fixtures.m`
- Create (generated, force-added): `optimal_control_examples/drone/quadrotor/data/quad_direct_M1.mat`
- Test: `optimal_control_examples/drone/quadrotor/tests/test_quads_direct.m`

**Interfaces:**
- Produces: `out = solve_quads_direct(model, x0, xT, geom, opts)`, `opts.N [100] .pairs .midpoints [true] .maxIter [3000] .print [0] .bump [0.05] .sense [+1] .X0 .U0`; `out` fields `.model .p .X [nx M x N+1] .U [nu M x N+1] .tN .J .success .ipoptStatus .iterations .maxDefect .lamDefect [nx M x N] .boundsMargin (struct .rotor .angle .accel) .boundsActive (logical, any margin < 1e-6) .minHnodes .geom .pairs .N .T .wall`. `rep = verify_quads_direct(out, opts)` with gates `.G1local .G1clear .G3 .G4 .bounds` and `.ok/.why`. `gen_quad_fixtures()` writes `data/quad_direct_M1.mat` (fields `out`, `rep`, `x0`, `xT`, `geom`).

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_quads_direct()
%% Purpose:
%
%   The rigid-body direct solve, M = 1, and the claims made for it:
%     1. the CasADi dynamics Function equals drone_field at 20 random
%        (x, u) to 1e-14, for pointmass AND rigid6dof (one copy of physics);
%     2. the M = 1 line solve succeeds, defects < 1e-9, the local residual
%        (oc.local_residual, trapezoid control reconstruction) < 1e-4 per
%        interval, the rotor and angle boxes are INACTIVE (reported), the
%        fine-grid |Theta| stays under thetaMax, and the cost exceeds the
%        hover cost 4 (m g / 4)^2 T (you cannot move for free);
%     3. the pointmass M = 4 circle solves with fine-grid clearance >= the
%        bar;
%     4. the fixture (if present) re-verifies.
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
root = fileparts(fileparts(fileparts(here)));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'), fullfile(fileparts(here), 'route1_direct'), fullfile(root, 'oclib'));
rngWas = rng(53);  restore = onCleanup(@() rng(rngWas));

for mdl = {'pointmass', 'rigid6dof'}
    p = drone_params(mdl{1});
    fdyn = quad_dynamics_function(p);
    worst = 0;
    for k = 1:20
        x = randn(p.nx, 1);  if p.nx == 12, x(7:8) = 0.4*x(7:8); end
        u = rand(p.nu, 1);
        [F, B] = drone_field(x, p);
        worst = max(worst, max(abs(full(fdyn(x, u)) - (F + B*u))));
    end
    ok = chk(ok, worst < 1e-14, sprintf('%s: CasADi dynamics = drone_field (%.1e)', p.model, worst));
end

[x0, xT, geom] = quad_layouts('line', 'rigid6dof');
o1 = solve_quads_direct('rigid6dof', x0, xT, geom);
rep = verify_quads_direct(o1, struct('mesh', false));
ok = chk(ok, o1.success && rep.ok, sprintf('rigid M = 1: %s, gates ok (%s)', o1.ipoptStatus, rep.why));
ok = chk(ok, o1.maxDefect < 1e-9, sprintf('defects %.1e', o1.maxDefect));
ok = chk(ok, rep.G1local.maxMiss < 1e-4, sprintf('G1 local residual max %.2e per interval', rep.G1local.maxMiss));
fprintf('  INFO  bounds margin: rotor %.3e N, angle %.3e rad (active = %d)\n', o1.boundsMargin.rotor, o1.boundsMargin.angle, o1.boundsActive);
ok = chk(ok, ~o1.boundsActive, 'rotor and angle boxes inactive on the M = 1 line (the indirect comparison is fair)');
ok = chk(ok, rep.bounds.maxAngleFine <= o1.p.thetaMax, sprintf('fine-grid |Theta| max %.3f <= thetaMax', rep.bounds.maxAngleFine));
Jhover = 4*(o1.p.m*o1.p.g/4)^2*geom.T;
ok = chk(ok, o1.J > Jhover, sprintf('J = %.6e > hover cost %.6e', o1.J, Jhover));

[x0, xT, geom] = quad_layouts('circle', 'pointmass');
o4 = solve_quads_direct('pointmass', x0, xT, geom);
rep4 = verify_quads_direct(o4, struct('mesh', false));
ok = chk(ok, o4.success && rep4.ok && rep4.G1clear.minHfine >= -1e-2*2*geom.r, ...
         sprintf('pointmass M = 4 circle: ok, fine-grid clearance %.2e', rep4.G1clear.minHfine));

fx = fullfile(here, 'data', 'quad_direct_M1.mat');
if isfile(fx)
    S = load(fx);
    repF = verify_quads_direct(S.out, struct('mesh', false));
    ok = chk(ok, repF.ok && abs(S.out.J - o1.J)/o1.J < 1e-6, sprintf('fixture re-verifies and matches this solve (%.1e rel)', abs(S.out.J - o1.J)/o1.J));
else
    fprintf('  INFO  no fixture yet (gen_quad_fixtures writes it)\n');
end

if ok, fprintf('TEST_QUADS_DIRECT: ALL PASS\n'); else, fprintf('TEST_QUADS_DIRECT: FAILURE\n'); end
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

Expected: `Unrecognized function or variable 'quad_dynamics_function'`.

- [ ] **Step 3: Implement the dynamics Function, the solver, the verifier**

`quad_dynamics_function.m`:

```matlab
function fdyn = quad_dynamics_function(p)
%% Purpose:
%
%   The library dynamics as a CasADi Function f(x, u) = F(x) + B(x) u,
%   built by calling drone_field ITSELF on symbolic inputs: one copy of
%   the physics. test_quads_direct pins it to drone_field numerically.
%
%  ASSUMPTIONS / NOTES:
%
% • Requires drone_field to run on casadi.MX (sin/cos/tan, products,
%   concatenation with doubles). If a future edit to drone_field breaks
%   that, fix drone_field (keeping it complex-step safe); never copy the
%   equations here.
%
%% Inputs:
%
%  p                        struct                  drone_params(model)
%
%% Outputs:
%
%  fdyn                     casadi.Function         f(x [nx], u [nu]) -> xdot
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'drone_dynamics'), fullfile(fileparts(here), 'route1_direct'));
casadi_path();
import casadi.*
xs = MX.sym('x', p.nx);
us = MX.sym('u', p.nu);
[F, B] = drone_field(xs, p);
fdyn = Function('fdyn', {xs, us}, {F + B*us});
end
```

`solve_quads_direct.m`:

```matlab
function out = solve_quads_direct(model, x0, xT, geom, opts)
%% Purpose:
%
%   M vehicles of one drone model (pointmass or rigid6dof) from rest to
%   rest with pairwise clearance, as a trapezoidal NLP in CasADi/IPOPT:
%
%       min  integral sum_i |u_i|^2 dt                (rotor thrusts, or
%                                                      thrust accelerations)
%       s.t. xdot = F(x) + B(x) u (the library's drone_field, one copy),
%            x(0) = x0, x(T) = xT,
%            |p_i - p_j| >= 2r at every node and Hermite midpoint
%            rigid6dof: 0 <= f <= fMax,  |phi|, |theta| <= thetaMax
%            pointmass: |a|^2 <= aMax^2
%
%   Initial guess: positions on the Route 1 swirl lines, velocities their
%   derivative, attitude level, rates zero, thrust at hover (pointmass:
%   a = g e3).
%
%  ASSUMPTIONS / NOTES:
%
% • Defects first in the Opti: multipliers by row range (never opti.dual).
% • boundsActive reports whether any box is within 1e-6 of active; the
%   indirect comparison is fair only when it is false.
%
%% Inputs:
%
%  model                    char                    'pointmass' | 'rigid6dof'
%
%  x0, xT                   [nx x M]                From quad_layouts
%
%  geom                     struct                  From quad_layouts
%
%  opts                     struct (optional)       .N [100] .pairs
%                                                   .midpoints [true]
%                                                   .maxIter [3000]
%                                                   .print [0] .bump
%                                                   [0.05] .sense [+1]
%                                                   .X0 .U0
%
%% Outputs:
%
%  out                      struct                  .model .p .X .U .tN .J
%                                                   .success .ipoptStatus
%                                                   .iterations .maxDefect
%                                                   .lamDefect
%                                                   .boundsMargin
%                                                   .boundsActive
%                                                   .minHnodes .geom
%                                                   .pairs .N .T .wall
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: one Crazyflie across the room:
     [a, b, g] = quad_layouts('line', 'rigid6dof');
     o = solve_quads_direct('rigid6dof', a, b, g);
     fprintf('M = 1 rigid: J = %.6e, %s, bounds active %d\n', o.J, o.ipoptStatus, o.boundsActive);
     return
end
if nargin < 5, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'), fullfile(fileparts(here), 'route1_direct'));
casadi_path();
import casadi.*
d = @(f, v) fieldd(opts, f, v);
p      = drone_params(model);
N      = d('N', 100);
T      = geom.T;
M      = geom.M;
pairs  = d('pairs', pair_list(M));
midpts = d('midpoints', true);
maxIt  = d('maxIter', 3000);
prnt   = d('print', 0);
bump   = d('bump', 0.05);
sense  = d('sense', 1);
nx = p.nx;  nu = p.nu;  nN = N + 1;  h = T/N;  tN = linspace(0, T, nN);
t0 = tic;

%% Variables, dynamics, defects (first)
opti = Opti();
X = opti.variable(nx*M, nN);
U = opti.variable(nu*M, nN);
fdyn = quad_dynamics_function(p);
Fmap = fdyn.map(nN);
Xd = MX.zeros(nx*M, nN);
for a = 1:M
    ia = nx*(a-1) + (1:nx);  ja = nu*(a-1) + (1:nu);
    Xd(ia, :) = Fmap(X(ia, :), U(ja, :));
end
opti.subject_to(X(:, 2:end) - X(:, 1:end-1) - (h/2)*(Xd(:, 1:end-1) + Xd(:, 2:end)) == 0);
opti.subject_to(X(:, 1) == x0(:));
opti.subject_to(X(:, end) == xT(:));

%% Cost
L = sum1(U.^2);
J = h*sum2(L) - (h/2)*(L(1) + L(end));
opti.minimize(J);

%% Boxes
switch p.model
    case 'rigid6dof'
        opti.subject_to(0 <= U(:) <= p.fMax);
        angRows = reshape(nx*(0:M-1) + [7; 8], [], 1);
        opti.subject_to(-p.thetaMax <= X(angRows, :) <= p.thetaMax);
    case 'pointmass'
        for a = 1:M
            ja = nu*(a-1) + (1:nu);
            opti.subject_to(sum1(U(ja, :).^2) <= p.aMax^2);
        end
end

%% Clearance on positions (rows 1:3 of each vehicle), nodes and Hermite midpoints
posRows = reshape(nx*(0:M-1) + [1; 2; 3], [], 1);
velRows = reshape(nx*(0:M-1) + [4; 5; 6], [], 1);
Pn = X(posRows, :);  Vn = X(velRows, :);
Pmid = (Pn(:, 1:end-1) + Pn(:, 2:end))/2 + (h/8)*(Vn(:, 1:end-1) - Vn(:, 2:end));
if midpts, Ps = [Pn, Pmid]; else, Ps = Pn; end
P = size(pairs, 2);
if P > 0
    ia = reshape(3*(pairs(1,:) - 1) + [1; 2; 3], [], 1);
    ib = reshape(3*(pairs(2,:) - 1) + [1; 2; 3], [], 1);
    Dp = Ps(ia, :) - Ps(ib, :);
    sq = Dp(1:3:end, :).^2 + Dp(2:3:end, :).^2 + Dp(3:3:end, :).^2;
    opti.subject_to(sq(:) >= (2*geom.r)^2);
end

%% Initial guess
X0 = d('X0', []);  U0 = d('U0', []);
if isempty(X0) || isempty(U0)
    [X0, U0] = guess(x0, xT, tN, p, M, bump, sense);
end
opti.set_initial(X, X0);  opti.set_initial(U, U0);

%% Solve
ip = struct('print_level', prnt, 'max_iter', maxIt, 'tol', 1e-10, 'acceptable_tol', 1e-8, 'linear_solver', 'mumps');
opti.solver('ipopt', struct('print_time', false), ip);
try
    sol = opti.solve();  solved = true;
catch
    sol = opti.debug;    solved = false;
end

%% Report
out = struct('model', p.model, 'p', p);
out.X = full(sol.value(X));  out.U = full(sol.value(U));  out.tN = tN;
out.J = full(sol.value(J));
out.ipoptStatus = opti.stats().return_status;  out.iterations = opti.stats().iter_count;
out.success = solved && any(strcmp(out.ipoptStatus, {'Solve_Succeeded', 'Solved_To_Acceptable_Level'}));
Xdn = zeros(nx*M, nN);
for a = 1:M
    ia = nx*(a-1) + (1:nx);  ja = nu*(a-1) + (1:nu);
    for k = 1:nN
        [F, B] = drone_field(out.X(ia, k), p);
        Xdn(ia, k) = F + B*out.U(ja, k);
    end
end
def = out.X(:, 2:end) - out.X(:, 1:end-1) - (h/2)*(Xdn(:, 1:end-1) + Xdn(:, 2:end));
out.maxDefect = max(abs(def(:)));
lamAll = full(sol.value(opti.lam_g));
out.lamDefect = reshape(lamAll(1:nx*M*N), nx*M, N);
switch p.model
    case 'rigid6dof'
        out.boundsMargin = struct('rotor', min(min(out.U(:)), p.fMax - max(out.U(:))), ...
                                  'angle', p.thetaMax - max(max(abs(out.X(angRows, :)))), 'accel', Inf);
    case 'pointmass'
        amag = zeros(M, nN);
        for a = 1:M, amag(a, :) = sqrt(sum(out.U(nu*(a-1) + (1:nu), :).^2, 1)); end
        out.boundsMargin = struct('rotor', Inf, 'angle', Inf, 'accel', p.aMax - max(amag(:)));
end
out.boundsActive = min([out.boundsMargin.rotor, out.boundsMargin.angle, out.boundsMargin.accel]) < 1e-6;
if P > 0
    hn = clearance_terms(out.X(posRows, :), pairs, geom);  out.minHnodes = min(hn(:));
else
    out.minHnodes = Inf;
end
out.geom = geom;  out.pairs = pairs;  out.N = N;  out.T = T;  out.bump = bump;  out.sense = sense;
out.midpoints = midpts;  out.wall = toc(t0);
end

% ------------------------------------------------------------------------
function [X0, U0] = guess(x0, xT, tN, p, M, bump, sense)
%% Purpose:
%
%   Swirl lines for the positions, their derivative for the velocities,
%   level attitude, zero rates, hover thrust (pointmass: a = g e3).
%
nN = numel(tN);  T = tN(end);  s = tN/T;
X0 = zeros(p.nx*M, nN);  U0 = zeros(p.nu*M, nN);
for a = 1:M
    ia = p.nx*(a-1);  ja = p.nu*(a-1);
    dv = xT(1:3, a) - x0(1:3, a);
    L = sqrt(sum(dv.^2));  if L > 0, dir = dv/L; else, dir = [1; 0; 0]; end
    perp = [-dir(2); dir(1); 0];
    if sqrt(sum(perp.^2)) < 1e-6, perp = [0; -dir(3); dir(2)]; end
    perp = sense*perp/sqrt(sum(perp.^2));
    X0(ia + (1:3), :) = x0(1:3, a) + dv*s + bump*perp*sin(pi*s);
    X0(ia + (4:6), :) = dv/T + bump*perp*(pi/T)*cos(pi*s);
    switch p.model
        case 'rigid6dof', U0(ja + (1:4), :) = (p.m*p.g/4);
        case 'pointmass', U0(ja + 3, :) = p.g;
    end
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

`verify_quads_direct.m`:

```matlab
function rep = verify_quads_direct(out, opts)
%% Purpose:
%
%   The gates on one quadrotor direct solve, with its own verdict:
%     G3        IPOPT status ok, re-evaluated defects < 1e-9
%     G1local   oc.local_residual: restart the true dynamics from each left
%               node with the trapezoid's piecewise-linear control and
%               measure the miss at the right node; max over intervals
%               < 1e-4 ("defect is not accuracy")
%     G1clear   fine-grid pairwise clearance of the Hermite-reconstructed
%               positions (velocities as the derivative) >= -1e-2 * 2r
%     bounds    fine-grid |Theta| <= thetaMax; margins reported
%     G4        (opts.mesh) re-solve at 2N; cost change < 1e-3
%
%% Inputs:
%
%  out                      struct                  From solve_quads_direct
%
%  opts                     struct (optional)       .mesh [false] .nSub
%                                                   [10] .tolLocal [1e-4]
%                                                   .tolDefect [1e-9]
%                                                   .tolMesh [1e-3]
%
%% Outputs:
%
%  rep                      struct                  .ok .why .G3 .G1local
%                                                   .G1clear .bounds .G4
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(here)));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'), fullfile(root, 'oclib'));
d = @(f, v) fieldd(opts, f, v);
doMesh = d('mesh', false);  nSub = d('nSub', 10);
tolLoc = d('tolLocal', 1e-4);  tolDef = d('tolDefect', 1e-9);  tolMesh = d('tolMesh', 1e-3);
p = out.p;  M = out.geom.M;  nx = p.nx;  nu = p.nu;
rep = struct('ok', false, 'why', '');

rep.G3 = struct('status', out.ipoptStatus, 'maxDefect', out.maxDefect, 'ok', out.success && out.maxDefect < tolDef);

%% G1 local
uOf = @(t) interp1(out.tN, out.U.', t, 'linear').';
rhs = @(t, z) stacked_field(z, uOf(t), p, M);
[dX, ~] = oc.local_residual(out.X, out.tN, rhs);
rep.G1local = struct('maxMiss', max(abs(dX(:))), 'perInterval', max(abs(dX), [], 1), 'ok', max(abs(dX(:))) < tolLoc);

%% G1 clearance and bounds on the fine grid
posRows = reshape(nx*(0:M-1) + [1; 2; 3], [], 1);
velRows = reshape(nx*(0:M-1) + [4; 5; 6], [], 1);
[Pf, ~] = hermite_resample(out.X(posRows, :), out.X(velRows, :), out.tN, nSub);
if size(out.pairs, 2) > 0
    hf = clearance_terms(Pf, out.pairs, out.geom);  minH = min(hf(:));
else
    minH = Inf;
end
rep.G1clear = struct('minHfine', minH, 'bar', -1e-2*2*out.geom.r, 'ok', minH >= -1e-2*2*out.geom.r);
maxAng = 0;
if strcmp(p.model, 'rigid6dof')
    angRows = reshape(nx*(0:M-1) + [7; 8], [], 1);
    rateRows = reshape(nx*(0:M-1) + [10; 11], [], 1);        % only an approximation to Thetadot; the check is on the interpolant
    Af = hermite_resample(out.X(angRows, :), out.X(rateRows, :), out.tN, nSub);
    maxAng = max(abs(Af(:)));
end
rep.bounds = struct('margin', out.boundsMargin, 'active', out.boundsActive, 'maxAngleFine', maxAng, ...
                    'ok', ~strcmp(p.model, 'rigid6dof') || maxAng <= p.thetaMax + 1e-9);

%% G4 mesh
rep.G4 = struct('ran', false, 'J2N', NaN, 'relChange', NaN, 'ok', true);
if doMesh && out.success
    X2 = hermite_resample(out.X, zeros(size(out.X)), out.tN, 2);     % linear in the state is enough for a warm start
    U2 = interp1(out.tN, out.U.', linspace(0, out.T, 2*out.N + 1)).';
    o2 = solve_quads_direct(p.model, reshape(out.X(:, 1), nx, M), reshape(out.X(:, end), nx, M), out.geom, ...
                            struct('N', 2*out.N, 'pairs', out.pairs, 'midpoints', out.midpoints, 'X0', X2, 'U0', U2));
    rep.G4.ran = true;  rep.G4.J2N = o2.J;  rep.G4.relChange = abs(o2.J - out.J)/abs(out.J);
    rep.G4.ok = o2.success && rep.G4.relChange < tolMesh;
end

if ~rep.G3.ok,          rep.why = sprintf('G3: IPOPT %s, max defect %.2e', out.ipoptStatus, out.maxDefect);
elseif ~rep.G1local.ok, rep.why = sprintf('G1 local residual %.2e', rep.G1local.maxMiss);
elseif ~rep.G1clear.ok, rep.why = sprintf('G1 fine-grid clearance %.2e', minH);
elseif ~rep.bounds.ok,  rep.why = sprintf('fine-grid |Theta| %.3f exceeds thetaMax', maxAng);
elseif ~rep.G4.ok,      rep.why = sprintf('G4 cost changed %.2e at 2N', rep.G4.relChange);
end
rep.ok = isempty(rep.why);
end

% ------------------------------------------------------------------------
function dz = stacked_field(z, u, p, M)
%% Purpose:
%
%   The true dynamics of M stacked vehicles with a supplied control.
%
dz = zeros(p.nx*M, 1);
for a = 1:M
    ia = p.nx*(a-1) + (1:p.nx);  ja = p.nu*(a-1) + (1:p.nu);
    [F, B] = drone_field(z(ia), p);
    dz(ia) = F + B*u(ja);
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

`gen_quad_fixtures.m`:

```matlab
function gen_quad_fixtures()
%% Purpose:
%
%   Write data/quad_direct_M1.mat: the rigid-body M = 1 line solve with its
%   verification (mesh gate run), the fixture run_quad_pmp seeds from. Run
%   once; the .mat is force-added.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none (writes data/quad_direct_M1.mat)
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(here);
[x0, xT, geom] = quad_layouts('line', 'rigid6dof');
out = solve_quads_direct('rigid6dof', x0, xT, geom);
rep = verify_quads_direct(out, struct('mesh', true));
assert(out.success && rep.ok, 'gen_quad_fixtures:notOk', 'the M = 1 solve did not pass its gates: %s', rep.why);
if ~isfolder(fullfile(here, 'data')), mkdir(fullfile(here, 'data')); end
save(fullfile(here, 'data', 'quad_direct_M1.mat'), 'out', 'rep', 'x0', 'xT', 'geom');
fprintf('gen_quad_fixtures: M = 1 rigid, J = %.6e, bounds active %d, mesh change %.2e\n', out.J, out.boundsActive, rep.G4.relChange);
end
```

- [ ] **Step 4: Run the test, generate the fixture, run again**

Run the test (fixture absent: INFO line), then `gen_quad_fixtures`, then the test again.
Expected: all `PASS` both times, the second with the fixture line. If `fdyn` fails to build because `drone_field` does not accept `MX`, fix `drone_field` (e.g. replace `zeros(3, 1)` with `0*x(1:3)` so the drift is symbolic) and re-run Stage 0's tests.

- [ ] **Step 5: Commit (force-add the fixture)**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/quadrotor/quad_dynamics_function.m optimal_control_examples/drone/quadrotor/solve_quads_direct.m optimal_control_examples/drone/quadrotor/verify_quads_direct.m optimal_control_examples/drone/quadrotor/gen_quad_fixtures.m optimal_control_examples/drone/quadrotor/tests/test_quads_direct.m
/Library/Developer/CommandLineTools/usr/bin/git add -f optimal_control_examples/drone/quadrotor/data/quad_direct_M1.mat
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone quadrotor: solve_quads_direct (CasADi/IPOPT trapezoid on the library's drone_field for pointmass and rigid6dof, M vehicles, clearance at nodes + midpoints, rotor/angle/accel boxes), verify_quads_direct (local residual, fine-grid clearance, angle, mesh), M = 1 fixture

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 3: `quad_pmp_rhs` and `quad_pmp_prop`

**Files:**
- Create: `optimal_control_examples/drone/quadrotor/quad_pmp_rhs.m`
- Create: `optimal_control_examples/drone/quadrotor/quad_pmp_prop.m`
- Test: `optimal_control_examples/drone/quadrotor/tests/test_quad_pmp_rhs.m`
- Test: `optimal_control_examples/drone/quadrotor/tests/test_quad_pmp_prop.m`

**Interfaces:**
- Produces: `[dy, f, H] = quad_pmp_rhs(y, p)` for `y = [x; lam]` [24 x 1]: `f = −B'lam/2`, `dy = [F + B f; −A(x, f)' lam]`, `H = f'f + lam'(F + B f)`. `[yEnd, PHI] = quad_pmp_prop(dt, y0, needSTM, p)`: the cart-pole propagator pattern (complex-step STM through `quad_pmp_rhs`, 24 + 576 equations, collapse contract).

- [ ] **Step 1: Write the failing tests**

```matlab
function ok = test_quad_pmp_rhs()
%% Purpose:
%
%   The quadrotor PMP field: f* = -B'lam/2 is the stationary point of
%   H = |f|^2 + lam'(F + B f) (dH/df = 0 by a complex step of H in f);
%   Legendre d2H/df2 = 2I; the costate rate is -dH/dx by a complex step of
%   H in x (through drone_field and the GENERATED Jacobian -- an
%   independent differentiation of the same physics); shapes.
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
rngWas = rng(61);  restore = onCleanup(@() rng(rngWas));
p = drone_params('rigid6dof');
worstU = 0;  worstX = 0;  hc = 1e-20;
for trial = 1:10
    x = randn(12, 1);  x(7:8) = 0.4*x(7:8);  x(10:12) = 0.5*x(10:12);
    lam = randn(12, 1);
    [dy, f, H] = quad_pmp_rhs([x; lam], p);
    ok = ok && isequal(size(dy), [24 1]) && isequal(size(f), [4 1]);
    Hof = @(xx, ff) ff.'*ff + lam.'*(fieldAt(xx, ff, p));
    ok = ok && abs(H - Hof(x, f)) < 1e-13*max(1, abs(H));
    for k = 1:4
        fc = f;  fc(k) = fc(k) + 1i*hc;
        worstU = max(worstU, abs(imag(Hof(x, fc))/hc));
    end
    for k = 1:12
        xc = x;  xc(k) = xc(k) + 1i*hc;
        worstX = max(worstX, abs(imag(Hof(xc, f))/hc + dy(12 + k))/max(1, abs(dy(12 + k))));
    end
end
ok = chk(ok, worstU < 1e-12, sprintf('dH/df = 0 at f* (complex step of H in f: %.1e)', worstU));
ok = chk(ok, worstX < 1e-10, sprintf('lamdot = -dH/dx (complex step of H in x: %.1e rel)', worstX));
fprintf('  INFO  Legendre d2H/df2 = 2I identically (H quadratic in f)\n');
ok = chk(ok, ok, 'shapes and H consistent');

if ok, fprintf('TEST_QUAD_PMP_RHS: ALL PASS\n'); else, fprintf('TEST_QUAD_PMP_RHS: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function xd = fieldAt(x, u, p)
%% Purpose:
%
%   F + B u.
%
[F, B] = drone_field(x, p);
xd = F + B*u;
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

```matlab
function ok = test_quad_pmp_prop()
%% Purpose:
%
%   The propagator contract: dt = 0 returns y0 and I; the STM matches a
%   complex step through the propagator without the STM; H is conserved
%   along a flight; malformed input refused by identifier.
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
p = drone_params('rigid6dof');
x0 = [-0.4; 0; 0.5; zeros(9, 1)];
lam0 = [0.02; 0; -0.3; 0.05; 0; -0.2; zeros(6, 1)];        % gentle: a small horizontal push, thrust near hover
y0 = [x0; lam0];
[yE, PH] = quad_pmp_prop(0, y0, true, p);
ok = chk(ok, isequal(yE, y0) && isequal(PH, eye(24)), 'dt = 0: y0 and I');
dt = 0.5;
[yE, PH] = quad_pmp_prop(dt, y0, true, p);
ok = chk(ok, all(isfinite(yE)) && all(isfinite(PH(:))), 'finite over 0.5 s');
worst = 0;  hc = 1e-20;
for k = 1:24
    yc = y0;  yc(k) = yc(k) + 1i*hc;
    yEc = quad_pmp_prop(dt, yc, false, p);
    worst = max(worst, max(abs(imag(yEc)/hc - PH(:, k)))/max(1, max(abs(PH(:, k)))));
end
ok = chk(ok, worst < 1e-7, sprintf('STM vs complex step through the propagator (%.1e rel)', worst));
[~, ~, H0] = quad_pmp_rhs(y0, p);  [~, ~, HE] = quad_pmp_rhs(yE, p);
ok = chk(ok, abs(HE - H0) < 1e-9*max(1, abs(H0)), sprintf('H conserved (%.1e rel)', abs(HE - H0)/max(1, abs(H0))));
try
    quad_pmp_prop(1, y0(1:23), false, p);
    ok = chk(ok, false, 'short y0 refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'quad_pmp_prop:input'), 'short y0 refused by quad_pmp_prop:input');
end

if ok, fprintf('TEST_QUAD_PMP_PROP: ALL PASS\n'); else, fprintf('TEST_QUAD_PMP_PROP: FAILURE\n'); end
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

- [ ] **Step 2: Run them to verify they fail** — `quad_pmp_rhs` undefined.

- [ ] **Step 3: Implement**

```matlab
function [dy, f, H] = quad_pmp_rhs(y, p)
%% Purpose:
%
%   The Pontryagin field of the quadrotor minimum-effort problem at fixed
%   final time: y = [x; lam], 24 states, cost integral |f|^2 dt,
%   xdot = F(x) + B(x) f (the library's drone_field).
%
%       H = f'f + lam'(F + B f),   dH/df = 2f + B'lam = 0  ->  f* = -B'lam/2
%       d2H/df2 = 2I > 0 (a minimum everywhere; no switching structure)
%       lamdot = -dH/dx = -A(x, f*)' lam,  A = drone_state_jac (generated)
%
%   Unconstrained rotors: the comparison with the direct solve is fair
%   only when the direct solution's boxes are inactive (reported there).
%
%  ASSUMPTIONS / NOTES:
%
% • Complex-step safe: A is the generated Jacobian, so a caller may
%   complex-step through this function (that is how the STM is built).
%
%% Inputs:
%
%  y                        [24 x 1]                [x; lam]
%
%  p                        struct                  drone_params('rigid6dof')
%
%% Outputs:
%
%  dy                       [24 x 1]                [xdot; lamdot]
%
%  f                        [4 x 1]                 The PMP rotor thrusts
%
%  H                        double                  The Hamiltonian
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the field at hover with a costate that asks for lift:
     pr = drone_params('rigid6dof');
     [dy, f] = quad_pmp_rhs([zeros(12, 1); 0; 0; -1; 0; 0; -0.5; zeros(6, 1)], pr);
     fprintf('f* = [%s] N, vdot = [%s]\n', sprintf('%.4f ', f), sprintf('%.3f ', dy(4:6)));
     return
end

x = y(1:12);  lam = y(13:24);
[F, B] = drone_field(x, p);
f = -(B.'*lam)/2;
A = drone_state_jac(x, f, p);
xd = F + B*f;
dy = [xd; -A.'*lam];
H = f.'*f + lam.'*xd;
end
```

```matlab
function [yEnd, PHI] = quad_pmp_prop(dt, y0, needSTM, p)
%% Purpose:
%
%   Propagate the quadrotor PMP flow with the 24x24 STM on request: the
%   oc.ms_bvp prob.prop contract, in the cart-pole propagator's shape. The
%   variational equations ride alongside the state (24 + 576 equations)
%   with the Jacobian by COMPLEX STEP through quad_pmp_rhs, which is safe
%   only because the state Jacobian inside the field is generated, not
%   complex-stepped.
%
%  ASSUMPTIONS / NOTES:
%
% • Tolerances 1e-12 / 1e-14. THROWS quad_pmp_prop:collapse on a truncated
%   grid or non-finite result (ode113 warns and stops early on blow-up);
%   quad_pmp_prop:input on malformed input. Real y0 required with the STM.
%
%% Inputs:
%
%  dt                       double                  Propagation time
%
%  y0                       [24 x 1]                [x; lam]
%
%  needSTM                  logical                 Return PHI
%
%  p                        struct                  drone_params('rigid6dof')
%
%% Outputs:
%
%  yEnd                     [24 x 1]                y(dt)
%
%  PHI                      [24 x 24] or []         d y(dt)/d y0
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: half a second of hover-ish flow:
     pr = drone_params('rigid6dof');
     [yE, PH] = quad_pmp_prop(0.5, [zeros(12, 1); 0; 0; -0.3; 0; 0; -0.2; zeros(6, 1)], true, pr);
     fprintf('z(0.5) = %.4f, cond(PHI) = %.3e\n', yE(3), cond(PH));
     return
end

if ~(isnumeric(y0) && numel(y0) == 24 && all(isfinite(y0(:))))
    error('quad_pmp_prop:input', 'y0 must be 24 finite values');
end
if needSTM && ~isreal(y0)
    error('quad_pmp_prop:input', 'y0 must be real when the STM is requested');
end
if ~(isnumeric(dt) && isreal(dt) && isscalar(dt) && isfinite(dt))
    error('quad_pmp_prop:input', 'dt must be a real finite scalar');
end
y0 = y0(:);
if dt == 0
    yEnd = y0;  PHI = [];  if needSTM, PHI = eye(24); end
    return
end
oo = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
if needSTM
    z0 = [y0; reshape(eye(24), 576, 1)];
    [T, Z] = ode113(@(t, z) rhs_with_stm(z, p), [0 dt], z0, oo);
    zEnd = Z(end, :).';  yEnd = zEnd(1:24);  PHI = reshape(zEnd(25:end), 24, 24);
else
    [T, Y] = ode113(@(t, y) quad_pmp_rhs(y, p), [0 dt], y0, oo);
    yEnd = Y(end, :).';  PHI = [];
end
if ~(T(end) == dt) || ~all(isfinite(yEnd)) || (needSTM && ~all(isfinite(PHI(:))))
    error('quad_pmp_prop:collapse', 'the PMP flow left the finite range over dt = %g (reached t = %g)', dt, T(end));
end
end

% ------------------------------------------------------------------------
function dz = rhs_with_stm(z, p)
%% Purpose:
%
%   zdot = [f(y); A(y) PHI], A by complex step through quad_pmp_rhs.
%
y = z(1:24);  PHI = reshape(z(25:end), 24, 24);
dy = quad_pmp_rhs(y, p);
A = zeros(24);  hc = 1e-20;
for col = 1:24
    yp = y;  yp(col) = yp(col) + 1i*hc;
    A(:, col) = imag(quad_pmp_rhs(yp, p))/hc;
end
dz = [dy; reshape(A*PHI, 576, 1)];
end
```

- [ ] **Step 4: Run both tests** — all `PASS`.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/quadrotor/quad_pmp_rhs.m optimal_control_examples/drone/quadrotor/quad_pmp_prop.m optimal_control_examples/drone/quadrotor/tests/test_quad_pmp_rhs.m optimal_control_examples/drone/quadrotor/tests/test_quad_pmp_prop.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone quadrotor: quad_pmp_rhs (f* = -B'lam/2, costate rate from the generated Jacobian, checked against complex steps of H) and quad_pmp_prop (complex-step STM, collapse contract)

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 4: `run_quad_pmp`

**Files:**
- Create: `optimal_control_examples/drone/quadrotor/run_quad_pmp.m`
- Test: `optimal_control_examples/drone/quadrotor/tests/test_quad_pmp.m`

**Interfaces:**
- Produces: `out = run_quad_pmp(opts)`, `opts.K [8] .tolR [1e-10] .plot [nargout == 0] .fixture [data/quad_direct_M1.mat] .direct (an out struct overriding the fixture)`; `out` fields `.ok .why .lam0 [12 x 1] .J .Jdirect .boundsActiveDirect .missTerminal .missEngine .missFlown .Hrel .adjRel .seed (.signCorr .flipped .ampRatio) .conj (.verdict .out) .notes .t .X .F .Lam .info`.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_quad_pmp()
%% Purpose:
%
%   The quadrotor indirect solve on the M = 1 fixture: N1-N5, X2 (when the
%   direct bounds are inactive), X3 (K = 8 vs 16), S2 printed and not
%   FAIL; a negated multiplier seed is flipped back; an uninformative seed
%   is refused.
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
root = fileparts(fileparts(fileparts(here)));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'), fullfile(root, 'oclib'));

o8 = run_quad_pmp(struct('K', 8, 'plot', false));
ok = chk(ok, o8.ok, sprintf('front door ok (%s)', o8.why));
ok = chk(ok, o8.info.converged && o8.info.normR < 1e-10, sprintf('N1 |R| = %.2e', o8.info.normR));
ok = chk(ok, o8.missTerminal < 1e-9, sprintf('N2 terminal miss %.2e', o8.missTerminal));
ok = chk(ok, o8.missFlown < 1e-6, sprintf('N3 flown miss %.2e', o8.missFlown));
ok = chk(ok, o8.Hrel < 1e-9, sprintf('N4 H spread %.2e', o8.Hrel));
ok = chk(ok, o8.adjRel < 1e-9, sprintf('N5 adjoint %.2e', o8.adjRel));
ok = chk(ok, o8.seed.signCorr > 0.9, sprintf('seed informative: corr %.3f, amplitude ratio %.3f', o8.seed.signCorr, o8.seed.ampRatio));
if o8.boundsActiveDirect
    fprintf('  INFO  X2 skipped: the direct solution has an active box (%s)\n', strjoin(o8.notes, '; '));
else
    ok = chk(ok, abs(o8.J - o8.Jdirect)/o8.Jdirect < 1e-2, sprintf('X2 J = %.6e vs direct %.6e', o8.J, o8.Jdirect));
end
fprintf('  INFO  S2 conjugate verdict: %s\n', o8.conj.verdict);
ok = chk(ok, ~strcmp(o8.conj.verdict, 'FAIL'), 'S2 not FAIL');
o16 = run_quad_pmp(struct('K', 16, 'plot', false));
ok = chk(ok, o16.ok && max(abs(o16.lam0 - o8.lam0)) < 1e-8, sprintf('X3 lam0 K = 8 vs 16: %.2e', max(abs(o16.lam0 - o8.lam0))));

S = load(fullfile(here, 'data', 'quad_direct_M1.mat'));
neg = S.out;  neg.lamDefect = -neg.lamDefect;
on = run_quad_pmp(struct('K', 8, 'plot', false, 'direct', neg));
ok = chk(ok, on.ok && on.seed.flipped && max(abs(on.lam0 - o8.lam0)) < 1e-8, 'negated seed flipped back to the same root');
bad = S.out;  bad.lamDefect = 1e-6*randn(size(bad.lamDefect));
try
    run_quad_pmp(struct('K', 8, 'plot', false, 'direct', bad));
    ok = chk(ok, false, 'uninformative seed refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'run_quad_pmp:sign'), 'uninformative seed refused by run_quad_pmp:sign');
end

if ok, fprintf('TEST_QUAD_PMP: ALL PASS\n'); else, fprintf('TEST_QUAD_PMP: FAILURE\n'); end
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

- [ ] **Step 2: Run it to verify it fails** — `run_quad_pmp` undefined.

- [ ] **Step 3: Implement**

```matlab
function out = run_quad_pmp(opts)
%% Purpose:
%
%   The quadrotor minimum-effort transfer solved INDIRECTLY: the 24-state
%   PMP boundary-value problem (quad_pmp_rhs) shot with oc.ms_bvp at
%   fixed T, seeded from the direct solve's defect multipliers through
%   oc.duals_to_costates (trapezoid stations), with the global sign
%   resolved by correlating the implied control f* = -B'lam/2 against the
%   direct control (the cart-pole's scalar-force pattern applied to four
%   rotors). Twelve unknowns lam(0) against twelve terminal conditions.
%
%   Gates (own verdict): N1 converged, N2 terminal miss (re-flight at
%   RelTol 2.5e-14), N3 flown miss (oc.fly_control with f(t) in closed
%   form from the costate arc), N4 H spread, N5 adjoint vs complex step,
%   X2 cost vs the direct solve ONLY when its boxes were inactive (else a
%   note), S2 conjugate test (FAIL fails the door).
%
%% Inputs:
%
%  opts                     struct (optional)       .K [8] .tolR [1e-10]
%                                                   .plot [nargout == 0]
%                                                   .fixture [data/
%                                                   quad_direct_M1.mat]
%                                                   .direct (a
%                                                   solve_quads_direct out
%                                                   struct, M = 1,
%                                                   overriding the fixture)
%
%% Outputs:
%
%  out                      struct                  .ok .why .lam0 .J
%                                                   .Jdirect
%                                                   .boundsActiveDirect
%                                                   .missTerminal
%                                                   .missEngine .missFlown
%                                                   .Hrel .adjRel .seed
%                                                   .conj .notes .t .X .F
%                                                   .Lam .info
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(here)));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'), fullfile(root, 'oclib'));
d = @(f, v) fieldd(opts, f, v);
K = d('K', 8);  tolR = d('tolR', 1e-10);  doPlot = d('plot', nargout == 0);
R = d('direct', []);
if isempty(R)
    S = load(d('fixture', fullfile(here, 'data', 'quad_direct_M1.mat')));
    R = S.out;
end
if R.geom.M ~= 1 || ~strcmp(R.model, 'rigid6dof')
    error('run_quad_pmp:fixture', 'the indirect solve is for ONE rigid6dof vehicle');
end
p = R.p;  T = R.T;  tN = R.tN;  x0 = R.X(:, 1);  xT = R.X(:, end);
notes = {};

%% Seed: multipliers -> costates, sign by the implied control
[lamS, tS] = oc.duals_to_costates(struct('scheme', 'trapezoid', 'mu', R.lamDefect, 'tNodes', tN, 'velRows', 4:6));
tS = tS(:).';
Xs = interp1(tN, R.X.', tS, 'pchip').';
fImp = zeros(4, numel(tS));
for k = 1:numel(tS)
    [~, Bk] = drone_field(Xs(:, k), p);
    fImp(:, k) = -(Bk.'*lamS(:, k))/2;
end
fDir = interp1(tN, R.U.', tS, 'pchip').';
% the hover thrust is common mode; correlate the DEVIATIONS from hover so a
% sign flip is visible against a near-constant background
fh = p.m*p.g/4;
a = fImp(:) - fh;  b = fDir(:) - fh;
den = sqrt(sum(a.^2)*sum(b.^2));
if ~(den > 0) || abs(sum(a.*b))/den < 0.5
    error('run_quad_pmp:sign', 'the seed carries no usable sign information (corr %.3f)', sum(a.*b)/max(den, eps));
end
signCorr = sum(a.*b)/den;
if signCorr < 0, lamS = -lamS;  a = -fImp(:) - fh; end
seedDiag = struct('signCorr', abs(signCorr), 'flipped', signCorr < 0, 'ampRatio', sqrt(sum(a.^2)/sum(b.^2)));

tGrid = linspace(0, T, K + 1);
Xg = interp1(tN, R.X.', tGrid, 'pchip').';
Lg = interp1(tS, lamS.', tGrid, 'pchip', 'extrap').';
seed = struct('tf', T, 'tGrid', tGrid, 'Y', [Xg; Lg]);
seed.Y(1:12, 1) = x0;

%% The problem, as three closures
prob = struct('ny', 24, 'freeIdx0', 13:24, ...
    'prop', @(dt, y0, needSTM) quad_pmp_prop(dt, y0, needSTM, p), ...
    'rhs',  @(y) quad_pmp_rhs(y, p), ...
    'terminal', @(y, needJ) terminalFcn(y, xT, needJ));
[~, info] = oc.ms_bvp(prob, seed, struct('fixedTf', true, 'tolR', tolR, 'maxIter', 60, 'polishMax', 5, 'keepSTMs', true));

%% Report
lam0 = info.Y(13:24, 1);
[t, Y] = ode113(@(tt, y) quad_pmp_rhs(y, p), linspace(0, T, 2001), [x0; lam0], odeset('RelTol', 2.5e-14, 'AbsTol', 1e-16));
reachedTf = numel(t) == 2001 && t(end) == T && all(isfinite(Y(:)));
X = Y(:, 1:12).';  Lam = Y(:, 13:24).';
Fr = zeros(4, numel(t));  H = zeros(1, numel(t));  adjRel = 0;  hc = 1e-20;
for k = 1:numel(t)
    [dy, Fr(:, k), H(k)] = quad_pmp_rhs(Y(k, :).', p);
    if mod(k, 200) == 1
        for q = 1:12
            yc = Y(k, :).';  yc(q) = yc(q) + 1i*hc;
            [~, ~, Hc] = quad_pmp_rhs(yc, p);
            adjRel = max(adjRel, abs(dy(12 + q) + imag(Hc)/hc)/max(1e-6, abs(imag(Hc)/hc)));
        end
    end
end
J = trapz(t, sum(Fr.^2, 1));
Hrel = (max(H) - min(H))/max(abs(H));
ppL = pchip(t, Lam);
zEnd = oc.fly_control(x0, [0 T], @(tt, x) flyRhs(x, ppval(ppL, min(max(tt, 0), T)), p), ...
                      struct('mode', 'span', 'solver', @ode113, 'RelTol', 1e-13, 'AbsTol', 1e-15));
missTerminal = max(abs(X(:, end) - xT));
missFlown = max(abs(zEnd - xT));
dtLast = info.tGrid(end) - info.tGrid(end-1);
yLast = quad_pmp_prop(dtLast, info.Y(:, end), false, p);
missEngine = max(abs(yLast(1:12) - xT));
conj = struct('verdict', 'NOT RUN', 'out', []);
if info.converged
    co = oc.ms_conjugate_test(info, struct('stateRows', 1:12, 'costateCols', 13:24, 'quotientDir', [], 'freeTime', false));
    conj.out = co;  conj.verdict = co.verdict;
end
if R.boundsActive
    notes{end+1} = 'direct solution has an active box: the indirect (unbounded) cost is a different problem; X2 skipped';
end
why = '';
if ~info.converged,            why = sprintf('engine not converged (|R| = %.2e)', info.normR);
elseif ~reachedTf,             why = 'reporting flight truncated or non-finite';
elseif ~(missTerminal < 1e-9), why = sprintf('terminal miss %.2e', missTerminal);
elseif ~(missFlown < 1e-6),    why = sprintf('flown miss %.2e', missFlown);
elseif ~(Hrel < 1e-9),         why = sprintf('H spread %.2e', Hrel);
elseif ~(adjRel < 1e-9),       why = sprintf('adjoint residual %.2e', adjRel);
elseif ~R.boundsActive && ~(abs(J - R.J)/R.J < 1e-2), why = sprintf('J = %.6e vs direct %.6e', J, R.J);
elseif strcmp(conj.verdict, 'FAIL'), why = 'conjugate point: FAIL';
end
out = struct('ok', isempty(why), 'why', why, 'lam0', lam0, 'J', J, 'Jdirect', R.J, 'boundsActiveDirect', R.boundsActive, ...
             'missTerminal', missTerminal, 'missEngine', missEngine, 'missFlown', missFlown, 'Hrel', Hrel, 'adjRel', adjRel, ...
             'seed', seedDiag, 'conj', conj, 'notes', {notes}, 't', t.', 'X', X, 'F', Fr, 'Lam', Lam, 'info', info);
if ~out.ok, warning('run_quad_pmp:notOk', 'NOT a solution: %s', why); end
if doPlot, plotQuad(out, R); end
if nargout == 0
    fprintf('quadrotor PMP-BVP (%s): J = %.6e (direct %.6e), terminal miss %.2e, conj %s\n', ...
            tern(out.ok, 'ok', ['NOT OK: ' why]), J, R.J, missTerminal, conj.verdict);
    clear out
end
end

% ------------------------------------------------------------------------
function [g, dgdy] = terminalFcn(y, xT, needJ)
%% Purpose:
%
%   x(T) = xT and its Jacobian.
%
g = y(1:12) - xT;
dgdy = [];  if needJ, dgdy = [eye(12), zeros(12)]; end
end

function dx = flyRhs(x, lam, p)
%% Purpose:
%
%   The TRUE dynamics with the PMP control from the costate arc and the
%   FLOWN state's B.
%
[F, B] = drone_field(x, p);
dx = F + B*(-(B.'*lam)/2);
end

function plotQuad(out, R)
%% Purpose:
%
%   Position, attitude and rotor thrusts, indirect vs direct.
%
figure('color', [1 1 1], 'Position', [100 100 1400 400]);
subplot(1, 3, 1);  plot(R.tN, R.X(1:3, :), '--', out.t, out.X(1:3, :), '-');  grid on;  title('position');
subplot(1, 3, 2);  plot(R.tN, R.X(7:9, :)*180/pi, '--', out.t, out.X(7:9, :)*180/pi, '-');  grid on;  title('Euler angles (deg)');
subplot(1, 3, 3);  plot(R.tN, R.U, '--', out.t, out.F, '-');  grid on;  title('rotor thrusts (N)');
end

function s = tern(c, a, b)
%% Purpose:
%
%   a if c, else b.
%
if c, s = a; else, s = b; end
end

function v = fieldd(s, f, v)
%% Purpose:
%
%   s.(f) if present and non-empty, else the default v.
%
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
```

- [ ] **Step 4: Run the test** — all `PASS`; the conjugate verdict is printed verbatim.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/quadrotor/run_quad_pmp.m optimal_control_examples/drone/quadrotor/tests/test_quad_pmp.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone quadrotor: run_quad_pmp front door (covector seed with sign by implied-control correlation, oc.ms_bvp fixed-tf on 12 costates, N1-N5/X2/X3/S2 gates, fairness note on active boxes), tested on the M = 1 fixture

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 5: Study script, movie, comparison table, runner, README

**Files:**
- Create: `optimal_control_examples/drone/quadrotor/quad_minenergy_study.m` (the cart-pole study's ten sections on the quadrotor: write it by transcribing `route2_indirect/drone_minenergy_study.m`'s structure with `quad_pmp_rhs`/`quad_pmp_prop`, `prob` as in `run_quad_pmp`, the oracle run inline in section 1 via `test_drone_physics`, and section 8 asserting against `run_quad_pmp`; the parameter block is `K = 8`, `selfCheck.strict = false`; it leaves `verdict` with fields `N1 N2 N3 N4 N5 X2 X3 libraryAgrees` and `study` with `conjVerdict`)
- Create: `optimal_control_examples/drone/quadrotor/movie_quads.m` (1280×720 MP4: each vehicle drawn as a "+" of two body-axis arms of length `l` scaled ×3, rotated by `R(Θ)`, over faint paths; `VideoWriter('MPEG-4')`, frames resized to 720×1280 as `movie_drones` does)
- Create: `optimal_control_examples/drone/quadrotor/quad_compare.m` (rows: `pointmass M = 4 direct`, `rigid6dof M = 1 direct`, `rigid6dof M = 1 indirect`, `rigid6dof M = 4 direct`; columns `model, M, route, J, Jhover, minHfine, boundsActive, localResidual, conj`; writes `results/quad_compare.csv`)
- Test: `optimal_control_examples/drone/quadrotor/tests/test_quad_study.m` (the Route 2 study test transcribed: runs the script, asserts the verdict fields above)
- Modify: `optimal_control_examples/drone/run_tests.m` (add `quadrotor` paths and the seven tests: `test_quad_layouts`, `test_quads_direct`, `test_quad_pmp_rhs`, `test_quad_pmp_prop`, `test_quad_pmp`, `test_quad_study`), `optimal_control_examples/drone/README.md` (Stage 3 row "built", the `quad_compare.csv` table)
- Create: `optimal_control_examples/drone/quadrotor/README.md`

- [ ] **Step 1: Write `quad_minenergy_study.m` and `test_quad_study.m`** following `route2_indirect/drone_minenergy_study.m` and its test one-for-one (sections 0-9; section 1 calls `test_drone_physics` and records its logical as `verdict.V1`; section 5 perturbs the costate seed by +5 % and requires the same `lam0` to 1e-8; section 6 flies once at RelTol 2.5e-14 and computes N2-N5 from that flight; section 7 prints Legendre `2I` and the conjugate verdict; section 8 asserts `lam0`, `J` and the verdict against `run_quad_pmp(struct('K', K, 'plot', false))`; section 9 plots position, angles, thrusts).

- [ ] **Step 2: Write `movie_quads.m`**

```matlab
function movie_quads(out, file, opts)
%% Purpose:
%
%   An MP4 of a quadrotor solve: each vehicle as a "+" of two body arms
%   (rotated by R(Theta), arms drawn 3x the true 46 mm so they read at
%   room scale) over faint paths. Frames 1280 x 720 (a multiple of 16).
%
%% Inputs:
%
%  out                      struct                  From solve_quads_direct
%                                                   (rigid6dof) or
%                                                   run_quad_pmp (.X, .t)
%
%  file                     char                    Output .mp4
%
%  opts                     struct (optional)       .fps [30] .nFrames [240]
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
nF = 240;  if isfield(opts, 'nFrames'), nF = opts.nFrames; end
if isfield(out, 'tN'), tt0 = out.tN; else, tt0 = out.t; end
p = drone_params('rigid6dof');
M = size(out.X, 1)/12;
tt = linspace(tt0(1), tt0(end), nF);
Xi = interp1(tt0, out.X.', tt, 'pchip').';
fig = figure('color', [1 1 1], 'Visible', 'off', 'Position', [100 100 1280 720]);
ax = axes(fig);  hold(ax, 'on');  grid(ax, 'on');  axis(ax, 'equal');  view(ax, 3);
xlim(ax, [-0.6 0.6]);  ylim(ax, [-0.6 0.6]);  zlim(ax, [0 1]);
cols = lines(M);  arm = 3*p.l;
for a = 1:M, ia = 12*(a-1) + (1:3);  plot3(ax, Xi(ia(1), :), Xi(ia(2), :), Xi(ia(3), :), '-', 'Color', [0.85 0.85 0.85]); end
hA = gobjects(M, 2);
for a = 1:M
    hA(a, 1) = plot3(ax, NaN, NaN, NaN, '-', 'Color', cols(a, :), 'LineWidth', 2);
    hA(a, 2) = plot3(ax, NaN, NaN, NaN, '-', 'Color', cols(a, :), 'LineWidth', 2);
end
vw = VideoWriter(file, 'MPEG-4');  vw.FrameRate = fps;  open(vw);
for f = 1:nF
    for a = 1:M
        ia = 12*(a-1);
        pos = Xi(ia + (1:3), f);  Rm = euler_zyx(Xi(ia + (7:9), f));
        e1 = Rm(:, 1)*arm;  e2 = Rm(:, 2)*arm;
        set(hA(a, 1), 'XData', pos(1) + [-e1(1) e1(1)], 'YData', pos(2) + [-e1(2) e1(2)], 'ZData', pos(3) + [-e1(3) e1(3)]);
        set(hA(a, 2), 'XData', pos(1) + [-e2(1) e2(1)], 'YData', pos(2) + [-e2(2) e2(2)], 'ZData', pos(3) + [-e2(3) e2(3)]);
    end
    title(ax, sprintf('t = %.2f s', tt(f)));
    fr = getframe(fig);
    writeVideo(vw, imresize(fr.cdata, [720 1280]));
end
close(vw);  close(fig);
end
```

- [ ] **Step 3: Write `quad_compare.m`**

```matlab
function tbl = quad_compare(opts)
%% Purpose:
%
%   The model comparison on the circle task: pointmass M = 4 direct,
%   rigid6dof M = 1 direct and indirect, rigid6dof M = 4 direct; written to
%   results/quad_compare.csv. The hover cost 4 (m g/4)^2 T per vehicle is
%   the floor every rigid-body row must exceed.
%
%% Inputs:
%
%  opts                     struct (optional)       .csv [results/
%                                                   quad_compare.csv]
%
%% Outputs:
%
%  tbl                      table                   One row per case
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
csvOut = fullfile(here, 'results', 'quad_compare.csv');  if isfield(opts, 'csv'), csvOut = opts.csv; end
pr = drone_params('rigid6dof');
rows = {};
[x0, xT, g] = quad_layouts('circle', 'pointmass');
o = solve_quads_direct('pointmass', x0, xT, g);  r = verify_quads_direct(o);
rows(end+1, :) = {'pointmass', 4, 'direct', o.J, NaN, r.G1clear.minHfine, o.boundsActive, r.G1local.maxMiss, 'n/a', r.ok};
S = load(fullfile(here, 'data', 'quad_direct_M1.mat'));
rows(end+1, :) = {'rigid6dof', 1, 'direct', S.out.J, 4*(pr.m*pr.g/4)^2*S.geom.T, Inf, S.out.boundsActive, S.rep.G1local.maxMiss, 'n/a', S.rep.ok};
oi = run_quad_pmp(struct('plot', false));
rows(end+1, :) = {'rigid6dof', 1, 'indirect', oi.J, 4*(pr.m*pr.g/4)^2*S.geom.T, Inf, false, NaN, oi.conj.verdict, oi.ok};
[x0, xT, g] = quad_layouts('circle', 'rigid6dof');
o4 = solve_quads_direct('rigid6dof', x0, xT, g);  r4 = verify_quads_direct(o4);
rows(end+1, :) = {'rigid6dof', 4, 'direct', o4.J, 4*4*(pr.m*pr.g/4)^2*g.T, r4.G1clear.minHfine, o4.boundsActive, r4.G1local.maxMiss, 'n/a', r4.ok};
tbl = cell2table(rows, 'VariableNames', {'model', 'M', 'route', 'J', 'Jhover', 'minHfine', 'boundsActive', 'localResidual', 'conj', 'ok'});
if ~isfolder(fileparts(csvOut)), mkdir(fileparts(csvOut)); end
writetable(tbl, csvOut);
disp(tbl);
end
```

Run it once; commit the CSV.

- [ ] **Step 4: Register, README, run everything**

Add to `drone/run_tests.m` (`quadrotor` and `quadrotor/tests` on the path; append the seven test names after Route 2's). Write `quadrotor/README.md` in the shape of `route2_indirect/README.md` (files table, the hover-cost floor, the fairness rule for the indirect comparison, the Euler-angle caveat). Update `drone/README.md` Stage 3 row and paste `quad_compare.csv`.

Run: `matlab -batch "cd('optimal_control_examples/drone'); exit(~run_tests())"` and the tree runner.
Expected: `23 of 23 passed`; tree exit 0.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/quadrotor/quad_minenergy_study.m optimal_control_examples/drone/quadrotor/movie_quads.m optimal_control_examples/drone/quadrotor/quad_compare.m optimal_control_examples/drone/quadrotor/results/quad_compare.csv optimal_control_examples/drone/quadrotor/README.md optimal_control_examples/drone/quadrotor/tests/test_quad_study.m optimal_control_examples/drone/run_tests.m optimal_control_examples/drone/README.md
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone quadrotor: study script, movie_quads, quad_compare.csv (pointmass / rigid direct / rigid indirect / rigid M = 4), tests registered (23), READMEs

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```
