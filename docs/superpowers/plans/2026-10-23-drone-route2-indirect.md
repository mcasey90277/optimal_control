# Drone Route 2 (Indirect PMP-BVP) Implementation Plan (Stage 2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Shoot the barrier form of the paper's multi-agent problem (its Hamiltonian system, eq. 7) with `oc.ms_bvp` at fixed final time, seeded from the Route 1 barrier solve, certify the root with `oc.ms_conjugate_test`, and ship a study script plus a three-way cost comparison (hard direct / barrier direct / indirect / the paper's OCFormer).

**Architecture:** `optimal_control_examples/drone/route2_indirect/`: `drone_pmp_rhs` (the `6M`-state Hamiltonian field, `v* = −λ`, analytic Jacobian from `barrier_terms`), `drone_pmp_prop` (ode113 + variational equations, the cart-pole collapse contract), `run_drones_pmp` (front door: two seeds, engine, gates N1-N6/X1-X3/S1-S2, `.ok/.why`, ε-ladder fallback), `drone_minenergy_study.m` (numbered study script), `drone_route_compare` (the CSV). Reads fixtures from `route1_direct/data/` and the plant from `drone_dynamics/`.

**Tech Stack:** MATLAB R2026a, `oclib` (`oc.ms_bvp`, `oc.duals_to_costates`, `oc.fly_control`, `oc.ms_conjugate_test`), `ode113`. CasADi only indirectly (Route 1's fixtures).

**Spec:** `docs/superpowers/specs/2026-10-23-drone-multiagent-design.md` (Stage 2 section). Stages 0 and 1 are DONE before this plan starts. RE-CHECK before execution: the Route 1 fixture's `lamDefect` sign convention and whether the M = 4 barrier root's conjugate verdict is PASS (if it is UNDETERMINED at K = 8, Task 3's gate reads "not FAIL" as written and the study reports it).

## Global Constraints

- House MATLAB style (header quartet, `nargin == 0` demos, no `i`/`j`, no `norm`, no `%#ok`, no emojis); paths relative to `mfilename('fullpath')`; `oclib` added as `fullfile(root, 'oclib')`; no `orbit_transfer` references.
- The PMP field uses the ANALYTIC Jacobian `A = [0, −I; −HU, 0]` from `barrier_terms`, never a complex step inside the propagator (the field's own derivative is exact and sparse).
- Verbatim gates: N1 `|R| < 1e-10`; N2 terminal miss `< 1e-9` at RelTol 2.5e-14; N3 flown miss `< 1e-6`; N4 H spread `< 1e-9` relative; N5 adjoint `< 1e-9` relative; X1 seed agreement `< 1e-2` relative RMS; X2 `|J_ind − J_dirB|/J_dirB < 1e-3`; X3 `λ(0)` from K = 8 vs K = 16 `< 1e-8`; S2 not FAIL on M = 4.
- Engine options: `fixedTf true`, `tolR 1e-10`, `maxIter 60`, `polishMax 5`, `keepSTMs true`, `K = 8` default.
- ε ladder (fallback only): `[1e-3, 1e-4, 1e-5]`; the rung used is reported in `out.epsRung`.
- Subagents never commit; commit messages end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`; do not stage the unrelated `missiles/` deletions.

## Review Focus

1. A 2D problem shot in 3D: the z rows of `x` and `λ` are identically zero along the flow (`∇U` has no z component when every agent has z = 0), so the `6M` system has `2M` trivially-satisfied equations and the Jacobian's z block is singular for the conjugate test. `run_drones_pmp` must reduce to the active rows (`activeRows = rows with d = 2 excluded`) for the engine and the conjugate test, and pad back on output. (Task 3 test: `freeIdx0` has `2M` entries for a 2D layout, `3M` for sphere.)
2. The seed from the direct multipliers may come back with the wrong global sign: resolve by correlation with `−V`, and refuse (by identifier) a correlation below 0.5 in magnitude. (Task 3 test with a deliberately negated `lamDefect`.)
3. The engine's root must be checked by an INDEPENDENT re-flight at a tighter tolerance than the propagator's, as the cart-pole demo does; the reporting flight must detect truncation. (Task 3.)
4. The conjugate test's UNDETERMINED verdict (an unresolved sample) must not be printed as PASS; the study script and `run_drones_pmp` print the verdict string verbatim. (Tasks 3-4.)
5. The barrier's `C²` kink at `h = l`: if any pair comes within `l` on the arc, `drone_pmp_prop`'s STM is only piecewise accurate; `run_drones_pmp` reports `minH` against `l` (N6) and the ladder note says so. (Task 3.)

---

### Task 1: `drone_pmp_rhs`

**Files:**
- Create: `optimal_control_examples/drone/route2_indirect/drone_pmp_rhs.m`
- Test: `optimal_control_examples/drone/route2_indirect/tests/test_drone_pmp_rhs.m`

**Interfaces:**
- Consumes: `barrier_terms`, `pair_list` (Stage 0).
- Produces: `[dy, v, H, A] = drone_pmp_rhs(y, pairs, geom)` for `y = [x; lam]` [6M x 1]: `dy = [−lam; −gU(x)]`, `v = −lam`, `H = ½|v|² + U + lam'v = U − ½|lam|²`, `A = [0, −I; −HU, 0]` (sparse, only when requested).

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_drone_pmp_rhs()
%% Purpose:
%
%   The Hamiltonian field of the barrier problem:
%     1. v* = -lam, and the costate rate is -dH/dx where H is formed HERE
%        as 1/2|v|^2 + U(x) + lam'v and differentiated by complex step
%        (U from barrier_terms is complex-step safe) -- the adjoint
%        equation checked against the Hamiltonian, not against itself;
%     2. the state rate is dH/dlam by the same complex step;
%     3. the analytic Jacobian A matches a complex step of the field;
%     4. Legendre: d2H/dv2 = I (printed, identically true);
%     5. shapes, and the M = 1 degenerate case (free flight: lam constant).
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
rngWas = rng(41);  restore = onCleanup(@() rng(rngWas));

[x0, ~, geom] = drone_layouts('sphere', 5);
M = 5;  pairs = pair_list(M);
x = x0(:) + 0.03*randn(3*M, 1);
lam = randn(3*M, 1);
y = [x; lam];
[dy, v, H, A] = drone_pmp_rhs(y, pairs, geom);
ok = chk(ok, isequal(size(dy), [6*M 1]) && isequal(size(A), [6*M 6*M]), 'shapes');
ok = chk(ok, isequal(v, -lam), 'v* = -lam');
Hof = @(xx, ll) 0.5*sum((-ll).^2) + barrier_terms(xx, pairs, geom) + ll.'*(-ll);
ok = chk(ok, abs(H - Hof(x, lam)) < 1e-14*max(1, abs(H)), 'H = 1/2|v|^2 + U + lam''v at v*');
hc = 1e-20;  worstX = 0;  worstL = 0;
for k = 1:3*M
    xc = x;  xc(k) = xc(k) + 1i*hc;
    worstX = max(worstX, abs(imag(Hof(xc, lam))/hc + dy(3*M + k))/max(1, abs(dy(3*M + k))));
    lc = lam;  lc(k) = lc(k) + 1i*hc;
    worstL = max(worstL, abs(imag(Hof(x, lc))/hc - dy(k))/max(1, abs(dy(k))));
end
ok = chk(ok, worstX < 1e-12, sprintf('lamdot = -dH/dx by complex step of H (%.1e rel)', worstX));
ok = chk(ok, worstL < 1e-12, sprintf('xdot = dH/dlam by complex step of H (%.1e rel)', worstL));
Acs = zeros(6*M);
for k = 1:6*M
    yc = y;  yc(k) = yc(k) + 1i*hc;
    Acs(:, k) = imag(drone_pmp_rhs(yc, pairs, geom))/hc;
end
ok = chk(ok, max(max(abs(full(A) - Acs))) < 1e-10*max(1, max(abs(Acs(:)))), ...
         sprintf('A vs complex step of the field (%.1e)', max(max(abs(full(A) - Acs)))));
ok = chk(ok, isequal(full(A(1:3*M, 3*M+1:end)), -eye(3*M)) && nnz(A(1:3*M, 1:3*M)) == 0 && nnz(A(3*M+1:end, 3*M+1:end)) == 0, ...
         'A = [0 -I; -HU 0] structure');
fprintf('  INFO  Legendre d2H/dv2 = I: identically, by H''s quadratic form in v\n');
[~, ~, g1] = drone_layouts('circle', 1);
dy1 = drone_pmp_rhs([randn(3,1); [1; 2; 3]], pair_list(1), g1);
ok = chk(ok, isequal(dy1, [-1; -2; -3; 0; 0; 0]), 'M = 1: free flight, lam constant');

if ok, fprintf('TEST_DRONE_PMP_RHS: ALL PASS\n'); else, fprintf('TEST_DRONE_PMP_RHS: FAILURE\n'); end
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

Expected: `Unrecognized function or variable 'drone_pmp_rhs'`.

- [ ] **Step 3: Implement**

```matlab
function [dy, v, H, A] = drone_pmp_rhs(y, pairs, geom)
%% Purpose:
%
%   The Pontryagin field of the barrier-penalised multi-agent problem at
%   fixed final time: y = [x; lam], 6M states.
%
%       H  = 1/2 |v|^2 + U(x) + lam' v
%       dH/dv = v + lam = 0   ->  v* = -lam      (d2H/dv2 = I: a minimum)
%       xdot   = -lam
%       lamdot = -dU/dx
%       A      = [0, -I; -d2U/dx2, 0]            (sparse, from barrier_terms)
%
%   H is constant along any trajectory of this autonomous flow (not zero:
%   T is fixed), so its spread measures integration quality, not
%   optimality (the cart-pole demo's caveat).
%
%  ASSUMPTIONS / NOTES:
%
% • Complex-step safe through barrier_terms (branch on real(h)); a test
%   differentiates H by complex step to check the adjoint equation.
% • A is requested only by the propagator's variational equations.
%
%% Inputs:
%
%  y                        [6M x 1]                [x; lam], agent-major
%
%  pairs                    [2 x P]                 From pair_list
%
%  geom                     struct                  From drone_layouts
%                                                   (.r .eps .l .ball)
%
%% Outputs:
%
%  dy                       [6M x 1]                [xdot; lamdot]
%
%  v                        [3M x 1]                The PMP control -lam
%
%  H                        double                  The Hamiltonian
%
%  A                        [6M x 6M] sparse        d(dy)/dy (nargout > 3)
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the field at the M = 4 circle's start with a unit costate:
     [a, ~, g] = drone_layouts('paper', 4);
     [dy, ~, H] = drone_pmp_rhs([a(:); ones(12, 1)], pair_list(4), g);
     fprintf('H = %.6f, |lamdot| = %.3e\n', H, sqrt(sum(dy(13:end).^2)));
     return
end

n = numel(y)/2;
x = y(1:n);  lam = y(n+1:end);
v = -lam;
if nargout > 3
    [U, gU, HU] = barrier_terms(x, pairs, geom);
    A = [sparse(n, n), -speye(n); -HU, sparse(n, n)];
else
    [U, gU] = barrier_terms(x, pairs, geom);
end
dy = [v; -gU];
H = 0.5*sum(v.^2) + U + lam.'*v;
end
```

- [ ] **Step 4: Run the test**

Expected: all `PASS`.

- [ ] **Step 5: Commit**

```bash
cd ~/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/route2_indirect/drone_pmp_rhs.m optimal_control_examples/drone/route2_indirect/tests/test_drone_pmp_rhs.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone route2: drone_pmp_rhs, the 6M-state Hamiltonian field of the barrier problem (v* = -lam, analytic sparse Jacobian), adjoint checked against a complex step of H

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 2: `drone_pmp_prop`

**Files:**
- Create: `optimal_control_examples/drone/route2_indirect/drone_pmp_prop.m`
- Test: `optimal_control_examples/drone/route2_indirect/tests/test_drone_pmp_prop.m`

**Interfaces:**
- Produces: `[yEnd, PHI] = drone_pmp_prop(dt, y0, needSTM, pairs, geom)`: the `oc.ms_bvp` `prob.prop` contract (ode113 at 1e-12/1e-14, variational equations with the analytic `A` when `needSTM`, `dt == 0` returns `y0` and `I`, THROWS `drone_pmp_prop:collapse` on a truncated grid or non-finite result, `drone_pmp_prop:input` on malformed input). `y0` may be complex when `needSTM` is false (so a test can complex-step through it).

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_drone_pmp_prop()
%% Purpose:
%
%   The propagator contract:
%     1. dt = 0 returns y0 and the identity;
%     2. the STM matches a complex step through the propagator WITHOUT the
%        STM (the ODE solver carries complex state; the field is
%        complex-step safe);
%     3. H is conserved along a flight to integrator tolerance;
%     4. malformed input is refused by identifier;
%     5. a collapse is reported as drone_pmp_prop:collapse (forced by a
%        pair starting INSIDE the quadratic branch with a huge eps, which
%        drives the costates to overflow within the span).
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
rngWas = rng(43);  restore = onCleanup(@() rng(rngWas));

[x0, xT, geom] = drone_layouts('paper', 4);
M = 4;  pairs = pair_list(M);
lam0 = -(xT(:) - x0(:))/10 + 0.02*randn(12, 1);       % roughly the straight-line costate
y0 = [x0(:); lam0];
[yE, PH] = drone_pmp_prop(0, y0, true, pairs, geom);
ok = chk(ok, isequal(yE, y0) && isequal(PH, eye(24)), 'dt = 0: y0 and I');

dt = 1.5;
[yE, PH] = drone_pmp_prop(dt, y0, true, pairs, geom);
ok = chk(ok, isequal(size(yE), [24 1]) && isequal(size(PH), [24 24]) && all(isfinite(PH(:))), 'shapes, finite');
worst = 0;  hc = 1e-20;
for k = 1:24
    yc = y0;  yc(k) = yc(k) + 1i*hc;
    yEc = drone_pmp_prop(dt, yc, false, pairs, geom);
    worst = max(worst, max(abs(imag(yEc)/hc - PH(:, k))));
end
ok = chk(ok, worst < 1e-7, sprintf('STM vs complex step through the propagator (%.1e; integrator-limited)', worst));
[~, ~, H0] = drone_pmp_rhs(y0, pairs, geom);
[~, ~, HE] = drone_pmp_rhs(yE, pairs, geom);
ok = chk(ok, abs(HE - H0) < 1e-10*max(1, abs(H0)), sprintf('H conserved over %.1f (%.1e rel)', dt, abs(HE - H0)/max(1, abs(H0))));

try
    drone_pmp_prop(1, [y0; 1], false, pairs, geom);
    ok = chk(ok, false, 'odd-length y0 refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'drone_pmp_prop:input'), 'odd-length y0 refused by drone_pmp_prop:input');
end
try
    drone_pmp_prop(NaN, y0, false, pairs, geom);
    ok = chk(ok, false, 'NaN dt refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'drone_pmp_prop:input'), 'NaN dt refused by drone_pmp_prop:input');
end

gBad = geom;  gBad.eps = 1e6;  gBad.r = 0.4;           % every pair deep in the quadratic branch
try
    drone_pmp_prop(50, y0, false, pairs, gBad);
    ok = chk(ok, false, 'collapse reported');
catch err
    ok = chk(ok, strcmp(err.identifier, 'drone_pmp_prop:collapse'), sprintf('collapse reported by identifier (%s)', err.identifier));
end

if ok, fprintf('TEST_DRONE_PMP_PROP: ALL PASS\n'); else, fprintf('TEST_DRONE_PMP_PROP: FAILURE\n'); end
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

Expected: `Unrecognized function or variable 'drone_pmp_prop'`.

- [ ] **Step 3: Implement**

```matlab
function [yEnd, PHI] = drone_pmp_prop(dt, y0, needSTM, pairs, geom)
%% Purpose:
%
%   Propagate the drone PMP flow for dt, with the 6M x 6M state-transition
%   matrix on request: the prob.prop contract of oc.ms_bvp. The
%   variational equations ride alongside the state with the ANALYTIC
%   Jacobian A = [0 -I; -HU 0] from drone_pmp_rhs (sparse; no complex step
%   inside the integrator).
%
%  ASSUMPTIONS / NOTES:
%
% • Tolerances 1e-12 / 1e-14: the STM feeds a Newton step.
% • THROWS drone_pmp_prop:collapse when ode113 returns a grid that does not
%   end at dt (it warns and stops early on a blow-up rather than throwing)
%   or a non-finite result; the engine turns a throw into a rejected
%   iterate. Everything else propagates with its own identifier.
% • y0 may be complex when needSTM is false (a caller complex-stepping
%   through the flow); with needSTM it must be real.
%
%% Inputs:
%
%  dt                       double                  Propagation time; 0
%                                                   returns y0 and I
%
%  y0                       [6M x 1]                [x; lam]
%
%  needSTM                  logical                 Return PHI as well
%
%  pairs                    [2 x P]                 From pair_list
%
%  geom                     struct                  From drone_layouts
%
%% Outputs:
%
%  yEnd                     [6M x 1]                y(dt)
%
%  PHI                      [6M x 6M] or []         d y(dt) / d y0
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: one second of the M = 4 flow from a straight-line costate:
     [a, b, g] = drone_layouts('paper', 4);
     [yE, PH] = drone_pmp_prop(1, [a(:); -(b(:) - a(:))/10], true, pair_list(4), g);
     fprintf('|y(1)| = %.4f, cond(PHI) = %.3e\n', sqrt(sum(yE.^2)), cond(PH));
     return
end

n2 = numel(y0);
if ~(isnumeric(y0) && mod(n2, 6) == 0 && n2 > 0 && all(isfinite(y0(:))))
    error('drone_pmp_prop:input', 'y0 must be 6M finite values');
end
if needSTM && ~isreal(y0)
    error('drone_pmp_prop:input', 'y0 must be real when the STM is requested');
end
if ~(isnumeric(dt) && isreal(dt) && isscalar(dt) && isfinite(dt))
    error('drone_pmp_prop:input', 'dt must be a real finite scalar');
end
y0 = y0(:);
if dt == 0
    yEnd = y0;
    PHI = [];  if needSTM, PHI = eye(n2); end
    return
end

oo = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
if needSTM
    z0 = [y0; reshape(eye(n2), n2*n2, 1)];
    [T, Z] = ode113(@(t, z) rhs_with_stm(z, n2, pairs, geom), [0 dt], z0, oo);
    zEnd = Z(end, :).';
    yEnd = zEnd(1:n2);
    PHI  = reshape(zEnd(n2+1:end), n2, n2);
else
    [T, Y] = ode113(@(t, y) drone_pmp_rhs(y, pairs, geom), [0 dt], y0, oo);
    yEnd = Y(end, :).';
    PHI  = [];
end
reachedEnd = (T(end) == dt);
if ~reachedEnd || ~all(isfinite(yEnd)) || (needSTM && ~all(isfinite(PHI(:))))
    error('drone_pmp_prop:collapse', ...
          'the PMP flow left the finite range over dt = %g (integrator reached t = %g)', dt, T(end));
end
end

% ------------------------------------------------------------------------
function dz = rhs_with_stm(z, n2, pairs, geom)
%% Purpose:
%
%   zdot = [f(y); A(y) PHI] with the analytic sparse A.
%
y = z(1:n2);
PHI = reshape(z(n2+1:end), n2, n2);
[dy, ~, ~, A] = drone_pmp_rhs(y, pairs, geom);
dz = [dy; reshape(A*PHI, n2*n2, 1)];
end
```

- [ ] **Step 4: Run the test**

Expected: all `PASS`. If the collapse case does not throw (the flow stays finite), raise `gBad.eps` to 1e9 in the test; the identifier is what is under test.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/route2_indirect/drone_pmp_prop.m optimal_control_examples/drone/route2_indirect/tests/test_drone_pmp_prop.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone route2: drone_pmp_prop (ode113 + analytic variational equations, the engine's collapse contract), STM checked by a complex step through the flow

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 3: `run_drones_pmp`, the front door

**Files:**
- Create: `optimal_control_examples/drone/route2_indirect/run_drones_pmp.m`
- Test: `optimal_control_examples/drone/route2_indirect/tests/test_drones_pmp.m`

**Interfaces:**
- Consumes: Route 1 fixtures (`out`/`outB` structs) or a fresh `solve_drones_direct` barrier solve; `oc.ms_bvp`, `oc.duals_to_costates`, `oc.fly_control`, `oc.ms_conjugate_test`.
- Produces: `out = run_drones_pmp(src, opts)`, `src` an M (loads `route1_direct/data/drones_direct_M<M>.mat` and uses `outB`) or a `solve_drones_direct` barrier out struct; `opts.K [8] .tolR [1e-10] .plot [nargout == 0] .epsLadder [1e-3 1e-4 1e-5] .conj [true]`. `out` fields: `.ok .why .lam0 [3M x 1] .J .JdirectB .JLB .missTerminal .missEngine .missFlown .Hrel .adjRel .minH .seed (.corr .flipped .relRMS) .conj (verdict string, .out) .epsRung .K .activeRows .t .X .Lam .V .info`.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_drones_pmp()
%% Purpose:
%
%   The indirect solve on the M = 4 barrier fixture, and the claims the
%   spec makes for it:
%     N1 converged;  N2 terminal miss < 1e-9;  N3 flown miss < 1e-6;
%     N4 H spread < 1e-9;  N5 adjoint < 1e-9;  N6 min h > l (reported);
%     X1 the two seeds agree (< 1e-2 rel RMS);  X2 J vs the direct barrier
%     cost (< 1e-3);  X3 lam0 from K = 8 and K = 16 (< 1e-8);
%     S2 the conjugate verdict is not FAIL (printed verbatim);
%   plus: a 2D layout uses 2M free costates (z rows excluded), a negated
%   multiplier seed is flipped back, and a seed with no sign information
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
addpath(here, fullfile(fileparts(here), 'drone_dynamics'), fullfile(fileparts(here), 'route1_direct'), fullfile(root, 'oclib'));

o8 = run_drones_pmp(4, struct('K', 8, 'plot', false));
ok = chk(ok, o8.ok, sprintf('front door ok (%s)', tern(isempty(o8.why), 'no reason recorded', o8.why)));
ok = chk(ok, o8.info.converged && o8.info.normR < 1e-10, sprintf('N1 |R| = %.2e', o8.info.normR));
ok = chk(ok, o8.missTerminal < 1e-9, sprintf('N2 terminal miss %.2e', o8.missTerminal));
ok = chk(ok, o8.missFlown < 1e-6, sprintf('N3 flown miss %.2e', o8.missFlown));
ok = chk(ok, o8.Hrel < 1e-9, sprintf('N4 H spread %.2e rel', o8.Hrel));
ok = chk(ok, o8.adjRel < 1e-9, sprintf('N5 adjoint %.2e rel', o8.adjRel));
fprintf('  INFO  N6 min h along the arc %.3e (l = %.0e)\n', o8.minH, 1e-4);
ok = chk(ok, o8.seed.relRMS < 1e-2, sprintf('X1 dual seed vs -V: %.2e rel RMS', o8.seed.relRMS));
ok = chk(ok, abs(o8.J - o8.JdirectB)/o8.JdirectB < 1e-3, sprintf('X2 J = %.6f vs direct barrier %.6f', o8.J, o8.JdirectB));
ok = chk(ok, numel(o8.activeRows) == 2*4*2 && numel(o8.lam0) == 12, '2D: 8 active costate rows (z excluded), lam0 padded to 12');
fprintf('  INFO  S2 conjugate verdict: %s\n', o8.conj.verdict);
ok = chk(ok, ~strcmp(o8.conj.verdict, 'FAIL'), 'S2 not FAIL');
ok = chk(ok, o8.epsRung == 1e-5, sprintf('solved at eps = %.0e without the ladder', o8.epsRung));

o16 = run_drones_pmp(4, struct('K', 16, 'plot', false));
ok = chk(ok, o16.ok && max(abs(o16.lam0 - o8.lam0)) < 1e-8, sprintf('X3 lam0 K = 8 vs 16: %.2e', max(abs(o16.lam0 - o8.lam0))));

S = load(fullfile(fileparts(here), 'route1_direct', 'data', 'drones_direct_M4.mat'));
neg = S.outB;  neg.lamDefect = -neg.lamDefect;
on = run_drones_pmp(neg, struct('K', 8, 'plot', false));
ok = chk(ok, on.ok && on.seed.flipped && max(abs(on.lam0 - o8.lam0)) < 1e-8, 'a negated multiplier seed is flipped back to the same root');
bad = S.outB;  bad.lamDefect = 1e-3*randn(size(bad.lamDefect));
try
    run_drones_pmp(bad, struct('K', 8, 'plot', false));
    ok = chk(ok, false, 'uninformative seed refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'run_drones_pmp:sign'), 'uninformative seed refused by run_drones_pmp:sign');
end

if ok, fprintf('TEST_DRONES_PMP: ALL PASS\n'); else, fprintf('TEST_DRONES_PMP: FAILURE\n'); end
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

Expected: `Unrecognized function or variable 'run_drones_pmp'`.

- [ ] **Step 3: Implement**

```matlab
function out = run_drones_pmp(src, opts)
%% Purpose:
%
%   The indirect solve of the barrier-penalised multi-agent drone problem:
%   form the Pontryagin boundary-value problem (drone_pmp_rhs, 6M states,
%   fixed T) and shoot it with oc.ms_bvp from a seed built out of a Route 1
%   barrier solve. Two seeds are built and compared (gate X1): lam = -V
%   from the direct solution (exact at the optimum, since v* = -lam) and
%   the defect multipliers through oc.duals_to_costates (the covector rule
%   on a new transcription), sign-resolved by correlation with -V. The
%   first seeds the engine.
%
%   2D layouts: the z rows of x and lam are identically zero along the
%   flow, so the engine and the conjugate test see only the ACTIVE rows
%   (out.activeRows); outputs are padded back to 3M.
%
%   Gates (the front door's own verdict, .ok/.why):
%     N1 converged, N2 terminal miss (independent re-flight at RelTol
%     2.5e-14), N3 flown miss (oc.fly_control with v = -lam(t) in closed
%     form), N4 H spread, N5 adjoint vs complex step of H, N6 min h vs l
%     (reported), X1 seed agreement, X2 cost vs the direct barrier cost,
%     S2 conjugate test (verdict reported verbatim; FAIL fails the door).
%
%   eps ladder (fallback only): if the engine fails from the eps = 1e-5
%   seed, walk eps through opts.epsLadder, re-seeding each rung from the
%   previous root; out.epsRung says where it ended.
%
%% Inputs:
%
%  src                      double | struct         M (loads route1_direct/
%                                                   data/drones_direct_M<M>
%                                                   .mat, uses outB) or a
%                                                   solve_drones_direct
%                                                   'barrier' out struct
%
%  opts                     struct (optional)       .K [8] .tolR [1e-10]
%                                                   .plot [nargout == 0]
%                                                   .epsLadder [1e-3 1e-4
%                                                   1e-5] .conj [true]
%
%% Outputs:
%
%  out                      struct                  .ok .why .lam0 .J
%                                                   .JdirectB .JLB
%                                                   .missTerminal
%                                                   .missEngine .missFlown
%                                                   .Hrel .adjRel .minH
%                                                   .seed .conj .epsRung
%                                                   .K .activeRows .t .X
%                                                   .Lam .V .info
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(here)));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'), fullfile(fileparts(here), 'route1_direct'), fullfile(root, 'oclib'));
d = @(f, v) fieldd(opts, f, v);
K      = d('K', 8);
tolR   = d('tolR', 1e-10);
doPlot = d('plot', nargout == 0);
ladder = d('epsLadder', [1e-3 1e-4 1e-5]);
doConj = d('conj', true);

%% The direct barrier solution this seeds from
if isnumeric(src)
    S = load(fullfile(fileparts(here), 'route1_direct', 'data', sprintf('drones_direct_M%d.mat', src)));
    R = S.outB;
else
    R = src;
end
if ~strcmp(R.mode, 'barrier')
    error('run_drones_pmp:mode', 'the seed must be a ''barrier'' direct solve (the hard problem has no smooth PMP field)');
end
geom = R.geom;  M = geom.M;  T = R.T;  pairs = R.pairs;  tN = R.tN;  N = R.N;
x0 = R.X(:, 1);  xT = R.X(:, end);
n = 3*M;
if geom.d == 2, act = setdiff(1:n, 3:3:n); else, act = 1:n; end
na = numel(act);
rowsY = [act, n + act];                                  % active rows of y = [x; lam]

%% Seed 1: lam = -V. Seed 2: the defect multipliers through the covector rule
lamV = -R.V;                                                          % [3M x N+1]
[lamS, tS] = oc.duals_to_costates(struct('scheme', 'trapezoid', 'mu', R.lamDefect, ...
                                         'tNodes', tN, 'velRows', 1:3));
tS = tS(:).';
lamVs = interp1(tN, lamV.', tS, 'pchip').';
num = sum(sum(lamS(act, :).*lamVs(act, :)));
den = sqrt(sum(sum(lamS(act, :).^2))*sum(sum(lamVs(act, :).^2)));
if ~(den > 0) || abs(num/den) < 0.5
    error('run_drones_pmp:sign', 'the multiplier seed carries no usable sign information (corr %.3f)', num/max(den, eps));
end
flipped = num < 0;
if flipped, lamS = -lamS; end
relRMS = sqrt(sum(sum((lamS(act, :) - lamVs(act, :)).^2))/sum(sum(lamVs(act, :).^2)));
seedDiag = struct('corr', abs(num/den), 'flipped', flipped, 'relRMS', relRMS);

%% The reduced problem for the engine
tGrid = linspace(0, T, K + 1);
Xg = interp1(tN, R.X.', tGrid, 'pchip').';
Lg = interp1(tN, lamV.', tGrid, 'pchip').';
Yfull = [Xg; Lg];  Yfull(1:n, 1) = x0;
seed = struct('tf', T, 'tGrid', tGrid, 'Y', Yfull(rowsY, :));
xTa = xT(act);
epsRung = geom.eps;
[info, gUsed] = shoot(seed, geom, pairs, act, n, na, xTa, K, tolR);
if ~info.converged
    % the ladder: walk eps down, re-seeding from each rung's root
    seedL = seed;
    for e = ladder
        gL = geom;  gL.eps = e;
        [infoL, gUsed] = shoot(seedL, gL, pairs, act, n, na, xTa, K, tolR);
        if infoL.converged
            seedL.Y = infoL.Y;  info = infoL;  epsRung = e;
        else
            break
        end
    end
end

%% Report: fly the answer and measure what the gates measure
lam0 = zeros(n, 1);  lam0(act) = info.Y(na + (1:na), 1);
y0 = [x0; lam0];
[t, Y] = ode113(@(tt, y) drone_pmp_rhs(y, pairs, gUsed), linspace(0, T, 2001), y0, ...
                odeset('RelTol', 2.5e-14, 'AbsTol', 1e-16));
reachedTf = (numel(t) == 2001) && (t(end) == T) && all(isfinite(Y(:)));
X = Y(:, 1:n).';  Lam = Y(:, n+1:end).';  V = -Lam;
H = zeros(1, numel(t));  adjRel = 0;  hc = 1e-20;
for k = 1:numel(t)
    [dy, ~, H(k)] = drone_pmp_rhs(Y(k, :).', pairs, gUsed);
    if mod(k, 200) == 1                               % N5 at 11 samples: lamdot vs -dH/dx by complex step
        gH = zeros(n, 1);
        for q = act
            yc = Y(k, :).';  yc(q) = yc(q) + 1i*hc;
            [~, ~, Hc] = drone_pmp_rhs(yc, pairs, gUsed);
            gH(q) = imag(Hc)/hc;
        end
        adjRel = max(adjRel, max(abs(dy(n + act) + gH(act)))/max(1, max(abs(gH(act)))));
    end
end
J = trapz(t, 0.5*sum(V.^2, 1));
Hrel = (max(H) - min(H))/max(abs(H));
hArc = clearance_terms(X, pairs, gUsed);
if isempty(hArc), minH = Inf; else, minH = min(hArc(:)); end
ppL = pchip(t, Lam);
zEnd = oc.fly_control(x0, [0 T], @(tt, x) -ppval(ppL, min(max(tt, 0), T)), ...
                      struct('mode', 'span', 'solver', @ode113, 'RelTol', 1e-13, 'AbsTol', 1e-15));
missTerminal = max(abs(X(act, end) - xT(act)));
missFlown    = max(abs(zEnd(act) - xT(act)));
dtLast = info.tGrid(end) - info.tGrid(end-1);
yLast  = drone_pmp_prop(dtLast, padY(info.Y(:, end), act, n), false, pairs, gUsed);
missEngine = max(abs(yLast(act) - xT(act)));

%% S2: the conjugate test on the reduced flow
conj = struct('verdict', 'NOT RUN', 'out', []);
if doConj && info.converged
    it = info;
    co = oc.ms_conjugate_test(it, struct('stateRows', 1:na, 'costateCols', na + (1:na), ...
                                        'quotientDir', [], 'freeTime', false));
    conj.out = co;  conj.verdict = co.verdict;
end

%% Verdict
why = '';
if ~info.converged,             why = sprintf('engine not converged (|R| = %.2e)', info.normR);
elseif ~reachedTf,              why = 'reporting flight truncated or non-finite';
elseif ~(missTerminal < 1e-9),  why = sprintf('terminal miss %.2e', missTerminal);
elseif ~(missFlown < 1e-6),     why = sprintf('flown miss %.2e', missFlown);
elseif ~(Hrel < 1e-9),          why = sprintf('H spread %.2e', Hrel);
elseif ~(adjRel < 1e-9),        why = sprintf('adjoint residual %.2e', adjRel);
elseif ~(abs(J - R.J)/R.J < 1e-3), why = sprintf('J = %.6f vs direct barrier %.6f', J, R.J);
elseif strcmp(conj.verdict, 'FAIL'), why = 'conjugate point: FAIL';
end
out = struct('ok', isempty(why), 'why', why, 'lam0', lam0, 'J', J, 'JdirectB', R.J, 'JLB', R.JLB, ...
             'missTerminal', missTerminal, 'missEngine', missEngine, 'missFlown', missFlown, ...
             'Hrel', Hrel, 'adjRel', adjRel, 'minH', minH, 'seed', seedDiag, 'conj', conj, ...
             'epsRung', epsRung, 'K', K, 'activeRows', act, 't', t.', 'X', X, 'Lam', Lam, 'V', V, 'info', info);
if ~out.ok, warning('run_drones_pmp:notOk', 'NOT a solution: %s', why); end
if doPlot, plotRoute2(out, R); end
if nargout == 0
    fprintf('drone PMP-BVP M = %d (%s): J = %.6f (direct barrier %.6f, hard-ref J_LB %.4f), terminal miss %.2e, conj %s\n', ...
            M, tern(out.ok, 'ok', ['NOT OK: ' why]), J, R.J, R.JLB, missTerminal, conj.verdict);
    clear out
end
end

% ------------------------------------------------------------------------
function [info, g] = shoot(seed, g, pairs, act, n, na, xTa, K, tolR)
%% Purpose:
%
%   One oc.ms_bvp call on the reduced (active-row) problem.
%
prob = struct('ny', 2*na, 'freeIdx0', na + (1:na), ...
    'prop', @(dt, y0, needSTM) reducedProp(dt, y0, needSTM, act, n, pairs, g), ...
    'rhs',  @(y) reducedRhs(y, act, n, pairs, g), ...
    'terminal', @(y, needJ) terminalFcn(y, xTa, na, needJ));
[~, info] = oc.ms_bvp(prob, seed, struct('fixedTf', true, 'tolR', tolR, 'maxIter', 60, ...
                                         'polishMax', 5, 'keepSTMs', true));
if numel(seed.tGrid) ~= K + 1
    warning('run_drones_pmp:K', 'seed grid has %d segments, not K = %d', numel(seed.tGrid) - 1, K);
end
end

function [yEnd, PHI] = reducedProp(dt, y0, needSTM, act, n, pairs, g)
%% Purpose:
%
%   Pad the active rows to the full 6M state, propagate, reduce back
%   (the STM's inactive rows and columns are dropped; they decouple
%   exactly because the z rows carry zero).
%
yf = padY(y0, act, n);
[yEf, PHIf] = drone_pmp_prop(dt, yf, needSTM, pairs, g);
rowsY = [act, n + act];
yEnd = yEf(rowsY);
PHI = [];  if needSTM, PHI = PHIf(rowsY, rowsY); end
end

function dy = reducedRhs(y, act, n, pairs, g)
%% Purpose:
%
%   The field on the active rows.
%
dyf = drone_pmp_rhs(padY(y, act, n), pairs, g);
dy = dyf([act, n + act]);
end

function yf = padY(y, act, n)
%% Purpose:
%
%   Active-row vector -> full 6M vector with zeros in the inactive rows.
%
na = numel(act);
yf = zeros(2*n, 1);
yf(act) = y(1:na);
yf(n + act) = y(na + (1:na));
end

function [gq, dgdy] = terminalFcn(y, xTa, na, needJ)
%% Purpose:
%
%   x(T) = xT on the active rows, and its Jacobian.
%
gq = y(1:na) - xTa;
dgdy = [];
if needJ, dgdy = [eye(na), zeros(na)]; end
end

function plotRoute2(out, R)
%% Purpose:
%
%   Indirect vs direct: paths and the control magnitude of agent 1.
%
figure('color', [1 1 1]);
subplot(1, 2, 1);  hold on
M = R.geom.M;
for a = 1:M
    ia = 3*(a-1) + (1:3);
    plot(R.X(ia(1), :), R.X(ia(2), :), '--', 'Color', [0.7 0.7 0.7]);
    plot(out.X(ia(1), :), out.X(ia(2), :), '-');
end
axis equal;  grid on;  title('paths: indirect (solid) vs direct barrier (dashed)');
subplot(1, 2, 2);
plot(R.tN, sqrt(sum(R.V(1:3, :).^2, 1)), '--', out.t, sqrt(sum(out.V(1:3, :).^2, 1)), '-');
grid on;  xlabel('t');  ylabel('|v_1|');  title('agent 1 speed');
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

- [ ] **Step 4: Run the test**

Run: `matlab -batch "cd('optimal_control_examples/drone/route2_indirect/tests'); test_drones_pmp"`
Expected: all `PASS` lines, with the two `INFO` lines printing N6 and the conjugate verdict verbatim. If N1 fails from the ε = 1e-5 seed and the ladder rescues it, the `epsRung == 1e-5` line fails: record the rung in the task report and open the question in the ledger rather than loosening the test.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/route2_indirect/run_drones_pmp.m optimal_control_examples/drone/route2_indirect/tests/test_drones_pmp.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone route2: run_drones_pmp front door (two seeds compared, active-row reduction for 2D, oc.ms_bvp fixed-tf, N1-N6/X1-X3/S2 gates with its own verdict, eps-ladder fallback), tested on the M = 4 fixture

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 4: `drone_minenergy_study.m` and its test

**Files:**
- Create: `optimal_control_examples/drone/route2_indirect/drone_minenergy_study.m`
- Test: `optimal_control_examples/drone/route2_indirect/tests/test_minenergy_study.m`

**Interfaces:**
- Produces: a script in the cart-pole study's shape (sections 0-9, numbered; every gate computed inline with its own PASS/FAIL; section 8 asserts the inline numbers against `run_drones_pmp`'s). Parameter block at the top: `M = 4`, `K = 8`, `selfCheck.strict = false`. Leaves `verdict` (struct of logicals per gate) and `study` (numbers) in the workspace.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_minenergy_study()
%% Purpose:
%
%   The study script runs end to end in a clean workspace, leaves a
%   verdict struct with every necessary gate true, the conjugate verdict
%   printed, and its inline numbers equal to run_drones_pmp's (section 8's
%   own assertion, re-checked here).
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
ws = runStudy(fullfile(here, 'drone_minenergy_study.m'));
ok = chk(ok, isfield(ws, 'verdict') && isfield(ws, 'study'), 'script leaves verdict and study');
if ~ok, fprintf('TEST_MINENERGY_STUDY: FAILURE\n'); return, end
gates = {'N1', 'N2', 'N3', 'N4', 'N5', 'X1', 'X2', 'X3'};
for g = gates
    ok = chk(ok, isfield(ws.verdict, g{1}) && ws.verdict.(g{1}), sprintf('verdict.%s true', g{1}));
end
ok = chk(ok, ischar(ws.study.conjVerdict) && ~strcmp(ws.study.conjVerdict, 'FAIL'), sprintf('S2 %s', ws.study.conjVerdict));
ok = chk(ok, ws.verdict.libraryAgrees, 'section 8: inline numbers agree with run_drones_pmp');

if ok, fprintf('TEST_MINENERGY_STUDY: ALL PASS\n'); else, fprintf('TEST_MINENERGY_STUDY: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ws = runStudy(file)
%% Purpose:
%
%   Run the script in this function's workspace and return it as a struct.
%
run(file);
vars = who;
ws = struct();
for k = 1:numel(vars)
    if ~strcmp(vars{k}, 'file'), ws.(vars{k}) = eval(vars{k}); end
end
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

Expected: `run` errors because the script file does not exist.

- [ ] **Step 3: Write the study script**

```matlab
%% DRONE_MINENERGY_STUDY  Multi-agent minimum-energy drones, indirectly, step by step.
%
%  Edit the parameter block, press Run. The steps are done here rather
%  than behind the front door (run_drones_pmp) so the machinery is
%  visible; the Jacobi determinant sweep stays in the library
%  (oc.ms_conjugate_test) because a second copy would be unverified, and
%  section 8 asserts the inline numbers against the front door.
%
%  Sections:
%    0  tolerances          one struct, feeding every printed line AND verdict
%    1  the plant           kinematic drones, the layout, J_LB
%    2  the problem         the barrier form, H, v* = -lam, A
%    3  the seed            lam = -V from the direct barrier solve; the
%                           multipliers through the covector rule (X1)
%    4  solve               oc.ms_bvp at fixed T on the active rows, STMs kept
%    5  independent check   a perturbed seed finds the same root (X3-like)
%    6  NECESSARY  (N1-N6)  converged, terminal, flown, H, adjoint, clearance
%    7  SUFFICIENCY (S1-S2) Legendre (identically), conjugate point
%    8  assert vs library   the inline numbers against run_drones_pmp
%    9  plot                paths vs direct, speeds, H(t), min h(t)
%
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------------------------------------------------------

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(here)));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'), fullfile(fileparts(here), 'route1_direct'), fullfile(root, 'oclib'));

%% Parameters
M = 4;  K = 8;
selfCheck = struct('strict', false);

%% 0. Tolerances
tol = struct('R', 1e-10, 'miss', 1e-9, 'flown', 1e-6, 'H', 1e-9, 'adj', 1e-9, 'seed', 1e-2, 'J', 1e-3, 'lam0', 1e-8);

%% 1. The plant and the layout
S = load(fullfile(fileparts(here), 'route1_direct', 'data', sprintf('drones_direct_M%d.mat', M)));
R = S.outB;  geom = R.geom;  T = R.T;  pairs = R.pairs;  n = 3*M;
if geom.d == 2, act = setdiff(1:n, 3:3:n); else, act = 1:n; end
na = numel(act);
JLB = lower_bound(reshape(R.X(:,1), 3, M), reshape(R.X(:,end), 3, M), T);
fprintf('1. M = %d kinematic drones, r = %.3f, T = %g, J_LB = %.4f, direct hard J = %.6f, direct barrier J = %.6f\n', ...
        M, geom.r, T, JLB, S.out.J, R.J);

%% 2. The problem: H = 1/2|v|^2 + U(x) + lam'v, v* = -lam, A = [0 -I; -HU 0]
y1 = [R.X(:, 1); -R.V(:, 1)];
[~, v1, H1, A1] = drone_pmp_rhs(y1, pairs, geom);
fprintf('2. at t = 0 with the direct seed: H = %.6f, |v*| = %.4f, nnz(A) = %d\n', H1, sqrt(sum(v1.^2)), nnz(A1));

%% 3. The seed, two ways
lamV = -R.V;
[lamS, tS] = oc.duals_to_costates(struct('scheme', 'trapezoid', 'mu', R.lamDefect, 'tNodes', R.tN, 'velRows', 1:3));
tS = tS(:).';
lamVs = interp1(R.tN, lamV.', tS, 'pchip').';
c = sum(sum(lamS(act,:).*lamVs(act,:)))/sqrt(sum(sum(lamS(act,:).^2))*sum(sum(lamVs(act,:).^2)));
if c < 0, lamS = -lamS; end
seedRelRMS = sqrt(sum(sum((lamS(act,:) - lamVs(act,:)).^2))/sum(sum(lamVs(act,:).^2)));
verdict.X1 = seedRelRMS < tol.seed;
fprintf('3. X1 dual seed vs -V: corr %.3f, rel RMS %.2e  [%s]\n', abs(c), seedRelRMS, pf(verdict.X1));

%% 4. Solve on the active rows
tGrid = linspace(0, T, K + 1);
Yfull = [interp1(R.tN, R.X.', tGrid, 'pchip').'; interp1(R.tN, lamV.', tGrid, 'pchip').'];
Yfull(1:n, 1) = R.X(:, 1);
rowsY = [act, n + act];
seed = struct('tf', T, 'tGrid', tGrid, 'Y', Yfull(rowsY, :));
pad = @(y) padRows(y, act, n);
prob = struct('ny', 2*na, 'freeIdx0', na + (1:na), ...
    'prop', @(dt, y0, needSTM) reduceProp(dt, y0, needSTM, act, n, pairs, geom), ...
    'rhs',  @(y) reduceVec(drone_pmp_rhs(pad(y), pairs, geom), act, n), ...
    'terminal', @(y, needJ) terminalRed(y, needJ, R.X(act, end), na));
[~, info] = oc.ms_bvp(prob, seed, struct('fixedTf', true, 'tolR', tol.R, 'maxIter', 60, 'polishMax', 5, 'keepSTMs', true));
lam0 = zeros(n, 1);  lam0(act) = info.Y(na + (1:na), 1);
verdict.N1 = info.converged && info.normR < tol.R;
fprintf('4. N1 engine: converged %d, |R| = %.2e, %d iterations  [%s]\n', info.converged, info.normR, info.iters, pf(verdict.N1));

%% 5. Independent check: a perturbed seed must find the same root
seedP = seed;  seedP.Y(na+1:end, :) = seedP.Y(na+1:end, :)*(1 + 0.05);
[~, infoP] = oc.ms_bvp(prob, seedP, struct('fixedTf', true, 'tolR', tol.R, 'maxIter', 60, 'polishMax', 5));
lam0P = zeros(n, 1);  lam0P(act) = infoP.Y(na + (1:na), 1);
verdict.X3 = infoP.converged && max(abs(lam0P - lam0)) < tol.lam0;
fprintf('5. X3 perturbed seed (+5%% costates): same lam0 to %.2e  [%s]\n', max(abs(lam0P - lam0)), pf(verdict.X3));

%% 6. NECESSARY: one flight, every number from it
[t, Y] = ode113(@(tt, y) drone_pmp_rhs(y, pairs, geom), linspace(0, T, 2001), [R.X(:,1); lam0], ...
                odeset('RelTol', 2.5e-14, 'AbsTol', 1e-16));
X = Y(:, 1:n).';  Lam = Y(:, n+1:end).';  V = -Lam;
H = zeros(1, numel(t));
for k = 1:numel(t), [~, ~, H(k)] = drone_pmp_rhs(Y(k,:).', pairs, geom); end
missT = max(abs(X(act, end) - R.X(act, end)));
ppL = pchip(t, Lam);
zEnd = oc.fly_control(R.X(:,1), [0 T], @(tt, x) -ppval(ppL, min(max(tt, 0), T)), struct('mode', 'span', 'solver', @ode113, 'RelTol', 1e-13, 'AbsTol', 1e-15));
missF = max(abs(zEnd(act) - R.X(act, end)));
Hrel = (max(H) - min(H))/max(abs(H));
adjRel = 0;  hc = 1e-20;
for k = 1:200:2001
    dy = drone_pmp_rhs(Y(k,:).', pairs, geom);
    for q = act
        yc = Y(k,:).';  yc(q) = yc(q) + 1i*hc;
        [~, ~, Hc] = drone_pmp_rhs(yc, pairs, geom);
        adjRel = max(adjRel, abs(dy(n + q) + imag(Hc)/hc)/max(1, abs(imag(Hc)/hc)));
    end
end
hArc = clearance_terms(X, pairs, geom);  minH = min(hArc(:));
J = trapz(t, 0.5*sum(V.^2, 1));
verdict.N2 = missT < tol.miss;  verdict.N3 = missF < tol.flown;  verdict.N4 = Hrel < tol.H;  verdict.N5 = adjRel < tol.adj;
verdict.X2 = abs(J - R.J)/R.J < tol.J;
fprintf('6. N2 terminal miss %.2e [%s]; N3 flown miss %.2e [%s]; N4 H spread %.2e [%s]; N5 adjoint %.2e [%s]\n', ...
        missT, pf(verdict.N2), missF, pf(verdict.N3), Hrel, pf(verdict.N4), adjRel, pf(verdict.N5));
fprintf('   N6 min h on the arc %.3e vs l = %.0e (reported); X2 J = %.6f vs direct barrier %.6f [%s]; hard reference %.6f\n', ...
        minH, geom.l, J, R.J, pf(verdict.X2), S.out.J);

%% 7. SUFFICIENCY
fprintf('7. S1 Legendre: d2H/dv2 = I identically (H is quadratic in v)  [PASS by form]\n');
co = oc.ms_conjugate_test(info, struct('stateRows', 1:na, 'costateCols', na + (1:na), 'quotientDir', [], 'freeTime', false));
fprintf('   S2 Jacobi: %s (first full-rank sample %d of %d)\n', co.verdict, co.firstFullRank, numel(co.t));
if selfCheck.strict && strcmp(co.verdict, 'FAIL'), error('study:conjugate', 'conjugate point refutes the root'); end

%% 8. Assert against the library
lib = run_drones_pmp(M, struct('K', K, 'plot', false));
verdict.libraryAgrees = lib.ok && max(abs(lib.lam0 - lam0)) < tol.lam0 && abs(lib.J - J) < 1e-9 && strcmp(lib.conj.verdict, co.verdict);
fprintf('8. run_drones_pmp: ok %d, lam0 agrees to %.2e, J agrees to %.2e, conj %s  [%s]\n', ...
        lib.ok, max(abs(lib.lam0 - lam0)), abs(lib.J - J), lib.conj.verdict, pf(verdict.libraryAgrees));

%% 9. Plot
figure('color', [1 1 1], 'Position', [100 100 1400 400]);
subplot(1, 3, 1);  hold on
for a = 1:M, ia = 3*(a-1) + (1:3);  plot(R.X(ia(1),:), R.X(ia(2),:), '--', 'Color', [0.7 0.7 0.7]);  plot(X(ia(1),:), X(ia(2),:), '-');  end
axis equal;  grid on;  title('paths: indirect (solid), direct barrier (dashed)');
subplot(1, 3, 2);  plot(t, H - H(1));  grid on;  xlabel('t');  title('H(t) - H(0)');
subplot(1, 3, 3);  plot(t, min(hArc, [], 1));  grid on;  xlabel('t');  ylabel('min h');  title('clearance along the arc');
study = struct('J', J, 'JdirectB', R.J, 'Jhard', S.out.J, 'JLB', JLB, 'lam0', lam0, 'missT', missT, 'missF', missF, ...
               'Hrel', Hrel, 'adjRel', adjRel, 'minH', minH, 'seedRelRMS', seedRelRMS, 'conjVerdict', co.verdict);

%% ------------------------------------------------------------------------
function s = pf(c)
%% Purpose:
%
%   PASS/FAIL text.
%
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function yf = padRows(y, act, n)
%% Purpose:
%
%   Active rows -> full 6M vector.
%
na = numel(act);  yf = zeros(2*n, 1);  yf(act) = y(1:na);  yf(n + act) = y(na + (1:na));
end

function v = reduceVec(vf, act, n)
%% Purpose:
%
%   Full 6M vector -> active rows.
%
v = vf([act, n + act]);
end

function [yEnd, PHI] = reduceProp(dt, y0, needSTM, act, n, pairs, geom)
%% Purpose:
%
%   Propagate on the full state, reduce back.
%
[yEf, PHIf] = drone_pmp_prop(dt, padRows(y0, act, n), needSTM, pairs, geom);
rowsY = [act, n + act];
yEnd = yEf(rowsY);  PHI = [];  if needSTM, PHI = PHIf(rowsY, rowsY); end
end

function [g, dgdy] = terminalRed(y, needJ, xTa, na)
%% Purpose:
%
%   x(T) = xT on the active rows, with its Jacobian (both outputs always
%   exist: the engine may ask for one or two).
%
g = y(1:na) - xTa;
dgdy = [];  if needJ, dgdy = [eye(na), zeros(na)]; end
end
```

- [ ] **Step 4: Run the test**

Expected: all `PASS`; the S2 line prints the verdict verbatim.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/route2_indirect/drone_minenergy_study.m optimal_control_examples/drone/route2_indirect/tests/test_minenergy_study.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone route2: drone_minenergy_study (numbered study script, every gate inline, asserted against the front door) with its test

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 5: `drone_route_compare`, runner registration, README

**Files:**
- Create: `optimal_control_examples/drone/route2_indirect/drone_route_compare.m`
- Create (generated, tracked): `optimal_control_examples/drone/route2_indirect/results/route_compare.csv`
- Create: `optimal_control_examples/drone/route2_indirect/README.md`
- Modify: `optimal_control_examples/drone/run_tests.m`, `optimal_control_examples/drone/README.md`

- [ ] **Step 1: Write the comparison writer**

```matlab
function tbl = drone_route_compare(opts)
%% Purpose:
%
%   The three-way cost comparison per M: J_LB, the hard direct optimum
%   (reference), the barrier direct cost, the indirect cost shot from it,
%   the paper's OCFormer and SympOCNet, and the conjugate verdict, written
%   to results/route_compare.csv. M = 4 and 8 read the fixtures; M = 16
%   solves both direct modes first (opts.Ms controls the list).
%
%% Inputs:
%
%  opts                     struct (optional)       .Ms [4 8 16] .K [8]
%                                                   .csv [results/
%                                                   route_compare.csv]
%
%% Outputs:
%
%  tbl                      table                   One row per M
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'), fullfile(fileparts(here), 'route1_direct'));
Ms = [4 8 16];  if isfield(opts, 'Ms'), Ms = opts.Ms; end
K = 8;          if isfield(opts, 'K'), K = opts.K; end
csvOut = fullfile(here, 'results', 'route_compare.csv');  if isfield(opts, 'csv'), csvOut = opts.csv; end
paper = containers.Map({4, 8, 16, 32}, {[0.2003 0.2001], [0.4036 0.4046], [0.8307 0.8345], [1.2531 1.4413]});
rows = {};
for M = Ms
    fx = fullfile(fileparts(here), 'route1_direct', 'data', sprintf('drones_direct_M%d.mat', M));
    if isfile(fx)
        S = load(fx);  oH = S.out;  oB = S.outB;
    else
        [x0, xT, geom] = drone_layouts('paper', M);
        oH = solve_drones_direct(x0, xT, geom, struct('mode', 'hard'));
        oB = solve_drones_direct(x0, xT, geom, struct('mode', 'barrier'));
    end
    oI = run_drones_pmp(oB, struct('K', K, 'plot', false));
    pp = paper(M);
    rows(end+1, :) = {M, oH.JLB, oH.J, oB.J, oI.J, pp(1), pp(2), oI.conj.verdict, oI.ok, oI.epsRung, oI.minH};
end
tbl = cell2table(rows, 'VariableNames', {'M', 'J_LB', 'J_hard', 'J_barrier_direct', 'J_indirect', 'OCFormer', 'SympOCNet', 'conj', 'indirectOk', 'epsRung', 'minH'});
if ~isfolder(fileparts(csvOut)), mkdir(fileparts(csvOut)); end
writetable(tbl, csvOut);
disp(tbl);
end
```

Run it once (`drone_route_compare` from the folder; M = 16 solves two direct NLPs and a 64-state shooting problem, allow several minutes). Commit the CSV.

- [ ] **Step 2: Register the tests and write the README**

Append to `drone/run_tests.m` `names` (and add `route2_indirect` + its `tests` to the `addpath`):

```matlab
         'test_drone_pmp_rhs'         % the Hamiltonian field vs a complex step of H
         'test_drone_pmp_prop'        % propagator, STM, collapse contract
         'test_drones_pmp'            % the solve, end to end, M = 4
         'test_minenergy_study'       % the study script (slowest)
```

`route2_indirect/README.md`:

```markdown
# route2_indirect

The barrier form of the paper's multi-agent problem (Zhang, Liu, Darbon &
Karniadakis, arXiv:2610.30388, eq. 7) solved INDIRECTLY: the Pontryagin
boundary-value problem shot with `oc.ms_bvp` at fixed T, seeded from
Route 1's barrier solve, certified with `oc.ms_conjugate_test`.

| file | what |
|---|---|
| `drone_pmp_rhs.m` | the 6M-state Hamiltonian field, `v* = -lam`, analytic sparse Jacobian |
| `drone_pmp_prop.m` | ode113 + variational equations, the engine's collapse contract |
| `run_drones_pmp.m` | front door: two seeds (lam = -V and the covector rule) compared, active-row reduction for 2D, gates N1-N6 / X1-X3 / S2, `.ok/.why`, eps-ladder fallback |
| `drone_minenergy_study.m` | the numbered study script: every gate inline, asserted against the front door |
| `drone_route_compare.m` | `results/route_compare.csv`: J_LB, hard, barrier-direct, indirect, OCFormer, SympOCNet, conjugate verdict per M |

The barrier is summed over pairs (the paper takes a max, which is only
piecewise smooth); the indirect cost therefore matches the barrier DIRECT
cost (X2), while the hard optimum is the reference. The conjugate verdict
is printed verbatim; UNDETERMINED is not PASS.
```

Update `drone/README.md`: Stage 2 row "built", paste the `route_compare.csv` columns `M, J_LB, J_hard, J_indirect, OCFormer, conj`.

- [ ] **Step 3: Run the folder runner and the tree runner**

Expected: `16 of 16 passed`; tree exit 0.

- [ ] **Step 4: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/route2_indirect/drone_route_compare.m optimal_control_examples/drone/route2_indirect/results/route_compare.csv optimal_control_examples/drone/route2_indirect/README.md optimal_control_examples/drone/run_tests.m optimal_control_examples/drone/README.md
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone route2: route_compare.csv (J_LB / hard / barrier-direct / indirect / OCFormer / conjugate verdict per M), tests registered (16), READMEs

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```
