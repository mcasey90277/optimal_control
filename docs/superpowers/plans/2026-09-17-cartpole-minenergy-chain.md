# Cart-Pole Sub-Project 1 Implementation Plan — the min-energy chain

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the pieces `ex3_cart_pole_pmp` already has into `cartpole_minenergy_chain.m` — the first of three stage-switched, resumable chain scripts spanning direct → indirect — and in doing so give the direct solve a real front door with two NLP backends on one transcription, analytic derivatives, an explicit seed probe, and a stage-9 export that the min-time and min-fuel chains will warm-start from.

**Architecture:** No new mathematics. The direct NLP that lives inside `gen_direct_ref.m` becomes `minenergy_direct_solve.m` (fmincon and CasADi/IPOPT behind one switch, sharing `trap_cost` and `trap_defects`); the seed block that lives inside `run_cartpole_pmp.m` becomes `cartpole_harvest_seed.m`, which `run_cartpole_pmp` then calls; a seed probe flies the harvested costates before anything shoots; the chain orchestrates one front door per stage and banks a file after each. The study script is untouched — stage 8 runs it headless as the indirect verification instrument.

**Tech Stack:** MATLAB R2026a, `fmincon` (sqp), CasADi 3.7.0 + IPOPT at `~/casadi-3.7.0` (optional backend; absence is a named skip, never a silent pass), `ode113`, `oclib/+oc` (`ms_bvp`, `ms_conjugate_test`, `duals_to_costates`, `fly_control`, `local_residual`).

**Spec:** `docs/superpowers/specs/2026-09-17-cartpole-chains-and-library-design.md` (sub-project 1), which revises `docs/superpowers/specs/2026-09-17-cartpole-three-objectives-design.md`.

## Global Constraints

- Plant constants exactly `m1 = 5`, `m2 = 1`, `L = 2`, `g = 9.8`, read from `cartpole_params()`; `t_f = 5` s fixed; `x(0) = [0;0;0;0]`; `x(t_f) = [0;pi;0;0]`; the fixture's relaxed bound `uMax = 2000` N; `N = 200`.
- House MATLAB style on every function: `%% Purpose: / %% Inputs: / %% Outputs: / %% Revision History:` closed by `%% ------------------------ Begin Code Sequence ---------------------------`; **no `%#ok` pragmas**; never `i`/`j` as loop or index variables; `nargin == 0` self-demo on library-style functions. Author `%  M. Casey`, `(c) 09/18/2026`, `%  Copyright Coorbital Inc.` Scripts carry the `transfer_study.m`/`build_70mN_library.m` header form instead of the quartet.
- `optimal_control_examples/` names no campaign on any executable line: no `orbit_transfer`, no CR3BP quantity, no pumpkyn. `oclib` and CasADi are allowed.
- `cartpole_common/`'s rule: objective-independent only. **Everything this plan creates knows the objective or the solve method, so everything lands in `ex3_cart_pole_pmp/`.** Nothing new goes into `cartpole_common/`.
- Nothing under `ex2_cart_pole_swing_up/` may be modified. `cartpole_minenergy_study.m` is not modified by this plan.
- Run MATLAB only as `/Applications/MATLAB_R2026a.app/bin/matlab -batch "..."`. Do **not** use the MATLAB MCP tool (the shared session in this repo is wedged). `timeout` does not exist on macOS: use `nohup ... > log 2>&1 &` and poll for anything long.
- Committed `.mat` files need `git add -f` (`*.mat` is gitignored). **This plan commits no new `.mat`**: chain outputs go to `results/`, the committed fixtures in `data/` are not regenerated.
- Gates are not widened to make a check pass. Measure, then set the gate with margin above the measured floor of the hardest case (spec L12), and record the measurement in the file.
- Every gate ships with a refusal test that feeds it a known-bad input and asserts a *named* error or a FAIL (spec L17), and every new check gets one deliberate-failure experiment recorded in the task report (spec L-dead). Nine checks on this demo have already turned out to be incapable of failing.
- Every saved artifact carries `meta` (plant, `tf`, `uMax`, `N`, backend, chain name, timestamp) and every reader asserts it before using the data (spec L16).
- A quantity that is the solver's own stopping test is not reported as an independent measurement; say which evaluation produced every residual (spec L14).

---

### Task 1: `minenergy_direct_solve` — the direct front door, `fmincon` backend, analytic derivatives

**Files:**
- Create: `optimal_control_examples/ex3_cart_pole_pmp/minenergy_direct_solve.m`
- Create: `optimal_control_examples/ex3_cart_pole_pmp/trap_cost.m`
- Create: `optimal_control_examples/ex3_cart_pole_pmp/trap_defects.m`
- Create: `optimal_control_examples/ex3_cart_pole_pmp/tests/test_minenergy_direct_solve.m`
- Modify: `optimal_control_examples/ex3_cart_pole_pmp/gen_direct_ref.m` (becomes a thin wrapper)

**Interfaces:**
- Consumes: `cartpole_params()`, `[F, G] = cartpole_field(x, p)`, `A = cartpole_state_jac(x, u, p)` (the 4×4 Jacobian of `F + G*u` with respect to `x`), the committed fixture `data/cartpole_direct_ref.mat` with fields `tN [1×201], X [4×201], U [1×201], muDefect [4×200], J, p, tf, uMax, N, solver`.
- Produces:
  - `out = minenergy_direct_solve(opts)` with `opts` fields `.N` [200], `.uMax` [2000], `.tf` [5], `.p` [`cartpole_params()`], `.backend` ['fmincon' | 'ipopt'], `.Z0` [the example's analytic guess], `.gradients` [true], `.display` ['iter'], `.maxIter` [5000]. Returns `out` with `.tN [1×nN]`, `.X [4×nN]`, `.U [1×nN]`, `.J`, `.muDefect [4×N]`, `.muBC0 [4×1]`, `.muBCf [4×1]`, `.Z`, `.p`, `.tf`, `.uMax`, `.N`, `.backend`, `.solver` (`.exitflag`, `.firstorderopt`, `.constrviolation`, `.iterations`, `.wall`), and `.meta` (see Global Constraints). **Throws by name** (`minenergy_direct_solve:notConverged`, `:infeasible`, `:notStationary`) instead of returning a bad point.
  - `J = trap_cost(u, h)` and `d = trap_defects(X, F, h)` — the ONE transcription. Pure arithmetic on their arguments, so they evaluate on doubles (fmincon) and on CasADi `MX` (Task 2) alike. `trap_defects` takes the already-evaluated field `F [4×nN]`, not the plant, precisely so it has no loop and no preallocation.

- [ ] **Step 1: Write the failing test**

Create `optimal_control_examples/ex3_cart_pole_pmp/tests/test_minenergy_direct_solve.m`:

```matlab
function ok = test_minenergy_direct_solve()
%% Purpose:
%
%   The direct front door must (1) carry a defect Jacobian that agrees with
%   COMPLEX-STEP differentiation of the defect function -- the generated
%   symbolic Jacobian against a numerical method that shares no line with
%   it; (2) reproduce the committed fixture's optimum, which was found by a
%   different code path (its own dynamics copy, finite-difference
%   gradients), so agreement is evidence about the optimum and not about
%   the code; (3) refuse, by name, a solve that did not converge.
%
%   The Jacobian check is the one with teeth: a wrong sign or a wrong
%   block placement in the sparse assembly leaves fmincon converging
%   slowly to the same optimum and no downstream test would ever see it.
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
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok   = true;
here = fileparts(fileparts(mfilename('fullpath')));
root = fileparts(fileparts(here));
addpath(here, fullfile(fileparts(here), 'cartpole_common'), fullfile(root, 'oclib'));
p  = cartpole_params();
R  = load(fullfile(here, 'data', 'cartpole_direct_ref.mat'));

%% (1) the defect Jacobian vs complex step, at a NON-optimal random point
%  (at the optimum many entries are small; a random point exercises every
%  block at full size)
rng(20260918);
N  = 12;  nN = N + 1;  h = R.tf/N;
Z  = [0.3*randn(4*nN, 1); 30*randn(nN, 1)];
[~, ceqA, ~, JcA] = minenergy_direct_solve_nlc(Z, nN, h, p);     % analytic, transposed
JcA = JcA.';                                                       % [nceq x nvar]
nvar = numel(Z);  ncon = numel(ceqA);
cols = unique([1:8, nN-1:nN+1, 2*nN, 3*nN+3, 4*nN+(1:3), 4*nN+nN-1, nvar]);  % every block type
worst = 0;
for c = cols
    Zc = complex(Z);  Zc(c) = Zc(c) + 1e-30i;
    [~, ceqC] = minenergy_direct_solve_nlc(Zc, nN, h, p);
    colCS = imag(ceqC)/1e-30;
    worst = max(worst, max(abs(full(JcA(:, c)) - colCS)) / max(1, max(abs(colCS))));
end
ok = chk(ok, worst < 1e-12, ...
         sprintf('defect Jacobian vs complex step: worst relative column error %.2e (%d columns)', worst, numel(cols)));
ok = chk(ok, isequal(size(JcA), [ncon, nvar]) && issparse(JcA), ...
         sprintf('Jacobian is sparse [%d x %d]', ncon, nvar));

%% (2) reproduces the committed fixture's optimum -- a different code path
out = minenergy_direct_solve(struct('display', 'off'));
ok = chk(ok, abs(out.J - R.J)/R.J < 1e-9, ...
         sprintf('J reproduces the fixture: %.9f vs %.9f (rel %.1e)', out.J, R.J, abs(out.J - R.J)/R.J));
ok = chk(ok, max(abs(out.X(:) - R.X(:))) < 1e-7 && max(abs(out.U - R.U)) < 1e-6, ...
         sprintf('states and control reproduce: max|dX| %.1e, max|dU| %.1e', ...
                 max(abs(out.X(:) - R.X(:))), max(abs(out.U - R.U))));
ok = chk(ok, max(abs(out.muDefect(:) - R.muDefect(:))) / max(abs(R.muDefect(:))) < 1e-6, ...
         sprintf('the defect multipliers reproduce (rel %.1e) -- these seed the costates', ...
                 max(abs(out.muDefect(:) - R.muDefect(:))) / max(abs(R.muDefect(:)))));
ok = chk(ok, out.solver.constrviolation < 1e-9 && out.solver.firstorderopt < 1e-6, ...
         sprintf('fmincon''s own diagnostics: violation %.1e, first-order %.1e', ...
                 out.solver.constrviolation, out.solver.firstorderopt));
ok = chk(ok, isequal(out.meta.p, p) && out.meta.tf == 5 && out.meta.uMax == 2000 && out.meta.N == 200 ...
             && strcmp(out.meta.backend, 'fmincon') && strcmp(out.meta.chain, 'minenergy'), ...
         'meta carries plant, tf, uMax, N, backend, chain');

%% (3) refusal by name: a solve that cannot converge must not return a point
try
    minenergy_direct_solve(struct('display', 'off', 'maxIter', 2));
    ok = chk(ok, false, 'REFUSAL: an unconverged solve returned instead of throwing');
catch err
    ok = chk(ok, startsWith(err.identifier, 'minenergy_direct_solve:'), ...
             sprintf('REFUSAL by name: %s', err.identifier));
end

if ok, fprintf('TEST_MINENERGY_DIRECT_SOLVE: ALL PASS\n');
else,  fprintf('TEST_MINENERGY_DIRECT_SOLVE: FAILURE (see lines above)\n');
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

Note the test calls `minenergy_direct_solve_nlc(Z, nN, h, p)` — the front door's
nonlinear-constraint function, which must therefore be a separate file (fmincon
needs it as a handle anyway, and the test needs to call it on complex `Z`).

- [ ] **Step 2: Run it to verify it fails**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); ok = test_minenergy_direct_solve()"
```

Expected: FAIL — `minenergy_direct_solve_nlc` does not exist.

- [ ] **Step 3: Write the transcription helpers and the front door**

`trap_cost.m`:

```matlab
function J = trap_cost(u, h)
%% Purpose:
%
%   Trapezoidal quadrature of the running cost u^2 on a uniform mesh. Pure
%   arithmetic on its argument: it evaluates on a double row (fmincon) and
%   on a CasADi MX row (IPOPT) alike -- this and trap_defects ARE the one
%   transcription both backends share.
%
%% Inputs:
%
%  u                        [1 x nN]                Control at the nodes
%
%  h                        double                  Mesh step (s)
%
%% Outputs:
%
%  J                        double / MX             h*(sum u^2 - (u1^2 + uN^2)/2)
%
%% Revision History:
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the exact integral of a linear control, 6 - 12 t on [0, 1]:
     tt = linspace(0, 1, 101);  uu = 6 - 12*tt;
     fprintf('trapezoid %.6f vs exact 12\n', trap_cost(uu, tt(2) - tt(1)));
     return
end

J = h*(sum(u.^2) - 0.5*(u(1)^2 + u(end)^2));
end
```

`trap_defects.m`:

```matlab
function d = trap_defects(X, F, h)
%% Purpose:
%
%   The trapezoidal defects, x_{k+1} - x_k - (h/2)(f_{k+1} + f_k), given the
%   FIELD ALREADY EVALUATED at every node. Taking F rather than the plant is
%   deliberate: no loop, no preallocation, so the same line runs on doubles
%   (fmincon) and on CasADi MX (IPOPT). The caller evaluates the field the
%   way its backend needs to.
%
%% Inputs:
%
%  X                        [4 x nN]                States at the nodes
%
%  F                        [4 x nN]                F(x_k) + G(x_k) u_k at
%                                                   the nodes
%
%  h                        double                  Mesh step (s)
%
%% Outputs:
%
%  d                        [4 x nN-1]              One defect column per
%                                                   interval
%
%% Revision History:
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a state that is exactly linear in time has zero defect:
     tt = linspace(0, 1, 6);  Xd = [tt; 2*tt; ones(1,6); 2*ones(1,6)];
     Fd = [ones(1,6); 2*ones(1,6); zeros(2,6)];
     fprintf('max |defect| for a linear state: %.1e\n', max(abs(trap_defects(Xd, Fd, 0.2)), [], 'all'));
     return
end

d = X(:,2:end) - X(:,1:end-1) - (h/2)*(F(:,2:end) + F(:,1:end-1));
end
```

`minenergy_direct_solve_nlc.m` — the constraint function fmincon calls, with the
analytic Jacobian assembled sparsely. Variable layout is state-first:
`Z = [q1(1:nN); q2(1:nN); q1d(1:nN); q2d(1:nN); u(1:nN)]`, so `X(i,k)` is `Z((i-1)*nN + k)`
and `u(k)` is `Z(4*nN + k)`. Constraint rows: `X(:,1) - x0` (rows 1:4),
`X(:,end) - xf` (rows 5:8), then `d(:)` column-major (row `8 + (k-1)*4 + i`).

```matlab
function [c, ceq, gradc, gradceq] = minenergy_direct_solve_nlc(Z, nN, h, p)
%% Purpose:
%
%   Boundary conditions (8 rows) then the trapezoidal defects (4 per
%   interval), in the order the harvest reads the multipliers by -- and,
%   when asked, their Jacobian assembled analytically from the generated
%   state Jacobian A = d(F + G u)/dx and the control column G.
%
%   For defect k = x_{k+1} - x_k - (h/2)(f_{k+1} + f_k):
%       d/dx_k     = -I - (h/2) A_k        d/dx_{k+1} = I - (h/2) A_{k+1}
%       d/du_k     =    - (h/2) G_k        d/du_{k+1} =   - (h/2) G_{k+1}
%   fmincon wants gradceq TRANSPOSED: [nvar x nceq].
%
%   Complex-step safe when gradients are not requested (the test
%   differentiates this function that way), because cartpole_field is.
%
%% Inputs:
%
%  Z                        [5 nN x 1]              Decision vector
%
%  nN                       double                  Number of nodes
%
%  h                        double                  Mesh step (s)
%
%  p                        struct                  Plant (cartpole_params)
%
%% Outputs:
%
%  c                        []                      No inequalities
%
%  ceq                      [8 + 4(nN-1) x 1]       Equality residuals
%
%  gradc                    []
%
%  gradceq                  sparse [5 nN x nceq]    Transposed Jacobian
%
%% Revision History:
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

x0 = [0; 0; 0; 0];  xf = [0; pi; 0; 0];
X = [Z(1:nN).'; Z(nN+1:2*nN).'; Z(2*nN+1:3*nN).'; Z(3*nN+1:4*nN).'];
u = Z(4*nN+1:5*nN).';
F = zeros(4, nN, 'like', Z);                     % 'like' keeps the complex step alive
for k = 1:nN
    [Fk, Gk] = cartpole_field(X(:,k), p);
    F(:,k) = Fk + Gk*u(k);
end
d   = trap_defects(X, F, h);
c   = [];
ceq = [X(:,1) - x0; X(:,end) - xf; d(:)];
gradc = [];
gradceq = [];
if nargout < 4, return; end

%% The Jacobian, as triplets. Column of X(i,k): (i-1)*nN + k; of u(k): 4*nN + k.
N     = nN - 1;
nvar  = 5*nN;
nceq  = 8 + 4*N;
colX  = @(i, k) (i-1)*nN + k;
colU  = @(k) 4*nN + k;
rowD  = @(i, k) 8 + (k-1)*4 + i;
% boundary rows: identity on X(:,1) and X(:,end)
rI = zeros(8,1);  cI = zeros(8,1);  vI = ones(8,1);
for i = 1:4
    rI(i)   = i;      cI(i)   = colX(i, 1);
    rI(4+i) = 4 + i;  cI(4+i) = colX(i, nN);
end
% defect rows: 4 x (4 + 4 + 1 + 1) entries per interval
nPer = 4*(4 + 4 + 1 + 1);
rD = zeros(nPer*N, 1);  cD = zeros(nPer*N, 1);  vD = zeros(nPer*N, 1);
m = 0;
for k = 1:N
    Ak  = cartpole_state_jac(X(:,k),   u(k),   p);
    Ak1 = cartpole_state_jac(X(:,k+1), u(k+1), p);
    [~, Gk]  = cartpole_field(X(:,k),   p);
    [~, Gk1] = cartpole_field(X(:,k+1), p);
    Bk  = -eye(4) - (h/2)*Ak;                    % d/dx_k
    Bk1 =  eye(4) - (h/2)*Ak1;                   % d/dx_{k+1}
    for i = 1:4
        for jx = 1:4
            m = m + 1;  rD(m) = rowD(i,k);  cD(m) = colX(jx, k);    vD(m) = Bk(i, jx);
            m = m + 1;  rD(m) = rowD(i,k);  cD(m) = colX(jx, k+1);  vD(m) = Bk1(i, jx);
        end
        m = m + 1;  rD(m) = rowD(i,k);  cD(m) = colU(k);    vD(m) = -(h/2)*Gk(i);
        m = m + 1;  rD(m) = rowD(i,k);  cD(m) = colU(k+1);  vD(m) = -(h/2)*Gk1(i);
    end
end
Jc = sparse([rI; rD], [cI; cD], [vI; vD], nceq, nvar);
gradceq = Jc.';
end
```

`minenergy_direct_solve.m`:

```matlab
function out = minenergy_direct_solve(opts)
%% Purpose:
%
%   THE direct front door for the minimum-effort swing-up: one trapezoidal
%   collocation NLP, min int u^2 dt at fixed t_f, solved by fmincon (sqp,
%   analytic gradient and sparse constraint Jacobian) or -- the same
%   transcription, evaluated symbolically -- by CasADi + IPOPT. The
%   transcription is the pair trap_cost / trap_defects; the field is
%   cartpole_field; the constraint function is minenergy_direct_solve_nlc.
%
%   This replaces the NLP that lived inside gen_direct_ref. It returns a
%   point only if it is feasible and first-order stationary by the
%   solver's own diagnostics; otherwise it THROWS by name. The committed
%   fixture data/cartpole_direct_ref.mat, found by the old code path (its
%   own dynamics copy, finite-difference gradients), is what
%   tests/test_minenergy_direct_solve holds this function to.
%
%% Inputs:
%
%  opts                     struct (optional)       .N [200], .uMax [2000],
%                                                   .tf [5], .p
%                                                   [cartpole_params()],
%                                                   .backend ['fmincon' |
%                                                   'ipopt'], .Z0 [the
%                                                   example's analytic
%                                                   guess], .gradients
%                                                   [true; fmincon only],
%                                                   .display ['iter'],
%                                                   .maxIter [5000],
%                                                   .casadiDir
%                                                   [~/casadi-3.7.0]
%
%% Outputs:
%
%  out                      struct                  .tN .X [4 x nN] .U
%                                                   [1 x nN] .J .muDefect
%                                                   [4 x N] .muBC0 .muBCf
%                                                   .Z .p .tf .uMax .N
%                                                   .backend .solver
%                                                   (.exitflag
%                                                   .firstorderopt
%                                                   .constrviolation
%                                                   .iterations .wall)
%                                                   .meta
%
%% Revision History:
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'cartpole_common'));
if nargin < 1, opts = struct(); end
d = @(f, v) fieldd(opts, f, v);
N       = d('N', 200);          nN = N + 1;
uMax    = d('uMax', 2000);
tf      = d('tf', 5);
p       = d('p', cartpole_params());
backend = d('backend', 'fmincon');
useGrad = d('gradients', true);
display = d('display', 'iter');
maxIter = d('maxIter', 5000);
casadiDir = d('casadiDir', fullfile(getenv('HOME'), 'casadi-3.7.0'));
tN = linspace(0, tf, nN);
h  = tN(2) - tN(1);

if nargin == 0 && nargout == 0
   %Demo: the fixture's problem, briefly:
     o = minenergy_direct_solve(struct('display', 'final'));
     fprintf('J = %.6f, max|u| = %.2f N, %d iterations\n', o.J, max(abs(o.U)), o.solver.iterations);
     return
end

%% The seed: the example's own analytic guess, unless the caller brought one
if isfield(opts, 'Z0') && ~isempty(opts.Z0)
    Z0 = opts.Z0(:);
    assert(numel(Z0) == 5*nN, 'minenergy_direct_solve:seed', 'Z0 has %d entries, need %d', numel(Z0), 5*nN);
else
    q1  = 0.3*sin(2*pi*tN/tf);
    q2  = pi*tN/tf;
    q1d = (pi/tf)*cos(pi*tN/tf);
    q2d = (pi/tf)*ones(1, nN);
    u0  = 40*(sin(2*pi*tN/tf) + 0.5*sin(4*pi*tN/tf));
    Z0  = [q1, q2, q1d, q2d, u0].';
end

t0 = tic;
switch backend
case 'fmincon'
    obj = @(Z) costWithGrad(Z, nN, h, useGrad);
    nlc = @(Z) minenergy_direct_solve_nlc(Z, nN, h, p);
    lb  = -inf(5*nN, 1);  ub = inf(5*nN, 1);
    lb(4*nN+1:end) = -uMax;  ub(4*nN+1:end) = uMax;
    fo = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', display, ...
        'MaxIterations', maxIter, 'MaxFunctionEvaluations', 1e7, ...
        'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-8, ...
        'SpecifyObjectiveGradient', useGrad, 'SpecifyConstraintGradient', useGrad);
    [Z, J, exitflag, output, lambda] = fmincon(obj, Z0, [], [], [], [], lb, ub, nlc, fo);
    mu = lambda.eqnonlin(:);
    solver = struct('exitflag', exitflag, 'firstorderopt', output.firstorderopt, ...
                    'constrviolation', output.constrviolation, 'iterations', output.iterations);
    % A positive exit flag alone is not evidence of a stationary point (flag
    % 2 is "step too small"): this point seeds a costate solve, so it must
    % be FEASIBLE and first-order optimal or it must not be returned at all.
    assert(exitflag > 0, 'minenergy_direct_solve:notConverged', 'fmincon exit flag %d', exitflag);
    assert(output.constrviolation < 1e-9, 'minenergy_direct_solve:infeasible', ...
           'constraint violation %.2e', output.constrviolation);
    assert(output.firstorderopt < 1e-4, 'minenergy_direct_solve:notStationary', ...
           'first-order optimality %.2e', output.firstorderopt);
case 'ipopt'
    [Z, J, mu, solver] = solveIpopt(Z0, nN, h, p, uMax, display, maxIter, casadiDir);   % Task 2
otherwise
    error('minenergy_direct_solve:backend', 'unknown backend "%s"', backend);
end
wall = toc(t0);
solver.wall = wall;

X = [Z(1:nN).'; Z(nN+1:2*nN).'; Z(2*nN+1:3*nN).'; Z(3*nN+1:4*nN).'];
U = Z(4*nN+1:5*nN).';
% the equality multipliers: 4 boundary rows at each end, then 4 per interval
muBC0    = mu(1:4);
muBCf    = mu(5:8);
muDefect = reshape(mu(9:end), 4, N);

meta = struct('p', p, 'tf', tf, 'uMax', uMax, 'N', N, 'backend', backend, ...
              'chain', 'minenergy', 'stage', 'direct', 'written', char(datetime('now')));
out = struct('tN', tN, 'X', X, 'U', U, 'J', J, 'muDefect', muDefect, 'muBC0', muBC0, ...
             'muBCf', muBCf, 'Z', Z, 'p', p, 'tf', tf, 'uMax', uMax, 'N', N, ...
             'backend', backend, 'solver', solver, 'meta', meta);
end

% ------------------------------------------------------------------------
function [J, g] = costWithGrad(Z, nN, h, useGrad)
%% Purpose:
%
%   trap_cost on the control block of Z, with its gradient: 2 h u_k inside,
%   h u_k at the two ends.
%
u = Z(4*nN+1:5*nN).';
J = trap_cost(u, h);
g = [];
if useGrad
    gu = 2*h*u;  gu(1) = h*u(1);  gu(end) = h*u(end);
    g  = [zeros(4*nN, 1); gu(:)];
end
end

function [Z, J, mu, solver] = solveIpopt(varargin)
%% Purpose:
%
%   Placeholder until Task 2 lands the CasADi backend; refuses by name so
%   the switch cannot silently fall through to fmincon.
%
error('minenergy_direct_solve:noCasadi', 'the ipopt backend arrives in Task 2');
end

function v = fieldd(s, f, v0)
%% Purpose:
%
%   s.(f) if present, else the default.
%
if isfield(s, f), v = s.(f); else, v = v0; end
end
```

`gen_direct_ref.m` becomes a wrapper that calls the front door and saves the
same field set the fixture has always had (`tN X U muDefect J p tf uMax N solver`)
— keep its header, replace its body with:

```matlab
here = fileparts(mfilename('fullpath'));
addpath(here);
out = minenergy_direct_solve(struct('display', 'iter'));
tN = out.tN;  X = out.X;  U = out.U;  muDefect = out.muDefect;  J = out.J;
p = out.p;  tf = out.tf;  uMax = out.uMax;  N = out.N;  solver = out.solver;
if ~isfolder(fullfile(here, 'data')), mkdir(fullfile(here, 'data')); end
save(fullfile(here, 'data', 'cartpole_direct_ref.mat'), ...
     'tN', 'X', 'U', 'muDefect', 'J', 'p', 'tf', 'uMax', 'N', 'solver');
fprintf('gen_direct_ref: J = %.6f, max|u| = %.2f N, %d iterations (see minenergy_direct_solve)\n', ...
        J, max(abs(U)), solver.iterations);
```

and delete its three local functions (`trap_cost`, `defects`, `ref_field`) — the
first two now live as files, the third was a duplicate of the plant. **Do not run
`gen_direct_ref`**; the committed fixture stays as the regression baseline.

- [ ] **Step 4: Run the test**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); exit(~test_minenergy_direct_solve())"
```

Expected: PASS on every check. Record in the report: the Jacobian's worst column
error, `J` to nine decimals against the fixture's `2779.381719`, the multiplier
agreement, the iteration count and wall time **with and without analytic
gradients** (set `opts.gradients` false for the comparison) — the gradient path
should be several times faster and land on the same optimum.

**Deliberate-failure experiment (required, record it):** in a scratch copy of
`minenergy_direct_solve_nlc.m` outside the repo, flip the sign of `Bk1`'s `(h/2)`
term; put the scratch folder first via `cd`; confirm the Jacobian check FAILS with
a column error of order 1; restore and confirm `git status` shows no repo change.

- [ ] **Step 5: Run the existing suite and commit**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); exit(~run_tests())"
```

Expected: 8 of 8 still pass (nothing they read has changed). Then:

```bash
cd /Users/msc/Desktop/optimal_control
git add optimal_control_examples/ex3_cart_pole_pmp/{minenergy_direct_solve.m,minenergy_direct_solve_nlc.m,trap_cost.m,trap_defects.m,gen_direct_ref.m,tests/test_minenergy_direct_solve.m}
git commit -m "cart-pole: a direct front door -- one transcription, analytic sparse Jacobian, refuses an unconverged point by name"
```

---

### Task 2: the CasADi/IPOPT backend on the same transcription, and X3

**Files:**
- Modify: `optimal_control_examples/ex3_cart_pole_pmp/minenergy_direct_solve.m` (replace the `solveIpopt` placeholder)
- Create: `optimal_control_examples/ex3_cart_pole_pmp/tests/test_minenergy_backends.m`

**Interfaces:**
- Consumes: Task 1's `trap_cost`, `trap_defects`, `cartpole_field`; CasADi at `~/casadi-3.7.0` (`import casadi.*`, `MX`, `Function`, `Opti`).
- Produces: `minenergy_direct_solve(struct('backend', 'ipopt'))` returning the same `out` shape, with `out.solver` carrying IPOPT's `.status`, `.iterations`, `.constrviolation` (measured from `ceq` at the returned point by the fmincon-path constraint function, so both backends are judged by the SAME residual), and `.firstorderopt` (IPOPT's `inf_du` at exit, from `sol.stats()`). `out.muDefect` is in **fmincon's sign convention** — see the sign rule below. Also `out.creg`, the labelled constraint-row registry, for the future `foc_check` promotion.
- The **sign rule** (spec L23): the front door applies a fixed constant `signIpopt` (±1) to `lam_g` so that the harvested multipliers match `fmincon`'s convention; the constant's value is *measured* by the test, not derived, and the test asserts the code's constant equals the measurement. Never `opti.dual()`; always `opti.lam_g` indexed by `creg` rows.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_minenergy_backends()
%% Purpose:
%
%   One transcription, two NLP solvers. Checks: (1) SYMBOL SAFETY -- the
%   CasADi Function built from cartpole_field, evaluated at a numeric
%   point, equals the fmincon path's numeric defect vector to round-off
%   (if the plant were not MX-safe the IPOPT path would silently solve a
%   different problem); (2) X3 -- both backends find the same optimum;
%   (3) the multiplier SIGN CONVENTION, measured: the same constraint
%   written x - c == 0 and c - x == 0 flips lam_g, and the constant the
%   front door applies makes IPOPT's harvested multipliers agree with
%   fmincon's -- the trap that cost the orbit work nine days;
%   (4) a missing CasADi is a NAMED refusal, not a fall-through.
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
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok   = true;
here = fileparts(fileparts(mfilename('fullpath')));
root = fileparts(fileparts(here));
addpath(here, fullfile(fileparts(here), 'cartpole_common'), fullfile(root, 'oclib'));
if isempty(which('casadi.MX')), addpath(fullfile(getenv('HOME'), 'casadi-3.7.0')); end
haveCasadi = ~isempty(which('casadi.MX'));
p = cartpole_params();

if ~haveCasadi
    fprintf('  SKIP  CasADi not found at ~/casadi-3.7.0: only the refusal check runs\n');
end

%% (1) symbol safety: the field as a CasADi Function vs the numeric loop
if haveCasadi
    import casadi.*
    xs = MX.sym('x', 4);  us = MX.sym('u');
    [Fs, Gs] = cartpole_field(xs, p);
    fdyn = Function('f', {xs, us}, {Fs + Gs*us});
    rng(20260918);  nN = 9;  h = 0.5;
    Xr = 0.4*randn(4, nN);  Ur = 20*randn(1, nN);
    Fnum = zeros(4, nN);
    for k = 1:nN, [Fk, Gk] = cartpole_field(Xr(:,k), p);  Fnum(:,k) = Fk + Gk*Ur(k); end
    Fm   = fdyn.map(nN);
    Fsym = full(Fm(Xr, Ur));
    dNum = trap_defects(Xr, Fnum, h);
    dSym = trap_defects(Xr, Fsym, h);
    ok = chk(ok, max(abs(dNum(:) - dSym(:))) < 1e-13, ...
             sprintf('symbol safety: defect vector numeric vs CasADi, max diff %.1e', max(abs(dNum(:) - dSym(:)))));
end

%% (2) X3: both backends, one transcription, one optimum
if haveCasadi
    oF = minenergy_direct_solve(struct('display', 'off', 'backend', 'fmincon'));
    oI = minenergy_direct_solve(struct('display', 'off', 'backend', 'ipopt'));
    dJ = abs(oI.J - oF.J)/oF.J;
    dU = sqrt(mean((oI.U - oF.U).^2))/sqrt(mean(oF.U.^2));
    dX = max(abs(oI.X(:) - oF.X(:)));
    ok = chk(ok, dJ < 1e-8, sprintf('X3 cost: fmincon %.9f, ipopt %.9f (rel %.1e)', oF.J, oI.J, dJ));
    ok = chk(ok, dU < 1e-6 && dX < 1e-6, sprintf('X3 control RMS %.1e, states max %.1e', dU, dX));
    ok = chk(ok, oI.solver.constrviolation < 1e-9, ...
             sprintf('ipopt point is feasible by the SAME residual fmincon is judged by: %.1e', oI.solver.constrviolation));
    ok = chk(ok, strcmp(oI.meta.backend, 'ipopt') && isfield(oI, 'creg') && any(strcmp({oI.creg.label}, 'defect')), ...
             'meta says ipopt and the constraint-row registry names the defect block');
    %% (3) the multipliers agree in fmincon's convention -- after the front
    %% door's fixed sign constant, NOT after a runtime sign fit
    dMu = max(abs(oI.muDefect(:) - oF.muDefect(:))) / max(abs(oF.muDefect(:)));
    ok = chk(ok, dMu < 1e-5, sprintf('multipliers agree in fmincon''s convention (rel %.1e)', dMu));
    % the convention itself, measured on a two-variable problem: min x^2 + y^2
    % s.t. x - 1 == 0 written both ways. fmincon: L = f + lam*ceq, so for
    % ceq = x - 1 the multiplier is -2 at the optimum.
    opti = Opti();  xv = opti.variable();  yv = opti.variable();
    opti.minimize(xv^2 + yv^2);  opti.subject_to(xv - 1 == 0);
    opti.solver('ipopt', struct('print_time', false, 'ipopt', struct('print_level', 0)));
    s1 = opti.solve();  lamA = full(s1.value(opti.lam_g));
    opti2 = Opti();  xv2 = opti2.variable();  yv2 = opti2.variable();
    opti2.minimize(xv2^2 + yv2^2);  opti2.subject_to(1 - xv2 == 0);
    opti2.solver('ipopt', struct('print_time', false, 'ipopt', struct('print_level', 0)));
    s2 = opti2.solve();  lamB = full(s2.value(opti2.lam_g));
    ok = chk(ok, abs(lamA + lamB) < 1e-6 && abs(abs(lamA) - 2) < 1e-6, ...
             sprintf('lam_g flips with constraint orientation (%+.4f vs %+.4f) -- the trap, demonstrated', lamA, lamB));
    signMeasured = sign(-2 / lamA);      % what maps this orientation onto fmincon's -2
    ok = chk(ok, signMeasured == minenergy_direct_solve_sign_ipopt(), ...
             sprintf('the front door''s fixed sign constant (%+d) equals the measured one (%+d)', ...
                     minenergy_direct_solve_sign_ipopt(), signMeasured));
end

%% (4) refusal by name when CasADi is absent. The front door's bootstrap
%% adds opts.casadiDir to the path when casadi.MX is not found; pointing it
%% at a directory that does not exist, in a session where CasADi is NOT
%% already on the path, must produce the named refusal and never fall
%% through to fmincon. (When CasADi is already on the path the bootstrap
%% is skipped, so this check is only meaningful after removing it.)
if haveCasadi
    cp = fileparts(fileparts(which('casadi.MX')));
    rmpath(cp);
end
try
    minenergy_direct_solve(struct('display', 'off', 'backend', 'ipopt', 'casadiDir', '/nonexistent/casadi'));
    ok = chk(ok, false, 'REFUSAL: ipopt backend ran with no CasADi reachable');
catch err
    ok = chk(ok, strcmp(err.identifier, 'minenergy_direct_solve:noCasadi'), ...
             sprintf('REFUSAL by name: %s', err.identifier));
end
if haveCasadi, addpath(cp); end

if ok, fprintf('TEST_MINENERGY_BACKENDS: ALL PASS\n');
else,  fprintf('TEST_MINENERGY_BACKENDS: FAILURE (see lines above)\n');
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

Note: the front door gains an `opts.casadiDir` field (default
`fullfile(getenv('HOME'), 'casadi-3.7.0')`) that its bootstrap uses; the refusal
check points it at a directory that does not exist. A refusal check that the
bootstrap could silently defeat would be a check that cannot fail.

- [ ] **Step 2: Run it to verify it fails**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); ok = test_minenergy_backends()"
```

Expected: FAIL — `minenergy_direct_solve_sign_ipopt` does not exist and the ipopt
backend throws the Task-1 placeholder error on the X3 step. First confirm CasADi
is actually present: `ls ~/casadi-3.7.0/+casadi | head -3`.

- [ ] **Step 3: Write the backend**

In `minenergy_direct_solve.m`, replace the `solveIpopt` placeholder with the real
backend, and add the sign constant as its own tiny file so the test can read it:

`minenergy_direct_solve_sign_ipopt.m`:

```matlab
function s = minenergy_direct_solve_sign_ipopt()
%% Purpose:
%
%   The FIXED sign that maps CasADi's opti.lam_g (for a constraint written
%   expr == 0 with expr = residual) onto fmincon's lambda.eqnonlin
%   convention (L = f + lambda' ceq). Measured by
%   tests/test_minenergy_backends on a two-variable problem and asserted
%   there; the front door applies this constant and never fits a sign at
%   run time -- a run-time fit would hide a mis-oriented constraint.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  s                        double                  +1 or -1
%
%% Revision History:
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

s = -1;    % <-- the implementer sets this from the test's measurement and
           %     records the measured value in this comment
end
```

The backend:

```matlab
function [Z, J, mu, solver, creg] = solveIpopt(Z0, nN, h, p, uMax, display, maxIter, casadiDir)
%% Purpose:
%
%   The SAME transcription -- trap_cost, trap_defects, cartpole_field --
%   evaluated on CasADi MX and handed to IPOPT. Constraint rows are
%   registered in creg by label so the multipliers are read back by ROW
%   RANGE from opti.lam_g, never through opti.dual() (whose orientation
%   canonicalisation cost the orbit work nine days), and the fixed sign
%   constant maps them onto fmincon's convention.
%
if isempty(which('casadi.MX')) && isfolder(casadiDir), addpath(casadiDir); end
assert(~isempty(which('casadi.MX')), 'minenergy_direct_solve:noCasadi', ...
       'CasADi not found on the path or at %s', casadiDir);
import casadi.*
N  = nN - 1;
x0 = [0; 0; 0; 0];  xf = [0; pi; 0; 0];

xs = MX.sym('x', 4);  us = MX.sym('u');
[Fs, Gs] = cartpole_field(xs, p);
fdyn = Function('f', {xs, us}, {Fs + Gs*us});

opti = Opti();
X = opti.variable(4, nN);
U = opti.variable(1, nN);
creg = struct('label', {}, 'rows', {});
Fall = fdyn.map(nN);  Fall = Fall(X, U);                        % [4 x nN]

r0 = size(opti.g, 1) + 1;  opti.subject_to(X(:,1) - x0 == 0);
creg(end+1) = struct('label', 'bc0', 'rows', r0:size(opti.g, 1));
r0 = size(opti.g, 1) + 1;  opti.subject_to(X(:,end) - xf == 0);
creg(end+1) = struct('label', 'bcf', 'rows', r0:size(opti.g, 1));
dd = trap_defects(X, Fall, h);
r0 = size(opti.g, 1) + 1;  opti.subject_to(reshape(dd, 4*N, 1) == 0);
creg(end+1) = struct('label', 'defect', 'rows', r0:size(opti.g, 1));
opti.subject_to(-uMax <= U <= uMax);
opti.minimize(trap_cost(U, h));

opti.set_initial(X, [Z0(1:nN).'; Z0(nN+1:2*nN).'; Z0(2*nN+1:3*nN).'; Z0(3*nN+1:4*nN).']);
opti.set_initial(U, Z0(4*nN+1:5*nN).');
po = struct();
po.print_time = false;
po.ipopt.max_iter = maxIter;
po.ipopt.tol = 1e-10;                   % tighter than the orbit family's 1e-7:
po.ipopt.constr_viol_tol = 1e-10;       % X3 compares J to 1e-8 relative
po.ipopt.print_level = tern(strcmp(display, 'off'), 0, 5);
po.ipopt.nlp_scaling_method = 'gradient-based';
opti.solver('ipopt', po);
sol = opti.solve();                     % IPOPT failure throws here, by CasADi

Xv = full(sol.value(X));  Uv = full(sol.value(U));
Z  = [Xv(1,:), Xv(2,:), Xv(3,:), Xv(4,:), Uv].';
J  = full(sol.value(opti.f));
lamAll = full(sol.value(opti.lam_g));
k0 = find(strcmp({creg.label}, 'bc0'), 1);
kf = find(strcmp({creg.label}, 'bcf'), 1);
kd = find(strcmp({creg.label}, 'defect'), 1);
mu = minenergy_direct_solve_sign_ipopt() * [lamAll(creg(k0).rows); lamAll(creg(kf).rows); lamAll(creg(kd).rows)];
% judged by the SAME residual as the fmincon path
[~, ceq] = minenergy_direct_solve_nlc(Z, nN, h, p);
st = sol.stats();
solver = struct('exitflag', double(st.success), 'firstorderopt', st.iterations.inf_du(end), ...
                'constrviolation', max(abs(ceq)), 'iterations', st.iter_count, 'status', st.return_status);
assert(st.success, 'minenergy_direct_solve:notConverged', 'ipopt: %s', st.return_status);
assert(solver.constrviolation < 1e-9, 'minenergy_direct_solve:infeasible', ...
       'constraint violation %.2e', solver.constrviolation);
end
```

Update the `case 'ipopt'` branch to receive `creg` and put it in `out` (`out.creg`;
the fmincon path sets `out.creg = struct('label', {}, 'rows', {})`). Add `tern` as
a local function (`if c, s = a; else, s = b; end`).

Details the implementer will have to settle by running, and must record:
`sol.stats()` field names in CasADi 3.7 (`iter_count`, `return_status`,
`success`, `iterations.inf_du`) — verify against the installed version and
adapt; whether `reshape(dd, 4*N, 1)` or `vec(dd)` is the accepted MX
column-stack; whether `-uMax <= U <= uMax` chains in the MATLAB interface (else
two `subject_to` calls). None of these changes the design.

- [ ] **Step 4: Run the test**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); exit(~test_minenergy_backends())"
```

Expected: PASS. Record the measured X3 numbers (`dJ`, `dU`, `dX`, `dMu`), IPOPT's
iteration count and wall time beside fmincon's, and the measured sign constant.

**If X3 fails at the stated tolerances**, that is a finding, not a tolerance to
widen: report which quantity disagrees and by how much, and whether tightening
IPOPT's `tol` closes it (a solver-tolerance effect) or not (a transcription or
sign defect). Do not commit a widened gate.

**Deliberate-failure experiment (required):** temporarily flip the sign constant
in a scratch copy of `minenergy_direct_solve_sign_ipopt.m` placed first via `cd`;
confirm check (3)'s multiplier agreement FAILS at relative error ~2 and the
constant-equals-measurement check FAILS; restore; confirm `git status` clean.

- [ ] **Step 5: Run the suite and commit**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); exit(~run_tests())"
cd /Users/msc/Desktop/optimal_control
git add optimal_control_examples/ex3_cart_pole_pmp/{minenergy_direct_solve.m,minenergy_direct_solve_sign_ipopt.m,tests/test_minenergy_backends.m}
git commit -m "cart-pole: the ipopt backend on the same transcription; lam_g by row range with a measured, fixed sign; X3 backend agreement"
```

---

### Task 3: `cartpole_harvest_seed`, the seed probe, and `run_cartpole_pmp` delegating

**Files:**
- Create: `optimal_control_examples/ex3_cart_pole_pmp/cartpole_harvest_seed.m`
- Create: `optimal_control_examples/ex3_cart_pole_pmp/cartpole_seed_probe.m`
- Create: `optimal_control_examples/ex3_cart_pole_pmp/tests/test_cartpole_harvest_seed.m`
- Modify: `optimal_control_examples/ex3_cart_pole_pmp/run_cartpole_pmp.m` (lines 77–109 become a call; gains `opts.seed`)

**Interfaces:**
- Consumes: `oc.duals_to_costates(struct('scheme','trapezoid','mu',muDefect,'tNodes',tN,'velRows',3:4))` returning `[lam [4×M], tStations [1×M]]`; `cartpole_field`; `cartpole_pmp_rhs(y, p)`; a direct result `R` with fields `tN, X, U, muDefect, tf` (the fixture, or Task 1's `out`).
- Produces:
  - `[seed, seedDiag] = cartpole_harvest_seed(R, K, p)` — `seed` has `.tf`, `.tGrid [1×K+1]`, `.Y [8×K+1]` (the `ms_bvp` seed shape), `seedDiag` has `.signCorr`, `.flipped`, `.ampRatio`. Throws `cartpole_harvest_seed:sign` exactly where `run_cartpole_pmp` threw `run_cartpole_pmp:sign`.
  - `probe = cartpole_seed_probe(seed, p)` — flies `[x0; seed.Y(5:8,1)]` under `cartpole_pmp_rhs` over `[0, tf]` at RelTol 1e-12, returns `.missTerminal`, `.reachedTf`, `.finite`, `.t`, `.Y`. **This is spec L15**: the measurement that says whether the harvested costates are a usable single-shot seed, made *before* anyone shoots.
  - `run_cartpole_pmp(struct('seed', seed))` uses the supplied seed and skips the harvest; `out.seed` is then the supplied diagnostics if given (`opts.seedDiag`) or a struct of NaNs marked `.supplied = true`.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_cartpole_harvest_seed()
%% Purpose:
%
%   The seed block moved out of run_cartpole_pmp into a front door. Checks:
%   (1) BITWISE agreement with the study script's own verbatim copy of the
%   same block on the same fixture -- the two copies are pinned to each
%   other, so a paraphrase in either shows; (2) run_cartpole_pmp through the
%   front door reproduces its committed regression lam0 exactly; (3) the
%   seed PROBE reports a finite terminal miss for the fixture's seed and
%   its number is recorded -- this is the honest "how good is the seed"
%   figure; (4) refusal by name on a fixture with no sign information.
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
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok   = true;
here = fileparts(fileparts(mfilename('fullpath')));
root = fileparts(fileparts(here));
addpath(here, fullfile(here, 'tests'), fullfile(fileparts(here), 'cartpole_common'), ...
        fullfile(fileparts(here), 'cartpole_common', 'tests'), fullfile(root, 'oclib'));
p = cartpole_params();
R = load(fullfile(here, 'data', 'cartpole_direct_ref.mat'));

%% (1) bitwise against the study script's verbatim copy
[seed, sd] = cartpole_harvest_seed(R, 8, p);
S = run_study_for_seed(fullfile(here, 'cartpole_minenergy_study.m'));
ok = chk(ok, isequal(sd.signCorr, S.seedDiag.signCorr) && isequal(sd.flipped, S.seedDiag.flipped) ...
             && isequal(sd.ampRatio, S.seedDiag.ampRatio), ...
         sprintf('seedDiag BITWISE equal to the study''s copy: corr %.15g, amp %.15g', sd.signCorr, sd.ampRatio));
ok = chk(ok, isequal(seed.Y, S.seed.Y) && isequal(seed.tGrid, S.seed.tGrid) && seed.tf == S.seed.tf, ...
         'seed.Y, .tGrid, .tf BITWISE equal to the study''s');

%% (2) the front door still reproduces its regression baseline
out = run_cartpole_pmp(struct('K', 8, 'plot', false));
ref = load(fullfile(here, 'data', 'cartpole_pmp_ref.mat'));
ok = chk(ok, max(abs(out.lam0 - ref.lam0)) < 1e-8, ...
         sprintf('run_cartpole_pmp via the front door: lam0 regression %.1e', max(abs(out.lam0 - ref.lam0))));
out2 = run_cartpole_pmp(struct('K', 8, 'plot', false, 'seed', seed, 'seedDiag', sd));
ok = chk(ok, isequal(out2.lam0, out.lam0), 'a supplied seed gives the IDENTICAL root (bitwise)');

%% (3) the seed probe: fly the harvested costates single-shot
pr = cartpole_seed_probe(seed, p);
ok = chk(ok, pr.reachedTf && pr.finite, sprintf('probe flight reaches t_f and is finite'));
ok = chk(ok, isfinite(pr.missTerminal) && pr.missTerminal > 0, ...
         sprintf('probe terminal miss recorded: %.3e (the seed''s single-shot quality; NOT a gate)', pr.missTerminal));

%% (4) refusal: a fixture whose control is identically zero carries no sign
Rz = R;  Rz.U = zeros(size(R.U));
try
    cartpole_harvest_seed(Rz, 8, p);
    ok = chk(ok, false, 'REFUSAL: a zero-control fixture was harvested');
catch err
    ok = chk(ok, strcmp(err.identifier, 'cartpole_harvest_seed:sign'), ...
             sprintf('REFUSAL by name: %s', err.identifier));
end

if ok, fprintf('TEST_CARTPOLE_HARVEST_SEED: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_HARVEST_SEED: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function S = run_study_for_seed(scriptPath)
%% Purpose:
%
%   Run the study script in a function workspace and hand back its seed
%   variables -- the verbatim copy this front door is pinned to.
%
run(scriptPath);
S = struct('seedDiag', seedDiag, 'seed', seed);
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
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); ok = test_cartpole_harvest_seed()"
```

Expected: FAIL — `cartpole_harvest_seed` does not exist.

- [ ] **Step 3: Write the front door, the probe, and the delegation**

`cartpole_harvest_seed.m` — the body is `run_cartpole_pmp.m` lines 84–109
**verbatim**, with `R` in place of the loaded struct and the two `assert`
identifiers retitled `cartpole_harvest_seed:sign`. Keep every comment; they record
measurements. The function returns `seed` and `seedDiag` exactly as those lines
build them. Header `Purpose:` states that this is the one home of the harvest
rule for this plant and that the study script keeps a verbatim copy pinned by
`test_cartpole_harvest_seed`.

`cartpole_seed_probe.m`:

```matlab
function probe = cartpole_seed_probe(seed, p)
%% Purpose:
%
%   Fly the harvested costates SINGLE-SHOT from the fixed departure state
%   under the PMP field and measure the terminal miss -- before anyone
%   shoots. The orbit work measured PMP flights from catalog costates
%   missing by 36,000-400,000 km while the recorded control landed within
%   10 km; the lesson was to run this probe first and, if it diverges, to
%   seed MULTIPLE shooting from the whole costate history rather than
%   single shooting from lam(0). This function is that probe for the
%   cart-pole; its number is reported, not gated.
%
%% Inputs:
%
%  seed                     struct                  .tf, .tGrid, .Y [8 x K+1]
%
%  p                        struct                  Plant
%
%% Outputs:
%
%  probe                    struct                  .missTerminal, .reachedTf,
%                                                   .finite, .t, .Y
%
%% Revision History:
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

xf = [0; pi; 0; 0];
y0 = [0; 0; 0; 0; seed.Y(5:8, 1)];
tOut = linspace(0, seed.tf, 501);
[t, Y] = ode113(@(tt, y) cartpole_pmp_rhs(y, p), tOut, y0, odeset('RelTol', 1e-12, 'AbsTol', 1e-14));
reachedTf = numel(t) == numel(tOut) && abs(t(end) - seed.tf) < 1e-12;
finite    = all(isfinite(Y(:)));
if reachedTf && finite, miss = max(abs(Y(end, 1:4).' - xf)); else, miss = Inf; end
probe = struct('missTerminal', miss, 'reachedTf', reachedTf, 'finite', finite, 't', t.', 'Y', Y.');
end
```

`run_cartpole_pmp.m`: replace lines 77–109 with

```matlab
%% Seed: the direct solve's multipliers become costates (cartpole_harvest_seed)
if isfield(opts, 'seed') && ~isempty(opts.seed)
    seed = opts.seed;
    assert(size(seed.Y, 2) == K + 1, 'run_cartpole_pmp:seedK', ...
           'the supplied seed has %d junctions but K = %d', size(seed.Y, 2) - 1, K);
    if isfield(opts, 'seedDiag'), seedDiag = opts.seedDiag;
    else, seedDiag = struct('signCorr', NaN, 'flipped', NaN, 'ampRatio', NaN, 'supplied', true);
    end
else
    [seed, seedDiag] = cartpole_harvest_seed(R, K, p);
end
```

and leave everything below it untouched (`seed.Y(1:4,1)` is already pinned inside
the harvest). Add `.seed` and `.seedDiag` to the header's Inputs block.

- [ ] **Step 4: Run the test and the suite**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); exit(~test_cartpole_harvest_seed())"
cd .. && /Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); exit(~run_tests())"
```

Expected: PASS, and the suite's 8 still pass — in particular `test_cartpole_pmp`'s
`lam0` regression and `test_minenergy_study`'s seedDiag pin. Record the probe's
terminal miss: it is the first honest figure for how far the direct seed's
costates fly on their own, and sub-project 2 quotes it.

**Deliberate-failure experiment (required):** in a scratch copy of
`cartpole_harvest_seed.m` placed first via `cd`, change `'pchip'` to `'linear'`
in the costate interpolation; confirm check (1) FAILS bitwise; restore.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add optimal_control_examples/ex3_cart_pole_pmp/{cartpole_harvest_seed.m,cartpole_seed_probe.m,run_cartpole_pmp.m,tests/test_cartpole_harvest_seed.m}
git commit -m "cart-pole: the harvest is a front door, pinned bitwise to the study's copy; a seed probe flies the costates before anyone shoots"
```

---

### Task 4: `cartpole_minenergy_chain.m` and its live-path test

**Files:**
- Create: `optimal_control_examples/ex3_cart_pole_pmp/cartpole_minenergy_chain.m`
- Create: `optimal_control_examples/ex3_cart_pole_pmp/tests/test_minenergy_chain.m`
- Modify: `.gitignore` (add `optimal_control_examples/*/results/`)

**Interfaces:**
- Consumes: every front door above; `test_cartpole_physics()`; `oc.local_residual(X, tGrid, rhs, opts)` and `oc.fly_control(x0, span, rhs, opts)`; `run_cartpole_pmp(struct('seed', ..., 'K', ..., 'plot', false))`; the study script via `run()` in a function workspace.
- Produces: the chain script (stages 0–9, 3–4 dark), banking under `outDir` (default `results/`): `minenergy_anchor.mat`, `minenergy_verify_direct.mat`, `minenergy_seed.mat`, `minenergy_shoot.mat`, `minenergy_verify_indirect.mat`, `cartpole_minenergy_chain.mat` (the export: `meta, J, lam0, t, X, U, Lam, seedDiag, probe, verdict, direct` summary), and `minenergy_states_control.png`. `chainOverrides` (`.outDir`, `.run.<stage>`, `.solve.backend`, `.solve.K`, `.prob.N`) is honoured with an `assert` on unknown fields. Every stage that is off asserts its input file exists and names the stage that makes it. Blockers are collected; the export stage refuses if any exist.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_minenergy_chain()
%% Purpose:
%
%   A chain script that has only ever run its cached path has never been
%   tested (the orbit work shipped one carrying three defects that way).
%   This runs the LIVE path end to end into a scratch outDir, then the
%   RESUME path against those banked files, then two refusals: a stage
%   whose input is missing, and a banked file whose metadata was tampered
%   with (a fixture at other constants must not be consumed silently).
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
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok   = true;
here = fileparts(fileparts(mfilename('fullpath')));
root = fileparts(fileparts(here));
addpath(here, fullfile(here, 'tests'), fullfile(fileparts(here), 'cartpole_common'), ...
        fullfile(fileparts(here), 'cartpole_common', 'tests'), fullfile(root, 'oclib'));
scratch = fullfile(tempdir, sprintf('minenergy_chain_%d', feature('getpid')));
if isfolder(scratch), rmdir(scratch, 's'); end
mkdir(scratch);

%% (1) LIVE path: every stage on, into scratch
ov = struct('outDir', scratch, 'run', struct('anchor', true, 'verifyDirect', true, 'harvest', true, ...
            'shoot', true, 'verifyIndirect', true, 'export', true, 'pictures', true));
C = run_chain(fullfile(here, 'cartpole_minenergy_chain.m'), ov);
want = {'minenergy_anchor.mat', 'minenergy_verify_direct.mat', 'minenergy_seed.mat', ...
        'minenergy_shoot.mat', 'minenergy_verify_indirect.mat', 'cartpole_minenergy_chain.mat', ...
        'minenergy_states_control.png'};
have = cellfun(@(f) isfile(fullfile(scratch, f)), want);
ok = chk(ok, all(have), sprintf('live path banked every stage file (%d of %d)', nnz(have), numel(want)));
E = load(fullfile(scratch, 'cartpole_minenergy_chain.mat'));
ok = chk(ok, isequal(E.meta.p, cartpole_params()) && E.meta.tf == 5 && strcmp(E.meta.chain, 'minenergy'), ...
         'export carries the plant, t_f and chain name');
ok = chk(ok, E.verdict.claim && E.verdict.necessary && E.verdict.sufficiency, ...
         'the study''s verdict, run headless in stage 8, claims');
ok = chk(ok, isempty(C.blockers), sprintf('no blockers (%d)', numel(C.blockers)));
ref = load(fullfile(here, 'data', 'cartpole_pmp_ref.mat'));
ok = chk(ok, max(abs(E.lam0 - ref.lam0)) < 1e-8, sprintf('export lam0 matches the regression baseline (%.1e)', max(abs(E.lam0 - ref.lam0))));
ok = chk(ok, abs(E.direct.J - E.J)/E.J < 2e-3 && E.direct.trueResidualMax < 1e-6, ...
         sprintf('direct and indirect costs agree (%.1e), true residual %.1e', abs(E.direct.J - E.J)/E.J, E.direct.trueResidualMax));
ok = chk(ok, E.direct.meshSpreadJ < 1e-3 && numel(E.direct.meshN) == 3, ...
         sprintf('mesh-density basin check ran at %d densities, J spread %.1e', numel(E.direct.meshN), E.direct.meshSpreadJ));
ok = chk(ok, isfinite(E.probe.missTerminal), sprintf('seed probe miss carried: %.2e', E.probe.missTerminal));

%% (2) RESUME path: solve stages off, banked files reused, same export
ov2 = ov;  ov2.run.anchor = false;  ov2.run.shoot = false;  ov2.run.harvest = false;
C2 = run_chain(fullfile(here, 'cartpole_minenergy_chain.m'), ov2);
E2 = load(fullfile(scratch, 'cartpole_minenergy_chain.mat'));
ok = chk(ok, isequal(E2.lam0, E.lam0) && isequal(E2.J, E.J) && isempty(C2.blockers), ...
         'resume path reproduces the export bitwise from the banked files');

%% (3) refusal: a stage whose input is missing names the stage that makes it
delete(fullfile(scratch, 'minenergy_seed.mat'));
ov3 = ov;  ov3.run.harvest = false;
try
    run_chain(fullfile(here, 'cartpole_minenergy_chain.m'), ov3);
    ok = chk(ok, false, 'REFUSAL: shoot ran with no banked seed');
catch err
    ok = chk(ok, contains(err.message, 'stage 6'), sprintf('REFUSAL names the maker: "%s"', err.message));
end

%% (4) refusal: tampered metadata in a banked file
A = load(fullfile(scratch, 'minenergy_anchor.mat'));
A.meta.p.m1 = 5.5;
save(fullfile(scratch, 'minenergy_anchor.mat'), '-struct', 'A');
ov4 = ov;  ov4.run.anchor = false;
try
    run_chain(fullfile(here, 'cartpole_minenergy_chain.m'), ov4);
    ok = chk(ok, false, 'REFUSAL: a banked file at other constants was consumed');
catch err
    ok = chk(ok, contains(err.identifier, 'meta'), sprintf('REFUSAL on metadata: %s', err.identifier));
end

rmdir(scratch, 's');
if ok, fprintf('TEST_MINENERGY_CHAIN: ALL PASS\n');
else,  fprintf('TEST_MINENERGY_CHAIN: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function C = run_chain(scriptPath, chainOverrides)
%% Purpose:
%
%   Run the chain SCRIPT in a function workspace with chainOverrides set,
%   and hand back what it leaves behind.
%
run(scriptPath);
C = struct('blockers', {blockers}, 'files', files);
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
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); ok = test_minenergy_chain()"
```

Expected: FAIL — the chain script does not exist.

- [ ] **Step 3: Write the chain**

`cartpole_minenergy_chain.m` (the chain `clearvars -except chainOverrides` and
reads the override struct exactly as `build_70mN_library.m` does):

```matlab
%% CARTPOLE_MINENERGY_CHAIN  The minimum-effort swing-up, direct through indirect
%
%   The companion of cartpole_minenergy_study: that script studies ONE solve
%   with every gate visible; this one runs the CHAIN that produces the
%   solve, with the same discipline -- parameter block first, every stage
%   named, every stage's output a file the next stage reads, a stage that
%   is off is skipped and its banked output reused. Nothing here is new
%   machinery: each stage calls one front door.
%
%     0  parameters and stage switches
%     1  prerequisites        the plant oracle passes; metadata defined
%     2  ANCHOR               the direct solve (minenergy_direct_solve)
%     3  walk                 DARK for min-energy (no homotopy needed)
%     4  refine               DARK for min-energy (no switches)
%     5  VERIFY, direct       first variation two ways, true residual,
%                             mesh-density basin check, bound activity, X3
%     6  HARVEST              duals -> costates; the SEED PROBE
%     7  SHOOT                run_cartpole_pmp on the banked seed
%     8  VERIFY, indirect     the study script, headless; its verdict banked
%     9  EXPORT + pictures    the .mat the min-time and min-fuel chains
%                             warm-start from
%
%   Outputs go to outDir (default results/). Inputs of a stage that is off
%   are read from outDir and must exist; the assert names the stage that
%   makes them. Every banked file carries meta, and every reader asserts it
%   before using the data. Blockers are collected, not thrown; the export
%   stage refuses if any exist.
%
%   A batch driver sets the struct `chainOverrides` (.outDir, .run.<stage>,
%   .solve.backend, .solve.K, .prob.N) before calling this script.
%
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------------------------------------------------------

clearvars -except chainOverrides
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'tests'), fullfile(fileparts(here), 'cartpole_common'), ...
        fullfile(fileparts(here), 'cartpole_common', 'tests'), ...
        fullfile(fileparts(fileparts(here)), 'oclib'));

%% ========================================================================
%  0. PARAMETERS AND STAGE SWITCHES
%% ========================================================================
plant = cartpole_params();
prob  = struct('tf', 5, 'uMax', 2000, 'N', 200, 'x0', [0;0;0;0], 'xf', [0;pi;0;0]);
solve = struct('backend', 'fmincon', 'K', 8, 'meshN', [100 200 400], ...
               'tolTrueResidual', 1e-6, 'tolMeshSpreadJ', 1e-3, 'tolX3J', 1e-8, ...
               'tolStationarity', 1e-6, 'tolAgreeJ', 2e-3);
run = struct('anchor', true, 'walk', false, 'refine', false, 'verifyDirect', true, ...
             'harvest', true, 'shoot', true, 'verifyIndirect', true, 'export', true, 'pictures', true);
outDir = fullfile(here, 'results');
if exist('chainOverrides', 'var')
    if isfield(chainOverrides, 'outDir'), outDir = chainOverrides.outDir; end
    for blk = {'run', 'solve', 'prob'}
        if isfield(chainOverrides, blk{1})
            for f = fieldnames(chainOverrides.(blk{1}))'
                assert(isfield(eval(blk{1}), f{1}), 'unknown %s field "%s"', blk{1}, f{1});
                eval(sprintf('%s.(f{1}) = chainOverrides.%s.(f{1});', blk{1}, blk{1}));
            end
        end
    end
end
assert(~run.walk && ~run.refine, 'stages 3 and 4 are dark for min-energy');
if ~isfolder(outDir), mkdir(outDir); end
files = struct('anchor',  fullfile(outDir, 'minenergy_anchor.mat'), ...
               'vdirect', fullfile(outDir, 'minenergy_verify_direct.mat'), ...
               'seed',    fullfile(outDir, 'minenergy_seed.mat'), ...
               'shoot',   fullfile(outDir, 'minenergy_shoot.mat'), ...
               'vind',    fullfile(outDir, 'minenergy_verify_indirect.mat'), ...
               'export',  fullfile(outDir, 'cartpole_minenergy_chain.mat'), ...
               'picture', fullfile(outDir, 'minenergy_states_control.png'));
blockers = {};

%% ========================================================================
%  1. PREREQUISITES -- the oracle, and the identity every file will carry
%% ========================================================================
assert(test_cartpole_physics(), 'chain:physics', 'the plant fails its oracle; nothing below means anything');
meta = struct('p', plant, 'tf', prob.tf, 'uMax', prob.uMax, 'N', prob.N, 'backend', solve.backend, ...
              'chain', 'minenergy', 'written', char(datetime('now')));
metaOk = @(L, what) assert(isfield(L, 'meta') && isequal(L.meta.p, plant) && L.meta.tf == prob.tf ...
                           && L.meta.uMax == prob.uMax && L.meta.N == prob.N && strcmp(L.meta.chain, 'minenergy'), ...
                           'chain:meta', '%s carries different metadata than this chain''s parameter block', what);
fprintf('1. plant oracle PASS; meta: m1 %g m2 %g L %g g %g, tf %g, uMax %g, N %d, backend %s\n', ...
        plant.m1, plant.m2, plant.L, plant.g, prob.tf, prob.uMax, prob.N, solve.backend);

%% ========================================================================
%  2. ANCHOR -- the direct solve
%% ========================================================================
if run.anchor
    A = minenergy_direct_solve(struct('N', prob.N, 'uMax', prob.uMax, 'tf', prob.tf, 'p', plant, ...
                                      'backend', solve.backend, 'display', 'final'));
    A.meta = meta;  A.meta.stage = 'anchor';
    save(files.anchor, '-struct', 'A');
else
    assert(isfile(files.anchor), 'anchor missing: %s (stage 2 makes it)', files.anchor);
    A = load(files.anchor);
end
metaOk(A, files.anchor);
fprintf('2. anchor: J = %.6f, max|u| = %.2f N, %d iterations (%s)\n', A.J, max(abs(A.U)), A.solver.iterations, A.backend);

%% ========================================================================
%  5. VERIFY, direct
%% ========================================================================
if run.verifyDirect
    nN = prob.N + 1;  h = A.tN(2) - A.tN(1);
    % (a) the first variation, TWO ways. With the solver's own multipliers
    % and OUR analytic gradients: does grad J + Jc' lambda vanish? Then the
    % same with complex-step gradients that share no line with the analytic
    % ones -- so a wrong analytic Jacobian that fmincon happened to tolerate
    % shows here.
    Z = A.Z;  mu = [A.muBC0; A.muBCf; A.muDefect(:)];
    [~, gJ] = chain_cost_grad(Z, nN, h);
    [~, ~, ~, JcT] = minenergy_direct_solve_nlc(Z, nN, h, plant);
    statAnalytic = max(abs(gJ + JcT*mu)) / max(1, max(abs(gJ)));
    gJcs = zeros(size(Z));  JcCS = zeros(numel(mu), numel(Z));
    for c = 1:numel(Z)
        Zc = complex(Z);  Zc(c) = Zc(c) + 1e-30i;
        gJcs(c) = imag(trap_cost(Zc(4*nN+1:5*nN).', h))/1e-30;
        [~, ceqC] = minenergy_direct_solve_nlc(Zc, nN, h, plant);
        JcCS(:, c) = imag(ceqC)/1e-30;
    end
    statCS = max(abs(gJcs + JcCS.'*mu)) / max(1, max(abs(gJcs)));
    % (b) the TRUE continuous residual: re-integrate the returned control
    % through the plant per interval (oc.local_residual), and end to end
    uOf = @(tt) interp1(A.tN, A.U, tt, 'pchip');
    rhs = @(tt, x) chain_fly_rhs(x, uOf(tt), plant);
    res = oc.local_residual(A.X, A.tN, rhs, struct('solver', @ode113, 'RelTol', 1e-12, 'AbsTol', 1e-14));
    trueResidualMax = max(abs(res(:)));
    zEnd = oc.fly_control(prob.x0, [0 prob.tf], rhs, struct('mode', 'span', 'solver', @ode113, 'RelTol', 1e-12, 'AbsTol', 1e-14));
    flownMiss = max(abs(zEnd - prob.xf));
    % (c) mesh density selects the basin (L2): cold solves at three densities
    meshJ = zeros(size(solve.meshN));
    for m = 1:numel(solve.meshN)
        o = minenergy_direct_solve(struct('N', solve.meshN(m), 'uMax', prob.uMax, 'tf', prob.tf, 'p', plant, ...
                                          'backend', solve.backend, 'display', 'off'));
        meshJ(m) = o.J;
    end
    meshSpreadJ = (max(meshJ) - min(meshJ)) / min(meshJ);
    % (d) bound activity
    boundFrac = max(abs(A.U)) / prob.uMax;
    % (e) X3: the other backend on the same transcription, if it is available
    other = tern(strcmp(solve.backend, 'fmincon'), 'ipopt', 'fmincon');
    x3 = struct('ran', false, 'dJ', NaN, 'why', '');
    try
        oX = minenergy_direct_solve(struct('N', prob.N, 'uMax', prob.uMax, 'tf', prob.tf, 'p', plant, ...
                                           'backend', other, 'display', 'off'));
        x3 = struct('ran', true, 'dJ', abs(oX.J - A.J)/A.J, 'why', '');
    catch err
        % ONLY an absent CasADi is a skip. A non-convergence or an infeasible
        % point on the other backend is a finding and must not read as "skipped".
        if strcmp(err.identifier, 'minenergy_direct_solve:noCasadi'), x3.why = err.message;
        else, rethrow(err);
        end
    end
    V = struct('statAnalytic', statAnalytic, 'statCS', statCS, 'trueResidualMax', trueResidualMax, ...
               'flownMiss', flownMiss, 'meshN', solve.meshN, 'meshJ', meshJ, 'meshSpreadJ', meshSpreadJ, ...
               'boundFrac', boundFrac, 'x3', x3, 'J', A.J, 'meta', meta);
    V.meta.stage = 'verifyDirect';
    save(files.vdirect, '-struct', 'V');
else
    assert(isfile(files.vdirect), 'direct verification missing: %s (stage 5 makes it)', files.vdirect);
    V = load(files.vdirect);
end
metaOk(V, files.vdirect);
fprintf('5. first variation: analytic %.2e, complex-step %.2e / %.0e\n', V.statAnalytic, V.statCS, solve.tolStationarity);
fprintf('   true residual (re-integrated, ode113 1e-12) %.2e / %.0e; flown miss %.2e\n', V.trueResidualMax, solve.tolTrueResidual, V.flownMiss);
fprintf('   mesh N = [%s]: J = [%s], spread %.2e / %.0e\n', num2str(V.meshN), num2str(V.meshJ, '%.6f '), V.meshSpreadJ, solve.tolMeshSpreadJ);
fprintf('   max|u| / uMax = %.3f (bound inactive, as the unconstrained comparison needs)\n', V.boundFrac);
if V.x3.ran, fprintf('   X3 backend agreement: dJ %.2e / %.0e\n', V.x3.dJ, solve.tolX3J);
else,        fprintf('   X3 SKIPPED: %s\n', V.x3.why);
end
if V.statCS > solve.tolStationarity, blockers{end+1} = sprintf('first variation (complex step) %.2e', V.statCS); end
if V.trueResidualMax > solve.tolTrueResidual, blockers{end+1} = sprintf('true residual %.2e', V.trueResidualMax); end
if V.meshSpreadJ > solve.tolMeshSpreadJ, blockers{end+1} = sprintf('mesh-density J spread %.2e', V.meshSpreadJ); end
if V.boundFrac > 0.9, blockers{end+1} = 'the relaxed bound is active'; end
if V.x3.ran && V.x3.dJ > solve.tolX3J, blockers{end+1} = sprintf('X3 backend disagreement %.2e', V.x3.dJ); end

%% ========================================================================
%  6. HARVEST -- duals -> costates, then the seed probe
%% ========================================================================
if run.harvest
    [seed, seedDiag] = cartpole_harvest_seed(A, solve.K, plant);
    probe = cartpole_seed_probe(seed, plant);
    Sd = struct('seed', seed, 'seedDiag', seedDiag, 'probe', probe, 'meta', meta);
    Sd.meta.stage = 'harvest';
    save(files.seed, '-struct', 'Sd');
else
    assert(isfile(files.seed), 'seed missing: %s (stage 6 makes it)', files.seed);
    Sd = load(files.seed);
end
metaOk(Sd, files.seed);
fprintf('6. harvest: sign corr %.4f (flipped %d), amplitude ratio %.4f\n', Sd.seedDiag.signCorr, Sd.seedDiag.flipped, Sd.seedDiag.ampRatio);
fprintf('   seed PROBE: single-shot flight of the harvested lam(0) misses by %.3e (reached t_f: %d)\n', ...
        Sd.probe.missTerminal, Sd.probe.reachedTf);

%% ========================================================================
%  7. SHOOT
%% ========================================================================
if run.shoot
    out = run_cartpole_pmp(struct('K', solve.K, 'plot', false, 'seed', Sd.seed, 'seedDiag', Sd.seedDiag));
    Sh = struct('out', out, 'meta', meta);  Sh.meta.stage = 'shoot';
    save(files.shoot, '-struct', 'Sh');
else
    assert(isfile(files.shoot), 'shoot missing: %s (stage 7 makes it)', files.shoot);
    Sh = load(files.shoot);
end
metaOk(Sh, files.shoot);
fprintf('7. shoot: %s, J = %.6f, terminal miss %.2e, flown miss %.2e\n', tern(Sh.out.ok, 'ok', ['NOT OK: ' Sh.out.why]), ...
        Sh.out.J, Sh.out.missTerminal, Sh.out.missFlown);
if ~Sh.out.ok, blockers{end+1} = ['shoot: ' Sh.out.why]; end

%% ========================================================================
%  8. VERIFY, indirect -- the study script, headless
%% ========================================================================
if run.verifyIndirect
    S = chain_run_study(fullfile(here, 'cartpole_minenergy_study.m'));
    close all
    Vi = struct('verdict', S.verdict, 'J', S.J, 'lam0', S.lam0, 'conj', S.conjOut, 'meta', meta);
    Vi.meta.stage = 'verifyIndirect';
    save(files.vind, '-struct', 'Vi');
else
    assert(isfile(files.vind), 'indirect verification missing: %s (stage 8 makes it)', files.vind);
    Vi = load(files.vind);
end
metaOk(Vi, files.vind);
dLam0 = max(abs(Vi.lam0 - Sh.out.lam0));
fprintf('8. study verdict: necessary %d, sufficiency %d, claim %d -- %s\n', Vi.verdict.necessary, ...
        Vi.verdict.sufficiency, Vi.verdict.claim, Vi.verdict.why);
fprintf('   the study''s own solve and stage 7''s agree on lam0 to %.1e\n', dLam0);
if ~Vi.verdict.claim, blockers{end+1} = ['study: ' Vi.verdict.why]; end
if dLam0 > 1e-6, blockers{end+1} = sprintf('study vs shoot lam0 %.1e', dLam0); end

%% ========================================================================
%  9. EXPORT + PICTURES -- refuses with blockers
%% ========================================================================
dJdi = abs(A.J - Sh.out.J)/Sh.out.J;
if dJdi > solve.tolAgreeJ, blockers{end+1} = sprintf('direct vs indirect J %.2e', dJdi); end
if ~isempty(blockers)
    fprintf('9. EXPORT REFUSED, %d blocker(s):\n', numel(blockers));
    fprintf('   - %s\n', blockers{:});
end
if run.export
    assert(isempty(blockers), 'chain:blockers', '%d blocker(s); nothing is exported', numel(blockers));
    direct = struct('J', A.J, 'trueResidualMax', V.trueResidualMax, 'flownMiss', V.flownMiss, ...
                    'statCS', V.statCS, 'meshN', V.meshN, 'meshJ', V.meshJ, 'meshSpreadJ', V.meshSpreadJ, ...
                    'boundFrac', V.boundFrac, 'x3', V.x3, 'backend', A.backend, 'tN', A.tN, 'X', A.X, 'U', A.U);
    Ex = struct('meta', meta, 'J', Sh.out.J, 'lam0', Sh.out.lam0, 't', Sh.out.t, 'X', Sh.out.X, 'U', Sh.out.U, ...
                'Lam', Sh.out.Lam, 'seedDiag', Sd.seedDiag, 'probe', rmfield(Sd.probe, {'t', 'Y'}), ...
                'verdict', Vi.verdict, 'direct', direct);
    Ex.meta.stage = 'export';
    save(files.export, '-struct', 'Ex');
    fprintf('9. exported %s: J = %.6f, lam0 = [%s]\n', files.export, Ex.J, sprintf('%.6f ', Ex.lam0));
end
if run.pictures
    fig = figure('color', [1 1 1], 'Visible', 'off');
    subplot(2,1,1); plot(Sh.out.t, Sh.out.X(2,:), 'k', A.tN, A.X(2,:), 'r--'); grid on
    ylabel('q_2 (rad)'); legend({'indirect (shoot)', 'direct (anchor)'}, 'Location', 'best');
    title(sprintf('cart-pole minimum effort: chain export, J = %.6f', Sh.out.J));
    subplot(2,1,2); plot(Sh.out.t, Sh.out.U, 'k', A.tN, A.U, 'r--'); grid on
    xlabel('t (s)'); ylabel('u (N)');
    exportgraphics(fig, files.picture, 'Resolution', 150);
    close(fig)
    fprintf('9. picture %s\n', files.picture);
end

% ------------------------------------------------------------------------
function [J, g] = chain_cost_grad(Z, nN, h)
%% Purpose:
%
%   trap_cost on Z's control block with its analytic gradient (the same
%   expression the front door hands fmincon).
%
u = Z(4*nN+1:5*nN).';
J = trap_cost(u, h);
gu = 2*h*u;  gu(1) = h*u(1);  gu(end) = h*u(end);
g = [zeros(4*nN, 1); gu(:)];
end

function dx = chain_fly_rhs(x, u, p)
%% Purpose:
%
%   The true dynamics with a supplied control.
%
[F, G] = cartpole_field(x, p);
dx = F + G*u;
end

function S = chain_run_study(scriptPath)
%% Purpose:
%
%   Run the study script in a function workspace; return what stage 8 banks.
%
run(scriptPath);
S = struct('verdict', verdict, 'J', J, 'lam0', lam0, 'conjOut', conjOut);
end

function s = tern(c, a, b)
%% Purpose:
%
%   a if c, else b.
%
if c, s = a; else, s = b; end
end
```

Two things the implementer settles by running: whether `clearvars -except
chainOverrides` inside a script invoked by `run()` from a function workspace
behaves (the test's `run_chain` passes `chainOverrides` as a local variable — if
`clearvars` in that context misbehaves, drop it and rely on the function
workspace's isolation, saying so); and the `eval`-based override merge, which the
implementer may replace with three explicit loops if `eval` on a struct name is
awkward — but the `assert` on an unknown field must stay.

Add to `.gitignore`:

```
optimal_control_examples/*/results/
```

- [ ] **Step 4: Run the test**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); exit(~test_minenergy_chain())"
```

Expected: PASS on all four sections. **This test runs the direct solve four
times (anchor plus three mesh densities) and the study once**; with analytic
gradients it should complete in a few minutes. Time it and record it. If it
exceeds ten minutes, reduce `solve.meshN` to `[100 200 300]` in the test's
override — not in the chain's defaults — and say so.

Record every number the chain prints: both stationarity residuals, the true
residual, the flown miss, the three mesh `J`s and their spread, the bound
fraction, X3's `dJ` (or its skip reason), the probe miss, the study's verdict,
and `dLam0`. **Then set the chain's `solve.tol*` values from those measurements
with margin** (spec L12) — the values in the code above are the plan's estimates
and must be replaced by measured floors; say in the report which moved and to what.

**Deliberate-failure experiments (required, two):** (i) tamper `A.U` in the banked
anchor by 1e-3 and confirm stage 5's true-residual blocker fires and the export is
refused; (ii) set `solve.backend` to `'ipopt'` via override and confirm X3 runs the
other way round and agrees.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add .gitignore optimal_control_examples/ex3_cart_pole_pmp/{cartpole_minenergy_chain.m,tests/test_minenergy_chain.m}
git commit -m "cart-pole: the min-energy chain -- stages 0-9, banked files, meta asserted at every reader, blockers refuse the export"
```

---

### Task 5: wire the tests in, and the docs

**Files:**
- Modify: `optimal_control_examples/ex3_cart_pole_pmp/run_tests.m` (four new names), `README.md` ("The files", "Fixture provenance", a new "Chain script" section), `CLAUDE.md` (one line), `optimal_control_examples/cartpole_common/README.md` (consumers: `cartpole_params` now has `minenergy_direct_solve` and the chain as callers)

**Interfaces:** consumes everything above; produces no code.

- [ ] **Step 1: `run_tests.m`**

Append to `names`, in this order (fast to slow):
`'test_minenergy_direct_solve'`, `'test_cartpole_harvest_seed'`,
`'test_minenergy_backends'`, `'test_minenergy_chain'` — the chain test last, since it
is the slowest. `run_all_tests` needs no change: it calls ex3's `run_tests`.

- [ ] **Step 2: README**

"The files": add `minenergy_direct_solve.m` (+ `_nlc`, `_sign_ipopt`), `trap_cost.m`,
`trap_defects.m`, `cartpole_harvest_seed.m`, `cartpole_seed_probe.m`,
`cartpole_minenergy_chain.m`; note `gen_direct_ref.m` is now a wrapper.

"Fixture provenance": the committed fixture was produced by the OLD
`gen_direct_ref` (its own dynamics copy, finite-difference gradients);
`test_minenergy_direct_solve` holds the new front door to it at the stated
tolerances; the fixture is deliberately **not** regenerated, so it remains a
regression baseline from a different code path.

New "Chain script" section: what a chain is for (the `build_70mN_library.m`
discipline, one paragraph), the stage list, how to run it, how to resume it,
what it exports and who reads it next (the min-time chain's anchor, per the
spec's L5), and the measured numbers from Task 4's report.

- [ ] **Step 3: CLAUDE.md and cartpole_common/README.md**

`CLAUDE.md`: the `ex3_cart_pole_pmp/` line gains "+ `cartpole_minenergy_chain.m`
(stages 0-9, direct through indirect, two NLP backends)". `cartpole_common/README.md`'s
Consumers section: `cartpole_params()` is now called by `minenergy_direct_solve`
and the chain as well as the study.

- [ ] **Step 4: Verify the whole tree and commit**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); exit(~run_all_tests())"
grep -rn --include='*.m' 'orbit_transfer\|pumpkyn\|cr3bp' . | grep -v ':[0-9]*: *%' && echo "CAMPAIGN LEAK -- stop" || echo "campaign-free: ok"
cd /Users/msc/Desktop/optimal_control
git add optimal_control_examples/ex3_cart_pole_pmp/{run_tests.m,README.md} optimal_control_examples/cartpole_common/README.md CLAUDE.md
git commit -m "cart-pole: the min-energy chain's tests wired in; docs for the front doors and the chain"
```

Expected: `run_all_tests` green with ex3 at 12 of 12; the grep finds nothing.

---

## Notes for the executor

- **Order is forced:** 1 → 2 → 3 → 4 → 5. Task 2 needs Task 1's transcription
  helpers; Task 4 needs everything.
- **Task 2 is the one that can be BLOCKED by the environment.** If CasADi's
  MATLAB interface rejects `cartpole_field` on `MX` (a `vertcat` of doubles with
  `MX`, say), report exactly what it rejected; the fallback is to build the
  `Function` from an explicit `MX` transcription of the same four lines inside
  `solveIpopt` — and then the symbol-safety check in the test is what pins that
  copy to the plant. Do not silently take the fallback.
- **The committed fixtures are not regenerated.** `gen_direct_ref` becomes a
  wrapper but is not run; `cartpole_pmp_ref.mat` is not touched. They are the
  regression baselines from the previous code path.
- **The study script is not modified.** Its verbatim seed block is what
  `cartpole_harvest_seed` is pinned to.
- **Every tolerance in the chain's parameter block is a plan estimate until Task
  4 measures it.** Replace estimates with measured floors plus margin; record both.
- **macOS has no `timeout`.** For the chain test use `nohup` and poll if it runs
  long.
