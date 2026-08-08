# Booster Landing (3-DOF PDG + TVLQR) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Solve the Falcon-9-class min-fuel landing burn two independent ways (Hermite-Simpson NLP and lossless-convexification), certify agreement + PMP structure, then close the loop with TVLQR and Monte-Carlo the landing.

**Architecture:** Campaign-style MATLAB: one front door (`run_booster_landing.m`) over `lib/` (params, dynamics, two guidance solvers, TVLQR, closed-loop sim, Monte Carlo), `certify/` (gates G1–G5), `viz/`, `tests/`, `doc/`. Phase 1 is vacuum/flat-Earth (convexification exactly valid); Phase 2 adds an opt-in drag branch to the collocation/tracking side only.

**Tech Stack:** MATLAB R2025b, CasADi 3.7 + IPOPT (at `~/casadi-3.7.0`), ode45, LaTeX note.

**Spec:** `docs/superpowers/specs/2026-08-08-booster-landing-design.md` — read it first.

## Global Constraints

- MATLAB is **R2025b only**: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "..."` (R2025a license is broken). Every run command below uses this binary.
- Standard run pattern (write it out fully in each step):
  `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; <command>"`
- CasADi: `addpath(fullfile(getenv('HOME'),'casadi-3.7.0'))` — done once in `setup_paths.m`, never per-file.
- **Pumpkyn house style** (invoke the `matlab-pumpkyn-style` skill when writing each lib file): full comment header (Purpose / INPUTS with sizes / OUTPUTS with sizes / REFERENCES), `%%` section comments, column-aligned `=` in parameter blocks, `if nargin==0` self-demo on library functions where it makes sense.
- **Never use `i` or `j` as loop variables** (use `k`, `kk`, `idx`, `krun`).
- **No `norm()`/`abs()` inside `pdg_dynamics`** — complex-step differentiation is a test gate; use `sqrt(sum(x.^2))` forms.
- Front-door contract: `run_booster_landing` with **no arguments** runs the complete nominal campaign from a fresh MATLAB (self-pathing, prints progress, writes `results/`). All knobs via an optional config struct, never by editing files.
- Determinism: every random draw behind `rng(seed)` with the seed in the params struct.
- Commit after every task, message prefix `booster_landing:`, footer `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`. Stage only booster_landing/docs files (the repo has unrelated dirty files — never `git add -A`).
- Results (`results/`) are generated products: git-ignored except `.gitkeep`.

---

### Task 1: Skeleton, `setup_paths`, `booster_params` + sanity test

**Files:**
- Create: `booster_landing/setup_paths.m`
- Create: `booster_landing/lib/booster_params.m`
- Create: `booster_landing/tests/test_params.m`
- Create: `booster_landing/results/.gitkeep`, `booster_landing/.gitignore`

**Interfaces:**
- Consumes: nothing.
- Produces: `P = booster_params()` — single struct every later task reads:
  `P.mdry P.m0 P.Tmax P.Tmin P.Isp P.g0 P.gvec(3x1) P.r0(3x1) P.v0(3x1) P.gs_deg P.theta_max_deg P.N P.tf_lo P.tf_hi P.pad_radius P.vtd_max P.seed P.drag.{on,rho0,H,Cd,A}` — and `setup_paths()` (adds campaign dirs + CasADi).

- [ ] **Step 1: Write the failing test**

`tests/test_params.m` (house pattern: self-bootstrapping, throws on failure):

```matlab
% TEST_PARAMS  Physical sanity of the booster parameter set.
%
% Checks the three facts the whole campaign leans on:
%   1. m0 > mdry (there is landing propellant);
%   2. min-throttle thrust-to-weight at DRY mass > 1 (hoverslam is forced:
%      the booster cannot hover, so the optimizer must produce a terminal
%      max-thrust arc arriving at v=0 exactly at touchdown);
%   3. bounds ordered: 0 < Tmin < Tmax, tf_lo < tf_hi, gs in (0,90) deg.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P = booster_params();
assert(P.m0 > P.mdry, 'no landing propellant');
twMin = P.Tmin / (P.mdry * P.g0);
assert(twMin > 1, 'min-throttle T/W at dry mass = %.3f, expected > 1', twMin);
assert(P.Tmin > 0 && P.Tmin < P.Tmax, 'thrust bounds disordered');
assert(P.tf_lo < P.tf_hi, 'tf bracket disordered');
assert(P.gs_deg > 0 && P.gs_deg < 90, 'glideslope out of range');
assert(P.drag.on == false, 'Phase 1 default must be vacuum');
fprintf('test_params PASS (min-throttle dry T/W = %.3f)\n', twMin);
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); addpath('tests'); test_params"`
Expected: FAIL — `setup_paths` / `booster_params` undefined.

- [ ] **Step 3: Implement `setup_paths.m` and `booster_params.m`**

`setup_paths.m`:

```matlab
function setup_paths()
% SETUP_PATHS  Paths for the booster_landing campaign (one campaign, one setup).
%
% Adds: this folder (front door), lib (dynamics + solvers + tracking),
% certify (gates), viz (plots + movie), tests, and CasADi 3.7.
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, 'lib'));
addpath(fullfile(here, 'certify'));
addpath(fullfile(here, 'viz'));
addpath(fullfile(here, 'tests'));
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
end
```

`lib/booster_params.m` — values verbatim from the spec (pumpkyn-style aligned block). `if nargin==0` self-demo prints the struct and the T/W numbers:

```matlab
function P = booster_params()
% BOOSTER_PARAMS  Falcon-9-class 3-DOF landing-burn parameter set.
%
% Single source of truth for physical constants, boundary conditions,
% solver grid and Monte-Carlo settings. Public F9 estimates; every number
% adjustable here and only here.
%
% INPUTS:  none
% OUTPUTS: P - parameter struct (fields documented inline below)
%
% REFERENCES:
%   [1] Blackmore, Acikmese, Scharf, "Minimum-Landing-Error Powered-Descent
%       Guidance for Mars Landing Using Convex Optimization," JGCD 2010.
%   [2] docs/superpowers/specs/2026-08-08-booster-landing-design.md

%% Vehicle (public Falcon 9 block-5 estimates):
P.mdry   = 25600;              % dry mass [kg]
P.m0     = 30000;              % mass at landing-burn start [kg]
P.Tmax   = 845e3;              % one Merlin 1D, sea level [N]
P.Tmin   = 0.40 * P.Tmax;      % ~40 percent min throttle [N]
P.Isp    = 282;                % sea-level Isp [s]
P.g0     = 9.80665;            % standard gravity [m/s^2]
P.gvec   = [0; 0; -P.g0];      % flat-Earth gravity, z up [m/s^2]

%% Boundary conditions (pad at origin, z up):
P.r0     = [500; 100; 2000];   % post-entry-burn position [m]
P.v0     = [-30; 0; -180];     % descending ~180 m/s [m/s]

%% Path constraints:
P.gs_deg        = 30;          % glideslope min elevation angle [deg]
P.theta_max_deg = Inf;         % thrust-pointing cone half-angle, Inf = off

%% Discretization / solver:
P.N      = 60;                 % Hermite-Simpson segments (collocation)
P.Nconv  = 120;                % trapezoid nodes for the convex solver
P.tf_lo  = 10;                 % free-final-time bracket [s]
P.tf_hi  = 50;

%% Atmosphere (Phase 2, OFF by default -- vacuum keeps convexification exact):
P.drag.on   = false;
P.drag.rho0 = 1.225;           % sea-level density [kg/m^3]
P.drag.H    = 8500;            % scale height [m]
P.drag.Cd   = 1.0;             % landing-leg config drag coefficient [-]
P.drag.A    = 10.75;           % reference area, 3.7 m diameter [m^2]

%% Success criteria / Monte Carlo:
P.pad_radius = 15;             % landing accuracy requirement [m]
P.vtd_max    = 2.0;            % max touchdown speed [m/s]
P.seed       = 42;             % rng seed, all random draws
end
```

- [ ] **Step 4: Run test to verify it passes**

Same command as Step 2. Expected: `test_params PASS (min-throttle dry T/W = 1.346)`.

- [ ] **Step 5: `.gitignore` + commit**

`booster_landing/.gitignore`:

```
results/*
!results/.gitkeep
*.asv
```

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/setup_paths.m booster_landing/lib/booster_params.m \
        booster_landing/tests/test_params.m booster_landing/.gitignore booster_landing/results/.gitkeep
git commit -m "booster_landing: skeleton, params, sanity test

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: `pdg_dynamics` with analytic Jacobians + complex-step test

**Files:**
- Create: `lib/pdg_dynamics.m`
- Create: `tests/test_dynamics_jac.m`

**Interfaces:**
- Consumes: `P` from `booster_params`.
- Produces: `[xdot, A, B] = pdg_dynamics(x, T, P)` — `x = [r; v; m]` (7x1), `T` thrust vector [N] (3x1), `xdot` 7x1, `A = ∂f/∂x` 7x7, `B = ∂f/∂T` 7x3. Honors `P.drag.on` (Phase 2 branch written NOW, tested now, default off — the opt-in-flag house pattern). Must be complex-step safe: no `norm`, no `abs`, no `max`.

- [ ] **Step 1: Write the failing test**

`tests/test_dynamics_jac.m`:

```matlab
% TEST_DYNAMICS_JAC  Complex-step vs analytic Jacobians of pdg_dynamics.
%
% Perturbs each of the 10 inputs (7 state + 3 thrust) with h=1e-30i and
% compares imag(f)/h against the analytic A and B columns. Run twice:
% vacuum (Phase 1 default) and with the drag branch forced on, at a state
% with nonzero velocity (the ||v||v drag term is not differentiable at
% v=0; we never linearize there).
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P  = booster_params();
x0 = [400; 80; 1500; -25; 5; -140; 28000];
T0 = [2e4; -1e4; 5e5];
h  = 1e-30;

for dragOn = [false true]
    P.drag.on = dragOn;
    [~, A, B] = pdg_dynamics(x0, T0, P);
    Acs = zeros(7,7);  Bcs = zeros(7,3);
    for k = 1:7
        xp = complex(x0);  xp(k) = xp(k) + 1i*h;
        Acs(:,k) = imag(pdg_dynamics(xp, complex(T0), P)) / h;
    end
    for k = 1:3
        Tp = complex(T0);  Tp(k) = Tp(k) + 1i*h;
        Bcs(:,k) = imag(pdg_dynamics(complex(x0), Tp, P)) / h;
    end
    errA = max(abs(A(:) - Acs(:))) / max(1, max(abs(Acs(:))));
    errB = max(abs(B(:) - Bcs(:))) / max(1, max(abs(Bcs(:))));
    assert(errA < 1e-12, 'A mismatch (drag=%d): %.2e', dragOn, errA);
    assert(errB < 1e-12, 'B mismatch (drag=%d): %.2e', dragOn, errB);
end
fprintf('test_dynamics_jac PASS (vacuum + drag)\n');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_dynamics_jac"`
Expected: FAIL — `pdg_dynamics` undefined.

- [ ] **Step 3: Implement `lib/pdg_dynamics.m`**

```matlab
function [xdot, A, B] = pdg_dynamics(x, T, P)
% PDG_DYNAMICS  3-DOF powered-descent dynamics with analytic Jacobians.
%
%   rdot = v
%   vdot = g + T/m + aD          aD = -(1/2) rho(z) Cd A |v| v / m  (opt-in)
%   mdot = -|T| / (Isp g0)
%
% Complex-step safe: magnitudes via sqrt(sum(.^2)), no norm/abs/max.
%
% INPUTS:
%   x - state [r(3); v(3); m] [7x1]  (SI: m, m/s, kg)
%   T - thrust vector [3x1] [N]
%   P - booster_params struct (uses gvec, Isp, g0, drag.*)
%
% OUTPUTS:
%   xdot - state derivative [7x1]
%   A    - d(xdot)/dx [7x7]   (only when requested)
%   B    - d(xdot)/dT [7x3]   (only when requested)
%
% REFERENCES:
%   [1] Acikmese & Ploen, "Convex Programming Approach to Powered Descent
%       Guidance for Mars Landing," JGCD 2007. (vacuum model)
v    = x(4:6);   m = x(7);
Tmag = sqrt(sum(T.^2));
aT   = T / m;
aD   = zeros(3,1);
if P.drag.on
    rho  = P.drag.rho0 * exp(-x(3) / P.drag.H);
    vmag = sqrt(sum(v.^2));
    kD   = 0.5 * rho * P.drag.Cd * P.drag.A;      % so aD = -kD |v| v / m
    aD   = -kD * vmag * v / m;
end
xdot = [v; P.gvec + aT + aD; -Tmag / (P.Isp * P.g0)];

if nargout > 1
    A          = zeros(7,7);
    A(1:3,4:6) = eye(3);
    A(4:6,7)   = -T / m^2;
    if P.drag.on
        % d(aD)/dv = -(kD/m)(|v| I + v v'/|v|); d/dz via drho/dz = -rho/H;
        % d/dm = +kD |v| v / m^2
        A(4:6,4:6) = -(kD/m) * (vmag*eye(3) + (v*v.')/vmag);
        A(4:6,3)   = A(4:6,3) + (kD*vmag/(P.drag.H*m)) * v;
        A(4:6,7)   = A(4:6,7) + kD * vmag * v / m^2;
    end
    B        = zeros(7,3);
    B(4:6,:) = eye(3) / m;
    B(7,:)   = -T.' / (Tmag * P.Isp * P.g0);
end
end
```

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: `test_dynamics_jac PASS (vacuum + drag)`.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/lib/pdg_dynamics.m booster_landing/tests/test_dynamics_jac.m
git commit -m "booster_landing: 3-DOF dynamics + analytic Jacobians, complex-step verified

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: Collocation guidance solver (`solve_pdg_colloc`)

**Files:**
- Create: `lib/solve_pdg_colloc.m`
- Create: `tests/test_colloc_smoke.m`

**Interfaces:**
- Consumes: `booster_params`, CasADi (`casadi.Opti`), `pdg_dynamics` *formulation* (the RHS is re-expressed symbolically inside the solver — CasADi variables can't flow through the plain-MATLAB `pdg_dynamics`; the certify layer (Task 5, gate G2) is what proves the two agree, by re-integrating the solution through `pdg_dynamics` with ode45).
- Produces: `sol = solve_pdg_colloc(P, opts)` with fields:
  - `sol.t` (1×(N+1)) node times, `sol.tf`, `sol.mf` final mass [kg]
  - `sol.X` (7×(N+1)) states at nodes; `sol.U` (3×(N+1)) node thrust [N]; `sol.Um` (3×N) midpoint thrust
  - `sol.lam_defect` (7×N) duals of the defect constraints, row-range extracted via `opti.lam_g` (house sign-bug lesson: never `opti.dual`)
  - `sol.stats` (IPOPT return status, iterations), `sol.P` (params snapshot)
  - `opts.N` (default `P.N`), `opts.init` (optional prior `sol` to warm-start), `opts.maxIter` (default 3000)

- [ ] **Step 1: Write the failing smoke test**

`tests/test_colloc_smoke.m`:

```matlab
% TEST_COLLOC_SMOKE  Coarse-grid (N=15) nonconvex NLP solve converges and
% obeys physics: solved status, mass in (mdry, m0), thrust annulus and
% glideslope satisfied at nodes, terminal state at the pad at rest.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
sol  = solve_pdg_colloc(P, struct('N', 15));
assert(sol.stats.success, 'IPOPT did not converge: %s', sol.stats.status);
assert(sol.mf > P.mdry && sol.mf < P.m0, 'final mass out of range');
Tmag = sqrt(sum(sol.U.^2, 1));
assert(all(Tmag >= P.Tmin - 1) && all(Tmag <= P.Tmax + 1), 'annulus violated');
rxy  = sqrt(sum(sol.X(1:2,:).^2, 1));
assert(all(rxy <= sol.X(3,:)/tand(P.gs_deg) + 1e-3), 'glideslope violated');
assert(max(abs(sol.X(1:6,end))) < 1e-3, 'terminal state not at rest on pad');
fprintf('test_colloc_smoke PASS  tf=%.2f s  mf=%.1f kg  fuel=%.1f kg\n', ...
        sol.tf, sol.mf, P.m0 - sol.mf);
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_colloc_smoke"`
Expected: FAIL — `solve_pdg_colloc` undefined.

- [ ] **Step 3: Implement `lib/solve_pdg_colloc.m`**

Formulation (write exactly this; header block per house standard):

```matlab
function sol = solve_pdg_colloc(P, opts)
% SOLVE_PDG_COLLOC  Min-fuel PDG as a free-tf Hermite-Simpson NLP (IPOPT).
%
% max m(tf) s.t. 3-DOF dynamics, thrust annulus Tmin<=|T|<=Tmax (nonconvex,
% handled directly), glideslope cone, optional pointing cone, r(tf)=v(tf)=0.
% Time normalized: t = tf*tau, tau in [0,1], tf a decision variable.
%
% INPUTS:
%   P    - booster_params struct
%   opts - (optional) .N segments [def P.N], .init prior sol (warm start),
%          .maxIter [def 3000]
% OUTPUTS:
%   sol  - .t .tf .mf .X(7xN+1) .U(3xN+1) .Um(3xN) .lam_defect(7xN)
%          .stats .P   (see plan Task 3 Interfaces)
%
% REFERENCES:
%   [1] Kelly, "An Introduction to Trajectory Optimization," SIAM Rev 2017.
if nargin < 2, opts = struct(); end
if ~isfield(opts,'N'),       opts.N = P.N;        end
if ~isfield(opts,'maxIter'), opts.maxIter = 3000; end
import casadi.*
N    = opts.N;
opti = casadi.Opti();

%% Decision variables:
X  = opti.variable(7, N+1);       % states at nodes
U  = opti.variable(3, N+1);       % thrust at nodes [N]
Um = opti.variable(3, N);         % thrust at HS midpoints [N]
tf = opti.variable();             % final time [s]

%% Dynamics as a CasADi expression (vacuum or drag per P.drag.on):
f = @(x,T) pdg_rhs_casadi(x, T, P);

%% Hermite-Simpson defects (register row range for dual extraction):
gDefStart = 1;                    % opti.g rows are appended in order
h = tf / N;
for k = 1:N
    xk = X(:,k);  xk1 = X(:,k+1);
    fk = f(xk, U(:,k));  fk1 = f(xk1, U(:,k+1));
    xm = 0.5*(xk + xk1) + (h/8)*(fk - fk1);
    fm = f(xm, Um(:,k));
    opti.subject_to(xk1 - xk - (h/6)*(fk + 4*fm + fk1) == 0);
end
nDefRows = 7 * N;

%% Path constraints (nodes and midpoints):
Uall = [U, Um];
for k = 1:size(Uall,2)
    T2 = sum(Uall(:,k).^2);
    opti.subject_to(P.Tmin^2 <= T2 <= P.Tmax^2);          % annulus
    if isfinite(P.theta_max_deg)
        opti.subject_to(Uall(3,k)^2 >= cosd(P.theta_max_deg)^2 * T2);
        opti.subject_to(Uall(3,k) >= 0);
    end
end
cotg = 1 / tand(P.gs_deg);
for k = 1:N+1
    opti.subject_to(X(1,k)^2 + X(2,k)^2 <= (cotg * X(3,k))^2);
    opti.subject_to(X(3,k) >= 0);
    opti.subject_to(X(7,k) >= P.mdry);
end

%% Boundary conditions + tf bounds:
opti.subject_to(X(1:3,1) == P.r0);
opti.subject_to(X(4:6,1) == P.v0);
opti.subject_to(X(7,1)   == P.m0);
opti.subject_to(X(1:6,end) == zeros(6,1));
opti.subject_to(P.tf_lo <= tf <= P.tf_hi);

%% Objective: min fuel == max final mass:
opti.minimize(-X(7,end));

%% Initial guess: straight line, gravity-cancelling thrust, or warm start:
if isfield(opts, 'init') && ~isempty(opts.init)
    s0 = opts.init;
    tauN = linspace(0,1,N+1);  tauM = (tauN(1:end-1)+tauN(2:end))/2;
    opti.set_initial(X,  interp1(s0.t/s0.tf, s0.X.',  tauN, 'pchip').');
    opti.set_initial(U,  interp1(s0.t/s0.tf, s0.U.',  tauN, 'pchip').');
    opti.set_initial(Um, interp1(s0.t/s0.tf, s0.U.',  tauM, 'pchip').');
    opti.set_initial(tf, s0.tf);
else
    tauN  = linspace(0,1,N+1);
    Xg    = [P.r0*(1-tauN); P.v0*(1-tauN); ...
             P.m0 + (P.mdry + 500 - P.m0)*tauN];
    Tg    = repmat(-0.9*P.m0*P.gvec, 1, N+1);      % ~hover thrust, straight up
    opti.set_initial(X, Xg);  opti.set_initial(U, Tg);
    opti.set_initial(Um, Tg(:,1:N));  opti.set_initial(tf, 30);
end

%% Solve:
sopts = struct('ipopt', struct('max_iter', opts.maxIter, ...
               'print_level', 3, 'tol', 1e-9));
opti.solver('ipopt', struct('print_time', false), sopts.ipopt);
try
    osol = opti.solve();
    ok   = true;
catch
    osol = opti.debug;                 % return best iterate for diagnosis
    ok   = false;
end

%% Package:
sol.tf   = full(osol.value(tf));
sol.t    = linspace(0, sol.tf, N+1);
sol.X    = full(osol.value(X));
sol.U    = full(osol.value(U));
sol.Um   = full(osol.value(Um));
sol.mf   = sol.X(7,end);
lam      = full(osol.value(opti.lam_g));
sol.lam_defect = reshape(lam(gDefStart : gDefStart+nDefRows-1), 7, N);
sol.stats = struct('success', ok, 'status', opti.stats.return_status, ...
                   'iter', opti.stats.iter_count);
sol.P    = P;
end

function xdot = pdg_rhs_casadi(x, T, P)
% PDG_RHS_CASADI  Same RHS as pdg_dynamics, CasADi-symbolic-safe.
% Gate G2 (certify_pdg) proves this matches pdg_dynamics by re-integration.
v    = x(4:6);  m = x(7);
Tmag = sqrt(sum(T.^2) + 1e-12);       % smooth at |T|=0 (never active: Tmin>0)
aD   = casadi.MX.zeros(3,1) * 0;      % vacuum default
if P.drag.on
    rho  = P.drag.rho0 * exp(-x(3)/P.drag.H);
    vmag = sqrt(sum(v.^2) + 1e-12);
    aD   = -0.5 * rho * P.drag.Cd * P.drag.A * vmag * v / m;
end
xdot = [v; P.gvec + T/m + aD; -Tmag/(P.Isp*P.g0)];
end
```

Implementation notes for the engineer:
- `opti.g` rows append in construction order; the defect block is constructed FIRST so its dual rows are `1 : 7N`. Keep that ordering or update `gDefStart`.
- `casadi.MX.zeros(3,1)*0` line: replace with a plain `aD = 0*v;` if MX zeros is awkward — anything symbolic-safe.
- If the cold start fails to converge at N=15, first try `opts.maxIter=6000`; if still failing, the documented fallback is warm-start from the convex solution (Task 4) — wire `opts.init` from it in the front door, but the smoke test must pass cold (spec requirement).

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: `test_colloc_smoke PASS tf=... mf=...`. Record the printed tf/mf/fuel numbers in the commit message — they are the first campaign truth.

- [ ] **Step 5: Nominal-grid solve (N=60) sanity, save reference**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; P=booster_params(); sol=solve_pdg_colloc(P); fprintf('N=60: tf=%.3f mf=%.2f\n', sol.tf, sol.mf); save(fullfile('results','pdg_colloc_nominal.mat'),'sol')"`
Expected: converged; tf within ~1 s and mf within ~5 kg of the N=15 answer (grid convergence sanity). Throttle profile should show max–min–max; eyeball `plot(sol.t, sqrt(sum(sol.U.^2,1)))` in Task 9's viz, not here.

- [ ] **Step 6: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/lib/solve_pdg_colloc.m booster_landing/tests/test_colloc_smoke.m
git commit -m "booster_landing: free-tf Hermite-Simpson min-fuel PDG solver (IPOPT)

N=15 smoke: tf=<fill from run> s, mf=<fill> kg. N=60 consistent.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: Convex guidance solver (`solve_pdg_convex`) + golden-section tf

**Files:**
- Create: `lib/solve_pdg_convex.m`
- Create: `tests/test_convex_lossless.m`

**Interfaces:**
- Consumes: `booster_params`, CasADi.
- Produces: `sol = solve_pdg_convex(P, opts)`:
  - inner fixed-tf convex solve + outer golden-section on tf over `[P.tf_lo, P.tf_hi]`
  - `sol.t` (1×Nconv), `sol.tf`, `sol.mf`
  - `sol.X` (7×Nconv) with `m = exp(z)` reconstructed; `sol.U` (3×Nconv) **thrust in N** (`T = m·u`, same units/meaning as colloc `sol.U`)
  - `sol.u` (3×Nconv) acceleration controls, `sol.sigma` (1×Nconv) slack — kept for the losslessness gate
  - `sol.lossless_gap` = `max |‖u‖ − σ|` over nodes
  - `sol.tf_curve` (Kx2, the (tf, mf) points the golden search evaluated — for the note's figure)
  - `opts.tolTf` (default 0.05 s)

- [ ] **Step 1: Write the failing test**

`tests/test_convex_lossless.m`:

```matlab
% TEST_CONVEX_LOSSLESS  Fixed-tf convexified PDG at a plausible tf: solved
% status, losslessness gap ||u||-sigma ~ 0 (the relaxation is TIGHT at the
% optimum -- the Acikmese/Ploen theorem, checked numerically), annulus in
% original variables, terminal conditions met.
%
% Uses solve_pdg_convex's single-tf mode (opts.tf fixed) to stay fast.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P   = booster_params();
sol = solve_pdg_convex(P, struct('tf', 25, 'Nconv', 60));
assert(sol.stats.success, 'convex solve failed: %s', sol.stats.status);
assert(sol.lossless_gap < 1e-4 * P.Tmax/P.m0, ...
       'relaxation not tight: gap=%.3e', sol.lossless_gap);
Tmag = sqrt(sum(sol.U.^2,1));
assert(all(Tmag >= P.Tmin*(1-1e-4)) && all(Tmag <= P.Tmax*(1+1e-4)), ...
       'annulus violated in original variables');
assert(max(abs(sol.X(1:6,end))) < 1e-2, 'terminal not met');
fprintf('test_convex_lossless PASS  gap=%.2e  mf=%.1f kg\n', ...
        sol.lossless_gap, sol.mf);
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_convex_lossless"`
Expected: FAIL — `solve_pdg_convex` undefined.

- [ ] **Step 3: Implement `lib/solve_pdg_convex.m`**

Structure: private `solve_fixed_tf(P, tf, Nc)` + outer golden section. Key formulation (Blackmore 2010 change of variables):

```matlab
function sol = solve_pdg_convex(P, opts)
% SOLVE_PDG_CONVEX  Min-fuel PDG via lossless convexification (fixed-tf
% convex subproblem + golden-section search on tf).
%
% Change of variables u=T/m, sigma=Gamma/m, z=ln m turns dynamics linear;
% the nonconvex annulus relaxes to ||u|| <= sigma with linear/quadratic
% Taylor bounds on sigma about the max-thrust depletion reference z0(t).
% The relaxation is provably tight (lossless) at the optimum; we CHECK
% that numerically (sol.lossless_gap) rather than assume it.
%
% Convex problem solved with IPOPT: local = global by convexity. A conic
% solver (ECOS/CVX) is a documented drop-in, deliberately not a dependency.
%
% INPUTS:
%   P    - booster_params struct
%   opts - (optional) .tf fixed final time (skips golden search),
%          .Nconv [def P.Nconv], .tolTf golden tolerance [def 0.05 s]
% OUTPUTS:
%   sol  - .t .tf .mf .X .U .u .sigma .lossless_gap .tf_curve .stats .P
%
% REFERENCES:
%   [1] Acikmese & Ploen, JGCD 2007.  [2] Blackmore et al., JGCD 2010.
if nargin < 2, opts = struct(); end
if ~isfield(opts,'Nconv'), opts.Nconv = P.Nconv; end
if ~isfield(opts,'tolTf'), opts.tolTf = 0.05;    end

if isfield(opts,'tf') && ~isempty(opts.tf)
    sol = solve_fixed_tf(P, opts.tf, opts.Nconv);
    sol.tf_curve = [opts.tf, sol.mf];
    return
end

%% Golden-section search on tf, maximizing mf (see note for unimodality):
phi = (sqrt(5)-1)/2;
a = P.tf_lo;  b = P.tf_hi;  curve = [];
c = b - phi*(b-a);  d = a + phi*(b-a);
sc = solve_fixed_tf(P, c, opts.Nconv);  sd = solve_fixed_tf(P, d, opts.Nconv);
curve = [curve; c, mf_or_neginf(sc); d, mf_or_neginf(sd)];
while (b - a) > opts.tolTf
    if mf_or_neginf(sc) > mf_or_neginf(sd)
        b = d;  d = c;  sd = sc;
        c = b - phi*(b-a);  sc = solve_fixed_tf(P, c, opts.Nconv);
        curve = [curve; c, mf_or_neginf(sc)];               %#ok<AGROW>
    else
        a = c;  c = d;  sc = sd;
        d = a + phi*(b-a);  sd = solve_fixed_tf(P, d, opts.Nconv);
        curve = [curve; d, mf_or_neginf(sd)];               %#ok<AGROW>
    end
end
if mf_or_neginf(sc) > mf_or_neginf(sd), sol = sc; else, sol = sd; end
sol.tf_curve = sortrows(curve, 1);
end

function v = mf_or_neginf(s)
% Infeasible tf (too short to stop) shows as solver failure -> -Inf.
if s.stats.success, v = s.mf; else, v = -Inf; end
end

function sol = solve_fixed_tf(P, tf, Nc)
% One convex subproblem at fixed tf, trapezoidal on the LINEAR dynamics.
import casadi.*
opti = casadi.Opti();
t   = linspace(0, tf, Nc);  h = t(2) - t(1);
R   = opti.variable(3, Nc);   V = opti.variable(3, Nc);
Z   = opti.variable(1, Nc);   Uu = opti.variable(3, Nc);
S   = opti.variable(1, Nc);   % sigma
al  = 1 / (P.Isp * P.g0);

%% Linear dynamics, trapezoid:
for k = 1:Nc-1
    opti.subject_to(R(:,k+1) == R(:,k) + (h/2)*(V(:,k)+V(:,k+1)));
    opti.subject_to(V(:,k+1) == V(:,k) + (h/2)*(2*P.gvec + Uu(:,k)+Uu(:,k+1)));
    opti.subject_to(Z(k+1)   == Z(k)   - (h/2)*al*(S(k)+S(k+1)));
end

%% Relaxed annulus + Taylor mass bounds about z0(t) (max-thrust depletion):
z0  = log(P.m0 - al*P.Tmax*t);            % reference depletion
zlb = log(P.m0 - al*P.Tmax*t);            % lower bound on z
zub = log(P.m0 - al*P.Tmin*t);            % upper bound on z
cotg = 1/tand(P.gs_deg);
for k = 1:Nc
    opti.subject_to(sum(Uu(:,k).^2) <= S(k)^2);
    opti.subject_to(S(k) >= 0);
    mu1 = P.Tmin*exp(-z0(k));  mu2 = P.Tmax*exp(-z0(k));
    dz  = Z(k) - z0(k);
    opti.subject_to(S(k) >= mu1*(1 - dz + dz^2/2));
    opti.subject_to(S(k) <= mu2*(1 - dz));
    opti.subject_to(zlb(k) <= Z(k) <= zub(k));
    opti.subject_to(R(1,k)^2 + R(2,k)^2 <= (cotg*R(3,k))^2);
    opti.subject_to(R(3,k) >= 0);
    if isfinite(P.theta_max_deg)
        opti.subject_to(Uu(3,k) >= cosd(P.theta_max_deg)*S(k));
    end
end

%% Boundary conditions, objective:
opti.subject_to(R(:,1) == P.r0);   opti.subject_to(V(:,1) == P.v0);
opti.subject_to(Z(1)   == log(P.m0));
opti.subject_to(R(:,end) == zeros(3,1));  opti.subject_to(V(:,end) == zeros(3,1));
opti.subject_to(Z(end) >= log(P.mdry));
opti.minimize(-Z(end));

opti.solver('ipopt', struct('print_time', false), ...
            struct('max_iter', 1000, 'print_level', 0, 'tol', 1e-10));
try
    osol = opti.solve();  ok = true;
catch
    osol = opti.debug;    ok = false;
end

%% Package in ORIGINAL variables (T = m u):
sol.t  = t;  sol.tf = tf;
m      = exp(full(osol.value(Z)));
sol.u  = full(osol.value(Uu));  sol.sigma = full(osol.value(S));
sol.X  = [full(osol.value(R)); full(osol.value(V)); m];
sol.U  = sol.u .* m;
sol.mf = m(end);
sol.lossless_gap = max(abs(sqrt(sum(sol.u.^2,1)) - sol.sigma));
sol.stats = struct('success', ok, 'status', opti.stats.return_status);
sol.P  = P;
end
```

Implementation notes:
- The bracket `[P.tf_lo, P.tf_hi]` must contain a feasible tf; if BOTH golden probes come back −Inf, error out with a clear message telling the user to widen the bracket, don't silently return garbage.
- `zlb` uses Tmax depletion (fastest mass loss = lowest mass = lower z bound), `zub` uses Tmin. Guard `P.m0 - al*P.Tmax*t > 0` — with tf ≤ 50 s and these params it holds (burn ≤ ~15.6 t), but assert it anyway at the top of `solve_fixed_tf`.
- Unimodality of mf(tf) is assumed for golden section (standard in this literature); `sol.tf_curve` is the evidence plot — if it ever shows multimodality, fall back to a fine scan. Say so in a comment.

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: `test_convex_lossless PASS gap=... mf=...`.

- [ ] **Step 5: Full golden-section run, compare to colloc by eye**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; P=booster_params(); sol=solve_pdg_convex(P); fprintf('convex: tf=%.3f mf=%.2f gap=%.2e\n', sol.tf, sol.mf, sol.lossless_gap); save(fullfile('results','pdg_convex_nominal.mat'),'sol')"`
Expected: tf and mf within ~1% of Task 3 Step 5's collocation numbers. If they disagree grossly, STOP and debug now (systematic-debugging skill) — do not proceed to certification with a known discrepancy.

- [ ] **Step 6: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/lib/solve_pdg_convex.m booster_landing/tests/test_convex_lossless.m
git commit -m "booster_landing: lossless-convexification PDG + golden-section tf

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: Certification gates (`certify_pdg`)

**Files:**
- Create: `certify/certify_pdg.m`
- Create: `tests/test_certify_nominal.m`

**Interfaces:**
- Consumes: `sol` structs from both solvers (Tasks 3–4), `pdg_dynamics`.
- Produces: `rep = certify_pdg(solC, solV, P)` — gates report struct:
  - `rep.G1_defect` (max HS defect re-evaluated at solution), `rep.G1_pass`
  - `rep.G2_resid` (ode45 re-integration: terminal position/velocity/mass error), `rep.G2_pass` (pos < 1 m, vel < 0.1 m/s, mass < 0.5 kg)
  - `rep.G3_dmf` (|mf colloc − convex|), `rep.G3_dtf`, `rep.G3_traj_Linf` (position, common time grid), `rep.G3_pass` (dmf < 1.0 kg, dtf < 0.2 s) — ADJUDICATED 2026-08-08: G3 measures the convex Taylor-bound model error (~0.7 kg genuine offset); gate raised 0.1→1.0 kg per user decision
  - `rep.G4_gap` (= solV.lossless_gap), `rep.G4_pass` (< 1e-4·Tmax/m0)
  - `rep.G5_structure` (throttle max–min–max: fraction of nodes on a bound ≥ 0.95, interior switch count ≤ 2), `rep.G5_primer_deg` (max angle between thrust and −λv), `rep.G5_pass` (structure ok AND primer < 1°)
  - `rep.all_pass`; `print_certify_report(rep)` nested pretty-printer (gate table, PASS/FAIL per row — the front door prints this)
- `solV = []` allowed (Phase 2 drag runs have no convex twin): G3/G4 report `skipped`, excluded from `all_pass`.

- [ ] **Step 1: Write the failing test**

`tests/test_certify_nominal.m`:

```matlab
% TEST_CERTIFY_NOMINAL  Full gate run on coarse solves (fast): all five
% gates must pass on the nominal vacuum problem. This is the campaign's
% core scientific claim -- two independent methods, one answer, PMP-shaped.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
solC = solve_pdg_colloc(P, struct('N', 30));
solV = solve_pdg_convex(P, struct('Nconv', 90));
rep  = certify_pdg(solC, solV, P);
print_certify_report(rep);
assert(rep.all_pass, 'certification failed -- see report above');
fprintf('test_certify_nominal PASS\n');
```

(Gate tolerances at coarse grids: if G3 dmf < 0.1 kg proves too tight at N=30 while fine grids pass, loosen the COARSE-test call site via an optional `tol` scale arg to `certify_pdg`, never the nominal tolerance itself.)

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_certify_nominal"`
Expected: FAIL — `certify_pdg` undefined.

- [ ] **Step 3: Implement `certify/certify_pdg.m`**

Key pieces (full header per standard; `print_certify_report` in the same file as a second function):

```matlab
function rep = certify_pdg(solC, solV, P, tolScale)
% CERTIFY_PDG  Gates G1-G5 for a PDG solution pair (colloc + convex).
% See plan Task 5 for gate definitions; report-only, throws never --
% callers decide what a FAIL means.
% INPUTS: solC colloc sol; solV convex sol or [] (G3/G4 skipped);
%         P params; tolScale (optional, def 1) scales agreement tols for
%         coarse-grid test calls.
% OUTPUTS: rep - gate report struct (fields in plan Task 5).
if nargin < 4, tolScale = 1; end

%% G1: HS defects re-evaluated in plain MATLAB (independent of CasADi):
N = size(solC.X,2) - 1;  h = solC.tf / N;  dmax = 0;
for k = 1:N
    fk  = pdg_dynamics(solC.X(:,k),   solC.U(:,k),   P);
    fk1 = pdg_dynamics(solC.X(:,k+1), solC.U(:,k+1), P);
    xm  = 0.5*(solC.X(:,k)+solC.X(:,k+1)) + (h/8)*(fk - fk1);
    fm  = pdg_dynamics(xm, solC.Um(:,k), P);
    d   = solC.X(:,k+1) - solC.X(:,k) - (h/6)*(fk + 4*fm + fk1);
    dmax = max(dmax, max(abs(d)));
end
rep.G1_defect = dmax;   rep.G1_pass = dmax < 1e-6;

%% G2: continuous residual -- fly the interpolated control with ode45
%% ("defect is not accuracy"): pchip on [nodes+midpoints] thrust.
tU = sort([solC.t, solC.t(1:end-1) + h/2]);
TU = zeros(3, numel(tU));
TU(:,1:2:end) = solC.U;  TU(:,2:2:end) = solC.Um;
ctrl = @(tt) interp1(tU.', TU.', min(tt, solC.tf), 'pchip').';
odef = @(tt, xx) pdg_dynamics(xx, ctrl(tt), P);
oo   = odeset('RelTol',1e-10,'AbsTol',1e-10);
[~, XX] = ode45(odef, [0 solC.tf], solC.X(:,1), oo);
ef = XX(end,:).' - solC.X(:,end);
rep.G2_pos = sqrt(sum(ef(1:3).^2));  rep.G2_vel = sqrt(sum(ef(4:6).^2));
rep.G2_dm  = abs(ef(7));
rep.G2_pass = rep.G2_pos < 1 && rep.G2_vel < 0.1 && rep.G2_dm < 0.5;

%% G3/G4: cross-method agreement + losslessness (skip if no convex twin):
if isempty(solV)
    rep.G3_pass = 'skipped';  rep.G4_pass = 'skipped';
else
    rep.G3_dmf = abs(solC.mf - solV.mf);
    rep.G3_dtf = abs(solC.tf - solV.tf);
    tq  = linspace(0, min(solC.tf, solV.tf), 200);
    rC  = interp1(solC.t.', solC.X(1:3,:).', tq.', 'pchip');
    rV  = interp1(solV.t.', solV.X(1:3,:).', tq.', 'pchip');
    rep.G3_traj_Linf = max(sqrt(sum((rC - rV).^2, 2)));
    rep.G3_pass = rep.G3_dmf < 0.1*tolScale && rep.G3_dtf < 0.2*tolScale;
    rep.G4_gap  = solV.lossless_gap;
    rep.G4_pass = rep.G4_gap < 1e-4 * P.Tmax / P.m0;
end

%% G5: PMP structure. Throttle bang-bang max-min-max + primer alignment.
Tmag = sqrt(sum(solC.U.^2, 1));
onLo = Tmag < P.Tmin * 1.001;   onHi = Tmag > P.Tmax * 0.999;
rep.G5_bound_frac = mean(onLo | onHi);
segs = diff([onHi(1), onHi]);                  % max<->min transitions
rep.G5_switches = sum(segs ~= 0);
structOk = rep.G5_bound_frac >= 0.95 && rep.G5_switches <= 2 ...
           && onHi(1) && onHi(end);            % max-first, max-last
% Primer: velocity costate from defect duals; thrust must align with -lam_v.
lamv = solC.lam_defect(4:6, :);                % 3xN, node-adjacent duals
Tdir = solC.U(:,1:end-1) ./ sqrt(sum(solC.U(:,1:end-1).^2,1));
pdir = -lamv ./ max(sqrt(sum(lamv.^2,1)), 1e-30);
cosang = sum(Tdir .* pdir, 1);
rep.G5_primer_deg = max(acosd(min(1, max(-1, cosang))));
rep.G5_pass = structOk && rep.G5_primer_deg < 1;

%% Verdict:
gates = [rep.G1_pass, rep.G2_pass, isequal(rep.G3_pass,true) || ...
         isequal(rep.G3_pass,'skipped'), isequal(rep.G4_pass,true) || ...
         isequal(rep.G4_pass,'skipped'), rep.G5_pass];
rep.all_pass = all(gates);
end
```

Notes:
- The primer-dual sign: if `G5_primer_deg` comes out near 180° instead of ~0°, the dual sign convention is flipped — flip `pdir` sign ONCE, document why in a comment referencing the `opti.lam_g` house lesson, and verify the switching structure still reads max–min–max. Do not add an `abs()` to make the test pass (that would hide a real sign error).
- `print_certify_report(rep)`: plain fprintf table, one row per gate, value + threshold + PASS/FAIL/skipped. Write it in the same file below `certify_pdg` — MATLAB allows multiple functions per file; expose it by making `certify/print_certify_report.m` a 3-line wrapper if needed elsewhere.

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: gate table printed, all PASS, `test_certify_nominal PASS`. If G5 primer fails on sign, apply the note above; if G3 fails, this is the moment the two methods genuinely disagree — systematic-debugging, likely suspects: Taylor mass-bound error (refine `Nconv`), or tf bracket edge.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/certify/certify_pdg.m booster_landing/tests/test_certify_nominal.m
git commit -m "booster_landing: certification gates G1-G5 (defect, residual, cross-method, lossless, PMP)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 6: TVLQR design (`tvlqr_design`)

**Files:**
- Create: `lib/tvlqr_design.m`
- Create: `tests/test_tvlqr_riccati.m`

**Interfaces:**
- Consumes: `sol` (colloc), `pdg_dynamics` analytic A,B.
- Produces: `ctrl = tvlqr_design(sol, P, opts)`:
  - `ctrl.tgrid` (1×M dense, M = 4·(N+1)), `ctrl.K` (3×7×M) gains, `ctrl.Pt` (7×7×M) Riccati solution
  - `ctrl.xnom` / `ctrl.Tnom` — pchip-interpolant handles `@(t)->7x1 / 3x1` of the nominal trajectory/control (what `sim_closed_loop` tracks)
  - `opts.Q` (default `diag([1e-4 1e-4 1e-4 1e-2 1e-2 1e-2 0])`), `opts.R` (default `1e-10*eye(3)`), `opts.Qf` (default `diag([1e-2 1e-2 1e-2 1 1 1 0])`) — SI units; mass column zero (mass is observable, not directly regulated)

- [ ] **Step 1: Write the failing test**

`tests/test_tvlqr_riccati.m`:

```matlab
% TEST_TVLQR_RICCATI  Riccati solution health along the nominal trajectory:
% P(t) symmetric positive semidefinite everywhere (PSD, not PD: the mass
% row/column is unweighted), gains finite, terminal condition P(tf)=Qf.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
sol  = solve_pdg_colloc(P, struct('N', 30));
assert(sol.stats.success);
ctrl = tvlqr_design(sol, P);
M    = numel(ctrl.tgrid);
for k = 1:M
    Pk = ctrl.Pt(:,:,k);
    assert(max(max(abs(Pk - Pk.'))) < 1e-6 * max(1,max(abs(Pk(:)))), ...
           'P not symmetric at k=%d', k);
    assert(min(eig((Pk+Pk.')/2)) > -1e-8 * max(1,max(abs(Pk(:)))), ...
           'P not PSD at k=%d', k);
end
assert(all(isfinite(ctrl.K(:))), 'non-finite gains');
Qf = diag([1e-2 1e-2 1e-2 1 1 1 0]);
assert(max(max(abs(ctrl.Pt(:,:,end) - Qf))) < 1e-9, 'terminal P ~= Qf');
fprintf('test_tvlqr_riccati PASS\n');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_tvlqr_riccati"`
Expected: FAIL — `tvlqr_design` undefined.

- [ ] **Step 3: Implement `lib/tvlqr_design.m`**

```matlab
function ctrl = tvlqr_design(sol, P, opts)
% TVLQR_DESIGN  Time-varying LQR about a PDG guidance trajectory.
%
% Linearizes pdg_dynamics along (x*(t), T*(t)) and integrates the
% differential Riccati equation backward:
%   -Pdot = A'P + PA - PBR^{-1}B'P + Q,   P(tf) = Qf
% Gains K(t) = R^{-1} B(t)' P(t), stored on a dense grid; controller is
%   T_cmd(t) = T*(t) - K(t) (x - x*(t))  (saturation applied by the sim).
%
% INPUTS:
%   sol  - collocation solution (Task 3 interface)
%   P    - booster_params
%   opts - (optional) .Q .R .Qf weight overrides (7x7, 3x3, 7x7)
% OUTPUTS:
%   ctrl - .tgrid(1xM) .K(3x7xM) .Pt(7x7xM) .xnom(@t->7x1) .Tnom(@t->3x1)
%
% REFERENCES:
%   [1] Anderson & Moore, "Optimal Control: Linear Quadratic Methods."
if nargin < 3, opts = struct(); end
if ~isfield(opts,'Q'),  opts.Q  = diag([1e-4 1e-4 1e-4 1e-2 1e-2 1e-2 0]); end
if ~isfield(opts,'R'),  opts.R  = 1e-10*eye(3);                            end
if ~isfield(opts,'Qf'), opts.Qf = diag([1e-2 1e-2 1e-2 1 1 1 0]);          end

%% Nominal interpolants (nodes+midpoints for thrust, nodes for state):
Nn = size(sol.X,2) - 1;  h = sol.tf/Nn;
tU = sort([sol.t, sol.t(1:end-1) + h/2]);
TU = zeros(3, numel(tU));  TU(:,1:2:end) = sol.U;  TU(:,2:2:end) = sol.Um;
ctrl.xnom = @(t) interp1(sol.t.', sol.X.', clampt(t, sol.tf), 'pchip').';
ctrl.Tnom = @(t) interp1(tU.',   TU.',    clampt(t, sol.tf), 'pchip').';

%% Backward Riccati on vec(P), dense output grid:
M     = 4*(Nn+1);
ctrl.tgrid = linspace(0, sol.tf, M);
Rinv  = inv(opts.R);
ric   = @(t, pv) ricrhs(t, pv, ctrl, P, opts.Q, Rinv);
oo    = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[~, PV] = ode45(ric, fliplr(ctrl.tgrid), opts.Qf(:), oo);   % integrates tf->0
PV    = flipud(PV);                                          % re-order 0->tf
ctrl.Pt = zeros(7,7,M);  ctrl.K = zeros(3,7,M);
for k = 1:M
    Pk = reshape(PV(k,:), 7, 7);  Pk = (Pk + Pk.')/2;
    [~, ~, Bk] = pdg_dynamics(ctrl.xnom(ctrl.tgrid(k)), ...
                              ctrl.Tnom(ctrl.tgrid(k)), P);
    ctrl.Pt(:,:,k) = Pk;
    ctrl.K(:,:,k)  = Rinv * (Bk.' * Pk);
end
end

function t = clampt(t, tf), t = min(max(t, 0), tf); end

function pdot = ricrhs(t, pv, ctrl, P, Q, Rinv)
Pk = reshape(pv, 7, 7);
[~, A, B] = pdg_dynamics(ctrl.xnom(t), ctrl.Tnom(t), P);
Pd   = -(A.'*Pk + Pk*A - Pk*B*Rinv*B.'*Pk + Q);
pdot = Pd(:);
end
```

Note: `ode45(ric, fliplr(tgrid), ...)` with a decreasing tspan integrates backward directly — verify PV row order corresponds to the flipped grid before `flipud` (check `PV(1,:)` reshapes to ≈ Qf).

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: `test_tvlqr_riccati PASS`.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/lib/tvlqr_design.m booster_landing/tests/test_tvlqr_riccati.m
git commit -m "booster_landing: TVLQR gains via backward Riccati along guidance

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 7: Closed-loop simulation (`sim_closed_loop`)

**Files:**
- Create: `lib/sim_closed_loop.m`
- Create: `tests/test_closed_loop_nominal.m`

**Interfaces:**
- Consumes: `sol`, `ctrl`, `pdg_dynamics`.
- Produces: `out = sim_closed_loop(sol, ctrl, P, dsp)`:
  - `dsp` dispersion struct (all optional, default 0): `.dr0(3x1) .dv0(3x1) .thrust_scale (1=nominal) .isp_scale .wind(3x1 const m/s, drag-frame — only felt when P.drag.on)`
  - `out.t`, `out.X` (Mx7), `out.Tcmd` (Mx3, post-saturation), `out.td` struct: `.r(3x1) .v(3x1) .m .miss` (horizontal distance to pad) `.vtd` (touchdown speed) — from the z=0 ode event
  - `out.sat_frac` — fraction of time on a thrust bound (diagnostic for the expected terminal-arc saturation)
- Controller law inside: `Traw = ctrl.Tnom(t) - K(t)*(x - ctrl.xnom(t))`, magnitude clamped to `[P.Tmin, P.Tmax]` (direction kept), `K(t)` linear-interpolated from `ctrl.K`. Truth dynamics apply `dsp.thrust_scale` and `dsp.isp_scale` (the plant differs from the model — that's the point).

- [ ] **Step 1: Write the failing test**

`tests/test_closed_loop_nominal.m`:

```matlab
% TEST_CLOSED_LOOP_NOMINAL  Zero-dispersion closed loop must reproduce the
% guidance: touchdown within 1 m and 0.1 m/s of the pad-at-rest target.
% Also: with a 50 m initial position offset, the tracker must still land
% inside the pad radius (the whole point of feedback).
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
sol  = solve_pdg_colloc(P, struct('N', 30));
ctrl = tvlqr_design(sol, P);

out0 = sim_closed_loop(sol, ctrl, P, struct());
assert(out0.td.miss < 1.0,  'nominal miss %.2f m', out0.td.miss);
assert(out0.td.vtd  < 0.1 + 1e-9, 'nominal touchdown speed %.3f', out0.td.vtd);

dsp  = struct('dr0', [50; -30; 0]);
out1 = sim_closed_loop(sol, ctrl, P, dsp);
assert(out1.td.miss < P.pad_radius, 'dispersed miss %.1f m', out1.td.miss);
assert(out1.td.vtd  < P.vtd_max,    'dispersed vtd %.2f m/s', out1.td.vtd);
fprintf('test_closed_loop_nominal PASS  (nom miss %.3f m, disp miss %.2f m)\n', ...
        out0.td.miss, out1.td.miss);
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_closed_loop_nominal"`
Expected: FAIL — `sim_closed_loop` undefined.

- [ ] **Step 3: Implement `lib/sim_closed_loop.m`**

```matlab
function out = sim_closed_loop(sol, ctrl, P, dsp)
% SIM_CLOSED_LOOP  Truth-model landing sim under TVLQR tracking.
%
% Plant = pdg_dynamics with dispersion multipliers (thrust/Isp bias, wind);
% controller = T*(t) - K(t) dx, magnitude-saturated to [Tmin, Tmax] with
% direction preserved. Integrates to the z=0 crossing (ode event).
%
% INPUTS:
%   sol,ctrl - Task 3 / Task 6 interfaces
%   P        - booster_params
%   dsp      - dispersions: .dr0 .dv0 [3x1], .thrust_scale .isp_scale
%              [scalars, def 1], .wind [3x1 m/s, needs P.drag.on] -- all
%              optional, defaults = nominal
% OUTPUTS:
%   out - .t(Mx1) .X(Mx7) .Tcmd(Mx3) .sat_frac
%         .td struct: .r .v .m .miss .vtd
d = struct('dr0',zeros(3,1), 'dv0',zeros(3,1), 'thrust_scale',1, ...
           'isp_scale',1, 'wind',zeros(3,1));
if nargin >= 4
    fn = fieldnames(dsp);
    for k = 1:numel(fn), d.(fn{k}) = dsp.(fn{k}); end
end

Pp = P;  Pp.Isp = P.Isp * d.isp_scale;      % plant params differ from model
x0 = [P.r0 + d.dr0; P.v0 + d.dv0; P.m0];

oo = odeset('RelTol',1e-8, 'AbsTol',1e-8, 'Events', @touchdown_event, ...
            'MaxStep', 0.25);
[tt, XX] = ode45(@(t,x) plant_rhs(t, x, sol, ctrl, Pp, d), ...
                 [0, 1.5*sol.tf], x0, oo);

out.t = tt;  out.X = XX;
Tc = zeros(numel(tt), 3);
for k = 1:numel(tt)
    Tc(k,:) = control_law(tt(k), XX(k,:).', ctrl, P).';
end
out.Tcmd = Tc;
Tmag = sqrt(sum(Tc.^2, 2));
out.sat_frac = mean(Tmag > P.Tmax*0.999 | Tmag < P.Tmin*1.001);
xe = XX(end,:).';
out.td = struct('r', xe(1:3), 'v', xe(4:6), 'm', xe(7), ...
                'miss', sqrt(sum(xe(1:2).^2)), 'vtd', sqrt(sum(xe(4:6).^2)));
end

function T = control_law(t, x, ctrl, P)
% TVLQR + magnitude saturation, direction preserved. K interp per column.
tq   = min(max(t, ctrl.tgrid(1)), ctrl.tgrid(end));
Kt   = zeros(3,7);
for r = 1:3
    Kt(r,:) = interp1(ctrl.tgrid.', squeeze(ctrl.K(r,:,:)).', tq, 'linear');
end
Traw = ctrl.Tnom(tq) - Kt * (x - ctrl.xnom(tq));
Tm   = sqrt(sum(Traw.^2));
T    = Traw * min(max(Tm, P.Tmin), P.Tmax) / max(Tm, 1e-9);
end

function xdot = plant_rhs(t, x, sol, ctrl, Pp, d)
T    = control_law(t, x, ctrl, Pp) * d.thrust_scale;
if Pp.drag.on                     % wind enters as airspeed shift
    xw = x;  xw(4:6) = x(4:6) - d.wind;
    xdot = pdg_dynamics(xw, T, Pp);
    xdot(1:3) = x(4:6);           % kinematics use ground velocity
else
    xdot = pdg_dynamics(x, T, Pp);
end
end

function [val, isterm, dir_] = touchdown_event(~, x)
val = x(3);  isterm = 1;  dir_ = -1;   % z falling through 0
end
```

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: PASS with nominal miss well under 1 m. If nominal miss is large: the guidance interpolants and the plant integrate slightly differently — check `ctrl.Tnom` uses nodes+midpoints (it does), then check the event isn't firing early (initial z must be positive). If the dispersed case misses: gains too soft — revisit `opts.Q/R` in Task 6, document the change there.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/lib/sim_closed_loop.m booster_landing/tests/test_closed_loop_nominal.m
git commit -m "booster_landing: closed-loop TVLQR landing sim with touchdown event

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 8: Monte Carlo (`run_monte_carlo`)

**Files:**
- Create: `lib/run_monte_carlo.m`
- Create: `tests/test_monte_carlo_small.m`

**Interfaces:**
- Consumes: `sol`, `ctrl`, `sim_closed_loop`.
- Produces: `mc = run_monte_carlo(sol, ctrl, P, opts)`:
  - `opts.Nrun` (default 200), `opts.sig` struct of 1σ values (defaults: `.r0 = [100;100;50]`, `.v0 = [10;10;10]`, `.thrust = 0.015`, `.isp = 0.01`, `.wind = [10;10;0]` — wind only drawn when `P.drag.on`)
  - `mc.land` (Nrun×2 touchdown xy), `mc.vtd` (Nrun×1), `mc.mprop` (Nrun×1 propellant remaining above dry), `mc.ok` (Nrun×1 logical: miss < pad_radius AND vtd < vtd_max AND m ≥ mdry), `mc.success_rate`
  - Deterministic: `rng(P.seed)` once at top; draws via `randn`.

- [ ] **Step 1: Write the failing test**

`tests/test_monte_carlo_small.m`:

```matlab
% TEST_MONTE_CARLO_SMALL  20-run MC sanity: executes, deterministic under
% the seed (two calls identical), success rate positive, outputs shaped.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
sol  = solve_pdg_colloc(P, struct('N', 30));
ctrl = tvlqr_design(sol, P);
mc1  = run_monte_carlo(sol, ctrl, P, struct('Nrun', 20));
mc2  = run_monte_carlo(sol, ctrl, P, struct('Nrun', 20));
assert(isequal(size(mc1.land), [20 2]) && numel(mc1.vtd) == 20);
assert(isequal(mc1.land, mc2.land), 'MC not deterministic under seed');
assert(mc1.success_rate > 0, 'zero successes at default dispersions');
fprintf('test_monte_carlo_small PASS  success=%.0f%%\n', 100*mc1.success_rate);
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_monte_carlo_small"`
Expected: FAIL — `run_monte_carlo` undefined.

- [ ] **Step 3: Implement `lib/run_monte_carlo.m`**

```matlab
function mc = run_monte_carlo(sol, ctrl, P, opts)
% RUN_MONTE_CARLO  Dispersed closed-loop landing campaign.
%
% Draws initial-state, thrust-bias, Isp-bias (and, with drag on, wind)
% dispersions and runs sim_closed_loop per sample. Deterministic under
% P.seed. Success = miss < pad_radius AND vtd < vtd_max AND m >= mdry.
%
% INPUTS:  sol, ctrl, P as usual; opts.Nrun [def 200], opts.sig 1-sigma
%          struct (see plan Task 8)
% OUTPUTS: mc - .land(Nx2) .vtd(Nx1) .mprop(Nx1) .ok(Nx1) .success_rate
if nargin < 4, opts = struct(); end
if ~isfield(opts,'Nrun'), opts.Nrun = 200; end
sig = struct('r0',[100;100;50], 'v0',[10;10;10], 'thrust',0.015, ...
             'isp',0.01, 'wind',[10;10;0]);
if isfield(opts,'sig')
    fn = fieldnames(opts.sig);
    for k = 1:numel(fn), sig.(fn{k}) = opts.sig.(fn{k}); end
end

rng(P.seed);
Nr = opts.Nrun;
mc.land = zeros(Nr,2);  mc.vtd = zeros(Nr,1);
mc.mprop = zeros(Nr,1); mc.ok  = false(Nr,1);
for krun = 1:Nr
    d = struct('dr0', sig.r0 .* randn(3,1), ...
               'dv0', sig.v0 .* randn(3,1), ...
               'thrust_scale', 1 + sig.thrust*randn(), ...
               'isp_scale',    1 + sig.isp*randn());
    if P.drag.on, d.wind = sig.wind .* randn(3,1); end
    out = sim_closed_loop(sol, ctrl, P, d);
    mc.land(krun,:) = out.td.r(1:2).';
    mc.vtd(krun)    = out.td.vtd;
    mc.mprop(krun)  = out.td.m - P.mdry;
    mc.ok(krun)     = out.td.miss < P.pad_radius && ...
                      out.td.vtd < P.vtd_max && out.td.m >= P.mdry;
    if mod(krun, 25) == 0, fprintf('  MC %d/%d\n', krun, Nr); end
end
mc.success_rate = mean(mc.ok);
end
```

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: PASS. Note the small-N success rate; the spec's ≥95% criterion is judged at Nrun=200 in the front door, not here.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/lib/run_monte_carlo.m booster_landing/tests/test_monte_carlo_small.m
git commit -m "booster_landing: dispersed Monte-Carlo landing campaign

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 9: Visualization (`plot_pdg_solution`, `plot_footprint`, `movie_landing`)

**Files:**
- Create: `viz/plot_pdg_solution.m`
- Create: `viz/plot_footprint.m`
- Create: `viz/movie_landing.m`
- Create: `tests/test_viz_smoke.m`

**Interfaces:**
- Consumes: `sol` (either solver), `mc`, `out` (closed-loop trace).
- Produces:
  - `fig = plot_pdg_solution(solC, solV, outfile)` — 2×2: (a) 3D trajectory + glideslope cone + thrust arrows every 5th node, both solutions overlaid; (b) throttle `‖T‖/Tmax` vs t (the max–min–max money plot), both solvers; (c) mass vs t; (d) speed vs t. `exportgraphics` to `outfile` PNG at 200 dpi if given.
  - `fig = plot_footprint(mc, P, outfile)` — landing scatter colored by success, pad-radius circle, 3σ ellipse (from `cov(mc.land)`), annotation with success rate and vtd stats.
  - `movie_landing(out, sol, P, outfile)` — house polished-graphics MP4: left panel booster altitude view (marker + thrust-vector arrow scaled by throttle, ground line, pad), right panel throttle trace with moving cursor; **frame size forced to 1280×720** (÷16 rule — the diagonal-streak H.264 lesson), 30 fps, `VideoWriter` MPEG-4. Loop over a fixed 12 s duration resampled from `out.t`.
- All three follow matlab-polished-graphics: fixed axis limits computed once (no autoscale jitter), title fixed-width, no emoji, colorblind-safe two-color scheme (`[0 0.447 0.741]` colloc / `[0.85 0.325 0.098]` convex).

- [ ] **Step 1: Write the failing smoke test**

`tests/test_viz_smoke.m`:

```matlab
% TEST_VIZ_SMOKE  Viz functions execute headless and write files.
% Not a beauty contest -- existence + nonzero size only. Movie smoke uses
% 2 seconds of frames to stay fast.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
solC = solve_pdg_colloc(P, struct('N', 30));
solV = solve_pdg_convex(P, struct('tf', solC.tf, 'Nconv', 60));
ctrl = tvlqr_design(solC, P);
out  = sim_closed_loop(solC, ctrl, P, struct());
mc   = run_monte_carlo(solC, ctrl, P, struct('Nrun', 8));

od = fullfile(tempdir, 'bl_viz_smoke');  if ~exist(od,'dir'), mkdir(od); end
plot_pdg_solution(solC, solV, fullfile(od, 'sol.png'));
plot_footprint(mc, P, fullfile(od, 'fp.png'));
movie_landing(out, solC, P, fullfile(od, 'mov.mp4'), struct('duration', 2));
for f = {'sol.png','fp.png','mov.mp4'}
    d = dir(fullfile(od, f{1}));
    assert(~isempty(d) && d.bytes > 1e3, 'missing/empty %s', f{1});
end
fprintf('test_viz_smoke PASS\n');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_viz_smoke"`
Expected: FAIL — functions undefined.

- [ ] **Step 3: Implement the three viz functions**

Follow the interface bullets above exactly; skeleton for the movie (the layout/streak-proofing parts that must not be improvised):

```matlab
function movie_landing(out, sol, P, outfile, opts)
% MOVIE_LANDING  Landing movie: booster + thrust vector | throttle trace.
% Frame size locked to 1280x720 (divisible by 16 -- H.264 shear guard).
% INPUTS: out (sim_closed_loop), sol (guidance), P, outfile .mp4,
%         opts.duration [s of playback, def 12], opts.fps [def 30]
if nargin < 5, opts = struct(); end
if ~isfield(opts,'duration'), opts.duration = 12; end
if ~isfield(opts,'fps'),      opts.fps = 30;      end

fig = figure('Visible','off','Position',[100 100 1280 720],'Color','w');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
axT = nexttile(tl);  axR = nexttile(tl);

%% Fixed limits computed ONCE (no autoscale jitter):
xmax = max([abs(out.X(:,1)); abs(out.X(:,2)); 50]);
zmax = max(out.X(:,3)) * 1.05;
Tmag = sqrt(sum(out.Tcmd.^2, 2));

nf = round(opts.duration * opts.fps);
ts = interp1(linspace(0,1,numel(out.t)), out.t, linspace(0,1,nf));
vw = VideoWriter(outfile, 'MPEG-4');  vw.FrameRate = opts.fps;  open(vw);
for k = 1:nf
    xk = interp1(out.t, out.X, ts(k)).';
    Tk = interp1(out.t, out.Tcmd, ts(k)).';
    % left: altitude view (downrange = sqrt(x^2+y^2) signed by x), booster
    % marker, thrust arrow DOWN-scaled by throttle, pad + ground line
    cla(axT);  hold(axT, 'on');
    plot(axT, sqrt(sum(out.X(:,1:2).^2,2)).*sign(out.X(:,1)+eps), ...
         out.X(:,3), '-', 'Color', [0.7 0.7 0.7]);
    dk = sqrt(sum(xk(1:2).^2)) * sign(xk(1)+eps);
    arr = -0.15 * zmax * Tk / P.Tmax;    % plume points opposite thrust
    quiver(axT, dk, xk(3), arr(1), arr(3), 0, 'r', 'LineWidth', 2, ...
           'MaxHeadSize', 0.5);
    plot(axT, dk, xk(3), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 9);
    plot(axT, [-P.pad_radius P.pad_radius], [0 0], 'g-', 'LineWidth', 4);
    xlim(axT, [-xmax xmax]*1.1);  ylim(axT, [-0.02*zmax zmax]);
    title(axT, sprintf('t = %6.2f s   alt = %7.1f m', ts(k), xk(3)));
    xlabel(axT, 'downrange [m]');  ylabel(axT, 'altitude [m]');
    % right: throttle trace + cursor
    cla(axR);  hold(axR, 'on');
    plot(axR, out.t, Tmag/P.Tmax, 'b-');
    yline(axR, P.Tmin/P.Tmax, 'k--');  yline(axR, 1, 'k--');
    xline(axR, ts(k), 'r-');
    ylim(axR, [0 1.1]);  xlim(axR, [0 out.t(end)]);
    xlabel(axR, 'time [s]');  ylabel(axR, 'throttle T/Tmax');
    title(axR, 'throttle');
    frame = getframe(fig);
    frame.cdata = imresize(frame.cdata, [720 1280]);   % divide-by-16 lock
    writeVideo(vw, frame);
end
close(vw);  close(fig);
end
```

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: `test_viz_smoke PASS`. Then open the smoke PNGs once (`open /tmp/.../sol.png` via the tempdir path printed) and eyeball: throttle plot must visibly show max–min–max.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/viz/plot_pdg_solution.m booster_landing/viz/plot_footprint.m \
        booster_landing/viz/movie_landing.m booster_landing/tests/test_viz_smoke.m
git commit -m "booster_landing: solution/footprint plots + landing movie (1280x720 lock)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 10: Front door (`run_booster_landing`) + flagship run

**Files:**
- Create: `run_booster_landing.m`
- Create: `tests/test_run_front_door.m`
- Create: `README.md`

**Interfaces:**
- Consumes: everything above.
- Produces: `R = run_booster_landing(cfg)` — the campaign. No args = full nominal Phase-1 campaign. `cfg` fields (all optional): `.P` (params override struct — merged onto `booster_params()`), `.doMovie` (default true), `.doMC` (default true), `.Nrun` (default 200), `.outdir` (default `results/`). Sequence, with `fprintf` stage banners:
  1. `setup_paths` (self — front door must work from a fresh MATLAB after `cd booster_landing`)
  2. solve colloc (N=60) → solve convex (golden tf) → `certify_pdg` → **print gate table**
  3. `tvlqr_design` → nominal `sim_closed_loop` → `run_monte_carlo`
  4. viz: `pdg_solution.png`, `footprint.png`, `landing.mp4` into outdir
  5. save `R` (all products: solC, solV, rep, ctrl, out0, mc, timestamps via `datetime`) to `outdir/booster_run.mat`
  6. final summary block: tf, mf, fuel used, gates verdict, MC success rate — the numbers a reader quotes.
- README.md: what the campaign is, how to run it (`matlab -batch "cd ...; run_booster_landing"`), what each folder holds, the spec/plan pointers.

- [ ] **Step 1: Write the failing test**

`tests/test_run_front_door.m`:

```matlab
% TEST_RUN_FRONT_DOOR  Front-door contract on a FAST config: runs end to
% end from only setup_paths, honors cfg overrides (the advertised-but-
% ignored-options bug is the canonical failure -- verify N actually
% reached the solver via solution size), writes its products.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

od  = fullfile(tempdir, 'bl_front_smoke');
cfg = struct('doMovie', false, 'Nrun', 6, 'outdir', od, ...
             'P', struct('N', 20, 'Nconv', 60));
R   = run_booster_landing(cfg);
assert(size(R.solC.X, 2) == 21, 'cfg.P.N did not reach the solver');
assert(R.rep.all_pass || ischar(R.rep.G3_pass), 'gates failed on fast config');
assert(isfile(fullfile(od, 'booster_run.mat')), 'products not written');
assert(isfile(fullfile(od, 'pdg_solution.png')), 'solution plot missing');
fprintf('test_run_front_door PASS\n');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_run_front_door"`
Expected: FAIL — `run_booster_landing` undefined.

- [ ] **Step 3: Implement `run_booster_landing.m` + README**

Front-door skeleton (ADJUSTABLE PARAMETERS block up top, per the standing principles; params override by field-merge):

```matlab
function R = run_booster_landing(cfg)
% RUN_BOOSTER_LANDING  Front door: full 3-DOF booster-landing campaign.
%
%   solve (colloc + convex) -> certify G1-G5 -> TVLQR -> closed loop ->
%   Monte Carlo -> plots + movie -> results/booster_run.mat
%
% Run me with no arguments for the nominal campaign:
%   /Applications/MATLAB_R2025b.app/bin/matlab -batch ...
%     "cd('~/Desktop/optimal_control/booster_landing'); run_booster_landing"
%
% INPUTS:
%   cfg - (optional) .P params overrides, .doMovie [true], .doMC [true],
%         .Nrun [200], .outdir ['results/']
% OUTPUTS:
%   R   - everything: .P .solC .solV .rep .ctrl .out0 .mc .when
%
% REFERENCES: spec at docs/superpowers/specs/2026-08-08-booster-landing-design.md
setup_paths;

%% ---------------- ADJUSTABLE PARAMETERS (defaults) ----------------
def = struct('doMovie', true, 'doMC', true, 'Nrun', 200, ...
             'outdir', fullfile(fileparts(mfilename('fullpath')), 'results'));
%% ------------------------------------------------------------------
if nargin < 1, cfg = struct(); end
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(cfg, fn{k}), cfg.(fn{k}) = def.(fn{k}); end
end
P = booster_params();
if isfield(cfg, 'P')
    pf = fieldnames(cfg.P);
    for k = 1:numel(pf), P.(pf{k}) = cfg.P.(pf{k}); end
end
if ~exist(cfg.outdir, 'dir'), mkdir(cfg.outdir); end

fprintf('=== [1/5] Guidance: collocation (N=%d) ===\n', P.N);
R.solC = solve_pdg_colloc(P);
fprintf('    tf=%.3f s  mf=%.2f kg\n', R.solC.tf, R.solC.mf);

fprintf('=== [2/5] Guidance: lossless convexification (golden tf) ===\n');
R.solV = solve_pdg_convex(P);
fprintf('    tf=%.3f s  mf=%.2f kg  gap=%.2e\n', ...
        R.solV.tf, R.solV.mf, R.solV.lossless_gap);

fprintf('=== [3/5] Certification ===\n');
R.rep = certify_pdg(R.solC, R.solV, P);
print_certify_report(R.rep);

fprintf('=== [4/5] Tracking: TVLQR + closed loop%s ===\n', ...
        ternary(cfg.doMC, ' + Monte Carlo', ''));
R.ctrl = tvlqr_design(R.solC, P);
R.out0 = sim_closed_loop(R.solC, R.ctrl, P, struct());
if cfg.doMC
    R.mc = run_monte_carlo(R.solC, R.ctrl, P, struct('Nrun', cfg.Nrun));
end

fprintf('=== [5/5] Products -> %s ===\n', cfg.outdir);
plot_pdg_solution(R.solC, R.solV, fullfile(cfg.outdir, 'pdg_solution.png'));
if cfg.doMC
    plot_footprint(R.mc, P, fullfile(cfg.outdir, 'footprint.png'));
end
if cfg.doMovie
    movie_landing(R.out0, R.solC, P, fullfile(cfg.outdir, 'landing.mp4'));
end
R.P = P;  R.when = datetime('now');
save(fullfile(cfg.outdir, 'booster_run.mat'), '-struct', 'R');

fprintf('\n==================== SUMMARY ====================\n');
fprintf('tf        %.3f s      fuel  %.1f kg (mf %.1f)\n', ...
        R.solC.tf, P.m0 - R.solC.mf, R.solC.mf);
fprintf('gates     %s\n', ternary(R.rep.all_pass, 'ALL PASS', 'FAILURES -- see table'));
fprintf('nom miss  %.2f m  vtd %.2f m/s\n', R.out0.td.miss, R.out0.td.vtd);
if cfg.doMC
    fprintf('MC        %d runs, success %.1f%%\n', cfg.Nrun, 100*R.mc.success_rate);
end
fprintf('=================================================\n');
end

function s = ternary(c, a, b), if c, s = a; else, s = b; end, end
```

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: `test_run_front_door PASS`.

- [ ] **Step 5: FLAGSHIP RUN — the real acceptance test**

Run (expect a few minutes): `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); run_booster_landing"`
Expected, per spec success criteria — verify each, honestly, before claiming done:
- both solvers converged; gate table ALL PASS; final masses agree < 1.0 kg (measured ~0.7 kg Taylor model offset — adjudicated, documented)
- throttle plot shows bang-bang min–max with terminal max-throttle arc (adjudicated: min–max IS the optimum here; open `results/pdg_solution.png` and LOOK)
- MC success ≥ 95% at 200 runs
- `results/landing.mp4` plays clean (open it — no diagonal streaks)
Record tf / mf / fuel / success-rate in the commit message. If any criterion fails, stop and fix (systematic-debugging) — do not commit a red flagship.

- [ ] **Step 6: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/run_booster_landing.m booster_landing/tests/test_run_front_door.m \
        booster_landing/README.md
git commit -m "booster_landing: front door + flagship run

Flagship: tf=<fill> s, fuel=<fill> kg, gates ALL PASS, MC <fill>%/200.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 11: Phase 2 — atmosphere + drag

**Files:**
- Modify: `run_booster_landing.m` (add `cfg.phase2` switch, default false)
- Create: `tests/test_phase2_drag.m`
- Create: `viz/plot_vacuum_vs_drag.m`

**Interfaces:**
- Consumes: everything; `pdg_dynamics` drag branch already exists and is Jacobian-tested (Task 2).
- Produces: `R = run_booster_landing(struct('phase2', true))` additionally runs:
  - drag-on collocation solve (`P.drag.on = true`), **warm-started** from the vacuum solution (`opts.init = R.solC`)
  - `certify_pdg(solD, [], P)` — G1/G2/G5 gates only (no convex twin; G3/G4 'skipped')
  - drag-aware TVLQR + closed loop + MC **with wind** (sig.wind active because `P.drag.on`)
  - `plot_vacuum_vs_drag(solC, solD, outfile)` — overlay: trajectory, throttle, mass; annotation box with Δfuel and Δtf (the "what drag-free guidance misses" number, headed for the note)
  - products: `results/phase2_*.png`, extended `booster_run.mat` fields `.solD .repD .mcD`

- [ ] **Step 1: Write the failing test**

`tests/test_phase2_drag.m`:

```matlab
% TEST_PHASE2_DRAG  Drag-on collocation solve (coarse, warm-started from
% vacuum) converges; fuel differs from vacuum by a NONZERO but sane amount
% (drag helps braking: expect LESS fuel, sanity band 0 < dfuel < 40% of
% vacuum fuel); certify gates G1/G2/G5 pass with G3/G4 skipped.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
solC = solve_pdg_colloc(P, struct('N', 30));
Pd   = P;  Pd.drag.on = true;
solD = solve_pdg_colloc(Pd, struct('N', 30, 'init', solC));
assert(solD.stats.success, 'drag solve failed');
fuelV = P.m0 - solC.mf;  fuelD = P.m0 - solD.mf;
assert(fuelD < fuelV, 'drag should reduce fuel (braking), got +%.1f kg', ...
       fuelD - fuelV);
assert(fuelV - fuelD < 0.4*fuelV, 'drag effect implausibly large');
rep = certify_pdg(solD, [], Pd);
assert(rep.G1_pass && rep.G2_pass && rep.G5_pass, 'drag gates failed');
assert(isequal(rep.G3_pass, 'skipped'), 'G3 should be skipped without twin');
fprintf('test_phase2_drag PASS  fuel vac=%.1f drag=%.1f kg\n', fuelV, fuelD);
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_phase2_drag"`
Expected: FAIL only if machinery breaks — Tasks 2/3/5 already built the branches. If it PASSES immediately, good: that's the opt-in-flag pattern paying off; proceed to Step 3.
(One likely genuine failure: G5 with drag — the PMP thrust-direction law changes when drag is present; if primer alignment fails ONLY in the drag case, relax G5 for drag runs to the bang-bang structure check alone, document why in `certify_pdg`'s header, and carry the primer question to the note.)

- [ ] **Step 3: Implement the front-door `phase2` branch + `plot_vacuum_vs_drag`**

In `run_booster_landing`, after stage 5, when `cfg.phase2`:

```matlab
if isfield(cfg,'phase2') && cfg.phase2
    fprintf('=== [P2] Drag-on re-solve (warm-started) ===\n');
    Pd = P;  Pd.drag.on = true;
    R.solD = solve_pdg_colloc(Pd, struct('init', R.solC));
    R.repD = certify_pdg(R.solD, [], Pd);
    print_certify_report(R.repD);
    ctrlD  = tvlqr_design(R.solD, Pd);
    R.mcD  = run_monte_carlo(R.solD, ctrlD, Pd, struct('Nrun', cfg.Nrun));
    plot_vacuum_vs_drag(R.solC, R.solD, fullfile(cfg.outdir, 'phase2_vac_vs_drag.png'));
    plot_footprint(R.mcD, Pd, fullfile(cfg.outdir, 'phase2_footprint.png'));
    fprintf('P2: fuel vac %.1f kg -> drag %.1f kg; MC(wind) %.1f%%\n', ...
            P.m0 - R.solC.mf, Pd.m0 - R.solD.mf, 100*R.mcD.success_rate);
end
```

`plot_vacuum_vs_drag`: 1×3 tiledlayout (trajectory downrange–altitude overlay; throttle overlay; mass overlay), Δfuel/Δtf annotation, house colors, exportgraphics.

- [ ] **Step 4: Run tests + Phase-2 flagship**

Run: the Step-2 test command (expect PASS), then
`/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); run_booster_landing(struct('phase2', true))"`
Expected: Phase 1 unchanged (gates still ALL PASS — this is the re-run-the-flagship-after-structural-change rule), Phase 2 converged, wind-MC success reported (≥95% not required by spec for Phase 2 — report the honest number), comparison PNG written.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/run_booster_landing.m booster_landing/tests/test_phase2_drag.m \
        booster_landing/viz/plot_vacuum_vs_drag.m
git commit -m "booster_landing: Phase 2 atmosphere -- drag re-solve, wind MC, vacuum comparison

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 12: Theory note (`doc/booster_landing_note.tex`)

**Files:**
- Create: `doc/booster_landing_note.tex`

**Interfaces:**
- Consumes: `results/` figures + the flagship numbers (real, from the runs — never invented).
- Produces: compiled PDF. Sections:
  1. **Problem** — 3-DOF PDG statement, F9-class parameters table, the hoverslam T/W argument
  2. **Two routes** — HS collocation NLP; lossless convexification (change of variables, the relaxation, statement of the Açıkmeşe–Ploen tightness theorem with proof sketch via the maximum principle; Taylor mass-bound construction)
  3. **Certification** — the five gates, table of flagship values (from `booster_run.mat`)
  4. **Tracking** — TVLQR derivation sketch, saturation-near-terminal-arc discussion (expected physics, shown honestly)
  5. **Monte Carlo** — dispersions table, footprint figure, success rate
  6. **Phase 2** — drag model, Δfuel result, what drag-free guidance misses
  7. **Future work** — 6-DOF, full return profile, SCvx, real conic solver
- References: Açıkmeşe & Ploen 2007; Blackmore, Açıkmeşe & Scharf 2010; Kelly 2017; Anderson & Moore. All real — no invented citations.

- [ ] **Step 1: Write the note**

Standard article class, house LaTeX conventions. Pull every quoted number from `results/booster_run.mat` (load it, read the fields) — if a number in the draft can't be traced to the .mat or a test output, delete it.

- [ ] **Step 2: Compile + clean aux**

Run: `cd /Users/msc/Desktop/optimal_control/booster_landing/doc && /Library/TeX/texbin/pdflatex -interaction=nonstopmode booster_landing_note.tex && /Library/TeX/texbin/pdflatex -interaction=nonstopmode booster_landing_note.tex && rm -f *.aux *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz`
Expected: clean compile, no undefined references.

- [ ] **Step 3: Verify figures/citations**

Run: `~/ai_council/venv/bin/python ~/Documents/myLatex/tools/paper_verify.py /Users/msc/Desktop/optimal_control/booster_landing/doc/booster_landing_note.tex --deep-figs`
Expected: green — no FAIL on cites/figures (house rule: mandatory before calling doc work done).

- [ ] **Step 4: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/doc/booster_landing_note.tex
git commit -m "booster_landing: theory note (PDG, convexification, gates, TVLQR, MC)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

(Then offer Mike the doc-review pass on the note — decided at spec time to review the note, not the spec.)

---

## Execution order & checkpoints

Tasks are sequential (each consumes the previous interfaces). Natural review checkpoints for Mike: after Task 5 (the scientific core — two methods, one answer), after Task 10 (flagship), after Task 12 (note).

## Self-review notes (done at write time)

- Spec coverage: every spec section maps to a task (problem/params→1–2, solvers→3–4, certify→5, tracking→6–7, MC→8, viz→9, front door→10, Phase 2→11, note→12, tests woven throughout, pumpkyn conventions in Global Constraints).
- Type consistency: `sol.{t,tf,mf,X,U}` identical meaning across both solvers (convex converts back to thrust-in-N); `ctrl.{tgrid,K,xnom,Tnom}` consumed by sim exactly as produced; `P` merge pattern uniform.
- Known judgment calls flagged inline: G5 primer sign (flip once, don't abs), coarse-grid G3 tolerance scaling, golden-section unimodality assumption, drag-case primer law.
