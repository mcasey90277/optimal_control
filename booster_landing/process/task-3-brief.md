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

