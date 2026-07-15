# GTO→ELFO min-time anchor (Route B) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Solve the GTO→ELFO hard all-burn (`s≡1`) minimum-time transfer to land a machine-tight, independently-verified `tfMin_ELFO` anchor.

**Architecture:** A new CasADi+IPOPT solver `casadi_mintime_freetf.m` mirrors the energy sibling `casadi_energy_freetf.m` (same two-primary Sundman clock + `cScale` free-`t_f` slack) but pins throttle `s≡1`, drops throttle from the control, and minimizes `t(τ_f)` directly — the objective *is* min-time. A driver `gen_elfo_mintime.m` warm-starts it from the lowest converged energy rung, saves the solution with an `s≡1` throttle row for drop-in compatibility, and re-verifies it solver-free via the existing `verify_elfo_seed.m`.

**Tech Stack:** MATLAB R2025b, CasADi 3.7.0 (`~/casadi-3.7.0`, prebuilt arm64), IPOPT/MUMPS.

## Global Constraints

- **MATLAB is R2025b only**: `/Applications/MATLAB_R2025b.app/bin/matlab` (R2025a license is broken).
- **Run headless**: `matlab -batch "cd('<elfo>'); setup_paths; <cmd>"` from the `elfo/` dir.
- **CasADi path**: solvers `addpath(getenv('HOME')/casadi-3.7.0)` via the `CASADI_PATH` env fallback already coded in the sibling.
- **elfo dir**: `/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo`.
- **Seed (fixed)**: `results/energy_elfo_f0990.mat` — the Route A floor, `tf` 6.2278 ND (27.61 d), edge 57.5%. Contains fields `X [9×nN], U [4×nN], sigma, rv0, rvf, tauf0, tf, moonZone, pSund, qSund`.
- **Every MATLAB function needs the commented header block** (purpose, inputs w/ sizes, outputs w/ sizes, references).
- **Never use `i`/`j` as loop indices** (imaginary unit).
- **Result `.mat` and `results/logs/` are gitignored** — do not `git add` them.
- **Acceptance for the anchor** (Task 2): `Solve_Succeeded`, `maxDefect < 1e-8` (target 1e-10+), rendezvous `< 1e-8`, `maxUnit < 1e-8`, `minR1 ≥ 0.95·r_perigee_GTO`, `tMonotone` true, and sanity band `tfMin_ELFO < 6.2278 ND` (≈ 5.5 ND expected). If `s≡1` cannot close the rendezvous at any `t_f`, STOP and report "ELFO min-time is not all-burn" (fallback out of scope).

---

### Task 1: `casadi_mintime_freetf.m` solver + smoke test

**Files:**
- Create: `NLP_lowThrust_GTO_tulip/elfo/casadi_mintime_freetf.m`
- Test: `NLP_lowThrust_GTO_tulip/elfo/test_mintime_freetf.m`

**Interfaces:**
- Consumes: the energy seed `results/energy_elfo_f0990.mat`; params from `minfuel_config()` + `cr3bp_lt_params()` (`p.Tmax, p.c, p.muStar, p.tStar`).
- Produces: `out = casadi_mintime_freetf(sigma, rv0, rvf, Tmax, cEx, muStar, X0, U0, tauf0, opts)` returning struct with `.X [9×nN] .U [3×nN] .tf(=tfMin) .cScale .mf .maxDefect .maxUnit .minR1 .tMonotone .lamDef [9×N] .lamAll .primerAlignDeg .lamMassEnd .success .ipoptStatus`. `opts` fields: `pSund[1.5] qSund[4] moonZone[0.15] cBox[[0.10 8]] c0[1] tfCapMult[4] maxIter[3000] warmTight[false]`.

- [ ] **Step 1: Write the failing smoke test**

Create `NLP_lowThrust_GTO_tulip/elfo/test_mintime_freetf.m`:

```matlab
% TEST_MINTIME_FREETF  Smoke test: casadi_mintime_freetf constructs, runs, and
% returns the all-burn (s==1) min-time LAYOUT on the ELFO energy seed. These
% assertions are convergence-INDEPENDENT (they hold after a short 60-iter solve).
% The numeric all-burn mass identity mf=1-(Tmax/c)*tf holds only once the defects
% are converged, so it is checked at the converged anchor in Task 2, not here.
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
S = load(fullfile(here,'results','energy_elfo_f0990.mat'));
o = struct('pSund',S.pSund,'qSund',S.qSund,'moonZone',S.moonZone, ...
           'maxIter',60,'warmTight',false);
out = casadi_mintime_freetf(S.sigma,S.rv0,S.rvf,p.Tmax,p.c,p.muStar, ...
                            S.X,S.U(1:3,:),S.tauf0,o);
nN = size(S.X,2);
assert(isequal(size(out.X),[9 nN]), 'X layout wrong: %dx%d', size(out.X,1),size(out.X,2));
% THE wiring proof: control is 3-row (throttle is NOT a decision variable => s==1
% is structural, thrust/mdot cannot carry a hidden s factor).
assert(isequal(size(out.U),[3 nN]), 'U must be 3-row steering, got %dx%d', size(out.U,1),size(out.U,2));
assert(all(isfield(out,{'tf','minR1','tMonotone','cScale','mf'})), 'missing min-time fields');
assert(isfinite(out.tf) && out.tf > 0, 'tf not a positive finite number: %g', out.tf);
assert(out.maxUnit < 1e-2, 'unit-steering not tracking: %.2e', out.maxUnit);
% informational (NOT asserted -- exact only at convergence, Task 2):
mf_pred = 1 - (p.Tmax/p.c)*out.tf;
fprintf('TEST_MINTIME_FREETF: PASS (tf=%.4f ND, maxUnit=%.2e; mf=%.4f vs all-burn pred %.4f)\n', ...
        out.tf, out.maxUnit, out.mf, mf_pred);
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo'); setup_paths; test_mintime_freetf"
```
Expected: FAIL — `Unrecognized function or variable 'casadi_mintime_freetf'`.

- [ ] **Step 3: Write the solver**

Create `NLP_lowThrust_GTO_tulip/elfo/casadi_mintime_freetf.m`:

```matlab
function out = casadi_mintime_freetf(sigma, rv0, rvf, Tmax, cEx, muStar, X0, U0, tauf0, opts)
% CASADI_MINTIME_FREETF  Hard all-burn (s==1) free-final-time min-TIME Sundman
% collocation with a two-primary clock (CasADi+IPOPT).
%
% Sibling of CASADI_ENERGY_FREETF: identical two-primary Sundman clock and cScale
% free-t_f slack, but throttle is HARD-PINNED s==1 (control is steering alpha
% only, ||alpha||=1) and the objective is min t(tau_f). No epsilon, no tfTarget,
% no throttle integrand -- the objective IS min-time, so IPOPT drives t_f down
% until the rendezvous BC can no longer be met all-burn. That floor is tfMin.
%
% State  x = [r(3); v(3); m; t; cScale]  (9).  Control u = alpha (3), ||alpha||=1.
% Always-burn: thrust accel = (Tmax/m)*alpha,  mdot = -Tmax/cEx.
%
% INPUTS:
%   sigma   - normalized independent-variable nodes [(N+1)x1], 0 -> 1
%   rv0,rvf - initial / target position-velocity (ND) [1x6]
%   Tmax    - max thrust (ND) [scalar]
%   cEx     - exhaust velocity (ND) [scalar]
%   muStar  - CR3BP mass ratio [scalar]
%   X0      - warm-start states [8x(N+1)] ([r;v;m;t]) or [9x(N+1)] ([...;cScale])
%   U0      - warm-start control [3x(N+1)] (alpha) or [4x(N+1)] ([alpha;s]); a
%             4th throttle row is dropped (s==1 is enforced here)
%   tauf0   - fixed total regularized length [scalar]
%   opts    - (optional) struct: .pSund[1.5] .qSund[4] .moonZone[0.15]
%             .cBox[[0.10 8]] .c0[1] .tfCapMult[4] .maxIter[3000] .warmTight[false]
%
% OUTPUTS:
%   out - struct: .X [9x(N+1)] .U [3x(N+1)] .tauf .tf(=tfMin) .cScale .mf
%         .maxDefect .maxUnit .minR1 .tMonotone (logical) .lamDef [9xN] .lamAll
%         .primerAlignDeg .lamMassEnd .success .ipoptStatus
%
% REFERENCES:
%   [1] casadi_energy_freetf.m (the energy/fuel sibling this mirrors).
%   [2] Betts, "Practical Methods for Optimal Control," SIAM (2010) -- sparse
%       free-final-time via a constant slack state.

if nargin < 10 || isempty(opts), opts = struct(); end
gd = @(f,d) getdef(opts,f,d);
pSund    = gd('pSund',    1.5);
qSund    = gd('qSund',    4);
moonZone = gd('moonZone', 0.15);
cBox     = gd('cBox',     [0.10 8]);
c0       = gd('c0',       1);
tfCapMult= gd('tfCapMult',4);
maxIter  = gd('maxIter',  3000);
warmTight= gd('warmTight',false);

cpath = getenv('CASADI_PATH');
if isempty(cpath), cpath = fullfile(getenv('HOME'), 'casadi-3.7.0'); end
addpath(cpath);
import casadi.*

sigma = sigma(:);  N = numel(sigma) - 1;  nN = N + 1;
dsig  = diff(sigma).';

% --- warm start: 9-row state, 3-row control ---------------------------------
assert(size(X0,1) >= 8, 'X0 must have >=8 rows ([r;v;m;t]); got %d', size(X0,1));
if size(X0,1) == 8
    X0 = [X0; c0*ones(1, size(X0,2))];
end
assert(size(U0,1) >= 3, 'U0 must have >=3 rows (alpha); got %d', size(U0,1));
U0 = U0(1:3,:);                               % drop any throttle row (s==1)
tf_ws = X0(8,end);
tfCap = tfCapMult * max(tf_ws, eps);

% --- symbolic always-burn two-primary Sundman dynamics dX/dtau --------------
x = MX.sym('x', 9);  u = MX.sym('u', 3);
r = x(1:3);  v = x(4:6);  m = x(7);  cs = x(9);   al = u;
dd = [r(1)+muStar; r(2); r(3)];               % vector from Earth
rr = [r(1)-1+muStar; r(2); r(3)];             % vector from Moon
d2 = dd.'*dd + 1e-12;   e2 = rr.'*rr + 1e-12;
r1 = sqrt(d2);  r2 = sqrt(e2);
d3 = d2^1.5;    r3 = e2^1.5;
gr = [r(1); r(2); 0] - (1-muStar)*dd/d3 - muStar*rr/r3;   % full gravity (muGain=1)
hv = [2*v(2); -2*v(1); 0];
accel = gr + hv + (Tmax/m)*al;                % s == 1
mdot  = -(Tmax/cEx);                          % s == 1
if moonZone > 0
    kappa = ( r1^(-qSund) + (r2/moonZone)^(-qSund) )^(-pSund/qSund);
else
    kappa = r1^pSund;
end
fdyn  = Function('f', {x,u}, {[ cs*kappa*[v; accel; mdot; 1]; 0 ]});  % 9x1
Fmap  = fdyn.map(nN);

% --- NLP --------------------------------------------------------------------
opti = Opti();
X    = opti.variable(9, nN);
U    = opti.variable(3, nN);
tauf = tauf0;                                 % FIXED (t_f floats via cScale)
F    = Fmap(X, U);

% trapezoidal defects in sigma: dX/dsigma = tauf * dX/dtau
D = X(:,2:end) - X(:,1:end-1) - tauf*(repmat(dsig,9,1)/2).*(F(:,1:end-1) + F(:,2:end));
opti.subject_to(D(:) == 0);

% unit steering direction
opti.subject_to((sum(U(1:3,:).^2, 1) - 1).' == 0);

% bounds (explicit two-sided); cScale in its box, t in [0, tfCap]
lbX = repmat([-3;-3;-3;-12;-12;-12;0.3;0;    cBox(1)], 1, nN);
ubX = repmat([ 3; 3; 3; 12; 12; 12;1.0;tfCap;cBox(2)], 1, nN);
opti.subject_to(X(:) >= lbX(:));   opti.subject_to(X(:) <= ubX(:));
lbU = repmat([-1.1;-1.1;-1.1], 1, nN);
ubU = repmat([ 1.1; 1.1; 1.1], 1, nN);
opti.subject_to(U(:) >= lbU(:));   opti.subject_to(U(:) <= ubU(:));

% boundary conditions; cScale free
opti.subject_to(X(1:6,1) == rv0(:));   opti.subject_to(X(7,1) == 1);   opti.subject_to(X(8,1) == 0);
opti.subject_to(X(1:6,nN) == rvf(:));

% objective: MIN-TIME
opti.minimize(X(8,nN));
opti.set_initial(X, X0);
opti.set_initial(U, U0);

p = struct;
p.print_time      = true;
p.ipopt.max_iter  = maxIter;
p.ipopt.tol       = 1e-7;
p.ipopt.constr_viol_tol = 1e-7;
p.ipopt.acceptable_tol  = 1e-5;
p.ipopt.acceptable_iter = 15;
p.ipopt.nlp_scaling_method = 'gradient-based';
p.ipopt.linear_solver = 'mumps';
p.ipopt.print_level = 5;
if warmTight
    p.ipopt.mu_strategy                 = 'monotone';
    p.ipopt.warm_start_init_point       = 'yes';
    p.ipopt.mu_init                     = 1e-4;
    p.ipopt.warm_start_bound_push       = 1e-9;
    p.ipopt.warm_start_bound_frac       = 1e-9;
    p.ipopt.warm_start_slack_bound_push = 1e-9;
    p.ipopt.warm_start_slack_bound_frac = 1e-9;
    p.ipopt.warm_start_mult_bound_push  = 1e-9;
else
    p.ipopt.mu_strategy           = 'monotone';
    p.ipopt.warm_start_init_point = 'yes';
    p.ipopt.mu_init               = 0.1;
end
opti.solver('ipopt', p);

success = true;  status = 'solved';
try
    sol = opti.solve();
    Xs = sol.value(X);  Us = sol.value(U);
    lamAll = full(sol.value(opti.lam_g));
    status = char(opti.return_status());
catch solveErr
    Xs = opti.debug.value(X);  Us = opti.debug.value(U);
    try
        lamAll = full(opti.debug.value(opti.lam_g));
    catch
        lamAll = [];
    end
    success = false;  status = solveErr.message;
end

% metrics
Fs = full(Fmap(Xs, Us));
Dd = Xs(:,2:end) - Xs(:,1:end-1) - tauf*(repmat(dsig,9,1)/2).*(Fs(:,1:end-1) + Fs(:,2:end));

% --- KKT multipliers -> discrete costates + PMP primer (all burn) -----------
lamDef = [];  primerAlignDeg = NaN;  lamMassEnd = NaN;
if numel(lamAll) >= 9*N
    lamDef = reshape(lamAll(1:9*N), 9, N);
    lamV   = lamDef(4:6, :);
    primer = -lamV ./ max(sqrt(sum(lamV.^2,1)), 1e-12);
    aMid   = 0.5*(Us(1:3,1:end-1) + Us(1:3,2:end));
    aMid   = aMid ./ max(sqrt(sum(aMid.^2,1)), 1e-12);
    cang   = sum(primer.*aMid, 1);            % every interval is a burn (s==1)
    if mean(cang) < 0, cang = -cang; end
    primerAlignDeg = mean(acosd(min(max(cang,-1),1)));
    lamMassEnd = lamDef(7,end);
end

% physicality diagnostics
r1N   = sqrt((Xs(1,:)+muStar).^2 + Xs(2,:).^2 + Xs(3,:).^2);
minR1 = min(r1N);
tMonotone = all(diff(Xs(8,:)) > 0);

out = struct('X', Xs, 'U', Us, 'tauf', tauf, 'tf', Xs(8,end), ...
             'cScale', Xs(9,end), 'mf', Xs(7,end), ...
             'maxDefect', max(abs(Dd(:))), ...
             'maxUnit', max(abs(sum(Us(1:3,:).^2,1) - 1)), ...
             'minR1', minR1, 'tMonotone', tMonotone, ...
             'lamDef', lamDef, 'lamAll', lamAll, ...
             'primerAlignDeg', primerAlignDeg, 'lamMassEnd', lamMassEnd, ...
             'success', success, 'ipoptStatus', status);
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo'); setup_paths; test_mintime_freetf"
```
Expected: PASS — prints `TEST_MINTIME_FREETF: PASS (tf=..., maxUnit=...; mf=... vs all-burn pred ...)`. Structural assertions only (`X` 9-row, `U` **3-row** = `s≡1` wired, `tf` finite/positive, `maxUnit < 1e-2`); the `mf`-vs-prediction number is printed for eyeballing but not asserted (it converges only at the anchor, Task 2).

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add NLP_lowThrust_GTO_tulip/elfo/casadi_mintime_freetf.m NLP_lowThrust_GTO_tulip/elfo/test_mintime_freetf.m
git commit -m "feat(elfo min-time): casadi_mintime_freetf hard all-burn solver + smoke test"
```

---

### Task 2: `gen_elfo_mintime.m` driver + verified anchor run

**Files:**
- Create: `NLP_lowThrust_GTO_tulip/elfo/gen_elfo_mintime.m`
- Uses (unchanged): `NLP_lowThrust_GTO_tulip/elfo/verify_elfo_seed.m` (reads `SEEDFILE`, `U(4,:)` throttle).

**Interfaces:**
- Consumes: `casadi_mintime_freetf` (Task 1); seed `results/energy_elfo_f0990.mat`; `minfuel_config`, `cr3bp_lt_params`.
- Produces: `outFile = gen_elfo_mintime(opts)` → `results/mintime_elfo.mat` with `X [9×nN]`, `U [4×nN]` (steering + an `s≡1` row for `verify_elfo_seed`/movie/export compatibility), `sigma, rv0, rvf, tauf0, tf(=tfMin), mf, cScale, maxDefect, minR1, pSund, qSund, moonZone`. `opts`: `seedFile[results/energy_elfo_f0990.mat] maxIter[3000] warmTight[false]`.

- [ ] **Step 1: Write the driver**

Create `NLP_lowThrust_GTO_tulip/elfo/gen_elfo_mintime.m`:

```matlab
function outFile = gen_elfo_mintime(opts)
% GEN_ELFO_MINTIME  Solve the GTO->ELFO minimum-time (hard all-burn, s==1)
% transfer to anchor the front. Loads the lowest converged energy rung, overrides
% throttle to s==1, minimizes t(tau_f) via casadi_mintime_freetf, saves the
% all-burn min-time solution (with an s==1 throttle row for downstream compat)
% and prints the acceptance diagnostics. Independent solver-free verification is
% a separate step (verify_elfo_seed on the saved file).
%
% INPUTS:
%   opts - (optional) struct: .seedFile[results/energy_elfo_f0990.mat]
%          .maxIter[3000] .warmTight[false]
%
% OUTPUTS:
%   outFile - results/mintime_elfo.mat: X[9xnN], U[4xnN] (alpha + s==1 row),
%             sigma, rv0, rvf, tauf0, tf(=tfMin), mf, cScale, maxDefect, minR1,
%             pSund, qSund, moonZone
%
% REFERENCES:
%   [1] casadi_mintime_freetf.m; [2] the Route B design spec
%       (docs/superpowers/specs/2026-07-15-elfo-mintime-route-b-design.md).

if nargin < 1, opts = struct(); end
gd = @(f,d) getdef(opts,f,d);
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
resDir = fullfile(here,'results');
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);

seedFile = gd('seedFile', fullfile(resDir,'energy_elfo_f0990.mat'));
S = load(seedFile);
fprintf('=== GEN_ELFO_MINTIME: seed %s (tf=%.4f ND, %.2f d, edge from energy) ===\n', ...
        seedFile, S.X(8,end), S.X(8,end)*p.tStar/86400);

% warm start: states from the seed, throttle overridden to s==1 (drop U row 4)
X0 = S.X;  U0 = S.U(1:3,:);
o = struct('pSund',S.pSund,'qSund',S.qSund,'moonZone',S.moonZone, ...
           'cBox',[0.10 8],'tfCapMult',4,'maxIter',gd('maxIter',3000), ...
           'warmTight',gd('warmTight',false));
out = casadi_mintime_freetf(S.sigma, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
                            X0, U0, S.tauf0, o);

% acceptance diagnostics
rperi = norm(S.rv0(1:3) - [-p.muStar 0 0]);   % GTO perigee radius (ND)
rferr = norm(out.X(1:6,end) - S.rvf(:));
fprintf('  ipopt: %s   success=%d\n', out.ipoptStatus, out.success);
fprintf('  tfMin = %.4f ND (%.2f d)   mf=%.4f (prop %.1f%%)   cScale=%.4f\n', ...
        out.tf, out.tf*p.tStar/86400, out.mf, 100*(1-out.mf), out.cScale);
fprintf('  maxDefect=%.2e  maxUnit=%.2e  rendezvous=%.2e\n', out.maxDefect, out.maxUnit, rferr);
fprintf('  minR1=%.4f (GTO perigee=%.4f)  tMonotone=%d  primerAlign=%.3f deg\n', ...
        out.minR1, rperi, out.tMonotone, out.primerAlignDeg);

% save with a 4-row U (s==1 row) for drop-in compat with verify_elfo_seed / movie
X = out.X;  U = [out.U; ones(1,size(out.U,2))];  sigma = S.sigma; %#ok<NASGU>
rv0 = S.rv0;  rvf = S.rvf;  tauf0 = S.tauf0; %#ok<NASGU>
pSund = S.pSund;  qSund = S.qSund;  moonZone = S.moonZone; %#ok<NASGU>
tf = out.tf;  mf = out.mf;  cScale = out.cScale; %#ok<NASGU>
maxDefect = out.maxDefect;  minR1 = out.minR1; %#ok<NASGU>
outFile = fullfile(resDir,'mintime_elfo.mat');
save(outFile,'X','U','sigma','rv0','rvf','tauf0','tf','mf','cScale', ...
     'maxDefect','minR1','pSund','qSund','moonZone');
fprintf('  saved %s\n', outFile);

% relabel: the mapped front in ELFO's own units
fprintf('\n  RELABEL: factor_ELFO = tf / tfMin_ELFO,  tfMin_ELFO = %.4f ND (%.2f d)\n', ...
        out.tf, out.tf*p.tStar/86400);
fprintf('GEN_ELFO_MINTIME DONE\n');
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
```

- [ ] **Step 2: Run the driver (the anchor solve)**

Run (this is the real ~2–10 min solve; watch for a MEX crash and just re-run if so):
```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo'); setup_paths; gen_elfo_mintime;"
```
Expected: prints `ipopt: Solve_Succeeded`, a `tfMin` line with `tf < 6.2278 ND` (≈ 5.5 ND / 24–25 d), `maxDefect` ≤ 1e-8, `rendezvous` ≤ 1e-8, `minR1 ≥ 0.95·GTO perigee`, `tMonotone=1`, and `saved .../mintime_elfo.mat`.

- [ ] **Step 3: Independently verify (solver-free) via the existing verifier**

Run:
```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo'); setup_paths; SEEDFILE=fullfile(pwd,'results','mintime_elfo.mat'); verify_elfo_seed"
```
Expected: prints `||X(:,end)-rvf|| ≤ 1e-8`, `INDEP defect (pure MATLAB, full gravity) ≤ 1e-8`, `maxUnit ≤ 1e-8`, and `VERIFY: PASS`. (The `s≡1` row makes the verifier recompute at full thrust — an independent MATLAB re-derivation of the CasADi solution.)

- [ ] **Step 4: Check the acceptance gate**

Confirm ALL hold from Steps 2–3 output:
- `Solve_Succeeded`, `success=1`; `verify_elfo_seed` prints `VERIFY: PASS`.
- `maxDefect ≤ 1e-8`, rendezvous `≤ 1e-8`, `maxUnit ≤ 1e-8`.
- `minR1 ≥ 0.95·rperi`, `tMonotone=1`, and `cScale` strictly interior to `[0.10, 8]` (no unphysical dive; time monotone; slack not pinned to a bound).
- **Converged all-burn mass identity:** `|mf − (1 − (Tmax/c)·tfMin)| ≤ 1e-6` — now exact because the defects are converged (this is the numeric `s≡1` check deferred from Task 1). The driver already prints `mf` and prop%; compute the prediction as `1 - (p.Tmax/p.c)*tf` and compare.
- `tfMin_ELFO < 6.2278 ND` and in the sanity band (~5–6.2 ND).

If instead the solve does NOT close the rendezvous (restoration failure / `rendezvous` stuck large / `success=0`): STOP. This is the pre-agreed signal that **ELFO min-time is not all-burn**; record it in `elfo/ELFO_RETARGET.md` and the ROADMAP as the finding, and do not fake an anchor. (Fallback throttle-free min-time is a separate effort, out of scope.)

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add NLP_lowThrust_GTO_tulip/elfo/gen_elfo_mintime.m
git commit -m "feat(elfo min-time): gen_elfo_mintime driver -> verified tfMin_ELFO anchor"
```

---

## Notes for the implementer

- The two `matlab -batch` **solve** runs need CasADi + R2025b; if that environment is unavailable, say so plainly — do NOT hand-edit numbers into `mintime_elfo.mat` or claim a pass.
- The MEX FATAL crash (~1 in 10 solves) is uncatchable; if a run dies with a non-MATLAB fatal, just re-run Step 2 (the driver is a single solve, no resume state needed).
- Do not `git add` anything under `results/` (`*.mat` and `results/logs/` are gitignored).
- After the anchor lands, updating `ELFO_RETARGET.md` / `ROADMAP.md` / memory with `tfMin_ELFO` and re-plotting the front in ELFO units are follow-ups, not part of this plan.
