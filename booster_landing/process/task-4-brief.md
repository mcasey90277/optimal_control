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

