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
% ADAPTATION FROM BRIEF (documented, same lesson as Task 3): the brief's
% verbatim raw-SI R,V (position O(1e3) m, velocity O(1e2) m/s) mixed with
% z,u,sigma (already O(1)-O(30)) gave IPOPT unreliable convergence -- a
% sweep of fixed-tf solves returned Infeasible_Problem_Detected at some tf
% and "success" with a BLOWN lossless gap (up to ~500) at others, even
% though the problem is convex and feasible throughout. Nondimensionalizing
% only R,V (by Lc=norm(r0), Vc=Lc/tf) fixed this; z,u,sigma are left in SI
% per the task context note. Every field of `sol` is unscaled back to SI
% before return, so the interface below is unchanged.
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
if ~isfinite(mf_or_neginf(sc)) && ~isfinite(mf_or_neginf(sd))
    error('solve_pdg_convex:infeasibleBracket', ...
        ['Both golden-section probes (tf=%.3f, tf=%.3f) failed to solve. ' ...
         'Widen the [P.tf_lo, P.tf_hi] bracket.'], c, d);
end
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
al = 1 / (P.Isp * P.g0);
assert(P.m0 - al*P.Tmax*tf > 0, ...
    'solve_pdg_convex:massDepletion', ...
    'max-thrust depletion reference goes non-positive at tf=%.3f (m0=%.1f, al*Tmax*tf=%.1f) -- shrink tf or widen mdry margin', ...
    tf, P.m0, al*P.Tmax*tf);
opti = casadi.Opti();
t   = linspace(0, tf, Nc);  h = t(2) - t(1);
Lc  = norm(P.r0);  Vc = Lc / tf;      % position/velocity scales (see note above)
Rh  = opti.variable(3, Nc);   Vh = opti.variable(3, Nc);   % R/Lc, V/Vc
Z   = opti.variable(1, Nc);   Uu = opti.variable(3, Nc);   % SI (already O(1)-ish)
S   = opti.variable(1, Nc);   % sigma, SI

%% Linear dynamics, trapezoid (R,V scaled; Z,Uu,S in SI -- see ADAPTATION note):
for k = 1:Nc-1
    opti.subject_to(Rh(:,k+1) == Rh(:,k) + (h/(2*tf))*(Vh(:,k)+Vh(:,k+1)));
    opti.subject_to(Vh(:,k+1) == Vh(:,k) + (h/(2*Vc))*(2*P.gvec + Uu(:,k)+Uu(:,k+1)));
    opti.subject_to(Z(k+1)    == Z(k)    - (h/2)*al*(S(k)+S(k+1)));
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
    opti.subject_to(Rh(1,k)^2 + Rh(2,k)^2 <= (cotg*Rh(3,k))^2);   % homogeneous in Lc
    opti.subject_to(Rh(3,k) >= 0);
    if isfinite(P.theta_max_deg)
        opti.subject_to(Uu(3,k) >= cosd(P.theta_max_deg)*S(k));
    end
end

%% Boundary conditions, objective:
opti.subject_to(Rh(:,1) == P.r0/Lc);   opti.subject_to(Vh(:,1) == P.v0/Vc);
opti.subject_to(Z(1)   == log(P.m0));
opti.subject_to(Rh(:,end) == zeros(3,1));  opti.subject_to(Vh(:,end) == zeros(3,1));
opti.subject_to(Z(end) >= log(P.mdry));
opti.minimize(-Z(end));

%% Initial guess (DEVIATION from brief's verbatim code, documented in the
%% task report): Opti defaults all variables to 0, which puts Z (~0) far
%% below its zlb/zub bounds (~10) and S (~0) far below its mu1/mu2 bounds
%% (tens of m/s^2) -- IPOPT's restoration phase can't recover from that and
%% returns Infeasible_Problem_Detected even though the problem is convex
%% and feasible. Seed with the straight-line position/velocity guess and
%% the max-thrust depletion reference (z0, and S at its matching upper
%% bound mu2, thrust pointed straight up) -- consistent with the Taylor
%% bounds already computed above, cheap, and it converges reliably.
tauN = t/tf;
opti.set_initial(Rh, (P.r0/Lc)*(1-tauN));
opti.set_initial(Vh, (P.v0/Vc)*(1-tauN));
opti.set_initial(Z, z0);
Sg = P.Tmax*exp(-z0);
opti.set_initial(S, Sg);
opti.set_initial(Uu, [zeros(2,Nc); Sg]);

opti.solver('ipopt', struct('print_time', false), ...
            struct('max_iter', 3000, 'print_level', 0, 'tol', 1e-8));
try
    osol = opti.solve();  ok = true;
catch
    osol = opti.debug;    ok = false;
end

%% Package in ORIGINAL variables (T = m u), unscaling Rh,Vh back to SI:
sol.t  = t;  sol.tf = tf;
m      = exp(full(osol.value(Z)));
sol.u  = full(osol.value(Uu));  sol.sigma = full(osol.value(S));
Rsi    = full(osol.value(Rh)) * Lc;
Vsi    = full(osol.value(Vh)) * Vc;
sol.X  = [Rsi; Vsi; m];
sol.U  = sol.u .* m;
sol.mf = m(end);
sol.lossless_gap = max(abs(sqrt(sum(sol.u.^2,1)) - sol.sigma));
sol.stats = struct('success', ok, 'status', opti.stats.return_status);
sol.P  = P;
end
