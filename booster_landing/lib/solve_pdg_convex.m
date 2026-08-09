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
% Convex problem solved with IPOPT: any point IPOPT converges to is global
% BY CONVEXITY OF THE FEASIBLE SET AND OBJECTIVE, verified after the fact
% via certify_pdg's G4 (relaxation tightness) and G3 (agreement with the
% independently-formulated nonconvex collocation solve) -- not by any
% globality property of IPOPT itself, which is a local NLP method, and not
% by the SOC being handed to a cone solver as a cone (the annulus is
% written here in the squared form sum(u.^2) <= sigma^2). Wording sharpened
% 2026-08-09 after an external review flagged the bare phrase "convex
% global solver" as claiming more than this implementation establishes. A
% conic solver (ECOS/CVX) is a documented drop-in, deliberately not a
% dependency.
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
%          .tf_curve is Kx3: [tf, mf_or_-Inf, validity_code] where
%          validity_code 3=valid+Solve_Succeeded, 2=valid+
%          Solved_To_Acceptable_Level-but-tight (ADAPTATION, task-7 fix
%          report round 3: gating is tightness-first now, not
%          status-string-first -- see mf_or_neginf's own note for why),
%          1=converged but relaxation not tight, 0=solver failed (see
%          mf_or_neginf below) -- so a downstream reader can see what the
%          golden search rejected and why. v (column 2) is mf for BOTH
%          code 2 and code 3 (both valid), -Inf otherwise.
%
% REFERENCES:
%   [1] Acikmese & Ploen, JGCD 2007.  [2] Blackmore et al., JGCD 2010.
if nargin < 2, opts = struct(); end
if ~isfield(opts,'Nconv'), opts.Nconv = P.Nconv; end
if ~isfield(opts,'tolTf'), opts.tolTf = 0.05;    end

if isfield(opts,'tf') && ~isempty(opts.tf)
    sol = solve_fixed_tf(P, opts.tf, opts.Nconv);
    % FIX (final-review, 2026-08-09): tf_curve col 2 must be mf_or_neginf's
    % OWN value (v0, -Inf when the iterate isn't valid -- see that
    % function's header for the code/tightness semantics), matching the
    % golden-search path below and the docstring's own "v (column 2) is mf
    % for BOTH code 2 and code 3, -Inf otherwise" -- was sol.mf
    % unconditionally, which reported a raw mf even for an invalid
    % (code<2) single-tf probe.
    [v0, code0] = mf_or_neginf(sol);
    sol.tf_curve = [opts.tf, v0, code0];
    return
end

%% Golden-section search on tf, maximizing mf (see note for unimodality):
phi = (sqrt(5)-1)/2;
a = P.tf_lo;  b = P.tf_hi;  curve = [];
c = b - phi*(b-a);  d = a + phi*(b-a);
[sc, vc, codec] = solve_fixed_tf_probe(P, c, opts.Nconv);
[sd, vd, coded] = solve_fixed_tf_probe(P, d, opts.Nconv);
curve = [curve; c, vc, codec; d, vd, coded];
if ~isfinite(vc) && ~isfinite(vd)
    error('solve_pdg_convex:infeasibleBracket', ...
        ['Both golden-section probes (tf=%.3f, tf=%.3f) failed to solve, ' ...
         'even after a retry each. Widen the [P.tf_lo, P.tf_hi] bracket.'], c, d);
end
while (b - a) > opts.tolTf
    if vc > vd
        b = d;  d = c;  sd = sc;  vd = vc;
        c = b - phi*(b-a);
        [sc, vc, codec] = solve_fixed_tf_probe(P, c, opts.Nconv);
        curve = [curve; c, vc, codec];                      %#ok<AGROW>
    else
        a = c;  c = d;  sc = sd;  vc = vd;
        d = a + phi*(b-a);
        [sd, vd, coded] = solve_fixed_tf_probe(P, d, opts.Nconv);
        curve = [curve; d, vd, coded];                      %#ok<AGROW>
    end
end
if vc > vd, sol = sc; else, sol = sd; end
sol.tf_curve = sortrows(curve, 1);
end

function [sol, v, code] = solve_fixed_tf_probe(P, tf, Nconv)
% SOLVE_FIXED_TF_PROBE  One golden-section probe, retried up to 3 attempts
% total if it doesn't land on a VALID iterate (code>=2 -- see mf_or_neginf).
%
% ADAPTATION (task-7 fix report round 3, documented, discovered while
% re-solving after the terminal-velocity BC change below): the golden
% search's probes are IPOPT solves at specific, golden-ratio-derived tf
% values and had no fallback if one landed on a bad iterate -- observed
% for real: the same call (tf=25.2786, Nconv=120) measured a TIGHT gap
% (6.5e-5) in one MATLAB process and an UNTIGHT gap (0.448) in another,
% identical inputs, separate runs -- genuine run-to-run nondeterminism in
% the fixed-tf IPOPT solve (most likely BLAS/IPOPT internal threading,
% not a real feasibility difference; the true feasibility wall, confirmed
% by a dedicated sweep, is a reproducible Infeasible_Problem_Detected
% starting near tf~27 s, where required propellant for that much loiter
% time hits P.mdry). This retry handles THAT case.
%
% It did NOT, by itself, fix a second and more important problem found
% right after (see mf_or_neginf's own ADAPTATION note for the actual
% root cause and fix): the golden search's final answer sat 52 kg below
% a direct fixed-tf query at collocation's own tf, because ONE golden-
% section probe deterministically converged to
% 'Solved_To_Acceptable_Level' with an excellent, tight gap, and the
% classification at the time rejected any non-'Solve_Succeeded' status
% outright regardless of tightness -- reproducing bit-for-bit on every
% retry, since it was a deterministic IPOPT outcome at that exact tf, not
% noise. Retrying alone cannot fix a deterministic misclassification;
% mf_or_neginf's tightness-first reordering is the actual fix. This
% retry stays for the genuinely nondeterministic case it targets.
%
% INPUTS:
%   P     - booster_params
%   tf    - probe final time [scalar, s]
%   Nconv - trapezoid nodes [scalar]
% OUTPUTS:
%   sol  - the accepted solve_fixed_tf result (best code among attempts)
%   v    - mf if code>=2, else -Inf (see mf_or_neginf)
%   code - mf_or_neginf's code for the accepted attempt
sol = solve_fixed_tf(P, tf, Nconv);
[v, code] = mf_or_neginf(sol);
attempt = 1;
while code < 2 && attempt < 3
    solK = solve_fixed_tf(P, tf, Nconv);
    [vK, codeK] = mf_or_neginf(solK);
    if codeK > code
        sol = solK;  v = vK;  code = codeK;
    end
    attempt = attempt + 1;
end
end

function [v, code] = mf_or_neginf(s)
% A probe counts toward the golden search ONLY if IPOPT converged to a
% feasible point (stats.success) AND the relaxation is actually tight
% there (lossless_gap under the same threshold test_convex_lossless
% checks). Without the tightness condition, a "successful" but untight
% iterate (observed at tf=17: status Solved_To_Acceptable_Level, gap
% ~1e-2, ~100x the optimum's ~1e-4) biases mf UPWARD -- exactly the
% direction that wins a maximization search -- and the golden section has
% no other defense against picking it.
%
% ADAPTATION (task-7 fix report round 3, documented): tightness is now
% checked BEFORE the exact status string, not after -- was: any status
% other than the literal 'Solve_Succeeded' (i.e. also
% 'Solved_To_Acceptable_Level') was rejected outright regardless of gap.
% That blanket rule turned out to be over-strict: re-solving after the
% terminal-velocity BC change (see booster_params.m), a genuine
% golden-section probe deterministically converged to
% 'Solved_To_Acceptable_Level' with an EXCELLENT gap (3.66e-5, same order
% as fully-"Solve_Succeeded" points elsewhere on the same curve) --
% reproducible bit-for-bit across three separate re-solves at that exact
% tf, so not the run-to-run nondeterminism documented in
% solve_fixed_tf_probe above. Rejecting it regardless of its (tight) gap
% censored a physically valid point sitting right next to the search's
% true optimum, and the golden section converged 52 kg short of it (see
% report). The original rationale -- guard against a bad iterate biasing
% mf upward -- only needs the TIGHTNESS test; whether IPOPT's internal
% tolerance tier says "fully succeeded" or "acceptable level" is not by
% itself evidence the relaxation failed. code 2 (tight, acceptable-level)
% is now valid (v=mf) alongside code 3 (tight, fully succeeded); the
% distinction is kept only as a diagnostic breadcrumb in tf_curve.
%   code: 3 = valid, fully converged; 2 = valid, acceptable-level only;
%         1 = converged but relaxation not tight (any status);
%         0 = solver failed/threw (opti.debug iterate, not a solution).
% Keyed to the GUIDANCE ceiling (task-7b) so the tightness classifier keeps
% the same meaning relative to the problem actually being solved. Slightly
% stricter than the old full-Tmax form; measured gaps (~5e-5) sit ~50x below
% either value, so no probe changes classification.
tightTol = 1e-4 * s.P.etaT * s.P.Tmax / s.P.m0;
if ~s.stats.success
    code = 0;
elseif s.lossless_gap >= tightTol
    code = 1;
elseif ~strcmp(s.stats.status, 'Solve_Succeeded')
    code = 2;
else
    code = 3;
end
if code >= 2, v = s.mf; else, v = -Inf; end
end

function sol = solve_fixed_tf(P, tf, Nc)
% One convex subproblem at fixed tf, trapezoidal on the LINEAR dynamics.
import casadi.*
al = 1 / (P.Isp * P.g0);
% GUIDANCE ceiling is DE-RATED (task-7b, P.etaT): every reference below that
% represents "the most thrust the GUIDANCE may use" -- the max-thrust
% depletion reference z0/zlb, its Taylor upper bound mu2, the matching
% initial guess Sg, and this depletion feasibility assert -- uses TmaxG.
% P.Tmin (engine floor) and the zub bound built from it are NOT de-rated.
TmaxG = P.etaT * P.Tmax;
assert(P.m0 - al*TmaxG*tf > 0, ...
    'solve_pdg_convex:massDepletion', ...
    'max-thrust depletion reference goes non-positive at tf=%.3f (m0=%.1f, al*TmaxG*tf=%.1f) -- shrink tf or widen mdry margin', ...
    tf, P.m0, al*TmaxG*tf);
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
z0  = log(P.m0 - al*TmaxG*t);             % reference depletion (de-rated)
zlb = log(P.m0 - al*TmaxG*t);             % lower bound on z  (de-rated)
zub = log(P.m0 - al*P.Tmin*t);            % upper bound on z
cotg = 1/tand(P.gs_deg);
for k = 1:Nc
    opti.subject_to(sum(Uu(:,k).^2) <= S(k)^2);
    opti.subject_to(S(k) >= 0);
    mu1 = P.Tmin*exp(-z0(k));  mu2 = TmaxG*exp(-z0(k));
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
% Terminal velocity ADJUDICATED 2026-08-08 (task-7 fix report round 3):
% Vh(:,end) == P.vf/Vc, was zeros(3,1) -- see P.vf's comment in
% booster_params.m for why (Tmin>weight makes v(tf)=0 singular).
opti.subject_to(Rh(:,1) == P.r0/Lc);   opti.subject_to(Vh(:,1) == P.v0/Vc);
opti.subject_to(Z(1)   == log(P.m0));
opti.subject_to(Rh(:,end) == zeros(3,1));  opti.subject_to(Vh(:,end) == P.vf/Vc);
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
Sg = TmaxG*exp(-z0);
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
