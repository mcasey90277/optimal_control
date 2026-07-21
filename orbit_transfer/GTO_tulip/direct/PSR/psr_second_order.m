function so = psr_second_order(solFile, opts)
% PSR_SECOND_ORDER  Second-order local-minimality certificate for a PSR solution.
%
% Upgrades the pipeline's first-order (PMP extremality) certificate toward
% LOCAL MINIMALITY by testing the second-order sufficient conditions (SSOSC)
% of the discrete NLP that PSR actually solves. At the KKT point, SSOSC holds
% iff the Hessian of the Lagrangian is positive definite on the CRITICAL CONE
% (the null space of the active constraints, restricted by strict
% complementarity). Instead of forming that (large) null space, this uses the
% equivalent, efficient KKT-INERTIA test:
%
%     reduced Hessian Z'HZ > 0   <=>   K = [ H  A' ; A  0 ]  has inertia
%                                       (n, m_active, 0)
%
% where H = Hessian of the Lagrangian (n x n), A = Jacobian of the ACTIVE
% constraints (m_active x n), n = #decision variables. K is sparse symmetric
% indefinite; its inertia comes from a sparse LDL^T factorization (MA57). If K
% has exactly n positive, m_active negative, and 0 zero eigenvalues, the
% reduced Hessian is positive definite and the solution is a STRICT local
% minimizer of the NLP.
%
% Why this handles the bang-bang degeneracy: the fuel objective is linear in
% the throttle, so naive control-space curvature is zero. But at a bang-bang
% solution the throttle sits at its bound on every arc with STRICT
% complementarity (switching function S != 0 there => nonzero bound
% multiplier), so those directions leave the critical cone. The reduced
% Hessian then measures the genuine curvature that remains -- the CR3BP
% dynamics (defect-constraint Hessians, multiplier-weighted) and the
% thrust-direction unit-sphere constraint -- over state-and-near-switch
% perturbations that stay feasible.
%
% SCOPE / WHAT THIS PROVES (be honest): when it CERTIFIES, this proves the
% DISCRETE (transcribed) problem is a strict local minimizer -- strong evidence
% at the mesh's resolution. It is NOT identical to continuous-OCP local
% minimality, which for a bang-bang extremal additionally requires the
% Maurer-Osmolovskii conditions (strengthened Legendre-Clebsch on any singular
% arc -- none here, all bang-bang; transversal switching Sdot != 0 at each
% switch -- reported as a complementary check; and no conjugate point of the
% associated switching-time problem).
%
% *** FINDING (2026-07-12, tested on the 1.15x refined solution) ***
% For an eps=0 MASS-optimal (bang-bang) solution this NLP-level test is
% typically NOT APPLICABLE, and that is STRUCTURAL, not a bug. The fuel
% objective is LINEAR in the throttle, so the solver leaves the throttle NEAR
% but not AT its bounds (s ~ 0.02 / 0.98, edge ~ 99%); the bound constraints
% are then inactive at KKT tolerance, the active set is ambiguous, and the
% reduced Hessian is degenerate (many flat + ambiguous-sign directions). The
% reconstruction here self-validates to ~1e-12 (it is correct); it simply
% reports "NOT APPLICABLE (bang-bang degeneracy)" in that case. The correct
% local-minimality certificate for bang-bang is the SWITCHING-TIME
% (Maurer-Osmolovskii) reduced Hessian: parameterize by the k switch times,
% integrate the arcs, and test the k x k Hessian of fuel-subject-to-terminal-
% constraints (well-posed, low-dimensional). That is the recommended next
% build. This function still applies cleanly to a REGULARIZED (eps>0) solution,
% where the throttle is a smooth strictly-convex interior control.
%
% SELF-VALIDATION: the NLP is reconstructed here to mirror
% lib/casadi_minfuel_sundman.m EXACTLY. Correctness is auto-checked before any
% verdict by (a) primal feasibility of the loaded solution and (b) the
% Lagrangian stationarity residual ||grad f + A' lam|| using the solution's
% own KKT multipliers (out.lamAll). A large residual means the reconstruction
% does not match the solver (or the seed is not a KKT point) -> the function
% ERRORS rather than issue a bogus certificate.
%
% INPUTS:
%   solFile - solution .mat in seed layout: out (X [8xnN], U [4xnN], and
%             REQUIRED out.lamAll = the NLP constraint multipliers), sigma,
%             tauf0, rv0, rvf, factor. (PSR direct/refined files and PSR_data
%             products all carry this.)
%   opts    - (optional) struct:
%             tolActive - |g_i| below which an inequality is ACTIVE [1e-7]
%             tolComp   - |lam_i| above which an active bound is STRICTLY
%                         complementary [1e-9]
%             tolEig    - magnitude below which a pivot eigenvalue is ZERO
%                         [1e-9 * median|H diag|]
%             statTol   - max Lagrangian-stationarity residual accepted in the
%                         self-check [1e-4]
%             verbose   - print progress [true]
%
% OUTPUTS:
%   so - struct:
%     .certLocalMin  - logical: SSOSC holds (strict local min of the NLP)
%     .n .mActive    - #variables, #active constraints used in K
%     .nPos .nNeg .nZero - inertia of K
%     .redHessMinEig - smallest eigenvalue of the reduced Hessian, recovered
%                      as (nNeg over target) diagnostics; NaN if not computed
%     .statResid     - Lagrangian stationarity residual (self-check)
%     .maxEqInfeas   - worst equality-constraint violation (self-check)
%     .nActiveBnd .nWeakActive - active inequality count; weakly-active
%                      (non-strictly-complementary) count (SSOSC assumes 0)
%     .switchTransversalMinAbs - min |dS/dtau| over switches (bang-bang
%                      transversality; > 0 required); NaN if unavailable
%     .verdict       - human-readable one-line verdict
%
% REFERENCES:
%   [1] Nocedal & Wright, Numerical Optimization, 2e, Ch. 12 (SSOSC), Ch. 16
%       (KKT inertia / reduced Hessian).
%   [2] Maurer & Osmolovskii, SIAM J. Control Optim. 42 (2004) (bang-bang SSC).
%   [3] lib/casadi_minfuel_sundman.m (the NLP this reconstructs).

if nargin < 2, opts = struct(); end
if ~isfield(opts,'eps'),       opts.eps       = 0;    end   % homotopy eps the
                            % solution was solved at (0 = bang-bang). MUST match,
                            % else the reconstructed Lagrangian won't be stationary.
if ~isfield(opts,'tolActive'), opts.tolActive = 1e-7; end
if ~isfield(opts,'tolComp'),   opts.tolComp   = 1e-9; end
if ~isfield(opts,'statTol'),   opts.statTol   = 1e-4; end
if ~isfield(opts,'verbose'),   opts.verbose   = true; end

vp = @(varargin) fprintf(varargin{:});
if ~opts.verbose, vp = @(varargin) []; end

% ---- load solution ----------------------------------------------------------
S = load(solFile);
X = S.out.X;  U = S.out.U;  sigma = S.sigma(:);  tauf0 = S.tauf0;
rv0 = S.rv0(:); rvf = S.rvf(:);
nN = size(X,2);  N = nN - 1;
n  = 12*nN;                                   % [X(:); U(:)]
tf = X(8,end);
assert(isfield(S.out,'lamAll') && numel(S.out.lamAll) >= 1, ...
    'psr_second_order:noDuals', 'out.lamAll (NLP multipliers) required');
lamAll = S.out.lamAll(:);
p = cr3bp_lt_params(0.025, 15, 2100);
Tmax = p.Tmax;  c = p.c;  muStar = p.muStar;  pSund = 1.5;
dsig = diff(sigma).';                          % 1 x N

% ---- reconstruct the EXACT NLP (mirror lib/casadi_minfuel_sundman.m) --------
cpath = getenv('CASADI_PATH');
if isempty(cpath), cpath = fullfile(getenv('HOME'), 'casadi-3.7.0'); end
addpath(cpath);  import casadi.*

x = MX.sym('x', 8);  u = MX.sym('u', 4);
r = x(1:3);  v = x(4:6);  m = x(7);  al = u(1:3);  s = u(4);
dd = [r(1)+muStar; r(2); r(3)];
rr = [r(1)-1+muStar; r(2); r(3)];
r1 = sqrt(dd.'*dd + 1e-12);
d3 = (dd.'*dd + 1e-12)^1.5;  r3 = (rr.'*rr + 1e-12)^1.5;
gr = [r(1); r(2); 0] - (1-muStar)*dd/d3 - muStar*rr/r3;
hv = [2*v(2); -2*v(1); 0];
accel = gr + hv + (s*Tmax/m)*al;
mdot  = -(Tmax/c)*s;
kappa = r1^pSund;
fdyn  = Function('f', {x,u}, {kappa*[v; accel; mdot; 1]});
Fmap  = fdyn.map(nN);
gint  = Function('g', {x,u}, {[kappa*s; kappa*s*(1-s)]});  % [q_fuel; q_smooth]
Gmap  = gint.map(nN);

opti = Opti();
Xv = opti.variable(8, nN);
Uv = opti.variable(4, nN);
tauf = tauf0;
Fx = Fmap(Xv, Uv);
D  = Xv(:,2:end) - Xv(:,1:end-1) - tauf*(repmat(dsig,8,1)/2).*(Fx(:,1:end-1) + Fx(:,2:end));
opti.subject_to(D(:) == 0);                                       % (1) defects  8N
opti.subject_to((sum(Uv(1:3,:).^2, 1) - 1).' == 0);              % (2) unit     nN
lbX = repmat([-3;-3;-3;-12;-12;-12;0.3;0], 1, nN);
ubX = repmat([ 3; 3; 3; 12; 12; 12;1.0; 2*tf], 1, nN);
opti.subject_to(Xv(:) >= lbX(:));                                % (3) X lower  8nN
opti.subject_to(Xv(:) <= ubX(:));                                % (4) X upper  8nN
lbU = repmat([-1.1;-1.1;-1.1;0], 1, nN);
ubU = repmat([ 1.1; 1.1; 1.1;1], 1, nN);
opti.subject_to(Uv(:) >= lbU(:));                                % (5) U lower  4nN
opti.subject_to(Uv(:) <= ubU(:));                                % (6) U upper  4nN
opti.subject_to(Xv(1:6,1) == rv0(1:6));                          % (7) rv0      6
opti.subject_to(Xv(7,1) == 1);                                   % (8) m0       1
opti.subject_to(Xv(8,1) == 0);                                   % (9) t0       1
opti.subject_to(Xv(1:6,nN) == rvf(1:6));                         % (10) rvf     6
opti.subject_to(Xv(8,nN) == tf);                                 % (11) tf      1
Gf = Gmap(Xv, Uv);
IntF = tauf * sum((dsig/2).*(Gf(1,1:end-1) + Gf(1,2:end)));      % fuel   Int[s]dt
IntS = tauf * sum((dsig/2).*(Gf(2,1:end-1) + Gf(2,2:end)));      % smooth Int[s(1-s)]dt
opti.minimize(IntF - opts.eps*IntS);                            % J(eps), matches solver

% ---- symbolic KKT pieces ----------------------------------------------------
xs = opti.x;  gs = opti.g;  fs = opti.f;  ls = opti.lam_g;
M  = size(gs,1);
Lag  = fs + ls.'*gs;
Hfun = Function('H',  {xs, ls}, {hessian(Lag, xs)});
Afun = Function('A',  {xs},     {jacobian(gs, xs)});
gfun = Function('gg', {xs},     {gs});
gLfun= Function('gL', {xs, ls}, {gradient(Lag, xs)});

wstar = [X(:); U(:)];
assert(numel(lamAll) == M, 'psr_second_order:dualCount', ...
    'lamAll has %d entries but the NLP has %d constraints -- reconstruction/mesh mismatch', ...
    numel(lamAll), M);

vp('psr_second_order: n=%d vars, M=%d constraints (nN=%d)\n', n, M, nN);

% ---- self-validation: feasibility + Lagrangian stationarity -----------------
% Classify rows from Opti's OWN canonical bounds (lbg/ubg), not a hardcoded
% order: equality <=> lbg==ubg; inequality is g in [lbg,ubg] (one side inf).
gnum  = full(gfun(wstar));
gLnum = full(gLfun(wstar, lamAll));
lbg = full(evalf(opti.lbg));  ubg = full(evalf(opti.ubg));   % -> double (incl. Inf)
isEqRow = (lbg == ubg) & isfinite(lbg);              % equality constraints
maxEqInfeas = max(abs(gnum(isEqRow) - lbg(isEqRow)));
ineqViol    = max([0; max(lbg - gnum, 0); max(gnum - ubg, 0)]);  % >=0
so.statResid = norm(gLnum, inf);  so.maxEqInfeas = maxEqInfeas;
vp('  self-check: max eq infeasibility %.2e, ineq violation %.2e, ||gradL||inf %.2e\n', ...
   maxEqInfeas, ineqViol, so.statResid);
if maxEqInfeas > 1e-5 || ineqViol > 1e-5 || so.statResid > opts.statTol
    error('psr_second_order:selfCheck', ...
        ['reconstruction/solution self-check FAILED (eqInfeas %.2e, ineqViol %.2e, ' ...
         '||gradL||inf %.2e) -- not a feasible KKT point of this NLP; refusing to certify'], ...
        maxEqInfeas, ineqViol, so.statResid);
end

% ---- active set + strict complementarity ------------------------------------
atLb    = ~isEqRow & isfinite(lbg) & (abs(gnum - lbg) < opts.tolActive);
atUb    = ~isEqRow & isfinite(ubg) & (abs(gnum - ubg) < opts.tolActive);
actIneq = atLb | atUb;                               % active inequalities
strictC = actIneq & (abs(lamAll) > opts.tolComp);    % strictly complementary
weakAct = actIneq & ~strictC;
so.nActiveBnd  = nnz(actIneq);
so.nWeakActive = nnz(weakAct);
vp('  active inequalities: %d (%d strictly complementary, %d weak)\n', ...
   nnz(actIneq), nnz(strictC), nnz(weakAct));

% active constraint set used to build K: all equalities + strictly-active ineq
actMask = isEqRow | strictC;
mActive = nnz(actMask);
so.n = n;  so.mActive = mActive;

% ---- assemble H, A_active; KKT inertia --------------------------------------
vp('  building Hessian + active Jacobian (CasADi)...\n');
Hnum = sparse(Hfun(wstar, lamAll));
Hnum = (Hnum + Hnum.')/2;                             % symmetrize (round-off)
Aall = sparse(Afun(wstar));
Aact = Aall(actMask, :);                              % m_active x n

K = [Hnum, Aact.'; Aact, sparse(mActive, mActive)];
vp('  KKT matrix %d x %d, LDL inertia...\n', size(K,1), size(K,1));
[~, Dblk, ~] = ldl(K, 'vector');
% Calibrate tolEig on the ACTUAL LDL pivot (block-eigenvalue) magnitudes, NOT
% on H's diagonal: the pivots are Schur complements on a different scale (up to
% ~cond(K)), so an H-diagonal threshold can land mid-distribution and mislabel
% genuine eigenvalues. One tol=0 pass returns all block eigenvalues; threshold
% at 1e-9 of the median nonzero magnitude.
[~, ~, ~, evK] = inertia_blockdiag(Dblk, 0);
nz = evK(evK ~= 0);
tolEig = 1e-9 * max(median(abs(nz)), 1);
so.tolEig = tolEig;
[nPos, nNeg, nZero] = inertia_blockdiag(Dblk, tolEig);
so.nPos = nPos;  so.nNeg = nNeg;  so.nZero = nZero;
% compact spectrum diagnostic (for a tolEig sweep without re-factorizing): the
% eigenvalues near zero are the only ambiguous ones; large ones are counted.
so.nPosLarge = sum(evK > 1);  so.nNegLarge = sum(evK < -1);
sm = evK(abs(evK) < 1);  [~, ord] = sort(abs(sm));
so.evNearZero = sm(ord(1:min(numel(sm), 3000)));   % |lambda|<1, closest-to-0 first

% spectrum diagnostic (tolEig sweep + negative-eigenvalue magnitudes over the
% FULL K spectrum): distinguishes conditioning noise (tiny ambiguous-sign
% eigenvalues) from genuine O(1) negative curvature.
if isfield(opts,'spectrumReport') && opts.spectrumReport
    fprintf('\n  -- spectrum diagnostic (n=%d, mActive=%d, |K|=%d) --\n', n, mActive, numel(evK));
    fprintf('  %-9s %-8s %-8s %-8s %-10s\n','tolEig','nPos','nNeg','nZero','excessNeg');
    for te = [1e-12 1e-10 1e-8 1e-6 1e-4 1e-2]
        nP = sum(evK > te);  nN = sum(evK < -te);  nZ = numel(evK) - nP - nN;
        fprintf('  %-9.0e %-8d %-8d %-8d %-10d\n', te, nP, nN, nZ, nN - mActive);
    end
    negs = abs(evK(evK < 0));  ed = [1 1e-2 1e-4 1e-6 1e-8 1e-10 0];
    fprintf('  negative eigenvalues by |magnitude|: total %d\n', numel(negs));
    for b = 1:numel(ed)-1
        fprintf('    [%.0e, %.0e): %d\n', ed(b), ed(b+1), sum(negs < ed(b) & negs >= ed(b+1)));
    end
end

% ---- verdict ----------------------------------------------------------------
% SSOSC  <=>  In(K) = (n, mActive, 0). The reduced Hessian's inertia is
% In(Z'HZ) = (nPos-r, nNeg-r, nZero-(mActive-r)) with r = rank(A_active); since
% r <= mActive, the reduced Hessian has AT LEAST (nNeg - mActive) negative
% eigenvalues regardless of rank. So:
so.excessNeg = max(0, nNeg - mActive);     % lower bound on negative red-Hess dirs
so.certLocalMin = (nPos == n) && (nNeg == mActive) && (nZero == 0) && (so.nWeakActive == 0);
so.redHessMinEig = NaN;

% ---- complementary bang-bang transversality check (Sdot != 0 at switches) ---
so.switchTransversalMinAbs = switch_transversality(S, p);

% Heuristic: is this an eps=0 bang-bang solution whose throttle sits NEAR (not
% at) its bounds? Then the active set is ambiguous and the raw NLP Hessian test
% is ILL-POSED (structural, not a solver defect) -- the mass-optimal objective
% is linear in the throttle, so the flat/near-flat throttle directions produce
% zero + ambiguous-sign eigenvalues. Detected when few throttle bounds are
% strictly active despite a near-bang edge fraction.
sthr = U(4,:);
nBangNodes = nnz(sthr < 0.05 | sthr > 0.95);
so.bangBangDegenerate = ~so.certLocalMin && (opts.eps == 0) && (so.nActiveBnd < 0.25*nBangNodes);

if so.certLocalMin
    so.verdict = sprintf(['LOCAL MIN (NLP SSOSC): reduced Hessian positive definite; ' ...
        'In(K)=(%d,%d,0)=(n,m_active,0)'], nPos, nNeg);
elseif so.bangBangDegenerate
    so.verdict = sprintf([...
        'NOT APPLICABLE (bang-bang degeneracy): the eps=0 mass-optimal NLP is ' ...
        'a degenerate KKT point for 2nd-order purposes -- throttle near (not at) ' ...
        'its bounds under a LINEAR objective => ambiguous active set (%d strictly ' ...
        'active bounds vs %d near-bang nodes), In(K)=(%d,%d,%d) with %d flat + ' ...
        '>=%d near-zero/neg reduced directions. This is STRUCTURAL for bang-bang ' ...
        'via a direct NLP, not a solver defect (reconstruction self-validated). ' ...
        'Use the switching-time (Maurer-Osmolovskii) 2nd-order test instead.'], ...
        so.nActiveBnd, nBangNodes, nPos, nNeg, nZero, nZero, so.excessNeg);
else
    so.verdict = sprintf(['INCONCLUSIVE: In(K)=(%d,%d,%d), expected (%d,%d,0); ' ...
        '%d flat, >=%d negative reduced directions -- inspect active set / rank'], ...
        nPos, nNeg, nZero, n, mActive, nZero, so.excessNeg);
end
vp('  %s\n', so.verdict);
if ~isnan(so.switchTransversalMinAbs)
    vp(['  bang-bang transversality (complementary ingredient): ' ...
        'min |dS/dtau| at switches = %.3e (> 0 good)\n'], so.switchTransversalMinAbs);
end
end

% =============================================================================
function [np, nn, nz, evAll] = inertia_blockdiag(D, tol)
% INERTIA_BLOCKDIAG  Inertia of the block-diagonal D from a sparse LDL (MA57):
% 1x1 and 2x2 blocks on the (super)diagonal. Returns (#pos, #neg, #zero) at
% threshold `tol`, plus evAll = ALL block eigenvalues (for a tol sweep).
nD = size(D,1);
dg = full(diag(D));  od = full(diag(D,1));
evAll = zeros(nD,1);  m = 0;  k = 1;
while k <= nD
    if k < nD && od(k) ~= 0                  % 2x2 block
        B = [dg(k) od(k); od(k) dg(k+1)];
        evAll(m+1:m+2) = eig(B);  m = m + 2;  k = k + 2;
    else                                     % 1x1 block
        m = m + 1;  evAll(m) = dg(k);  k = k + 1;
    end
end
evAll = evAll(1:m);
np = sum(evAll >  tol);
nn = sum(evAll < -tol);
nz = sum(abs(evAll) <= tol);
end

% =============================================================================
function md = switch_transversality(S, p)
% SWITCH_TRANSVERSALITY  min |dS/dtau| at the switching-function zero crossings,
% the strong bang-bang second-order ingredient (transversal crossing). Uses the
% dual-map costates (mode 'd') to form S(tau); returns NaN if unavailable.
md = NaN;
try
    tmp = [tempname '.mat'];
    out = S.out;  sigma = S.sigma;  tauf0 = S.tauf0; %#ok<NASGU>
    rv0 = S.rv0;  rvf = S.rvf;  factor = S.factor;   %#ok<NASGU>
    save(tmp, 'out', 'sigma', 'tauf0', 'rv0', 'rvf', 'factor');
    [~, sd, info] = sms_seed_duals(tmp, 40, 1e-4, 'd');
    delete(tmp);
    lam = info.Y16(9:16, :);  X = info.X;  tauN = info.tauN;
    Sn  = 1 - sqrt(sum(lam(4:6,:).^2,1)).*sd.c./X(7,:) - lam(7,:);
    cr  = find(diff(sign(Sn)) ~= 0);
    slopes = zeros(1, numel(cr));
    for q = 1:numel(cr)
        k = cr(q);
        slopes(q) = abs((Sn(k+1) - Sn(k)) / (tauN(k+1) - tauN(k)));
    end
    if ~isempty(slopes), md = min(slopes); end
catch
end
end
