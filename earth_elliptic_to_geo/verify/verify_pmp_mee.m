function ver = verify_pmp_mee(out, par, sigma, opts)
% VERIFY_PMP_MEE  First-order PMP consistency from the MEE/sigma-domain NLP's
% KKT duals -- the MEE analog of verify_pmp_2body.m, built for this
% transcription's own dual->costate map (casadi_lt_mee.m has no cScale-style
% per-node slack state -- DeltaL is one free scalar column -- so the
% Campaign-B cScale-anomaly class does not apply here; still verified, not
% assumed, per node-level gates below).
%
% PIPELINE:
%   1. mee_dual_to_costate: interval defect duals out.lamDef -> nodal costate
%      lam, via the step-weighted adjacent-interval average (mandatory
%      correctness fix, DESIGN_dual_map.md), one-sided at the endpoints.
%   2. mee_primer_switch: forms the primer vector primerVec = (Tmax/m)*Ldot*
%      p_el - G*K_L*e3 (p_el = B(X)'*lam_el, the MEE primer; G the sigma-
%      domain Hamiltonian bracket; K_L the Ldot-on-control coupling
%      coefficient -- see that file's header for the full derivation) and the
%      switching function S = C1*Ldot0 - lam_t*K_L*beta_3 (C1 the thr-
%      coefficient of G), both derived by EXACT stationarity of the
%      sigma-domain Hamiltonian H_sigma = (DeltaL/Ldot(X,U))*G(X,U,lam) w.r.t.
%      beta (Lagrange multiplier for ||beta||=1) and thr respectively -- no
%      approximation, the Ldot-on-control coupling is kept, not dropped.
%   3. Global sign resolution: IPOPT/CasADi's opti.dual Lagrangian-sign
%      convention (L = f +/- lam'*g) is a standard, non-anomalous ambiguity
%      of a SINGLE overall sign on the whole costate vector; resolved here by
%      the same empirical best-fit-cosine rule verify_pmp_2body.m uses
%      (mean cosine of thrust vs. predicted direction on burn arcs > 0),
%      applied to lam BEFORE recomputing primerVec/S (not by negating S
%      directly -- S's C1 term has a "+1" constant piece from the physical
%      running-cost weight that must NOT flip with the dual-sign ambiguity,
%      so S has to be recomputed from the flipped lam, not sign-flipped
%      itself).
%   4. Gates (Campaign-B T1 acceptance, binding): primer misalignment < 1 deg
%      (MEDIAN over burn nodes) and switching-sign agreement >= 99% of ALL
%      nodes (S<0 <=> thr=1, S>0 <=> thr=0). Also reports the tangential
%      residual (I-beta*beta')*primerVec at burn nodes (its norm
%      distribution) -- NOTE: despite the similar name, this is NOT
%      DESIGN_dual_map.md's multi-group Lagrangian-residual T1 test (that
%      test assembles the FULL NLP Lagrangian gradient from ALL dual groups
%      -- defect, cone, terminal/equality, bound multipliers -- and is not
%      built here; this transcription has no cone constraint on beta to
%      begin with, since ||beta||=1 is enforced by construction, not as an
%      NLP constraint with its own multiplier). This quantity is
%      mathematically the same information as the primer misalignment angle
%      recast in norm form (both come from decomposing primerVec against the
%      burn direction beta), reported separately by explicit request. The
%      full multi-group Lagrangian-residual T1 remains future work. Also
%      reports a coupling-strength diagnostic K_L/Ldot0 (should be small for
%      "low thrust"; NOT assumed anywhere in the derivation, just reported).
%
% INPUTS:
%   out   - casadi_lt_mee result struct [fields .X .U .dL .lamDef]
%   par   - kepler_lt_params struct [.Tmax .c .mu]
%   sigma - node grid used for the solve, monotonic [(N+1)x1]
%   opts  - struct: .eps (homotopy parameter the solution was solved at;
%           default 0 -- the switching-function derivation is exact ONLY at
%           eps=0 since the eps>0 running cost is quadratic in thr, breaking
%           the linear-fractional/monotone-in-thr argument; asserts eps==0)
%
% OUTPUTS:
%   ver - struct:
%     .primerMedianDeg .primerMeanDeg  - primer misalignment on burns [deg]
%     .overallSignPct                  - pct of ALL nodes with correct S sign
%     .burnSignPct .coastSignPct       - same, split by arc (diagnostic)
%     .tangentialResidNormRelMedian .tangentialResidNormRelMax
%                                       - ||tangential primerVec component|| /
%                                         ||primerVec|| on burns -- primer-
%                                         angle information in norm form, NOT
%                                         DESIGN_dual_map.md's full multi-
%                                         group Lagrangian-residual T1 (see
%                                         pipeline step 4 note above; that
%                                         test remains future work)
%     .KLoverLdot0Median .KLoverLdot0Max
%                                       - |K_L/Ldot0| coupling-strength stats
%     .maxSwitchAlignErr               - max |t of nearest S=0 crossing - t of
%                                         nearest solver thr=0.5 switch| (H1
%                                         pinpoint-zero check, physical time)
%     .singularArcNodes                - count of nodes inside a run of >=3
%                                         consecutive |S|<tol nodes (H2 check;
%                                         expect 0)
%     .lamMendRel                      - |lam_m(end)|/max|lam_m| (free-mass
%                                         transversality; bonus, not gated)
%     .pass                            - primerMedianDeg<1 && overallSignPct>=99
%     .sigma .L .t .thr .burn .S .primerDeg .KLoverLdot0N .tCross .tSwitch
%                                       - per-node arrays for fig_switching.m
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/verify_pmp_2body.m (Cartesian sibling; same
%       empirical-sign / gate-reporting pattern, DIFFERENT dual->costate map
%       since this transcription has no cScale).
%   [2] earth_elliptic_to_geo/mee_dual_to_costate.m, mee_primer_switch.m (the
%       two pieces this file orchestrates).
%   [3] earth_elliptic_to_geo/DESIGN_dual_map.md (Campaign-B's own T1
%       acceptance test is the full multi-group Lagrangian residual, NOT the
%       tangentialResid quantity reported here -- see pipeline step 4 note).
%   [4] Haberkorn, Martinon, Gergaud, JGCD 27(6), 2004, Fig. 16 (H1/H2).
if nargin < 4, opts = struct(); end
d = @(f, v) optdef(opts, f, v);
eps0 = d('eps', 0);
assert(eps0 == 0, ['verify_pmp_mee: switching-function derivation assumes ' ...
    'eps=0 (fuel); the eps>0 homotopy running cost is quadratic in thr and ' ...
    'breaks the linear-fractional-in-thr monotonicity argument S is built on.']);

sigma = sigma(:);
X = out.X;  U = out.U;  dL = out.dL;
Nn = size(X, 2);
assert(size(U, 2) == Nn && numel(sigma) == Nn, 'verify_pmp_mee: size mismatch X/U/sigma');
burn = U(4, :) > 0.5;

% --- step 1: dual -> nodal costate ------------------------------------------
lam = mee_dual_to_costate(out.lamDef, sigma);

% --- step 2 (pass 1, raw sign) + step 3: resolve the global sign ------------
[primerVec0, ~, ~] = mee_primer_switch(X, U, lam, sigma, dL, par);
pvn0 = sqrt(sum(primerVec0.^2, 1));
cosb0 = zeros(1, Nn);
for k = 1:Nn
    cosb0(k) = dot(U(1:3, k), -primerVec0(:, k)) / max(pvn0(k), 1e-30);
end
if mean(cosb0(burn)) < 0
    lam = -lam;
end

% --- step 2 (pass 2, correctly-signed lam) ----------------------------------
[primerVec, S, info] = mee_primer_switch(X, U, lam, sigma, dL, par);
pvn = sqrt(sum(primerVec.^2, 1));
cosb = zeros(1, Nn);
for k = 1:Nn
    cosb(k) = dot(U(1:3, k), -primerVec(:, k)) / max(pvn(k), 1e-30);
end
primerDegAll = real(acosd(max(-1, min(1, cosb))));
primerDegAll(~burn) = nan;

% --- step 4: gates -----------------------------------------------------------
ver = struct();
ver.primerMedianDeg = median(primerDegAll(burn));
ver.primerMeanDeg   = mean(primerDegAll(burn));

predBurn = S < 0;
ver.overallSignPct = 100 * mean(predBurn == burn);
ver.burnSignPct    = 100 * mean(S(burn)  < 0);
ver.coastSignPct   = 100 * mean(S(~burn) > 0);

% Tangential residual: (I - beta*beta')*primerVec at burn nodes. This is the
% primer-angle information recast in norm form -- NOT DESIGN_dual_map.md's
% full multi-group Lagrangian-residual T1 (defect + cone + terminal + bound
% duals); that test remains future work (see pipeline step 4 note above).
nB = nnz(burn);  Rtan = zeros(1, nB);  pvB = pvn(burn);
Ub = U(1:3, burn);  pvVecB = primerVec(:, burn);
for q = 1:nB
    b = Ub(:, q);
    Rtan(q) = norm(pvVecB(:, q) - dot(pvVecB(:, q), b) * b);
end
RtanRel = Rtan ./ max(pvB, 1e-30);
ver.tangentialResidNormRelMedian = median(RtanRel);
ver.tangentialResidNormRelMax    = max(RtanRel);

ver.KLoverLdot0Median = median(abs(info.KLoverLdot0));
ver.KLoverLdot0Max    = max(abs(info.KLoverLdot0));

% H1: pinpoint-zero check -- align S=0 crossings (in time) against the
% solver's own thr=0.5 crossings
t = X(7, :);
crossI = find(diff(sign(S)) ~= 0);
tCross = zeros(1, numel(crossI));
for q = 1:numel(crossI)
    k = crossI(q);
    tCross(q) = t(k) + (0 - S(k)) * (t(k+1) - t(k)) / (S(k+1) - S(k));
end
swI = find(diff(double(burn)) ~= 0);
tSwitch = 0.5 * (t(swI) + t(swI + 1));
if isempty(tSwitch) || isempty(tCross)
    ver.maxSwitchAlignErr = nan;
else
    err = zeros(1, numel(tSwitch));
    for q = 1:numel(tSwitch)
        err(q) = min(abs(tCross - tSwitch(q)));
    end
    ver.maxSwitchAlignErr = max(err);
end

% H2: no singular arc -- flag any run of >=3 consecutive near-zero-S nodes
tolS = 1e-6 * max(abs(S));
nearZero = abs(S) < max(tolS, 1e-12);
runLen = 0;  singCount = 0;
for k = 1:Nn
    if nearZero(k), runLen = runLen + 1; else, runLen = 0; end
    if runLen >= 3, singCount = singCount + 1; end
end
ver.singularArcNodes = singCount;

ver.lamMendRel = abs(lam(6, end)) / max(abs(lam(6, :)));

ver.pass = ver.primerMedianDeg < 1.0 && ver.overallSignPct >= 99;

% --- arrays for fig_switching.m ---------------------------------------------
ver.sigma = sigma(:).';
ver.L     = pi + ver.sigma * dL;
ver.t     = t;
ver.thr   = U(4, :);
ver.burn  = burn;
ver.S     = S;
ver.primerDeg   = primerDegAll;
ver.KLoverLdot0N = info.KLoverLdot0;
ver.tCross  = tCross;
ver.tSwitch = tSwitch;

fprintf(['verify_pmp_mee: primer median %.3f deg (mean %.3f) | sign agree ' ...
         '%.2f%% (burn %.1f%% coast %.1f%%) | tangentialResid median %.2e ' ...
         'max %.2e | K_L/Ldot0 median %.2e max %.2e | switchAlignErr %.2e | ' ...
         'singularArcNodes %d | lamM rel %.1e | pass=%d\n'], ...
    ver.primerMedianDeg, ver.primerMeanDeg, ver.overallSignPct, ...
    ver.burnSignPct, ver.coastSignPct, ver.tangentialResidNormRelMedian, ...
    ver.tangentialResidNormRelMax, ver.KLoverLdot0Median, ver.KLoverLdot0Max, ...
    ver.maxSwitchAlignErr, ver.singularArcNodes, ver.lamMendRel, ver.pass);
end
