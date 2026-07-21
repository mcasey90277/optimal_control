function cert = certify_minfuel_pmp(solFile, makePlot)
% CERTIFY_MINFUEL_PMP  Tier-1 PMP certification of the Sundman min-fuel solution.
%
% STATUS (2026-07-08): PARTIAL -- returns a REVIEW verdict, not PASS. This
% continuous-adjoint recovery DEFEATS the ~1e11 conditioning wall in the solve
% (hard-constrained saddle/lsqlin -> adjoint recursion residual ~1e-14), but it
% cannot fully certify: the recursion B_k here is only a 2nd-order (trapezoidal)
% discretization of the CONTINUOUS adjoint, and over the ~40-rev arc that
% approximation error -- amplified by the perigee sensitivity -- makes ~1/4 of
% the NLP burn-node thrust directions mutually unfittable by any single costate.
% The correct certificate for a DISCRETE optimum uses the EXACT discrete adjoint
% (the NLP's own KKT duals), not a reconstructed continuous costate. Kept as a
% consistency diagnostic + scaffold (report/figure/switch-alignment reusable).
% See TIER1_PMP_CERTIFICATION_SCOPE.md (FINDING) for the full account.
%
% Independent verification that the converged Sundman-regularized minimum-fuel
% GTO -> tulip solution is a genuine Pontryagin extremal. The position-velocity
% costate block is a LINEAR homogeneous system driven only by the state; the
% optimal thrust direction is the primer alpha = -lam_v/||lam_v||, so at burn
% nodes the DIRECTION of lam_v is known exactly. We recover the costate history
% by requiring the discrete adjoint recursion AND the primer-direction data to
% hold simultaneously, solved as a single sparse least-squares system over all
% node costates.
%
% Why a global solve, not forward shooting: over the ~40-rev spiral the
% homogeneous costate map amplifies by ~1e11 (reciprocal-eigenvalue/symplectic
% growth) even WITH the Sundman kappa scaling, so propagating lam(0) forward and
% fitting is hopelessly ill-conditioned. The block-bidiagonal global solve is
% BVP-stable for exactly the reason multiple shooting beats single shooting: it
% never forward-multiplies the ill-conditioned transition map. This mirrors the
% primer-vector verification of ../../min_energy_tutorial/primer_check.m, adapted to
% a long multi-rev arc.
%
% Recovery pipeline:
%   1. recover lam_r,lam_v (up to a global positive scale) from the adjoint
%      recursion + primer directions (sparse LS, one normalization row);
%   2. lam_m is then closed-form: dlam_m/dt = -||lam_v|| s Tmax/m^2 and the
%      free-final-mass transversality lam_m(tau_f)=0 fix it up to the same scale;
%   3. the scale c is pinned in closed form by S = 1 - c*W = 0 at the switch
%      nodes (least squares over the 25 switches).
% Then S(t) = 1 - ||lam_v|| c/m - lam_m is checked for the bang-bang sign law.
%
% INPUTS:
%   solFile  - (optional) path to a certified-solution .mat OR a struct with
%              fields {out,sigma,tauf0,pSund,eps,rv0,rvf}. Default:
%              'sundman_minfuel_certified.mat' beside this file. [char|struct]
%   makePlot - (optional) draw the S(t) vs throttle figure [logical, default true]
%
% OUTPUTS:
%   cert - certification struct:
%       .recursionResid  max ||Lam_{k+1}-B_k Lam_k|| after the LS solve [scalar]
%       .primerDirErr    max ||alpha + lam_v/||lam_v|||| over burn nodes [scalar]
%       .scale           recovered positive costate scale c [scalar]
%       .scaleSpread     std/mean of per-switch scale estimates (consistency) [scalar]
%       .signMatchFrac   fraction of strict (full/coast) nodes whose S sign
%                        matches the throttle law (S<0<->burn) [scalar]
%       .nSwitchNLP      throttle switches in the NLP solution [scalar]
%       .nSwitchS        sign changes of S [scalar]
%       .switchesMatched NLP switches with an S zero-crossing within tolNodes
%       .transversality  lam_m(tau_f) (0 by construction; reported for record)
%       .S, .throttle, .tDays   per-node arrays for plotting
%       .passed          all checks within tolerance [logical]
%
% REFERENCES:
%   [1] Lawden, "Optimal Trajectories for Space Navigation," 1963 (primer vector).
%   [2] Bertrand & Epenoy, OCAM 23(4), 2002 (min-fuel switching function).
%   [3] cf. ../../min_energy_tutorial/primer_check.m (verify-along-trajectory idea).

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(solFile)
    solFile = fullfile(here, 'sundman_minfuel_certified.mat');
end
if nargin < 2 || isempty(makePlot), makePlot = true; end

% ---- load solution ------------------------------------------------------
if isstruct(solFile), S = solFile; else, S = load(solFile); end
out   = S.out;
sigma = S.sigma(:);
if isfield(S,'tauf0'), tauf = S.tauf0; else, tauf = out.tauf; end
if isfield(S,'pSund'), pSund = S.pSund; else, pSund = 1.5; end

addpath(here);
p      = cr3bp_lt_params(0.025, 15, 2100);   % matches the certified run
muStar = p.muStar;  Tmax = p.Tmax;  c = p.c;  tStar = p.tStar;

% ---- node arrays --------------------------------------------------------
X = out.X;  U = out.U;
r     = X(1:3,:);   m = X(7,:);   tphys = X(8,:);
alpha = U(1:3,:);   s = U(4,:);
nN    = size(X,2);  Nseg = nN - 1;
dsig  = diff(sigma).';                    % 1 x Nseg
burn  = s > 0.5;
swIdx = find(diff(double(burn)) ~= 0);    % throttle switch nodes

% ---- per-node gravity gradient G, dynamics matrix M, Sundman kappa ------
% M = [ 0  -G ; -I  -Hc' ],  lam_r' = -G lam_v,  lam_v' = -lam_r - Hc' lam_v.
Hct = [0 -2 0; 2 0 0; 0 0 0];  I3 = eye(3);  I6 = eye(6);
kappa = zeros(1,nN);  Mst = zeros(6,6,nN);
for k = 1:nN
    rk = r(:,k);
    dd = [rk(1)+muStar; rk(2); rk(3)];   rr = [rk(1)-1+muStar; rk(2); rk(3)];
    nd = sqrt(dd.'*dd);  nr = sqrt(rr.'*rr);
    G  = diag([1 1 0]) - (1-muStar)*(I3/nd^3 - 3*(dd*dd.')/nd^5) ...
                       -    muStar *(I3/nr^3 - 3*(rr*rr.')/nr^5);
    Mst(:,:,k) = [zeros(3), -G; -I3, -Hct];
    kappa(k)   = nd^pSund;
end

% ---- one-step adjoint maps B_k (trapezoidal in sigma) -------------------
% Lam_{k+1} = B_k Lam_k,  B_k = (I - aR M_{k+1})^{-1}(I + aL M_k).
B = cell(1,Nseg);
for k = 1:Nseg
    aL = tauf*dsig(k)/2*kappa(k);
    aR = tauf*dsig(k)/2*kappa(k+1);
    B{k} = (I6 - aR*Mst(:,:,k+1)) \ (I6 + aL*Mst(:,:,k));
end

% ---- recover costates: HARD-constrained least squares (saddle point) ----
% Solve   min_y ||P y||^2   s.t.   R y = 0 (adjoint recursion, HARD),
%                                   n' y = 1 (scale + sign normalization),
% where y = [Lam_1; ...; Lam_nN], P projects lam_v onto the plane perpendicular
% to the primer direction at burn nodes ((I - alpha alpha') lam_v = 0), and R is
% the block-bidiagonal recursion Lam_{k+1} - B_k Lam_k = 0.
%
% Enforcing R as a HARD equality constraint (not an LS penalty) is the crux: the
% KKT saddle-point solve never forward-multiplies the ~1e11-amplifying transition
% map, so it is immune to that dynamic range (the reason collocation beats
% shooting). A penalized A\b instead trades the huge recursion residual against
% the O(1) primer residual and collapses lam_v to zero -- the failure of the
% earlier min-norm formulation.
idx = @(k) (k-1)*6 + (1:6);
nCol = 6*nN;

% R: recursion, row-normalized (RHS is 0, so scaling leaves the constraint set)
iiR = []; jjR = []; vvR = [];  row = 0;
for k = 1:Nseg
    rows = row + (1:6);
    [gj,gi] = meshgrid(idx(k+1), rows);
    iiR = [iiR; gi(:)]; jjR = [jjR; gj(:)]; vvR = [vvR; reshape(I6,[],1)];
    [gj,gi] = meshgrid(idx(k), rows);
    iiR = [iiR; gi(:)]; jjR = [jjR; gj(:)]; vvR = [vvR; reshape(-B{k},[],1)];
    row = row + 6;
end
nR = row;
R  = sparse(iiR, jjR, vvR, nR, nCol);
rnorm = sqrt(sum(R.^2, 2));  rnorm(rnorm == 0) = 1;
R = spdiags(1./rnorm, 0, nR, nR) * R;

% SIGNED primer fit: the projector (I-alpha alpha')lam_v=0 is sign-BLIND (lam_v
% parallel OR anti-parallel both give zero residual), so it cannot pin the sign
% and admits spurious sign-flipped arcs. PMP requires lam_v = -rho*alpha with
% rho >= 0 (thrust opposes the primer). Enforce that with a magnitude variable
% rho_k >= 0 per burn node and the residual  lam_v,k + rho_k alpha_k = 0, solved
% as a bound-constrained sparse LS (lsqlin). The rho_k >= 0 bound is exactly what
% breaks the sign symmetry: a wrong-sign lam_v forces rho_k=0 and pays the full
% ||lam_v,k||^2 penalty, so the optimizer avoids it.
bnodes = find(burn);  nB = numel(bnodes);
nx = nCol + nB;                              % [y ; rho]
[~, kb] = max(s);  kbLoc = find(bnodes == kb, 1);

% objective rows C x = d:  lam_v,k + rho_k alpha_k = 0
iiC = []; jjC = []; vvC = [];  crow = 0;
for t = 1:nB
    k = bnodes(t);  rows = crow + (1:3);
    colsLv = (k-1)*6 + (4:6);
    [gj,gi] = meshgrid(colsLv, rows);
    iiC = [iiC; gi(:)]; jjC = [jjC; gj(:)]; vvC = [vvC; reshape(I3,[],1)];
    iiC = [iiC; rows(:)]; jjC = [jjC; (nCol+t)*ones(3,1)]; vvC = [vvC; alpha(:,k)];
    crow = crow + 3;
end
C = sparse(iiC, jjC, vvC, crow, nx);  dvec = zeros(crow,1);

% equality: recursion (hard) + scale normalization rho(kb)=1
Aeq = [ R, sparse(nR, nB);
        sparse(1, nCol), sparse(1, kbLoc, 1, 1, nB) ];
beq = [ zeros(nR,1); 1 ];

lb = -inf(nx,1);  lb(nCol+1:end) = 0;        % rho >= 0
ub =  inf(nx,1);
loo = optimoptions('lsqlin','Algorithm','interior-point', ...
                   'Display','off','MaxIterations',400, ...
                   'OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10);
x = lsqlin(C, dvec, [], [], Aeq, beq, lb, ub, [], loo);
Lam  = reshape(x(1:nCol), 6, nN);
lamv = Lam(4:6,:);  nlamv = sqrt(sum(lamv.^2,1));

% recursion residual against the ORIGINAL (unnormalized) map -- should be ~0
recR = 0;
for k = 1:Nseg
    recR = max(recR, norm(Lam(:,k+1) - B{k}*Lam(:,k)));
end

% ---- lam_m and the global scale (closed form) ---------------------------
% lam_m(t) = c*(Q(tf)-Q(t)),  Q = cumulative int of ||lam_v_base|| s Tmax/m^2.
qint = kappa .* (nlamv .* s .* Tmax ./ m.^2);      % integrand in sigma (base scale)
Q = zeros(1,nN);
for k = 1:Nseg
    Q(k+1) = Q(k) + tauf*dsig(k)/2*(qint(k) + qint(k+1));
end
W = nlamv .* c ./ m + (Q(end) - Q);                % S = 1 - scale*W
Wsw = W(swIdx);
scale = sum(Wsw) / sum(Wsw.^2);                    % LS of 1 = scale*Wsw
perSwScale = 1 ./ Wsw;                             % per-switch scale estimates

lam_m = scale * (Q(end) - Q);
Sv    = 1 - scale * W;

% ---- metrics ------------------------------------------------------------
strict = (s > 0.95) | (s < 0.05);
sgnok  = ((Sv < 0) & burn) | ((Sv >= 0) & ~burn);
signMatchFrac = mean(sgnok(strict));

sSwIdx = find(diff(double(Sv < 0)) ~= 0);
tolNodes = 3;  matched = 0;
for i = 1:numel(swIdx)
    if any(abs(sSwIdx - swIdx(i)) <= tolNodes), matched = matched + 1; end
end

primerDirErr = 0;
if any(burn)
    d = alpha(:,burn) + lamv(:,burn)./max(nlamv(burn),1e-30);
    primerDirErr = max(sqrt(sum(d.^2,1)));
end

cert = struct();
cert.recursionResid  = recR;
cert.primerDirErr    = primerDirErr;
cert.scale           = scale;
cert.scaleSpread     = std(perSwScale)/mean(perSwScale);
cert.signMatchFrac   = signMatchFrac;
cert.nSwitchNLP      = numel(swIdx);
cert.nSwitchS        = numel(sSwIdx);
cert.switchesMatched = matched;
cert.transversality  = abs(lam_m(end));
cert.S               = Sv;
cert.throttle        = s;
cert.tDays           = tphys * tStar / 86400;
cert.lamv            = lamv;        % debug
cert.nlamv           = nlamv;       % debug
cert.bnodes          = find(burn);  % debug
cert.passed = (primerDirErr < 1e-2) && (signMatchFrac > 0.99) && ...
              (matched == numel(swIdx)) && (cert.scaleSpread < 0.05);

% ---- report -------------------------------------------------------------
fprintf('\n=== Tier-1 PMP certification ===\n');
fprintf('  adjoint recursion resid : %.3e\n', recR);
fprintf('  primer direction error  : %.3e   (max over burn nodes)\n', primerDirErr);
fprintf('  costate scale c         : %.4g   (per-switch spread %.2f%%)\n', scale, 100*cert.scaleSpread);
fprintf('  switching sign match    : %.2f%%  of full/coast nodes\n', 100*signMatchFrac);
fprintf('  switches NLP / S        : %d / %d   matched %d/%d (<=%d nodes)\n', ...
        cert.nSwitchNLP, cert.nSwitchS, matched, cert.nSwitchNLP, tolNodes);
fprintf('  VERDICT: %s\n\n', ternary(cert.passed,'PASS (PMP conditions verified)','REVIEW (see metrics)'));

% ---- figure -------------------------------------------------------------
if makePlot
    figure('Color','w','Name','Tier-1 PMP certification');
    yyaxis left
    plot(cert.tDays, Sv, 'LineWidth', 1.2); hold on; yline(0, ':k');
    ylabel('switching function  S = 1 - ||\lambda_v|| c/m - \lambda_m');
    yyaxis right
    plot(cert.tDays, s, 'LineWidth', 1.0);
    ylabel('throttle  s'); ylim([-0.05 1.05]);
    for i = 1:numel(swIdx)
        xline(cert.tDays(swIdx(i)), '-', 'Alpha', 0.12, 'HandleVisibility','off');
    end
    xlabel('time (days)');
    title(sprintf('S sign-changes at %d/%d throttle switches   (primer err %.1e)', ...
                  matched, cert.nSwitchNLP, primerDirErr));
    grid on
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
