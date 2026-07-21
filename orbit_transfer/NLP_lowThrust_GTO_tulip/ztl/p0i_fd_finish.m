% P0I_FD_FINISH  Preflight: finish the 75 mN energy solve with GN on the
% CENTRAL-DIFFERENCE Jacobian (the CS Jacobian is corrupted).
%
% p0h_diag_floor found the smoking gun at the P0h floor (||R|| = 1.6e-2):
% the complex-step Jacobian through ode113 disagrees with central
% differences at O(1) (rel err 7.77, cond 9.5e10 vs 1.4e8) -- the solver has
% been steering with a corrupted J. This rerun swaps in the FD Jacobian
% (14 integrations per iteration, ~1 min) in the same equilibrated
% truncated-SVD GN + Armijo loop. Convergence here = the ladder anchor AND
% definitive proof that derivative quality (Zhang ingredient (a)) was the
% campaign-wide rate limiter; Z0's variational STM then replaces FD for
% speed and exactness.
%
% Requires: results/p0h_gn_finish.mat (the floor iterate).
% Output:   results/p0i_fd_finish.mat

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');

A = load(fullfile(resDir, 'p0h_gn_finish.mat'));
a = A.anchor;
Tmax = (a.Tmax_mN/25)*a.P.Tmax25;
resFun = @(lam) shoot_residual_energy(lam, a.tf, a.rv0, 1, a.rvf, Tmax, a.P.c, a.P.muStar);

fdJ = @(lam) local_fd_jacobian(resFun, lam);

lam = a.lam0(:);
R = resFun(lam);  rn = norm(R);
fprintf('=== P0i: FD-Jacobian GN finish @ %g mN from ||R|| = %.3e ===\n', a.Tmax_mN, rn);

maxIter = 80;  relTrunc = 1e-12;  stallCount = 0;
for it = 1:maxIter
    J = fdJ(lam);
    cs = max(abs(J), [], 1);  cs(cs == 0) = 1;
    Js = J ./ cs;
    [U, S, V] = svd(Js);  sv = diag(S);
    keep = sv >= relTrunc*sv(1);
    dz = -(V(:, keep) * ((U(:, keep)'*R) ./ sv(keep))) ./ cs(:);

    alpha = 1;  accepted = false;
    for ls = 1:14
        lamT = lam + alpha*dz;
        Rt = resFun(lamT);
        if all(isfinite(Rt)) && norm(Rt) < rn*(1 - 1e-4*alpha)
            accepted = true;  break
        end
        alpha = alpha/2;
    end
    if accepted
        prog = 1 - norm(Rt)/rn;
        lam = lamT;  R = Rt;  rn = norm(R);
        if prog < 1e-3, stallCount = stallCount + 1; else, stallCount = 0; end
    else
        relTrunc = relTrunc*100;
        stallCount = stallCount + 1;
        fprintf('  it %2d: rejected at all alpha -> relTrunc = %.0e\n', it, relTrunc);
    end
    fprintf('  it %2d: ||R|| = %.6e  alpha = %.3g  cond(Js) = %.2e\n', ...
            it, rn, alpha*accepted, cond(Js));
    if rn < 1e-8, break; end
    if stallCount >= 8, fprintf('  8 stall iterations -- stopping.\n'); break; end
end

% --- accounting -------------------------------------------------------------
optsInt = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, yI] = ode113(@lt_pmp_eom_energy, [0 a.tf], [a.rv0(:); 1; lam], ...
                 optsInt, Tmax, a.P.c, a.P.muStar);
lamvMag = sqrt(sum(yI(:,11:13).^2, 2));
u  = min(max(Tmax*(lamvMag./yI(:,7) + yI(:,14)/a.P.c), 0), 1);
mF = yI(end,7);
anchor = struct('Tmax_mN', a.Tmax_mN, 'tf', a.tf, 'lam0', lam, 'resNorm', rn, ...
    'iters', it, 'mProp_kg', a.P.m0kg*(1-mF), ...
    'dV_kms', a.P.c*log(1/mF)*a.P.lStar/a.P.tStar, 'uMin', min(u), 'uMax', max(u), ...
    'fracSatHi', mean(u > 0.999), 'fracSatLo', mean(u < 1e-3), ...
    'rv0', a.rv0, 'rvf', a.rvf, 'P', a.P);
save(fullfile(resDir, 'p0i_fd_finish.mat'), 'anchor');

fprintf(['\nP0i RESULT @ %g mN: ||R|| = %.3e after %d iters  prop = %.4f kg  ' ...
         'dV = %.4f km/s\n  throttle: min %.3f  max %.3f  sat-hi %.1f%%  sat-lo %.1f%%\n'], ...
    a.Tmax_mN, rn, it, anchor.mProp_kg, anchor.dV_kms, ...
    anchor.uMin, anchor.uMax, 100*anchor.fracSatHi, 100*anchor.fracSatLo);
if rn < 1e-8
    fprintf('GATE P0i: PASS -- LADDER ANCHOR CONVERGED at %g mN (FD J; CS J was the killer).\n', a.Tmax_mN);
else
    fprintf('GATE P0i: FAIL at %.3e -- floor survives accurate J; formulation-level next.\n', rn);
end

% ---------------------------------------------------------------------------
function J = local_fd_jacobian(resFun, lam)
% Central-difference Jacobian [7x7]; step scaled per component.
n = numel(lam);
J = zeros(n);
for k = 1:n
    h = max(1e-6*abs(lam(k)), 1e-8);
    ep = zeros(n,1);  ep(k) = h;
    J(:,k) = (resFun(lam+ep) - resFun(lam-ep))/(2*h);
end
end
