% P0H_GN_FINISH  Preflight: finish the 75 mN energy solve with an
% equilibrated truncated-SVD Gauss-Newton step (the anti-crawl stepper).
%
% P0g's warm restarts descended monotonically (75 mN: 1.42 -> 0.379; 50 mN:
% 6.05 -> 0.99) but LM's damping crawls at ~0.2%/iter -- the same signature
% ifs_solve2's equilibrated truncated-SVD GN step broke in the IFS work.
% Here the system is only 7 unknowns (J is 7x7), so the tsvd step is
% trivial: column-equilibrate J, SVD, pseudo-inverse step with adaptive
% relative truncation, Armijo backtracking on ||R||, LM-style fallback if a
% step is rejected at all alphas.
%
% Start: P0g's best iterate (75 mN, ||R|| = 0.379). Gate: ||R|| <= 1e-8.
%
% Requires: results/p0g_warm_restart.mat.
% Output:   results/p0h_gn_finish.mat

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');

G = load(fullfile(resDir, 'p0g_warm_restart.mat'));
[rv0, rvf, P] = ztl_endpoints();
tfL = G.tfL;

[~, kb] = min([G.runs.resNorm]);
r0   = G.runs(kb);
Tmax = (r0.Tmax_mN/25)*P.Tmax25;
fprintf('=== P0h: tsvd-GN finish @ %g mN from ||R|| = %.3e ===\n', ...
        r0.Tmax_mN, r0.resNorm);

warning('off', 'MATLAB:ode113:IntegrationTolNotMet');
resFun = @(lam) shoot_residual_energy(lam, tfL, rv0, 1, rvf, Tmax, P.c, P.muStar);

lam = r0.lamSol(:);
[R, J] = resFun(lam);
rn = norm(R);
maxIter  = 400;
relTrunc = 1e-10;                       % relative SVD truncation (adaptive)
hist = nan(maxIter, 3);
stallCount = 0;

for it = 1:maxIter
    % column equilibration (unknown scales differ by ~1e2-1e3)
    cs = max(abs(J), [], 1);  cs(cs == 0) = 1;
    Js = J ./ cs;

    [U, S, V] = svd(Js);
    sv = diag(S);
    keep = sv >= relTrunc*sv(1);
    dz = -(V(:, keep) * ((U(:, keep)'*R) ./ sv(keep))) ./ cs(:);

    % Armijo backtracking on ||R||
    alpha = 1;  accepted = false;
    for ls = 1:12
        lamT = lam + alpha*dz;
        Rt = resFun(lamT);
        if all(isfinite(Rt)) && norm(Rt) < rn*(1 - 1e-4*alpha)
            accepted = true;  break
        end
        alpha = alpha/2;
    end

    if accepted
        lam = lamT;  R = Rt;  rnNew = norm(R);
        prog = 1 - rnNew/rn;  rn = rnNew;
        [~, J] = resFun(lam);                       % J at the new point
        if prog < 1e-3, stallCount = stallCount + 1; else, stallCount = 0; end
    else
        % no acceptable step: tighten truncation (drop weakest directions)
        relTrunc = relTrunc*100;
        stallCount = stallCount + 1;
        fprintf('  it %3d: step rejected at all alpha -> relTrunc = %.1e\n', it, relTrunc);
    end

    hist(it, :) = [rn, alpha*accepted, cond(Js)];
    if mod(it, 10) == 0 || rn < 1e-8
        fprintf('  it %3d: ||R|| = %.6e  alpha = %.3g  cond(Js) = %.2e  trunc = %.0e\n', ...
                it, rn, alpha*accepted, cond(Js), relTrunc);
    end
    if rn < 1e-8, break; end
    if stallCount >= 20
        fprintf('  20 consecutive stall/reject iterations -- stopping.\n');
        break
    end
    if relTrunc > 1e-2 && stallCount > 5
        relTrunc = 1e-10;                            % relax again after escape
    end
end

% --- accounting -------------------------------------------------------------
optsInt = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, yI] = ode113(@lt_pmp_eom_energy, [0 tfL], [rv0(:); 1; lam], ...
                 optsInt, Tmax, P.c, P.muStar);
lamvMag = sqrt(sum(yI(:,11:13).^2, 2));
u  = min(max(Tmax*(lamvMag./yI(:,7) + yI(:,14)/P.c), 0), 1);
mF = yI(end,7);
anchor = struct('Tmax_mN', r0.Tmax_mN, 'tf', tfL, 'lam0', lam, 'resNorm', rn, ...
    'iters', it, 'mProp_kg', P.m0kg*(1-mF), ...
    'dV_kms', P.c*log(1/mF)*P.lStar/P.tStar, 'uMin', min(u), 'uMax', max(u), ...
    'fracSatHi', mean(u > 0.999), 'fracSatLo', mean(u < 1e-3), ...
    'rv0', rv0, 'rvf', rvf, 'P', P, 'hist', hist);
save(fullfile(resDir, 'p0h_gn_finish.mat'), 'anchor');
fprintf('saved %s\n', fullfile(resDir, 'p0h_gn_finish.mat'));

fprintf(['\nP0h RESULT @ %g mN: ||R|| = %.3e after %d iters  prop = %.4f kg  ' ...
         'dV = %.4f km/s\n  throttle: min %.3f  max %.3f  sat-hi %.1f%%  sat-lo %.1f%%\n'], ...
    anchor.Tmax_mN, rn, it, anchor.mProp_kg, anchor.dV_kms, ...
    anchor.uMin, anchor.uMax, 100*anchor.fracSatHi, 100*anchor.fracSatLo);
if rn < 1e-8
    fprintf('GATE P0h: PASS -- LADDER ANCHOR CONVERGED at %g mN.\n', anchor.Tmax_mN);
else
    fprintf('GATE P0h: FAIL -- tsvd-GN also floors (at %.3e); next = direct-side seed.\n', rn);
end
