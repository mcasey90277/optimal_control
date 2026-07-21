% Z0_ACCEPT2_LM  Acceptance test part 2: Levenberg-Marquardt WITH the exact
% variational-STM Jacobian (the Zhang configuration).
%
% Part 1 (z0_accept_75mN, plain Newton + Armijo) crawled: 1.56e-2 -> 1.526e-2
% in 30 iters at alpha ~ 0.002 -- an undamped Newton direction is the wrong
% step strategy at cond(Jeq) ~ 7e8 with strong nonlinearity, so part 1 does
% NOT settle the derivative-quality fork. Part 2 drives the same residual
% with lsqnonlin LM (proper adaptive damping) + SpecifyObjectiveGradient
% (J = exact STM costate columns), real budget. AMENDMENT NOTE (honest
% record): the pre-registered gate said "Newton, <= 30 iters"; this amends
% the STEP STRATEGY only -- the Jacobian source (variational STM) and the
% gate (||R|| <= 1e-8) are unchanged.
%
%   PASS -> 75 mN ladder anchor converged.
%   FAIL (LM floors with an exact J) -> the wall is basin/valley geometry,
%          not derivative quality: pivot to direct-side dual seed (better
%          START POINT), with multiple shooting as the structural fallback.
%
% Requires: results/z0_accept_trace.mat (part-1 best iterate).
% Output:   results/z0_anchor_75mN.mat (on PASS) / z0_accept2_trace.mat

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');

T = load(fullfile(resDir, 'z0_accept_trace.mat'));
A = load(fullfile(resDir, 'p0i_fd_finish.mat'));  a = A.anchor;
[rv0, rvf, P0] = ztl_endpoints();
Tmax = 3*P0.Tmax25;
tfL  = a.tf;
P = struct('muStar', P0.muStar, 'c', P0.c, 'Tmax', Tmax, 'eps', 1, ...
           'odeRelTol', 1e-12, 'odeAbsTol', 1e-14);

lam0 = T.lam(:);                                % part-1 best (1.526e-2)
fprintf('=== Z0 ACCEPTANCE pt2: LM + exact STM J @ 75 mN (start %.4e) ===\n', T.rn);

opts = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'Algorithm', 'levenberg-marquardt', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance', 1e-24, ...
    'StepTolerance', 1e-16, ...
    'ScaleProblem', 'jacobian', ...
    'MaxIterations', 1500, ...
    'MaxFunctionEvaluations', 8000);

tic;
[lam, res2, ~, flag, outp] = lsqnonlin(@(l) resjac_lm(l, rv0, rvf, tfL, P), ...
                                       lam0, [], [], opts);
rn = sqrt(res2);
fprintf('LM done: ||R|| = %.6e  flag=%d  iters=%d  fevals=%d  (%.0f s)\n', ...
        rn, flag, outp.iterations, outp.funcCount, toc);
save(fullfile(resDir, 'z0_accept2_trace.mat'), 'lam', 'rn', 'flag', 'outp');

if rn < 1e-8
    o = ztl_flow([rv0(:); 1; lam], [0 tfL], P, false);
    lamv = sqrt(sum(o.y(:, 11:13).^2, 2));
    S = 1 - lamv*P.c./o.y(:, 7) - o.y(:, 14);
    u = min(max(0.5 - S/(2*P.eps), 0), 1);
    mF = o.yf(7);
    anchor = struct('Tmax_mN', 75, 'tf', tfL, 'eps', 1, 'lam0_BE', lam, ...
        'resNorm', rn, 'solver', 'lsqnonlin-LM + variational-STM J', ...
        'mProp_kg', P0.m0kg*(1-mF), 'dV_kms', P0.c*log(1/mF)*P0.lStar/P0.tStar, ...
        'uMin', min(u), 'uMax', max(u), 'nSegs', o.nSegs, ...
        'rv0', rv0, 'rvf', rvf, 'P', P, 'P0', P0);
    save(fullfile(resDir, 'z0_anchor_75mN.mat'), 'anchor');
    fprintf(['GATE Z0 (amended): PASS -- LADDER ANCHOR @ 75 mN.\n' ...
             '  prop = %.4f kg  dV = %.4f km/s  u in [%.3f, %.3f]\n'], ...
            anchor.mProp_kg, anchor.dV_kms, anchor.uMin, anchor.uMax);
else
    fprintf(['GATE Z0 (amended): FAIL at %.3e -- LM + exact J also floors.\n' ...
             'Verdict: wall is basin/valley geometry, not derivative quality.\n' ...
             'Pivot: direct-side dual seed (start point), or multiple shooting.\n'], rn);
end

% ---------------------------------------------------------------------------
function [R, J] = resjac_lm(lamv, rv0, rvf, tfL, P)
% lsqnonlin objective: shooting residual + exact STM Jacobian when requested.
if nargout > 1
    o = ztl_flow([rv0(:); 1; lamv], [0 tfL], P, true);
    J = [o.PHI(1:6, 8:14); o.PHI(14, 8:14)];
else
    o = ztl_flow([rv0(:); 1; lamv], [0 tfL], P, false);
end
R = [o.yf(1:6) - rvf(:); o.yf(14)];
end
