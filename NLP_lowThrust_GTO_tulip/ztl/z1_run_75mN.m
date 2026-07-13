% Z1_RUN_75MN  THE Z1 test: converge the 75 mN energy anchor by MULTIPLE
% shooting, where single shooting (Z0 acceptance) floored at 4.95e-3.
%
% Seeds by chopping the banked Z0 single-shooting iterate
% (z0_accept2_trace.mat, ||R|| = 4.95e-3) into M dynamically-consistent arcs,
% then runs ztl_ms_solve (augmented-QR LM on the exact block Jacobian).
%
% GATE Z1: ||R|| <= 1e-8.
%   PASS -> the 75 mN LADDER ANCHOR is converged (first fully converged
%           indirect multi-rev solve of the campaign). lam0 = z(1:7) is a
%           bona-fide single-shooting solution (continuity = 0 at convergence).
%           Z3 marches thrust down from it. Saves results/z1_anchor_75mN.mat.
%   FAIL -> record the floor + trace; try more nodes (M) or perigee-aware
%           node placement before conceding.
%
% Optional: pass M on the command line via `Mnodes` in the base workspace;
% default 26. Requires results/z0_accept2_trace.mat + p0i_fd_finish.mat.

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');

if evalin('base', 'exist(''Mnodes'',''var'')'), M = evalin('base', 'Mnodes'); else, M = 52; end
if evalin('base', 'exist(''solver'',''var'')'), solver = evalin('base', 'solver'); else, solver = 'lm'; end

T = load(fullfile(resDir, 'z0_accept2_trace.mat'));
A = load(fullfile(resDir, 'p0i_fd_finish.mat'));  a = A.anchor;
[rv0, rvf, P0] = ztl_endpoints();
Tmax = 3*P0.Tmax25;  tfL = a.tf;
P = struct('muStar', P0.muStar, 'c', P0.c, 'Tmax', Tmax, 'eps', 1, ...
           'odeRelTol', 1e-13, 'odeAbsTol', 1e-15);
lam0 = T.lam(:);

fprintf('=== Z1: multiple shooting @ 75 mN, eps=1, M=%d arcs, solver=%s ===\n', M, solver);
[z0, prob, si] = ztl_ms_seed(lam0, rv0, rvf, tfL, P, M);
fprintf('seed: maxCont=%.2e (on flow), termErr=%.4e (= SS floor)\n', ...
        si.maxContSeed, si.termErrSeed);

tStart = tic;
switch solver
    case 'tr'
        [z, out] = ztl_ms_solve_tr(z0, prob, struct('tolR', 1e-9, 'maxIter', 200));
    case 'wlm'
        [z, out] = ztl_ms_solve_wlm(z0, prob, struct('tolR', 1e-9, 'maxIter', 200));
    case 'newton'
        [z, out] = ztl_ms_solve_newton(z0, prob, struct('tolR', 1e-9, 'maxIter', 120));
    case 'trr'
        [z, out] = ztl_ms_solve_trr(z0, prob, struct('tolR', 1e-10, 'maxIter', 400));
    case 'lm'
        [z, out] = ztl_ms_solve(z0, prob, struct('tolR', 1e-9, 'maxIter', 300, 'mu0', 1e-6));
end
fprintf('MS solve: ||R||=%.4e  flag=%d  iters=%d  (%.0f s)\n', ...
        out.resNorm, out.flag, out.iters, toc(tStart));

save(fullfile(resDir, 'z1_trace.mat'), 'z', 'out', 'prob', 'M', 'lam0');

% --- verdict + anchor product ------------------------------------------------
lamStar = z(1:7);
o = ztl_flow([rv0(:); 1; lamStar], [0 tfL], P, false);
lamv = sqrt(sum(o.y(:, 11:13).^2, 2));
S = 1 - lamv*P.c./o.y(:, 7) - o.y(:, 14);
u = min(max(0.5 - S/(2*P.eps), 0), 1);
mF = o.yf(7);
% independent single-shooting residual at the extracted lam0
ssRes = norm([o.yf(1:6) - rvf(:); o.yf(14)]);

if out.resNorm < 1e-8
    anchor = struct('Tmax_mN', 75, 'tf', tfL, 'eps', 1, 'lam0_BE', lamStar, ...
        'resNorm', out.resNorm, 'ssResCheck', ssRes, 'M', M, 'iters', out.iters, ...
        'solver', 'multiple-shooting augmented-QR LM + variational STM', ...
        'mProp_kg', P0.m0kg*(1-mF), 'dV_kms', P0.c*log(1/mF)*P0.lStar/P0.tStar, ...
        'uMin', min(u), 'uMax', max(u), 'rv0', rv0, 'rvf', rvf, 'P', P, 'P0', P0);
    save(fullfile(resDir, 'z1_anchor_75mN.mat'), 'anchor');
    fprintf(['\nGATE Z1: PASS -- 75 mN LADDER ANCHOR CONVERGED (M=%d, %d iters).\n' ...
             '  ||R||_MS = %.3e   single-shooting check = %.3e\n' ...
             '  prop = %.4f kg   dV = %.4f km/s   u in [%.3f, %.3f]\n' ...
             '  saved %s\n'], M, out.iters, out.resNorm, ssRes, ...
            anchor.mProp_kg, anchor.dV_kms, anchor.uMin, anchor.uMax, ...
            fullfile(resDir, 'z1_anchor_75mN.mat'));
else
    fprintf(['\nGATE Z1: FAIL -- floors at %.3e (flag %d) with M=%d.\n' ...
             '  Try more nodes (set Mnodes) or perigee-aware placement.\n'], ...
            out.resNorm, out.flag, M);
end
