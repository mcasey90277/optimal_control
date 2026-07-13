% Z1_RESUME_TR  Warm-restart the trust-region MS solve from the saved iterate
% and grind toward 1e-8. Resumable: re-run to continue further.
%
% The TR run (z1_run_75mN, solver 'tr') converges monotonically but at ~1%/it
% (nonlinearity keeps Delta small); it hit the iteration cap at 1.4e-5 still
% descending. This continues from results/z1_trace.mat with a large budget and
% a warm Delta, and re-saves the iterate each call so progress is never lost.
%
% Sets/uses base var `moreIter` (default 600).

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');

if evalin('base', 'exist(''moreIter'',''var'')'), moreIter = evalin('base','moreIter'); else, moreIter = 600; end

S = load(fullfile(resDir, 'z1_trace.mat'));      % z, out, prob, M, lam0
z0 = S.z;  prob = S.prob;
Delta0 = 0.02;
if isfield(S,'out') && isfield(S.out,'Delta') && isfinite(S.out.Delta), Delta0 = S.out.Delta; end

[~, ~, ri0] = ztl_ms_residual(z0, prob, false);
fprintf('=== Z1 RESUME (TR) M=%d: start ||R||=%.4e (termErr=%.2e maxCont=%.2e), Delta0=%.2e ===\n', ...
        prob.M, norm(ztl_ms_residual(z0, prob, false)), ri0.termErr, ri0.maxCont, Delta0);

[z, out] = ztl_ms_solve_tr(z0, prob, ...
    struct('tolR', 1e-9, 'maxIter', moreIter, 'Delta0', Delta0));

S.z = z;  S.out = out;                            % overwrite trace with progress
save(fullfile(resDir, 'z1_trace.mat'), '-struct', 'S');
fprintf('RESUME done: ||R||=%.4e  flag=%d  iters=%d\n', out.resNorm, out.flag, out.iters);

% --- if converged, save the anchor product ---------------------------------
if out.resNorm < 1e-8
    [rv0, rvf, P0] = ztl_endpoints();
    P = prob.P;  tfL = prob.tNodes(end);
    lamStar = z(1:7);
    o = ztl_flow([rv0(:); 1; lamStar], [0 tfL], P, false);
    lamv = sqrt(sum(o.y(:,11:13).^2, 2));
    Sv = 1 - lamv*P.c./o.y(:,7) - o.y(:,14);
    u = min(max(0.5 - Sv/(2*P.eps), 0), 1);
    mF = o.yf(7);
    ssRes = norm([o.yf(1:6)-rvf(:); o.yf(14)]);
    anchor = struct('Tmax_mN', 75, 'tf', tfL, 'eps', 1, 'lam0_BE', lamStar, ...
        'resNorm', out.resNorm, 'ssResCheck', ssRes, 'M', prob.M, ...
        'solver', 'multiple-shooting SVD trust-region + variational STM', ...
        'mProp_kg', P0.m0kg*(1-mF), 'dV_kms', P0.c*log(1/mF)*P0.lStar/P0.tStar, ...
        'uMin', min(u), 'uMax', max(u), 'rv0', rv0, 'rvf', rvf, 'P', P, 'P0', P0);
    save(fullfile(resDir, 'z1_anchor_75mN.mat'), 'anchor');
    fprintf(['GATE Z1: PASS -- 75 mN LADDER ANCHOR CONVERGED.\n' ...
             '  ||R||=%.3e  SS-check=%.3e  prop=%.4f kg  dV=%.4f km/s  u in [%.3f,%.3f]\n'], ...
            out.resNorm, ssRes, anchor.mProp_kg, anchor.dV_kms, anchor.uMin, anchor.uMax);
else
    fprintf('GATE Z1: not yet (||R||=%.3e); re-run z1_resume_tr to continue.\n', out.resNorm);
end
