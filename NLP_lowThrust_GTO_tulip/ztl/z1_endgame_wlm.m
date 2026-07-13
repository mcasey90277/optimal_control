% Z1_ENDGAME_WLM  Endgame refinement from the trust-region iterate: a FIXED
% consistently-weighted LM to break the column-SVD conditioning noise floor
% (~5e-7 continuity) that the trust region asymptotes on.
%
% The TR (column-scaled, objective-safe) floors at the linear-solve noise of
% cond(Jc)~1e9. Two-sided equilibration would fix the conditioning but changes
% the objective (ascent) FAR from the solution. From the TR iterate (~1e-6,
% terminal already ~1e-9) we are NEAR the solution, so a fixed-weight solve
% (Dr = row equilibration, held constant -> a legitimate objective since
% ||Dr R||->0 iff ||R||->0) has a short, safe path and should reach 1e-8/1e-9.
%
% Loads results/z1_trace.mat; runs ztl_ms_solve_wlm; saves the anchor on
% success.

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');

S = load(fullfile(resDir, 'z1_trace.mat'));
z0 = S.z;  prob = S.prob;
[~, ~, ri0] = ztl_ms_residual(z0, prob, false);
fprintf('=== Z1 ENDGAME (fixed-W LM) M=%d: start ||R||=%.4e (term=%.2e cont=%.2e) ===\n', ...
        prob.M, norm(ztl_ms_residual(z0, prob, false)), ri0.termErr, ri0.maxCont);

[z, out] = ztl_ms_solve_wlm(z0, prob, struct('tolR', 1e-9, 'maxIter', 150, 'mu0', 1e-8));
fprintf('ENDGAME done: ||R||=%.4e  flag=%d  iters=%d\n', out.resNorm, out.flag, out.iters);

S.z = z;  S.out = out;
save(fullfile(resDir, 'z1_trace.mat'), '-struct', 'S');

if out.resNorm < 1e-8
    [rv0, rvf, P0] = ztl_endpoints();  P = prob.P;  tfL = prob.tNodes(end);
    lamStar = z(1:7);
    o = ztl_flow([rv0(:); 1; lamStar], [0 tfL], P, false);
    lamv = sqrt(sum(o.y(:,11:13).^2, 2));
    Sv = 1 - lamv*P.c./o.y(:,7) - o.y(:,14);  u = min(max(0.5 - Sv/(2*P.eps), 0), 1);
    mF = o.yf(7);  ssRes = norm([o.yf(1:6)-rvf(:); o.yf(14)]);
    anchor = struct('Tmax_mN', 75, 'tf', tfL, 'eps', 1, 'lam0_BE', lamStar, ...
        'resNorm', out.resNorm, 'ssResCheck', ssRes, 'M', prob.M, ...
        'solver', 'MS trust-region (geodesic) + fixed-W LM endgame', ...
        'mProp_kg', P0.m0kg*(1-mF), 'dV_kms', P0.c*log(1/mF)*P0.lStar/P0.tStar, ...
        'uMin', min(u), 'uMax', max(u), 'rv0', rv0, 'rvf', rvf, 'P', P, 'P0', P0);
    save(fullfile(resDir, 'z1_anchor_75mN.mat'), 'anchor');
    fprintf(['GATE Z1: PASS -- 75 mN LADDER ANCHOR CONVERGED.\n' ...
             '  ||R||=%.3e  SS-check=%.3e  prop=%.4f kg  dV=%.4f km/s  u in [%.3f,%.3f]\n'], ...
            out.resNorm, ssRes, anchor.mProp_kg, anchor.dV_kms, anchor.uMin, anchor.uMax);
else
    fprintf('GATE Z1: endgame reached %.3e (flag %d).\n', out.resNorm, out.flag);
end
