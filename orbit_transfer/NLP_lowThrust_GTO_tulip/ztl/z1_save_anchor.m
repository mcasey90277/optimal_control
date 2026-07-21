% Z1_SAVE_ANCHOR  Save the Z1 anchor product from the current z1_trace iterate.
% The MS trust-region (perigee + geodesic) converged the 75 mN energy anchor to
% ||R|| ~ 1.4e-6 (terminal BC to ~2e-9); below that a conditioning floor blocks
% (see ZTL_RESULTS Z1). This is a legitimate near-extremal and a usable Z3
% ladder anchor -- save it as results/z1_anchor_75mN.mat.

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');

S = load(fullfile(resDir, 'z1_trace.mat'));
z = S.z;  prob = S.prob;
[R, ~, ri] = ztl_ms_residual(z, prob, true);
rn = norm(R);

[rv0, rvf, P0] = ztl_endpoints();  P = prob.P;  tfL = prob.tNodes(end);
lamStar = z(1:7);
o = ztl_flow([rv0(:); 1; lamStar], [0 tfL], P, false);
lamv = sqrt(sum(o.y(:,11:13).^2, 2));
Sv = 1 - lamv*P.c./o.y(:,7) - o.y(:,14);
u = min(max(0.5 - Sv/(2*P.eps), 0), 1);
mF = o.yf(7);
ssRes = norm([o.yf(1:6)-rvf(:); o.yf(14)]);

anchor = struct('Tmax_mN', 75, 'tf', tfL, 'eps', 1, 'lam0_BE', lamStar, ...
    'resNorm', rn, 'termErr', ri.termErr, 'maxCont', ri.maxCont, ...
    'ssResCheck', ssRes, 'M', prob.M, ...
    'solver', 'MS SVD trust-region + perigee nodes + geodesic acceleration', ...
    'converged', rn < 1e-8, 'mProp_kg', P0.m0kg*(1-mF), ...
    'dV_kms', P0.c*log(1/mF)*P0.lStar/P0.tStar, ...
    'uMin', min(u), 'uMax', max(u), 'nodeStates', {reshape([[rv0(:);1;lamStar], ...
        reshape(z(8:end), 14, [])], 14, [])}, ...
    'z', z, 'prob_tNodes', prob.tNodes, 'rv0', rv0, 'rvf', rvf, 'P', P, 'P0', P0);
save(fullfile(resDir, 'z1_anchor_75mN.mat'), 'anchor');

fprintf('=== Z1 ANCHOR (75 mN, eps=1) saved ===\n');
fprintf('  ||R|| = %.3e   termErr = %.2e   maxCont = %.2e\n', rn, ri.termErr, ri.maxCont);
fprintf('  prop = %.4f kg   dV = %.4f km/s   throttle u in [%.3f, %.3f] (interior)\n', ...
        anchor.mProp_kg, anchor.dV_kms, anchor.uMin, anchor.uMax);
fprintf('  saved %s\n', fullfile(resDir, 'z1_anchor_75mN.mat'));
