% Z1_RESUME_SUN  Warm-restart the Sundman MS trust-region solve toward 1e-8.
% Resumable. Base var moreIter [500].

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');
warning('off','MATLAB:ode89:IntegrationTolNotMet');  warning('off','MATLAB:ode113:IntegrationTolNotMet');

if evalin('base','exist(''moreIter'',''var'')'), moreIter = evalin('base','moreIter'); else, moreIter = 500; end
S = load(fullfile(resDir,'z1_sun_trace.mat'));   % z, out, prob, M, pSund, lam0
z0 = S.z;  prob = S.prob;
if ~isfield(prob,'resFun'), prob.resFun = @ztl_ms_residual_sun; end
Delta0 = 0.5;
if isfield(S,'out') && isfield(S.out,'Delta') && isfinite(S.out.Delta) && S.out.Delta>0, Delta0 = S.out.Delta; end

[~, ~, ri0] = ztl_ms_residual_sun(z0, prob, false);
fprintf('=== Z1-SUN RESUME M=%d pSund=%.2f: start ||R||=%.4e (term=%.2e cont=%.2e) ===\n', ...
        prob.M, S.pSund, norm(ztl_ms_residual_sun(z0,prob,false)), ri0.termErr, ri0.maxCont);

[z, out] = ztl_ms_solve_tr(z0, prob, struct('tolR',1e-9,'maxIter',moreIter,'Delta0',Delta0));
S.z = z;  S.out = out;
save(fullfile(resDir,'z1_sun_trace.mat'), '-struct', 'S');
fprintf('RESUME done: ||R||=%.4e  flag=%d  iters=%d\n', out.resNorm, out.flag, out.iters);

if out.resNorm < 1e-8
    [rv0, rvf, P0] = ztl_endpoints();  P = prob.P;
    lamStar = z(1:7);  tauF = z(end);
    o = ztl_flow_sun([rv0(:);1;lamStar;0], tauF, [0 1], P, false);  mF = o.Yf(7);
    anchor = struct('Tmax_mN',75,'tf',prob.tf,'eps',1,'pSund',S.pSund,'lam0_BE',lamStar, ...
        'tauF',tauF,'resNorm',out.resNorm,'M',prob.M, ...
        'solver','Sundman MS trust-region + geodesic', ...
        'mProp_kg',P0.m0kg*(1-mF),'dV_kms',P0.c*log(1/mF)*P0.lStar/P0.tStar, ...
        'rv0',rv0,'rvf',rvf,'P',P,'P0',P0);
    save(fullfile(resDir,'z1_sun_anchor_75mN.mat'),'anchor');
    fprintf(['GATE Z1-SUN: PASS -- CONVERGED. ||R||=%.3e  prop=%.4f kg  dV=%.4f km/s\n'], ...
            out.resNorm, anchor.mProp_kg, anchor.dV_kms);
else
    fprintf('GATE Z1-SUN: reached %.3e (flag %d); re-run to continue.\n', out.resNorm, out.flag);
end
