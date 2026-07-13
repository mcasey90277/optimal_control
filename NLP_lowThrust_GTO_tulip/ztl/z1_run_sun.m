% Z1_RUN_SUN  Sundman-regularized multiple shooting for the 75 mN energy anchor
% -- the structural fix for the conditioning noise floor (cond(J) 9.5e8 -> 5.3e7
% at pSund=1.5, SUN_BUILD.md / ZTL_RESULTS). Seeds from the CLEAN Z0 single-
% shooting lam0 (the physical MS anchor's lam0 has a large single-shooting
% residual -- load-bearing gaps), builds the Sundman MS seed, and runs the same
% trust-region + geodesic solver via prob.resFun = @ztl_ms_residual_sun.
%
% GATE Z1-Sun: ||R|| <= 1e-8. Base vars: Mnodes[104], pSund[1.5].

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');
warning('off','MATLAB:ode89:IntegrationTolNotMet');  warning('off','MATLAB:ode113:IntegrationTolNotMet');

if evalin('base','exist(''Mnodes'',''var'')'), M = evalin('base','Mnodes'); else, M = 104; end
if evalin('base','exist(''pSund'',''var'')'), pSund = evalin('base','pSund'); else, pSund = 1.5; end

T = load(fullfile(resDir,'z0_accept2_trace.mat'));   % clean Z0 lam0 (SS res 4.95e-3)
[rv0, rvf, P0] = ztl_endpoints();  Tmax = 3*P0.Tmax25;  tfL = 1.15*6.29081541876621;
lam0 = T.lam(:);
P = struct('muStar',P0.muStar,'c',P0.c,'Tmax',Tmax,'eps',1,'pSund',pSund, ...
           'odeRelTol',1e-13,'odeAbsTol',1e-15);

fprintf('=== Z1-SUN: Sundman MS @ 75 mN, eps=1, M=%d, pSund=%.2f ===\n', M, pSund);
[z0, prob, si] = ztl_ms_seed_sun(lam0, rv0, rvf, tfL, P, M);
prob.resFun = @ztl_ms_residual_sun;
fprintf('seed: maxCont=%.2e (on flow), termErr=%.4e, tauF=%.3f\n', ...
        si.maxContSeed, si.termErrSeed, si.tauF);

[z, out] = ztl_ms_solve_tr(z0, prob, struct('tolR',1e-9,'maxIter',300));
fprintf('SUN solve: ||R||=%.4e  flag=%d  iters=%d\n', out.resNorm, out.flag, out.iters);
save(fullfile(resDir,'z1_sun_trace.mat'), 'z','out','prob','M','pSund','lam0');

lamStar = z(1:7);  tauF = z(end);
o = ztl_flow_sun([rv0(:);1;lamStar;0], tauF, [0 1], P, false);
mF = o.Yf(7);
if out.resNorm < 1e-8
    anchor = struct('Tmax_mN',75,'tf',tfL,'eps',1,'pSund',pSund,'lam0_BE',lamStar, ...
        'tauF',tauF,'resNorm',out.resNorm,'M',M, ...
        'solver','Sundman MS trust-region + geodesic', ...
        'mProp_kg',P0.m0kg*(1-mF),'dV_kms',P0.c*log(1/mF)*P0.lStar/P0.tStar, ...
        'rv0',rv0,'rvf',rvf,'P',P,'P0',P0);
    save(fullfile(resDir,'z1_sun_anchor_75mN.mat'),'anchor');
    fprintf(['GATE Z1-SUN: PASS -- 75 mN ANCHOR CONVERGED (Sundman).\n' ...
             '  ||R||=%.3e  prop=%.4f kg  dV=%.4f km/s\n  saved %s\n'], ...
            out.resNorm, anchor.mProp_kg, anchor.dV_kms, fullfile(resDir,'z1_sun_anchor_75mN.mat'));
else
    fprintf('GATE Z1-SUN: reached %.3e (flag %d).\n', out.resNorm, out.flag);
end
