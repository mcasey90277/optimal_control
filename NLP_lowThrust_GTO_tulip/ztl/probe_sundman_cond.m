% PROBE_SUNDMAN_COND  The decisive S5 test: cond of the COLUMN-EQUILIBRATED
% Sundman MS Jacobian vs the physical one (the solver column-scales J, so this
% is the conditioning that sets the noise floor). Also S3 (block J vs FD) and
% S4 (seed consistency).

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
warning('off','MATLAB:ode89:IntegrationTolNotMet'); warning('off','MATLAB:ode113:IntegrationTolNotMet');
A = load('results/z1_anchor_75mN.mat');  a = A.anchor;
[rv0, rvf, P0] = ztl_endpoints();  Tmax = 3*P0.Tmax25;  tfL = a.tf;  lam0 = a.lam0_BE(:);
colcond = @(J) cond(J ./ (sqrt(sum(J.^2,1)).' + (sqrt(sum(J.^2,1)).'==0)).');

% --- physical reference (perigee1 nodes, M=104) ---
Pp = struct('muStar',P0.muStar,'c',P0.c,'Tmax',Tmax,'eps',1,'odeRelTol',1e-13,'odeAbsTol',1e-15);
[zp, pp] = ztl_ms_seed(lam0, rv0, rvf, tfL, Pp, 104, 'perigee1');
[~, Jp] = ztl_ms_residual(zp, pp, true);
fprintf('PHYSICAL (perigee1, M=104): cond(col-scaled J) = %.2e\n', colcond(Jp));

% --- S3 block-Jacobian correctness (Sundman, small M) ---
Ps = struct('muStar',P0.muStar,'c',P0.c,'Tmax',Tmax,'eps',1,'pSund',2,'odeRelTol',1e-12,'odeAbsTol',1e-14);
[z6, p6, si6] = ztl_ms_seed_sun(lam0, rv0, rvf, tfL, Ps, 6);
[~, J6] = ztl_ms_residual_sun(z6, p6, true);
n6 = numel(z6);  Jfd = zeros(n6);
for k = 1:n6
    h = 1e-7*max(1,abs(z6(k)));  e=zeros(n6,1); e(k)=h;
    Jfd(:,k) = (ztl_ms_residual_sun(z6+e,p6,false) - ztl_ms_residual_sun(z6-e,p6,false))/(2*h);
end
fprintf('S3 block J vs FD (Sundman, M=6): rel err = %.2e\n', norm(J6-Jfd)/norm(Jfd));
fprintf('S4 seed (M=6): maxCont=%.2e  termErr=%.2e  tauF=%.3f\n', si6.maxContSeed, si6.termErrSeed, si6.tauF);

% --- S5 the decisive cond comparison (Sundman, M=104) ---
for pS = [1.5 2]
    Ps = struct('muStar',P0.muStar,'c',P0.c,'Tmax',Tmax,'eps',1,'pSund',pS,'odeRelTol',1e-13,'odeAbsTol',1e-15);
    [zs, ps, sis] = ztl_ms_seed_sun(lam0, rv0, rvf, tfL, Ps, 104);
    [~, Js] = ztl_ms_residual_sun(zs, ps, true);
    fprintf('SUNDMAN pS=%.1f (M=104): cond(col-scaled J) = %.2e   seed maxCont=%.2e termErr=%.2e\n', ...
            pS, colcond(Js), sis.maxContSeed, sis.termErrSeed);
end
