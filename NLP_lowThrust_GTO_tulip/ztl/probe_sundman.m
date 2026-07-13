% PROBE_SUNDMAN  Core hypothesis test BEFORE building the MS layer: does
% Sundman regularization lower the per-arc STM amplification (=> cond)?
% Also gates ztl_A_sun vs FD (S2) and finds tauF (t(tau_f)=tf).

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
warning('off', 'MATLAB:ode89:IntegrationTolNotMet');
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');

A = load('results/z1_anchor_75mN.mat');  a = A.anchor;
[rv0, rvf, P0] = ztl_endpoints();  Tmax = 3*P0.Tmax25;  tfL = a.tf;
lam0 = a.lam0_BE(:);
for pS = [1.5 2 3]
    P = struct('muStar',P0.muStar,'c',P0.c,'Tmax',Tmax,'eps',1, ...
               'pSund',pS,'odeRelTol',1e-13,'odeAbsTol',1e-15);

    % --- S2: ztl_A_sun vs FD-of-field (once, pS=1.5) ---
    if pS == 1.5
        Y0 = [rv0(:); 1; lam0; 0];
        Acs = ztl_A_sun(Y0, P, 'medium');  Afd = zeros(15);
        for k = 1:15
            h = 1e-7*max(1,abs(Y0(k)));  e = zeros(15,1);  e(k)=h;
            Afd(:,k) = (ztl_eom_sun(Y0+e,P,'medium') - ztl_eom_sun(Y0-e,P,'medium'))/(2*h);
        end
        fprintf('S2 ztl_A_sun vs FD: rel err = %.2e\n', norm(Acs-Afd)/norm(Afd));
    end

    % --- find tauF: integrate dY/dtau until t = tf ---
    Y0 = [rv0(:); 1; lam0; 0];
    oo = odeset('RelTol',1e-13,'AbsTol',1e-15,'Events',@(tau,Y) tf_event(tau,Y,tfL));
    [~,~,tauE] = ode89(@(tau,Y) ztl_eom_sun(Y,P,'medium'), [0 1e6], Y0, oo);
    tauF = tauE(1);

    % --- Sundman per-arc ||Phi|| (uniform sigma) at M=104 ---
    M = 104;  sN = linspace(0,1,M+1);  Yk = Y0;  ampS = zeros(1,M);
    for k = 1:M
        o = ztl_flow_sun(Yk, tauF, [sN(k) sN(k+1)], P, true);
        ampS(k) = norm(o.PHI);  Yk = o.Yf;
    end
    fprintf('pS=%.1f: tauF=%.3f  Sundman worst||Phi||=%.2e  med=%.2e  max/med=%.0f\n', ...
            pS, tauF, max(ampS), median(ampS), max(ampS)/median(ampS));
end

% --- physical per-arc ||Phi|| (uniform t) at M=104 for reference ---
Pp = struct('muStar',P0.muStar,'c',P0.c,'Tmax',Tmax,'eps',1,'odeRelTol',1e-13,'odeAbsTol',1e-15);
M = 104;  tN = linspace(0,tfL,M+1);  Yk = [rv0(:);1;lam0];  ampP = zeros(1,M);
for k = 1:M
    o = ztl_flow(Yk, [tN(k) tN(k+1)], Pp, true);  ampP(k)=norm(o.PHI);  Yk=o.Yf;
end
fprintf('PHYSICAL (uniform t): worst||Phi||=%.2e  med=%.2e  max/med=%.0f\n', ...
        max(ampP), median(ampP), max(ampP)/median(ampP));

function [v,isterm,dir] = tf_event(~, Y, tfL)
v = Y(15) - tfL;  isterm = 1;  dir = 1;
end
