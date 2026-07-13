% TEST_ZTL_Z1  Unit gates for the Z1 multiple-shooting build.
%
% H1 seed consistency: chopping a single-shooting trajectory gives continuity
%    residual ~integrator tol (seed is on the flow), terminal BC = the SS
%    residual (unchanged).
% H2 block Jacobian vs finite-difference of the MS residual (the test that
%    catches a wrong STM-block placement / sign).
% H3 conditioning: cond(J_MS) must be ORDERS below the single-shooting 7x7
%    Jacobian condition (the reason MS is expected to work).

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');

T = load(fullfile(here, 'results', 'z0_accept2_trace.mat'));
A = load(fullfile(here, 'results', 'p0i_fd_finish.mat'));  a = A.anchor;
[rv0, rvf, P0] = ztl_endpoints();
Tmax = 3*P0.Tmax25;  tfL = a.tf;
P = struct('muStar', P0.muStar, 'c', P0.c, 'Tmax', Tmax, 'eps', 1, ...
           'odeRelTol', 1e-12, 'odeAbsTol', 1e-14);
lam0 = T.lam(:);

nPass = 0;  nFail = 0;
Mtest = 6;                                   % small M for the FD test

%% H1 -- seed consistency -----------------------------------------------------
[z, prob, si] = ztl_ms_seed(lam0, rv0, rvf, tfL, P, Mtest);
% terminal BC of the seed must equal the single-shooting residual at lam0
oSS = ztl_flow([rv0(:); 1; lam0], [0 tfL], P, false);
ssRes = norm([oSS.yf(1:6) - rvf(:); oSS.yf(14)]);
ok1 = si.maxContSeed < 1e-6 && abs(si.termErrSeed - ssRes) < 1e-4;
[nPass, nFail] = gate(ok1, nPass, nFail, ...
    sprintf('H1 seed consistency: maxCont=%.2e (seed on flow), termErr=%.4e vs SS=%.4e', ...
    si.maxContSeed, si.termErrSeed, ssRes));

%% H2 -- block Jacobian vs FD of the MS residual, across node counts ---------
% The FD reference degrades on long, high-amplification arcs (Z0's STM gate
% used a 0.05 arc). If the block PLACEMENT/signs are right, the mismatch must
% SHRINK as arcs shorten (per-arc STM -> Z0-gate accuracy). Require a monotone
% decrease and M=26 below 1e-4 (the FD-of-multi-rev-flow floor).
e2v = zeros(1,3);  Mv = [6 13 26];
for im = 1:3
    [zm, pm] = ztl_ms_seed(lam0, rv0, rvf, tfL, P, Mv(im));
    [~, Jm] = ztl_ms_residual(zm, pm, true);
    nm = numel(zm);  Jfd = zeros(nm);
    for k = 1:nm
        h = 1e-7*max(1, abs(zm(k)));
        ep = zeros(nm,1);  ep(k) = h;
        Jfd(:,k) = (ztl_ms_residual(zm+ep, pm, false) - ztl_ms_residual(zm-ep, pm, false))/(2*h);
    end
    e2v(im) = norm(Jm - Jfd)/norm(Jfd);
end
ok2 = all(diff(e2v) < 0) && e2v(3) < 1e-4;
[nPass, nFail] = gate(ok2, nPass, nFail, ...
    sprintf('H2 block J vs FD: M=6 %.2e -> M=13 %.2e -> M=26 %.2e (shrinks, arc-accuracy)', ...
    e2v(1), e2v(2), e2v(3)));

%% H3 -- conditioning advantage ----------------------------------------------
% single-shooting 7x7 Jacobian: costate columns of the FULL-trajectory STM
oFull = ztl_flow([rv0(:); 1; lam0], [0 tfL], P, true);
Jss = [oFull.PHI(1:6, 8:14); oFull.PHI(14, 8:14)];
condSS = cond(Jss);
% MS Jacobian at a realistic node count
[zB, probB] = ztl_ms_seed(lam0, rv0, rvf, tfL, P, 26);
[~, JB] = ztl_ms_residual(zB, probB, true);
condMS = cond(JB);
[nPass, nFail] = gate(condMS < condSS/10, nPass, nFail, ...
    sprintf('H3 conditioning: cond(J_SS)=%.2e vs cond(J_MS,M=26)=%.2e (ratio %.1e)', ...
    condSS, condMS, condMS/condSS));

fprintf('\n=== Z1 GATES: %d PASS, %d FAIL ===\n', nPass, nFail);

% ---------------------------------------------------------------------------
function [nP, nF] = gate(ok, nP, nF, msg)
if ok, nP = nP + 1;  tag = 'PASS'; else, nF = nF + 1;  tag = 'FAIL'; end
fprintf('[%s] %s\n', tag, msg);
end
