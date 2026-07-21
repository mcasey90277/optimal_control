% TEST_ZTL_Z0  Unit gates G1-G5 for the Z0 build (Z0_BUILD.md par.7).
%
% G1 legacy equivalence at eps=1 (EOM + automaton + costate map, vs
%    lt_pmp_eom_energy as oracle)
% G2 ztl_A vs central-difference OF THE FIELD (field differencing is valid)
% G3 variational STM vs central-difference OF THE FLOW on a short arc
% G4 event integrity at eps=0.5 (reproducibility, S at events, u continuity)
% G5 saltation-corrected STM across one eps=0 switch vs FD of the flow

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');

A = load(fullfile(here, 'results', 'p0i_fd_finish.mat'));
a = A.anchor;
[rv0, rvf, P0] = ztl_endpoints();
Tmax = 3*P0.Tmax25;                            % 75 mN
tfL  = a.tf;

P = struct('muStar', P0.muStar, 'c', P0.c, 'Tmax', Tmax, 'eps', 1, ...
           'odeRelTol', 1e-13, 'odeAbsTol', 1e-15);
kMap  = 2*Tmax/P0.c;                           % legacy -> BE costate map
lamBE = kMap*a.lam0(:);
y0BE  = [rv0(:); 1; lamBE];

nPass = 0;  nFail = 0;

%% G1 -- legacy equivalence at eps=1 -----------------------------------------
% throttle agreement at t=0
lamL = a.lam0(:);
uLeg = min(max(Tmax*(sqrt(sum(lamL(4:6).^2))/1 + lamL(7)/P0.c), 0), 1);
[~, aux0] = ztl_eom(y0BE, P, 'medium');
uBE = min(max(real(aux0.u), 0), 1);
e1a = abs(uLeg - uBE);

% terminal-state agreement over the full 13-rev arc. Oracle uses the SAME
% integrator (ode89) as ztl_flow: diag_g1_integrator_floor measured the
% cross-integrator floor ode113-vs-ode89 at 9.6e-7 over 13 revs, while
% like-for-like agreement is 3.1e-8 (AbsTol-on-rescaled-costates floor).
o1 = ztl_flow(y0BE, [0 tfL], P, false);
optsInt = odeset('RelTol', 1e-13, 'AbsTol', 1e-15);
[~, yL] = ode89(@(t,y) lt_pmp_eom_energy(t, y, Tmax, P0.c, P0.muStar), ...
                [0 tfL], [rv0(:); 1; lamL], optsInt);
e1b = max(abs(o1.yf(1:7) - yL(end, 1:7).'));
[nPass, nFail] = gate(e1a < 1e-12 && e1b < 1e-7, nPass, nFail, ...
    sprintf('G1 legacy equivalence (ode89 oracle): |du0|=%.2e, term-state err=%.2e (segs=%d)', e1a, e1b, o1.nSegs));

%% G2 -- ztl_A vs FD of the field ---------------------------------------------
e2 = 0;
for rgc = {'on', 'off', 'medium'}
    rg = rgc{1};
    Acs = ztl_A(y0BE, P, rg);
    Afd = zeros(14);
    for k = 1:14
        h = 1e-7*max(1, abs(y0BE(k)));
        ep = zeros(14,1);  ep(k) = h;
        Afd(:,k) = (ztl_eom(y0BE+ep, P, rg) - ztl_eom(y0BE-ep, P, rg))/(2*h);
    end
    e2 = max(e2, norm(Acs - Afd)/norm(Afd));
end
[nPass, nFail] = gate(e2 < 1e-6, nPass, nFail, ...
    sprintf('G2 A exactness (worst regime): rel err = %.2e', e2));

%% G3 -- variational STM vs FD of the flow (short arc, no events) --------------
tShort = 0.05;
oS = ztl_flow(y0BE, [0 tShort], P, true);
PhiFD = zeros(14);
for k = 1:14
    h = 1e-7*max(1, abs(y0BE(k)));
    ep = zeros(14,1);  ep(k) = h;
    op = ztl_flow(y0BE+ep, [0 tShort], P, false);
    om = ztl_flow(y0BE-ep, [0 tShort], P, false);
    PhiFD(:,k) = (op.yf - om.yf)/(2*h);
end
e3 = norm(oS.PHI - PhiFD)/norm(PhiFD);
[nPass, nFail] = gate(e3 < 1e-6, nPass, nFail, ...
    sprintf('G3 variational STM (0.05 arc): rel err = %.2e (segs=%d)', e3, oS.nSegs));

%% G4 -- event integrity at eps=0.5 --------------------------------------------
P4 = P;  P4.eps = 0.5;
o4a = ztl_flow(y0BE, [0 tfL], P4, false);
P4b = P4;  P4b.odeRelTol = 1e-12;  P4b.odeAbsTol = 1e-14;
o4b = ztl_flow(y0BE, [0 tfL], P4b, false);
if isempty(o4a.events)
    [nPass, nFail] = gate(false, nPass, nFail, 'G4: NO events at eps=0.5 (unexpected)');
else
    nCmp = min(numel(o4a.events), numel(o4b.events));
    dtEv = max(abs([o4a.events(1:nCmp).t] - [o4b.events(1:nCmp).t]));
    dS   = max(abs(abs([o4a.events.S]) - P4.eps));
    % u continuity across each event
    duMax = 0;
    for kev = 1:numel(o4a.events)
        yev = o4a.events(kev).yEv;
        [~, auxA] = ztl_eom(yev, P4, o4a.events(kev).from);
        [~, auxB] = ztl_eom(yev, P4, o4a.events(kev).to);
        uA = min(max(real(auxA.u),0),1);  uB = min(max(real(auxB.u),0),1);
        duMax = max(duMax, abs(uA - uB));
    end
    ok4 = numel(o4a.events) == numel(o4b.events) && dtEv < 5e-8 && ...
          dS < 1e-9 && duMax < 1e-8;
    [nPass, nFail] = gate(ok4, nPass, nFail, ...
        sprintf('G4 events eps=0.5: n=%d/%d, dt_repro=%.1e, |S|-eps=%.1e, du=%.1e', ...
        numel(o4a.events), numel(o4b.events), dtEv, dS, duMax));
end

%% G5 -- saltation across one eps=0 switch --------------------------------------
P5 = P;  P5.eps = 0;
probe = ztl_flow(y0BE, [0 2.5], P5, false);          % find the first switch
if isempty(probe.events)
    [nPass, nFail] = gate(false, nPass, nFail, 'G5: no eps=0 switch found in [0,2.5]');
else
    t5 = probe.events(1).t*1.02 + 0.01;              % span containing 1 switch
    o5 = ztl_flow(y0BE, [0 t5], P5, true);
    Phi5FD = zeros(14);
    for k = 1:14
        h = 1e-7*max(1, abs(y0BE(k)));
        ep = zeros(14,1);  ep(k) = h;
        op = ztl_flow(y0BE+ep, [0 t5], P5, false);
        om = ztl_flow(y0BE-ep, [0 t5], P5, false);
        Phi5FD(:,k) = (op.yf - om.yf)/(2*h);
    end
    e5 = norm(o5.PHI - Phi5FD)/norm(Phi5FD);
    [nPass, nFail] = gate(e5 < 1e-5 && numel(o5.events) == 1, nPass, nFail, ...
        sprintf('G5 saltation STM (1 switch @ t=%.3f, |Sdot|=%.2e): rel err = %.2e', ...
        probe.events(1).t, abs(probe.events(1).Sdot), e5));
end

fprintf('\n=== Z0 GATES: %d PASS, %d FAIL ===\n', nPass, nFail);

% ---------------------------------------------------------------------------
function [nP, nF] = gate(ok, nP, nF, msg)
if ok, nP = nP + 1;  tag = 'PASS'; else, nF = nF + 1;  tag = 'FAIL'; end
fprintf('[%s] %s\n', tag, msg);
end
