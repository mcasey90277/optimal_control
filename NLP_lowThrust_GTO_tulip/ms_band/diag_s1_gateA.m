% DIAG_S1_GATEA  Attribute the Gate A(ii) mismatch + calibrate check (iii-b).
%
% Question 1: is the 8.5e-8 sigma-vs-time state mismatch an EOM error or
% integrator-noise amplification through the perigee passes?
%   (a) short-span comparison (first 0.5 ND of time): amplification is
%       small there, so a real EOM inconsistency would survive; noise ~1e-12.
%   (b) time-domain SELF-noise: same integration at AbsTol 1e-15 vs 3e-15;
%       their mismatch is the amplification floor no cross-domain
%       comparison can beat.
%   (c) mismatch profile vs time (where does it grow?).
% Question 2: detection power of the offset-lamT conservation check on a
% SHORT span (full span diverges into the Earth singularity): measure the
% kap*Hval drift over sigma in [0,3] with the correct EOM. The flipped
% -sign EOM is checked by sms_eom_flipsign_scratch.m (separate scratch file).
setup_paths;
ref  = run_gto_tulip_indirect(false);
lam7 = ref.zSol(1:7);
tf   = ref.zSol(8);
prob = sms_problem(1.00, 0.3);
epsA = 0.3;

y0T = [prob.rv0; prob.m0; lam7];
y0S = [prob.rv0; prob.m0; 0; lam7; 0];
[~, Ht0] = sms_eom(0, y0S, prob.Tmax, prob.c, prob.muStar, epsA, prob.pSund);
y0S(16) = -Ht0;

fT = @(t, y) lt_pmp_eom_minfuel(t, y, prob.Tmax, prob.c, prob.muStar, epsA);
fS = @(s, y) sms_eom(s, y, prob.Tmax, prob.c, prob.muStar, epsA, prob.pSund);

% (a) short-span cross-domain comparison: integrate sigma-domain to t = 0.5
optE  = odeset(prob.odeOpts, 'Events', @(s, y) diag_event(y, 0.5));
solSs = ode113(fS, [0 20], y0S, optE);
solTs = ode113(fT, [0 0.5], y0T, prob.odeOpts);
tSs   = solSs.y(8, :);
keep  = tSs > 0.01 & tSs < 0.49;
YTs   = deval(solTs, tSs(keep));
fprintf('(a) short span (t<=0.5 ND): max|r,v,m err| = %.3e  costates %.3e\n', ...
        max(max(abs(solSs.y(1:7, keep) - YTs(1:7, :)))), ...
        max(max(abs(solSs.y(9:15, keep) - YTs(8:14, :)))));

% (b) time-domain self-noise: AbsTol 1e-15 vs 3e-15 over full span
solT1 = ode113(fT, [0 tf], y0T, prob.odeOpts);
solT2 = ode113(fT, [0 tf], y0T, odeset(prob.odeOpts, 'AbsTol', 3e-15));
tq    = linspace(0.02*tf, 0.98*tf, 400);
selfN = max(max(abs(deval(solT1, tq) - deval(solT2, tq))));
fprintf('(b) time-domain self-noise (AbsTol 1e-15 vs 3e-15): %.3e\n', selfN);

% (c) cross-domain mismatch profile vs time
optE2 = odeset(prob.odeOpts, 'Events', @(s, y) diag_event(y, tf));
solS  = ode113(fS, [0 300], y0S, optE2);
tS    = solS.y(8, :);
kk    = find(tS > 0.02*tf & tS < 0.98*tf);
YT    = deval(solT1, tS(kk));
errP  = max(abs(solS.y(1:7, kk) - YT(1:7, :)), [], 1);
edges = linspace(0.02*tf, 0.98*tf, 13);
for b = 1:12
    inB = tS(kk) >= edges(b) & tS(kk) < edges(b+1);
    if any(inB)
        fprintf('(c) t in [%5.2f %5.2f]: max err %.3e\n', ...
                edges(b), edges(b+1), max(errP(inB)));
    end
end

% (2) offset-lamT conservation on a short span, correct EOM
y0b = y0S;  y0b(16) = y0S(16) + 0.5;
solB = ode113(fS, [0 3], y0b, prob.odeOpts);
Hs   = zeros(1, size(solB.y, 2));
for k = 1:numel(Hs)
    [~, Htk] = sms_eom(0, solB.y(:, k), prob.Tmax, prob.c, prob.muStar, ...
                       epsA, prob.pSund);
    r1k = sqrt(sum((solB.y(1:3, k) - [-prob.muStar; 0; 0]).^2));
    Hs(k) = r1k^prob.pSund*(Htk + solB.y(16, k));
end
fprintf('(2) short-span offset-lamT: kapHval(0) = %.6e  drift = %.3e\n', ...
        Hs(1), max(abs(Hs - Hs(1))));

% (2b) detection power: same short span with the dkapdr sign FLIPPED
fF   = @(s, y) eom_flipped(y, prob.Tmax, prob.c, prob.muStar, epsA, prob.pSund);
solF = ode113(fF, [0 3], y0b, prob.odeOpts);
HsF  = zeros(1, size(solF.y, 2));
for k = 1:numel(HsF)
    [~, Htk] = sms_eom(0, solF.y(:, k), prob.Tmax, prob.c, prob.muStar, ...
                       epsA, prob.pSund);
    r1k = sqrt(sum((solF.y(1:3, k) - [-prob.muStar; 0; 0]).^2));
    HsF(k) = r1k^prob.pSund*(Htk + solF.y(16, k));
end
fprintf('(2b) FLIPPED dkapdr sign, same span: drift = %.3e (should be >> gate)\n', ...
        max(abs(HsF - HsF(1))));

% -------------------------------------------------------------------------
function dY = eom_flipped(Y, Tmax, c, muStar, epsSmooth, pSund)
% EOM_FLIPPED  sms_eom with the dkapdr*Hval adjoint term sign flipped
% (detection-power calibration only; wraps the true EOM and re-adds 2x the
% kappa-gradient term).
%
% INPUTS:
%   Y         - augmented state [16x1]
%   Tmax, c, muStar, epsSmooth, pSund - as SMS_EOM
%
% OUTPUTS:
%   dY - d Y / d sigma with lamRdot = kap*(-G'lamV) + dkapdr*Hval [16x1]
[dY, Ht] = sms_eom(0, Y, Tmax, c, muStar, epsSmooth, pSund);
dd  = [Y(1) + muStar; Y(2); Y(3)];
d1  = sqrt(sum(dd.^2));
dkapdr = pSund * d1^(pSund - 2) * dd;
dY(9:11) = dY(9:11) + 2*dkapdr*(Ht + Y(16));
end

% -------------------------------------------------------------------------
function [val, isterm, dir] = diag_event(y, tStop)
% DIAG_EVENT  Stop when the carried time state reaches tStop.
%
% INPUTS:
%   y     - augmented state [16x1]
%   tStop - stop time (ND) [scalar]
%
% OUTPUTS:
%   val    - y(8) - tStop [scalar]
%   isterm - 1
%   dir    - +1
val    = y(8) - tStop;
isterm = 1;
dir    = 1;
end
