function out = run_gto_tulip_minfuel(lamSeed, tfFactor, makePlot)
% RUN_GTO_TULIP_MINFUEL  Min-fuel GTO -> tulip transfer, indirect side.
%
% Solves the fixed-tf minimum-fuel TPBVP by smoothing continuation
% (SOLVE_MINFUEL_INDIRECT) from a SUPPLIED costate seed, then reports the
% burn/coast structure. The seed matters: min-time costates -- raw or
% rescaled -- are structurally wrong for min-fuel (the "+1" in the
% switching function anchors an absolute costate scale, and the true
% min-fuel lambda_r/lambda_v ratio differs; shooting stalls at
% ||R|| ~ 0.4-1 from any such seed). Obtain the seed from the direct
% solution's KKT multipliers via the covector mapping -- see
% GTO_tulip/NLP_lowThrust_GTO_Tulip_minfuel.m, which runs
% the full direct-then-indirect pipeline and calls this driver's
% machinery.
%
% INPUTS:
%   lamSeed  - initial costate seed [7x1] (Lagrange form: lambda_m(tf)=0
%              convention), e.g. [-mu_1(1:6); 1-mu_1(7)] from the NLP
%   tfFactor - (optional) tf as a multiple of the min time [default 1.2]
%   makePlot - (optional) true for switching/throttle plots [default false]
%
% OUTPUTS:
%   out - results struct:
%           .lamSol      converged initial costates [7x1]
%           .tf_ND, .tf_days   fixed transfer time
%           .resNorm     terminal residual 2-norm
%           .mProp_kg    propellant used (min-time reference: 2.9247 kg)
%           .dV_kms      total delta-V
%           .burnFrac    fraction of the arc with u > 0.5
%           .nCoasts     number of coast arcs
%           .tau, .rv, .S, .u   trajectory, switching fn, throttle
%
% REFERENCES:
%   [1] Bertrand & Epenoy, OCAM 23(4), 2002.

if nargin < 2 || isempty(tfFactor), tfFactor = 1.2;  end
if nargin < 3 || isempty(makePlot), makePlot = false; end

muStar = 0.012150585609624;
lStar  = 389703.264829278;
tStar  = 382981.289129055;

muEarth = 6.67384e-20*(1 - muStar)*(5.9736E24 + 7.35E22);
sma = (6378+350 + 6378+35786)/2;
ecc = (35786 - 350)/(2*sma);
[r0, v0] = pumpkyn.cr3bp.orb2eci(muEarth, [sma, ecc, 0, -25*pi/180, 0, 0], 2);
rv0 = pumpkyn.cr3bp.fromPCI(0, [r0, v0], muStar, tStar, lStar, 1);

[~, x0Tulip] = pumpkyn.cr3bp.getTulip((5/6)*2*pi, 7, -1, 1e-12);
[~, rvTgt] = pumpkyn.cr3bp.prop((5/6)*2*pi, x0Tulip, muStar);
[~, idx_f] = max(rvTgt(:,5));
rvf = rvTgt(idx_f, :);

m0kg = 15;
g0   = 9.80665*tStar^2/(1000*lStar);
Tmax = (0.025/m0kg)*tStar^2/(lStar*1000);
c    = (2100/tStar)*g0;

tfMinRef = 6.2906939607;
tf = tfFactor*tfMinRef;

fprintf('min-fuel indirect: tf = %.6f ND = %.4f days (%.2fx min-time)\n', ...
        tf, tf*tStar/86400, tfFactor);
[lamSol, resNorm] = solve_minfuel_indirect(rv0, 1, rvf, tf, lamSeed, ...
                        Tmax, c, muStar, [0.03, 0.01, 3e-3, 1e-3]);

epsFinal = 1e-3;
optsInt  = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[tau, rv] = ode113(@lt_pmp_eom_minfuel, [0 tf], [rv0(:); 1; lamSol], ...
                   optsInt, Tmax, c, muStar, epsFinal);

lamvMag = sqrt(sum(rv(:,11:13).^2, 2));
S  = 1 - lamvMag.*c./rv(:,7) - rv(:,14);
u  = (1 - tanh(S./(2*epsFinal)))/2;

mProp    = m0kg*(1 - rv(end,7));
dVtot    = c*log(1/rv(end,7))*lStar/tStar;
burnFrac = trapz(tau, double(u > 0.5))/tf;
nCoasts  = sum(diff([false; S(:) > 0]) == 1);

out = struct('lamSol', lamSol, 'tf_ND', tf, 'tf_days', tf*tStar/86400, ...
             'resNorm', resNorm, 'mProp_kg', mProp, 'dV_kms', dVtot, ...
             'burnFrac', burnFrac, 'nCoasts', nCoasts, ...
             'tau', tau, 'rv', rv, 'S', S, 'u', u);

fprintf('||R|| = %.3g\n', resNorm);
fprintf('propellant = %.4f kg (min-time: 2.9247), dV = %.4f km/s\n', ...
        mProp, dVtot);
fprintf('burn fraction = %.1f%%, coast arcs = %d\n', 100*burnFrac, nCoasts);

if makePlot
    figure('Color', 'w');
    subplot(2,1,1);
    plot(tau*tStar/86400, S, 'k', 'LineWidth', 1); hold on;
    yline(0, ':');
    xlabel('time [days]'); ylabel('S(t)'); grid on;
    title('min-fuel switching function');
    subplot(2,1,2);
    plot(tau*tStar/86400, u, 'k', 'LineWidth', 1);
    xlabel('time [days]'); ylabel('throttle u'); ylim([-0.05 1.05]); grid on;
end
end
