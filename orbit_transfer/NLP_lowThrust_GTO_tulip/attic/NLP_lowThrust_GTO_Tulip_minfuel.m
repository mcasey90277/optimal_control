function out = NLP_lowThrust_GTO_Tulip_minfuel(N, tfFactor, legStart, makePlot)
% NLP_LOWTHRUST_GTO_TULIP_MINFUEL  Min-fuel tulip-arrival leg, direct-first.
%
% MIN-FUEL variant, posed on the FINAL LEG of the GTO -> tulip transfer as
% a mid-transfer replan: start from the min-time arc's state at
% tau = legStart (position, velocity, AND remaining mass), fix the leg
% time at tfFactor times the leg's minimum time, and minimize propellant.
% By the principle of optimality the tail of the full min-time solution IS
% the leg's min-time solution, so the leg's minimum time and its optimal
% costates are known exactly -- no new min-time solve.
%
% The arrival state is the TULIP POINT ONE COAST DOWNSTREAM of the
% original rvf: burn the min-time tail (arriving at rvf at tfMinLeg), then
% ride the tulip ballistically for the remaining time. That makes the
% burn+coast warm start feasible by construction up to the interior-point
% throttle clip (s in [0.02, 0.98] vs the propagated s = 1/0 states;
% measured max defect ~1e-4) -- pinning the original rvf instead leaves an
% O(0.3)-ND defect cliff at the clipped final node, which feasibility
% restoration cannot heal (measured: flag -2). With a near-feasible start
% the NLP's only job is optimization: beat the warm start's
% min-time-style propellant by redistributing thrust.
%
% Why a leg and not the full transfer: the full 25-mN min-fuel spiral has
% ~40 revolutions and ~80 burn/coast switches -- research-grade for both
% methods (measured: the NLP grinds in feasibility mode without
% optimizing; covector-seeded shooting hits integrator breakdown at a
% perigee dive). The leg (default: the last ~1.3 ND, a few revolutions,
% a handful of switches) exercises every piece of min-fuel machinery at
% tractable difficulty.
%
% Architecture ("direct finds, indirect certifies"):
%   1. solve the NLP from a BURN+COAST warm start (min-time tail + tulip
%      coast; dynamically consistent except one healable spike at the
%      pinned final node);
%   2. reconstruct a shooting seed from the NLP arc (COSTATE_SEED_FROM_NLP:
%      linear costate ODE + LSQ fit to primer directions; fmincon's own
%      eqnonlin multipliers at an lbfgs stall are ~100x off in scale and
%      printed only as a diagnostic);
%   3. attempt the indirect polish (smoothing continuation). HONEST
%      STATUS: the polish does NOT converge -- the min-fuel shooting
%      basin on a multi-rev CR3BP leg is smaller than the best seed here
%      (residual progression 1.55/0.83/0.33/0.14 across four seed
%      strategies). Closing it needs switch-time multi-arc shooting or
%      multiple shooting; see the theory note S6 and tutorial Phase H2.
%
% INPUTS:
%   N        - (optional) trapezoidal segments [default 3000]
%   tfFactor - (optional) leg time as multiple of the leg's min time
%              [default 1.3]
%   legStart - (optional) start of the leg along the min-time arc (ND)
%              [default 4.0 -- several perigee passes, so the optimum has
%              real multi-switch structure; 0 = the full research-grade
%              case]
%   makePlot - (optional) true for throttle-profile overlay [default true]
%
% OUTPUTS:
%   out - results struct: .mProp_kg, .dV_kms, .burnFrac (NLP values),
%         .maxDefect, .tauMesh, .X, .U, .exitflag, .lamSeed,
%         .mPropMinTime_kg (leg's min-time propellant, the number to beat),
%         .indirect (polished indirect result struct, [] if not converged)
%
% REFERENCES:
%   [1] Betts, SIAM 2010, Ch. 4 (transcription; covector mapping).
%   [2] Bertrand & Epenoy, OCAM 23(4), 2002 (indirect-side smoothing).

if nargin < 1 || isempty(N),        N = 3000;        end
if nargin < 2 || isempty(tfFactor), tfFactor = 1.3;  end
if nargin < 3 || isempty(legStart), legStart = 4.0;  end
if nargin < 4 || isempty(makePlot), makePlot = true;  end

setup_paths();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'lowThrust_GTO_tulip'));

muStar = 0.012150585609624;
lStar  = 389703.264829278;
tStar  = 382981.289129055;
m0kg   = 15;
g0     = 9.80665*tStar^2/(1000*lStar);
Tmax   = (0.025/m0kg)*tStar^2/(lStar*1000);
c      = (2100/tStar)*g0;

% --- full min-time arc (endpoints + leg state + burn-phase guess) ----------
muEarth = 6.67384e-20*(1 - muStar)*(5.9736E24 + 7.35E22);
sma = (6378+350 + 6378+35786)/2;
ecc = (35786 - 350)/(2*sma);
[r0d, v0d] = pumpkyn.cr3bp.orb2eci(muEarth, [sma, ecc, 0, -25*pi/180, 0, 0], 2);
rv0Full = pumpkyn.cr3bp.fromPCI(0, [r0d, v0d], muStar, tStar, lStar, 1);
% (the ORIGINAL tulip target needs no explicit computation here: the
% min-time arc below terminates on it, and the leg's target is the
% phase-shifted coast terminus computed later)

tfMinFull = 6.2906939607;
zMinTime  = [190.4760481; -79.7060409; -0.4298691037; 0.3011592775; ...
              0.5866700046; -0.007117348902; 4.329378839];
optsInt = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[tauF, yF] = ode113(@lt_pmp_eom, [0 tfMinFull], [rv0Full(:); 1; zMinTime], ...
                    optsInt, Tmax, c, muStar);
[tauF, keep] = unique(tauF, 'stable');
yF = yF(keep, :);

% Leg quantities (principle of optimality: the arc tail IS the leg's
% min-time solution, with costates yLeg(8:14))
xLegStart = interp1(tauF, yF, legStart, 'pchip');
rv0   = xLegStart(1:6);
m0    = xLegStart(7);
tfMinLeg = tfMinFull - legStart;
tf    = tfFactor*tfMinLeg;
mPropMinTime = m0kg*Tmax*tfMinLeg/c;     % full burn over the min-time leg

fprintf('leg: start tau = %.3f, m0 = %.4f, tfMinLeg = %.4f, tf = %.4f ND (%.4f d)\n', ...
        legStart, m0, tfMinLeg, tf, tf*tStar/86400);
fprintf('warm-start (burn-then-coast) propellant = %.4f kg (the number to beat)\n', ...
        mPropMinTime);

% --- burn + coast warm start ------------------------------------------------
maskB = tauF >= legStart;
tauB  = tauF(maskB) - legStart;
yB    = yF(maskB, :);
if tauB(1) > 0, tauB = [0; tauB]; yB = [xLegStart; yB]; end %#ok<AGROW>

coastFun = @(t, x) lt_dynamics_throttle(x, [0;0;0;0], Tmax, c, muStar);
[tauC, yC] = ode113(coastFun, [0, tf - tfMinLeg], yB(end,1:7).', optsInt);
[tauC, keep] = unique(tauC, 'stable');
yC = yC(keep, :);

% arrival state = coast terminus (exactly feasible warm start)
rvf    = yC(end, 1:6);
tauAll = [tauB; tfMinLeg + tauC(2:end)];
xAll   = [yB(:,1:7); yC(2:end,:)];
uAll   = [0.98*ones(numel(tauB),1); 0.02*ones(numel(tauC)-1,1)];
lamVB  = yB(:,11:13);
alphB  = -lamVB./sqrt(sum(lamVB.^2, 2));
vC     = yC(2:end, 4:6);
alphC  = vC./sqrt(sum(vC.^2, 2));
alphAll = [alphB; alphC];

sigFull = tauAll./tauAll(end);
sigma   = interp1(linspace(0,1,numel(sigFull)).', sigFull, ...
                  linspace(0,1,N+1).');
sigma(1) = 0; sigma(end) = 1;
assert(all(diff(sigma) > 0), 'minfuel:mesh', 'mesh not increasing');

tMesh = sigma.*tf;
Xg    = interp1(tauAll, xAll,    tMesh, 'pchip').';
sG    = min(max(interp1(tauAll, uAll, tMesh, 'pchip').', 0.02), 0.98);
aG    = interp1(tauAll, alphAll, tMesh, 'pchip').';
aG    = aG./sqrt(sum(aG.^2, 1));
Z0    = [Xg(:); reshape([aG.*sG; sG], [], 1)];

[~, ceq0] = nlp_constraints_minfuel(Z0, sigma, tf, Tmax, c, muStar);
fprintf('burn+coast warm start: max defect = %.3g\n', max(abs(ceq0(1:7*N))));

% pin the leg's initial mass (not 1) -- handled by solve via rv0/m0 bounds:
[~, nlp] = solve_minfuel_nlp(Z0, sigma, tf, rv0, m0, rvf, Tmax, c, muStar);
fprintf('NLP: flag %d, mf = %.6f, max defect = %.3g\n', ...
        nlp.exitflag, nlp.mf, nlp.maxDefect);

mProp    = m0kg*(m0 - nlp.mf);
dVtot    = c*log(m0/nlp.mf)*lStar/tStar;
sSol     = nlp.U(4, :);
burnFrac = trapz(tMesh, double(sSol > 0.5).')/tf;

fprintf('\n=== min-fuel NLP (flag %d) ===\n', nlp.exitflag);
fprintf('leg tf (fixed) = %.6f ND = %.4f days\n', tf, tf*tStar/86400);
fprintf('propellant = %.4f kg (min-time leg: %.4f), dV = %.4f km/s\n', ...
        mProp, mPropMinTime, dVtot);
fprintf('burn fraction = %.1f%%, max defect = %.3g\n', 100*burnFrac, nlp.maxDefect);

% --- covector mapping (diagnostic) -----------------------------------------
% The KKT multipliers ARE the discrete costates in principle, but at an
% lbfgs step-tolerance stall fmincon's multiplier estimates are
% feasibility-accurate and optimality-loose (measured: ~100x too small,
% throttle channel saturated at coast). Printed for the record; the seed
% below anchors the scale physically instead.
muDef   = reshape(full(nlp.eqnonlin(1:7*N)), 7, N);   % sparse from fmincon
lamCov  = [-muDef(1:6,1); 1 - muDef(7,1)];
fprintf('covector lambda(0) (diagnostic) = [%s]\n', sprintf('%.6g ', lamCov));

% --- shooting seed: costate reconstruction from the NLP arc -----------------
% The costate ODE is linear along the known trajectory; fit lambda(0) to
% the NLP's primer directions + lambda_m(tf) = 0 by linear least squares
% (see COSTATE_SEED_FROM_NLP). This is the covector mapping done properly.
lamSeed = costate_seed_from_nlp(tMesh, nlp.X, nlp.U, Tmax, c, muStar);
fprintf('reconstructed seed lambda(0) = [%s]\n', sprintf('%.6g ', lamSeed));

% --- indirect polish --------------------------------------------------------
indirect = [];
try
    [lamSol, resNorm] = solve_minfuel_indirect(rv0, m0, rvf, tf, lamSeed, ...
                            Tmax, c, muStar, [0.03, 0.01, 3e-3, 1e-3]);
    if resNorm < 1e-3
        epsF = 1e-3;
        [tauI, rvI] = ode113(@lt_pmp_eom_minfuel, [0 tf], ...
                             [rv0(:); m0; lamSol], optsInt, ...
                             Tmax, c, muStar, epsF);
        lamvMag = sqrt(sum(rvI(:,11:13).^2, 2));
        SI  = 1 - lamvMag.*c./rvI(:,7) - rvI(:,14);
        uI  = (1 - tanh(SI./(2*epsF)))/2;
        indirect = struct('lamSol', lamSol, 'resNorm', resNorm, ...
            'mProp_kg', m0kg*(m0 - rvI(end,7)), ...
            'dV_kms', c*log(m0/rvI(end,7))*lStar/tStar, ...
            'burnFrac', trapz(tauI, double(uI > 0.5))/tf, ...
            'nCoasts', sum(diff([false; SI(:) > 0]) == 1), ...
            'tau', tauI, 'rv', rvI, 'S', SI, 'u', uI);
        fprintf('\n=== indirect polish ===\n');
        fprintf('||R|| = %.3g, propellant = %.4f kg, dV = %.4f km/s\n', ...
                resNorm, indirect.mProp_kg, indirect.dV_kms);
        fprintf('burn fraction = %.1f%%, coast arcs = %d\n', ...
                100*indirect.burnFrac, indirect.nCoasts);
    else
        fprintf('indirect polish did not converge (||R|| = %.3g)\n', resNorm);
    end
catch polishErr
    fprintf('indirect polish errored: %s\n', polishErr.message);
end

out = struct('mProp_kg', mProp, 'dV_kms', dVtot, 'burnFrac', burnFrac, ...
             'maxDefect', nlp.maxDefect, 'tauMesh', tMesh.', ...
             'X', nlp.X, 'U', nlp.U, 'exitflag', nlp.exitflag, ...
             'lamSeed', lamSeed, 'mPropMinTime_kg', mPropMinTime, ...
             'indirect', indirect);

if makePlot
    figure('Color', 'w');
    plot(tMesh*tStar/86400, sSol, 'k', 'LineWidth', 1); hold on;
    if ~isempty(indirect)
        plot(indirect.tau*tStar/86400, indirect.u, 'r:', 'LineWidth', 1);
        legend('direct NLP', 'indirect (smoothed)', 'Location', 'best');
    end
    xlabel('time [days]'); ylabel('throttle s');
    ylim([-0.05 1.05]); grid on;
    title('min-fuel leg throttle: bang-bang structure');
end
end
