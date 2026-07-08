function out = NLP_lowThrust_GTO_Tulip_energy(N, tfFactor, makePlot, doIndirect)
% NLP_LOWTHRUST_GTO_TULIP_ENERGY  Min-ENERGY GTO -> tulip, FULL spiral.
%
% The point of this driver: min-energy has a CONTINUOUS control, so unlike
% min-fuel it can be posed on the ENTIRE GTO -> tulip transfer (~40 revs),
% not just an arrival leg. Fixed transfer time tf = tfFactor x the min
% time. Direct-first ("direct finds, indirect certifies"):
%   1. solve the collocation NLP (SOLVE_ENERGY_NLP) from a warm start built
%      by time-stretching the converged min-time arc onto [0, tf];
%   2. attempt the indirect polish (SOLVE_ENERGY_INDIRECT) seeded from the
%      min-time costates -- the honest stress test of whether single
%      shooting survives the ~40-perigee sensitivity for min-energy.
%
% Reuses the min-fuel transcription verbatim (LT_DYNAMICS_THROTTLE,
% NLP_CONSTRAINTS_MINFUEL: defects + throttle cone); only the objective
% changes (int 1/2 s^2 dt).
%
% INPUTS:
%   N        - (optional) trapezoidal segments [default 4000]
%   tfFactor - (optional) tf as a multiple of the full min time
%              [default 1.1 -- mild throttling, easy warm start]
%   makePlot - (optional) throttle + trajectory plots [default true]
%
% OUTPUTS:
%   out - struct: .energy, .mf, .mProp_kg, .dV_kms, .maxDefect, .exitflag,
%         .tauMesh, .X, .U, .sigma, .tf_ND, .tf_days, .indirect
%
% REFERENCES:
%   [1] Betts, SIAM 2010, Ch. 4.
%   [2] Caillau, Gergaud, Noailles, JOTA 2003 (min-energy transfer).

if nargin < 1 || isempty(N),         N = 4000;         end
if nargin < 2 || isempty(tfFactor),  tfFactor = 1.1;   end
if nargin < 3 || isempty(makePlot),  makePlot = true;  end
if nargin < 4 || isempty(doIndirect), doIndirect = true; end

setup_paths();
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'lowThrust_GTO_tulip'));

muStar = 0.012150585609624;
lStar  = 389703.264829278;
tStar  = 382981.289129055;
m0kg   = 15;
g0     = 9.80665*tStar^2/(1000*lStar);
Tmax   = (0.025/m0kg)*tStar^2/(lStar*1000);
c      = (2100/tStar)*g0;

% --- endpoints (same GTO and tulip target as the min-time/min-fuel work) ---
muEarth = 6.67384e-20*(1 - muStar)*(5.9736E24 + 7.35E22);
sma = (6378+350 + 6378+35786)/2;
ecc = (35786 - 350)/(2*sma);
[r0d, v0d] = pumpkyn.cr3bp.orb2eci(muEarth, [sma, ecc, 0, -25*pi/180, 0, 0], 2);
rv0 = pumpkyn.cr3bp.fromPCI(0, [r0d, v0d], muStar, tStar, lStar, 1);
[~, x0Tulip] = pumpkyn.cr3bp.getTulip((5/6)*2*pi, 7, -1, 1e-12);
[~, rvTgt]   = pumpkyn.cr3bp.prop((5/6)*2*pi, x0Tulip, muStar);
[~, idx_f]   = max(rvTgt(:,5));
rvf = rvTgt(idx_f, :);

% --- warm start: min-time burn (feasible) + ballistic coast ----------------
% A time-stretched min-time arc violates the dynamics across all 40 revs and
% fmincon cannot restore feasibility (measured: it grinds at ~0.04). Instead
% build a DYNAMICALLY CONSISTENT guess: burn the min-time arc to the tulip
% (exactly feasible), then coast ballistically for the extra time, and pin
% the target at the coast terminus (phase-shifted, one coast downstream).
% This is the same construction that made the min-fuel leg feasible, applied
% from tau = 0 to the whole transfer.
tfMinFull = 6.2906939607;
zMinTime  = [190.4760481; -79.7060409; -0.4298691037; 0.3011592775; ...
              0.5866700046; -0.007117348902; 4.329378839];
optsInt = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
tf = tfFactor*tfMinFull;
fprintf('min-energy FULL spiral: tf = %.4f ND = %.3f d (%.2fx min-time), N = %d\n', ...
        tf, tf*tStar/86400, tfFactor, N);

[tauB, yB] = ode113(@lt_pmp_eom, [0 tfMinFull], [rv0(:); 1; zMinTime], ...
                    optsInt, Tmax, c, muStar);
[tauB, kB] = unique(tauB, 'stable');  yB = yB(kB, :);

coastFun   = @(t, x) lt_dynamics_throttle(x, [0;0;0;0], Tmax, c, muStar);
[tauC, yC] = ode113(coastFun, [0, tf - tfMinFull], yB(end,1:7).', optsInt);
[tauC, kC] = unique(tauC, 'stable');  yC = yC(kC, :);
rvf = yC(end, 1:6);                       % target := coast terminus

tauAll  = [tauB; tfMinFull + tauC(2:end)];
xAll    = [yB(:,1:7); yC(2:end,:)];
sAll    = [0.98*ones(numel(tauB),1); 0.02*ones(numel(tauC)-1,1)];
lamVB   = yB(:,11:13);  alphB = -lamVB./sqrt(sum(lamVB.^2, 2));
vC      = yC(2:end,4:6); alphC = vC./sqrt(sum(vC.^2, 2));
alphAll = [alphB; alphC];

% density-matched mesh from the concatenated adaptive grid
sigAd = unique(tauAll)/tauAll(end);
sigma = interp1(linspace(0,1,numel(sigAd)).', sigAd, linspace(0,1,N+1).');
sigma(1) = 0; sigma(end) = 1;
assert(all(diff(sigma) > 0), 'energy:mesh', 'mesh not increasing');

tMesh = sigma.*tf;
Xg = interp1(tauAll, xAll,    tMesh, 'pchip').';
sG = min(max(interp1(tauAll, sAll, tMesh, 'pchip').', 0.02), 0.98);
aG = interp1(tauAll, alphAll, tMesh, 'pchip').';  aG = aG./sqrt(sum(aG.^2, 1));
Xg(1:6, 1)   = rv0(:);  Xg(7,1) = 1;
Xg(1:6, end) = rvf(:);
Z0 = [Xg(:); reshape([aG.*sG; sG], [], 1)];

[~, ceq0] = nlp_constraints_minfuel(Z0, sigma, tf, Tmax, c, muStar);
fprintf('burn+coast warm start: max defect = %.3g\n', max(abs(ceq0(1:7*N))));

% --- direct solve ----------------------------------------------------------
[~, nlp] = solve_energy_nlp(Z0, sigma, tf, rv0, 1, rvf, Tmax, c, muStar);
mProp = m0kg*(1 - nlp.mf);
dVtot = c*log(1/nlp.mf)*lStar/tStar;
fprintf('\n=== min-energy DIRECT (flag %d) ===\n', nlp.exitflag);
fprintf('energy J = %.6g, max defect = %.3g\n', nlp.energy, nlp.maxDefect);
fprintf('propellant = %.4f kg, dV = %.4f km/s, mf = %.6f\n', mProp, dVtot, nlp.mf);
fprintf('throttle range s in [%.3f, %.3f]  (min-energy ramps, not bang-bang)\n', ...
        min(nlp.U(4,:)), max(nlp.U(4,:)));

% --- indirect polish (stress test on the full spiral) ----------------------
indirect = [];
if doIndirect
try
    % seed from the DIRECT solution via the covector mapping (the proper
    % way): reconstruct lambda(0) from the NLP arc's primer directions +
    % lambda_m(tf)=0, scaled by the smooth stationarity u = S_e.
    lamSeed = costate_seed_from_nlp_energy(tMesh, nlp.X, nlp.U, Tmax, c, muStar);
    fprintf('reconstructed energy seed lambda(0) = [%s]\n', sprintf('%.5g ', lamSeed));
    [lamSol, resNorm, flag] = solve_energy_indirect(rv0, 1, rvf, tf, lamSeed, ...
                                Tmax, c, muStar);
    if resNorm < 1e-3
        [tauI, rvI] = ode113(@lt_pmp_eom_energy, [0 tf], [rv0(:); 1; lamSol], ...
                             optsInt, Tmax, c, muStar);
        indirect = struct('lamSol', lamSol, 'resNorm', resNorm, 'flag', flag, ...
                          'tau', tauI, 'rv', rvI);
        fprintf('indirect CONVERGED: ||R|| = %.3g\n', resNorm);
    else
        indirect = struct('lamSol', lamSol, 'resNorm', resNorm, 'flag', flag);
        fprintf('indirect did not converge on the full spiral (||R|| = %.3g)\n', resNorm);
    end
catch meErr
    fprintf('indirect attempt errored: %s\n', meErr.message);
end
end

out = struct('energy', nlp.energy, 'mf', nlp.mf, 'mProp_kg', mProp, ...
             'dV_kms', dVtot, 'maxDefect', nlp.maxDefect, ...
             'exitflag', nlp.exitflag, 'tauMesh', tMesh.', 'X', nlp.X, ...
             'U', nlp.U, 'sigma', sigma.', 'tf_ND', tf, ...
             'tf_days', tf*tStar/86400, 'indirect', indirect);

if makePlot
    figure('Color','w','Position',[100 100 900 700]);
    subplot(2,1,1);
    plot3(nlp.X(1,:), nlp.X(2,:), nlp.X(3,:), 'b', 'LineWidth', 0.8); hold on;
    plot3(-muStar,0,0,'o','MarkerFaceColor',[0.1 0.35 0.8],'MarkerSize',9);
    plot3(1-muStar,0,0,'o','MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',7);
    grid on; axis equal; view(-37,22);
    xlabel('x'); ylabel('y'); zlabel('z');
    title('min-energy full spiral (direct NLP)');
    subplot(2,1,2);
    plot(tMesh*tStar/86400, nlp.U(4,:), 'k', 'LineWidth', 1);
    xlabel('time [days]'); ylabel('throttle s'); ylim([-0.05 1.05]); grid on;
    title('min-energy throttle: continuous saturated ramp');
end
end
