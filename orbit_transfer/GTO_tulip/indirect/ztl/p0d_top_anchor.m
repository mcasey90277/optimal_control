% P0D_TOP_ANCHOR  Preflight retry: converge the Z3 top anchor at 200 mN using
% pumpkyn's analytic-STM min-time solver + the old energy shooter.
%
% Supersedes P0b/P0c after their diagnosis (see ZTL_RESULTS.md P0 section):
% the old CS+dogleg min-time ladder never actually converged above 1.2x (the
% duplicate tfMin entries in the 2025 table were stalls), so P0c's energy
% probe was seeded from garbage. This probe replaces the machinery, not the
% idea:
%
%   Stage 1  min-time at 200 mN with pumpkyn.cr3bp.tfMin (analytic STM
%            Jacobian -- Zhang ingredient (a), already built): direct jump
%            from the validated 25 mN costates, thrust-march fallback.
%   Stage 2  min-energy at 200 mN with solve_energy_indirect, seeded from
%            the Stage-1 costates (beta-rescaled onto S_e), at the FIXED
%            ladder time tf = 1.15 x tfMin(25 mN) = 7.2344 ND. tfMin is
%            nonincreasing in thrust, so this tf keeps >= 15% margin at
%            EVERY rung -- the fixed-tf ladder needs no min-time table.
%
% Output: results/p0d_top_anchor.mat + verdicts.

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

[rv0, rvf, P] = ztl_endpoints();

zGuess25 = [190.476497248065; -79.7064866984696; -0.430399154713168; ...
             0.301159446575878; 0.586671892449694; -0.00711582435720301; ...
             4.32931089137559; 6.29081541876621];   % validated 25 mN min-time
tfLadder = 1.15*zGuess25(8);                        % fixed ladder tf = 7.2344

% --- residual evaluator (same conditions as pumpkyn's shootingResidual) ----
mintimeRes = @(z, Tmax) local_mintime_res(z, rv0, rvf, Tmax, P.c, P.muStar);

% Bad mid-iteration costates crash trajectories into the Earth singularity;
% the resulting huge residual steers the solver away, so the per-eval ODE
% tolerance warnings are pure noise (a stuck direct-jump attempt flooded
% 2.6M of them). March warm instead of jumping, and silence the warning.
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');

%% Stage 1: min-time thrust march 25 -> 200 mN (analytic STM) ---------------
% ztl_mintime_solve = pumpkyn's residual/J at REAL budgets (its own wrapper
% caps MaxFunctionEvaluations at 100, which stalled every rung; P0d run 2).
fprintf('=== P0d Stage 1: min-time march via ztl_mintime_solve ===\n');
TmaxTop = 8*P.Tmax25;
fMarch = [1, 1.25, 1.56, 1.95, 2.44, 3.05, 3.81, 4.77, 5.96, 8];
z = zGuess25;
for k = 1:numel(fMarch)
    Tk = fMarch(k)*P.Tmax25;
    tic;
    [z, rn, mtOut] = ztl_mintime_solve(rv0, rvf, z, Tk, P.c, P.muStar, 1500);
    fprintf('  f=%5.2f (%6.1f mN): tf_min=%.6f  ||R||=%.3e  flag=%d  switches=%d  (%.0f s)\n', ...
            fMarch(k), 25*fMarch(k), z(8), rn, mtOut.flag, mtOut.nSwitch, toc);
end
zTop = z;  rnTop = rn;  tfMinTop = z(8);
okS1 = rnTop < 1e-8;
fprintf('Stage 1 %s: tf_min(200 mN) = %.6f, ||R|| = %.3e\n', ...
    ternary(okS1,'PASS','FAIL'), tfMinTop, rnTop);
fprintf('fixed ladder tf = %.4f -> margin over tfMin(200) = %.1f%%\n', ...
    tfLadder, 100*(tfLadder/tfMinTop - 1));

%% Stage 2: min-energy @ 200 mN, tf = 7.2344, seeded from Stage 1 -----------
fprintf('\n=== P0d Stage 2: min-energy @ 200 mN, tf = %.4f ===\n', tfLadder);
lamMT = zTop(1:7);
Se0raw = TmaxTop*(norm(lamMT(4:6))/1 + lamMT(7)/P.c);
betaGrid = [0.9, 0.6, 1.5, 0.3]/Se0raw;

best = struct('resNorm', inf);
for kb = 1:numel(betaGrid)
    beta = betaGrid(kb);
    fprintf('--- seed scale beta = %.3e (S_e(0) = %.2f) ---\n', beta, beta*Se0raw);
    [lamSol, rnE, flag] = solve_energy_indirect(rv0, 1, rvf, tfLadder, ...
                              beta*lamMT, TmaxTop, P.c, P.muStar);
    if rnE < best.resNorm
        best = struct('resNorm', rnE, 'lamSol', lamSol, 'flag', flag, 'beta', beta);
    end
    if rnE < 1e-8, break; end
end

optsInt = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, yI] = ode113(@lt_pmp_eom_energy, [0 tfLadder], [rv0(:); 1; best.lamSol], ...
                 optsInt, TmaxTop, P.c, P.muStar);
lamvMag = sqrt(sum(yI(:,11:13).^2, 2));
Se = TmaxTop*(lamvMag./yI(:,7) + yI(:,14)/P.c);
u  = min(max(Se, 0), 1);
mF = yI(end,7);

anchor = struct('Tmax_mN', 200, 'tf', tfLadder, 'tfMinTop', tfMinTop, ...
    'zMinTime', zTop, 'rnMinTime', rnTop, ...
    'lam0Energy', best.lamSol, 'rnEnergy', best.resNorm, 'flag', best.flag, ...
    'beta', best.beta, 'mProp_kg', P.m0kg*(1-mF), ...
    'dV_kms', P.c*log(1/mF)*P.lStar/P.tStar, ...
    'uMin', min(u), 'uMax', max(u), 'fracSatHi', mean(u > 0.999), ...
    'fracSatLo', mean(u < 1e-3), 'rv0', rv0, 'rvf', rvf, 'P', P);
save(fullfile(resDir, 'p0d_top_anchor.mat'), 'anchor');
fprintf('saved %s\n', fullfile(resDir, 'p0d_top_anchor.mat'));

fprintf(['Stage 2 RESULT: ||R|| = %.3e (flag %d, beta %.3e)  prop = %.4f kg  ' ...
         'dV = %.4f km/s\n  throttle: min %.3f  max %.3f  sat-hi %.1f%%  sat-lo %.1f%%\n'], ...
    best.resNorm, best.flag, best.beta, anchor.mProp_kg, anchor.dV_kms, ...
    anchor.uMin, anchor.uMax, 100*anchor.fracSatHi, 100*anchor.fracSatLo);
if best.resNorm < 1e-8
    fprintf('GATE P0d: PASS -- Z3 top anchor (min-energy @ 200 mN) CONVERGED.\n');
else
    fprintf('GATE P0d: NOT CONVERGED (Stage 1 %s) -- see log.\n', ternary(okS1,'ok','failed too'));
end

% ---------------------------------------------------------------------------
function R = local_mintime_res(z, rv0, rvf, Tmax, c, muStar)
% Min-time residual [8x1] at z = [lambda0(7); tf]: rendezvous + lam_m + H.
y0 = [rv0(:); 1; z(1:7)];
[~, Y] = pumpkyn.cr3bp.tfMinProp(z(8), y0, Tmax, c, muStar);
yf = Y(end, 1:14).';
[~, Hf] = pumpkyn.cr3bp.tfMinEoM(z(8), yf, Tmax, c, muStar);
R = [yf(1:6) - rvf(:); yf(14); Hf];
end

function v = ternary(cond, a, b)
if cond, v = a; else, v = b; end
end
