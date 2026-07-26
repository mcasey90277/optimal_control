%% SOLVE_MINFUEL  Run the min-fuel GTO -> south-pole tulip transfer solve.
%
% Solves the minimum-fuel low-thrust transfer in the Earth-Moon CR3BP as a
% direct NLP (Sundman-regularized trapezoidal collocation, CasADi + IPOPT) by
% (1) loading the collocation-feasible energy seed, (2) re-timing it to the
% transfer time you pick below, (3) solving min-ENERGY at that time, then
% (4) sweeping the energy->fuel homotopy down to sharp bang-bang min-fuel.
%
% Edit the PARAMETERS block, then run the whole file. The solution lands in the
% workspace as `best` (and `tbl`, the per-homotopy-step table).
%
% Reference default (tf_factor = 1.15): ~25 switches, defect ~1e-14,
% propellant ~2.264 kg, dV ~3.3696 km/s. Expect ~10-15 min at 4000 nodes.
%
% Requires: CasADi on $CASADI_PATH (default ~/casadi-3.7.0).

clear; clc;

%% ======================= PARAMETERS (edit me) ==========================

% --- transfer time -------------------------------------------------------
% Set as a multiple of the minimum-time transfer t_f^min (= 27.8845 days).
% Reliable band ~1.12 - 1.25 (and 1.85); the hard transition band ~1.01-1.11
% is not solvable this way. 1.15 reproduces the certified reference solution.
tf_factor = 1.15;

% --- spacecraft (PINNED to the stored seed; see note) --------------------
% The energy seed was built for this exact configuration. Changing these
% invalidates the seed and the solve will almost certainly fail -- leave as is
% unless you also regenerate minfuel_from_energy_seed.mat.
thrust_N = 0.025;    % max thrust [N]  (25 mN)
m0_kg    = 15;       % wet mass   [kg]
Isp_s    = 2100;     % specific impulse [s]

% --- solver / output -----------------------------------------------------
maxIter    = 1500;   % IPOPT iterations per homotopy step
saveResult = true;   % write the solution to results/
doPlot     = true;   % throttle history + trajectory plots

%% =======================================================================
%  Nothing below here needs editing for a normal run.
%% =======================================================================

here = fileparts(mfilename('fullpath'));
addpath(here);
cfg = minfuel_config();                  % single source of truth (anchors, schedules)
pSund = cfg.pSund;                        % Sundman power dt/dtau = r1^pSund

% physics constants and the transfer time
p  = cr3bp_lt_params(thrust_N, m0_kg, Isp_s);
tf = tf_factor * cfg.tfMin;               % target transfer TIME (ND)
fprintf('\n=== min-fuel solve: tf_factor=%.3f  ->  t_f=%.4f ND (%.2f days) ===\n', ...
        tf_factor, tf, tf*p.tStar/86400);

% --- (1) load the collocation-feasible energy seed -----------------------
% Carries the true endpoints (rv0, rvf) and a feasible trajectory at 1.15x.
S   = load(fullfile(here, 'minfuel_from_energy_seed.mat'));
rv0 = S.rv0;   rvf = S.rvf;

% --- (2) map into Sundman coordinates, RE-TIMED to the chosen tf ---------
% sundman_seed_map builds tSeed = sgNorm*tf, so passing the target tf linearly
% re-times the seed. Positions (hence kappa) are unchanged; the energy solve
% below adapts the re-timed guess to a feasible optimum at this tf.
[sigma, X0, U0, tauf0] = sundman_seed_map(S.nlp.X, S.nlp.U, tf, S.sigma, ...
                                          pSund, p.muStar, rv0, rvf);

% --- (3) min-ENERGY solve at this tf (loose warm start = a genuine move) --
fprintf('\n-- solving min-energy (eps = 1) to seed the homotopy --\n');
oE = casadi_minfuel_sundman(sigma, tf, rv0, rvf, p.Tmax, p.c, p.muStar, ...
                            X0, U0, tauf0, pSund, maxIter, 1, false);
fprintf('   energy: success=%d  maxDefect=%.2e  edge=%.1f%%\n', ...
        oE.success, oE.maxDefect, 100*oE.edge);
if ~(oE.success && oE.maxDefect < 1e-6)
    error('solve_minfuel:energyFailed', ...
        ['min-energy did not converge at tf_factor=%.3f. The coast is likely too ' ...
         'large (outside the solvable band ~1.12-1.25, 1.85). Pick a tf_factor in ' ...
         'that band, or build an energy-backbone seed at this tf first.'], tf_factor);
end

% --- (4) energy -> fuel homotopy: sharpen to bang-bang min-fuel ----------
fprintf('\n-- sweeping energy -> fuel (eps: 1 -> 0) --\n');
epsSched = [1 0.6 0.35 0.2 0.12 0.07 0.04 0.03 0.022 0.016 0.012 0.009 ...
            0.006 0.004 0.0025 0.0015 0.001 0.0005 0];
if saveResult
    if ~exist(cfg.dirs.root, 'dir'), mkdir(cfg.dirs.root); end
    saveFile = fullfile(cfg.dirs.root, sprintf('solve_minfuel_f%04d.mat', round(1000*tf_factor)));
else
    saveFile = '';
end
[best, tbl] = sundman_homotopy(p, rv0, rvf, sigma, oE.X, oE.U, tauf0, pSund, ...
                               epsSched, maxIter, saveFile);

%% ---------------------------- report -----------------------------------
dV_kms   = p.c * log(1/best.mf) * p.lStar / p.tStar;
prop_kg  = p.m0kg * (1 - best.mf);
certOK   = isfield(best,'certified') && best.certified;

fprintf('\n=========================== RESULT ============================\n');
fprintf('  transfer time   : %.4f ND  (%.2f days,  %.3f x min-time)\n', tf, tf*p.tStar/86400, tf_factor);
fprintf('  delta-V         : %.4f km/s\n', dV_kms);
fprintf('  propellant      : %.4f kg  (final mass fraction %.4f)\n', prop_kg, best.mf);
fprintf('  throttle switches: %d   (bang-bang fraction %.1f%%)\n', best.switches, 100*best.edge);
fprintf('  max defect      : %.2e   primer alignment: %.3f deg\n', best.maxDefect, best.primerAlignDeg);
fprintf('  certified (tight): %d\n', certOK);
if saveResult, fprintf('  saved to        : %s\n', saveFile); end
fprintf('===============================================================\n');
if ~certOK
    warning('solve_minfuel:uncertified', ...
        'no homotopy step converged tight -- treat this solution as UNCERTIFIED.');
end

%% ---------------------------- plots ------------------------------------
if doPlot
    tDays = best.X(8,:) * p.tStar / 86400;     % physical time per node [days]
    s     = best.U(4,:);                        % throttle history

    figure('Name', sprintf('min-fuel tf_factor=%.3f', tf_factor), 'Color', 'w');

    subplot(1,2,1);
    plot(tDays, s, '-', 'LineWidth', 1.2); grid on;
    xlabel('time [days]'); ylabel('throttle s'); ylim([-0.05 1.05]);
    title(sprintf('throttle  (%d switches, %.4f km/s)', best.switches, dV_kms));

    subplot(1,2,2);
    plot(best.X(1,:), best.X(2,:), '-', 'LineWidth', 0.8); hold on; grid on; axis equal;
    plot(-p.muStar,   0, 'b.', 'MarkerSize', 22);          % Earth
    plot(1-p.muStar,  0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', [.5 .5 .5]);  % Moon
    plot(rv0(1), rv0(2), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');  % start (GTO)
    plot(rvf(1), rvf(2), 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % target (tulip)
    xlabel('x (rotating, ND)'); ylabel('y (rotating, ND)');
    legend({'trajectory','Earth','Moon','start','target'}, 'Location', 'best');
    title('GTO \rightarrow south-pole tulip');

    if saveResult
        pngFile = fullfile(cfg.dirs.root, sprintf('solve_minfuel_f%04d.png', round(1000*tf_factor)));
        exportgraphics(gcf, pngFile, 'Resolution', 150);
        fprintf('  figure saved to : %s\n', pngFile);
    end
end
