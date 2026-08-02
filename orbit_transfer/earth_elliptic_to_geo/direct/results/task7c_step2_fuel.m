% TASK7C_STEP2_FUEL  1 N fuel (eps=0) solve at c_tf=1.5, npr=25 (~1740
% nodes), warm-started from the certified 2.5 N fuel solution via the C-law
% dL rescale -- mirrors run_ladder.m's own fuel-warm-start chain exactly,
% but driven standalone (not through run_ladder, whose hard assert on fuel
% certification would kill the whole call rather than letting us apply the
% small-N-first fallback per the Task 7c recipe).
here = '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo';
cd(here); addpath(here);
resDir = fullfile(here, 'results');

S1anchor = load(fullfile(resDir, 'MEE_mintime_T10.mat'));
tfMinAnchor1N = S1anchor.out.tfmin;

Sfuel25 = load(fullfile(resDir, 'MEE_M2_2p5N.mat'));
prevFuel = Sfuel25.res;

thrustN = 1; prevThrust = 2.5;
dLGuessFuel = prevFuel.fuel.dL * (prevThrust / thrustN);
fprintf('C-law fuel dL rescale: dL(2.5N)=%.6f -> dLGuess(1N)=%.6f rad -> revsGuess=%.4f\n', ...
        prevFuel.fuel.dL, dLGuessFuel, dLGuessFuel/(2*pi));
fprintf('tfMinAnchor(1N) = %.6f ND (from certified step-1b anchor)\n', tfMinAnchor1N);

fuelCfg = struct('thrustN', thrustN, 'ctf', 1.5, 'tfMinAnchor', tfMinAnchor1N, ...
    'tag', mee_fuel_tag(thrustN), 'seedThr', 0.4, 'betaMode', 'tangential', ...
    'nodesPerRev', 25, 'maxIter', 1500, 'm0kg', 1500, 'ispS', 2000);
fuelCfg.warmStart = struct('sigma', prevFuel.sigma, 'X', prevFuel.fuel.X, ...
    'U', prevFuel.fuel.U, 'dL', dLGuessFuel);

t0 = tic;
res = run_transfer_mee(fuelCfg);
fprintf('STEP2 (npr=25 direct) done in %.1f s: certified=%d defect=%.3e m_f=%.4f sw=%d revs=%.4f edge=%.2f%%\n', ...
        toc(t0), res.report.certified, res.report.defect, res.report.m_f_kg, ...
        res.report.switches, res.report.revs, 100*res.report.edge);
save(fullfile(resDir, 'task7c_step2_out.mat'), 'res');
