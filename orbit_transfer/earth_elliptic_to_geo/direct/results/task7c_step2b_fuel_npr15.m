% TASK7C_STEP2B_FUEL_NPR15  1 N fuel solve, small-N-first (npr=15, ~1044
% nodes) warm from the 2.5 N fuel solution -- the npr=25 direct attempt
% (task7c_step2_fuel.m) got to essentially machine precision at eps=0.001
% (defect=5.12e-14) but then REGRESSED hard at the final eps=0 step
% (defect=6.31e-02, landing back in the same bad discrete-switch basin as
% the very first direct attempt) -- not certified. Retrying at a smaller,
% cheaper node density per the Task 7c recipe's small-N-first fallback.
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

fuelCfg = struct('thrustN', thrustN, 'ctf', 1.5, 'tfMinAnchor', tfMinAnchor1N, ...
    'tag', 'MEE_M2_1N_npr15', 'seedThr', 0.4, 'betaMode', 'tangential', ...
    'nodesPerRev', 15, 'maxIter', 1500, 'm0kg', 1500, 'ispS', 2000);
fuelCfg.warmStart = struct('sigma', prevFuel.sigma, 'X', prevFuel.fuel.X, ...
    'U', prevFuel.fuel.U, 'dL', dLGuessFuel);

t0 = tic;
res = run_transfer_mee(fuelCfg);
fprintf('STEP2b (npr=15) done in %.1f s: certified=%d defect=%.3e m_f=%.4f sw=%d revs=%.4f edge=%.2f%%\n', ...
        toc(t0), res.report.certified, res.report.defect, res.report.m_f_kg, ...
        res.report.switches, res.report.revs, 100*res.report.edge);
save(fullfile(resDir, 'task7c_step2b_out.mat'), 'res');
