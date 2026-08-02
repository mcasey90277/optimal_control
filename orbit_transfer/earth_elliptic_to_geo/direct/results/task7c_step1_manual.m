% TASK7C_STEP1_MANUAL  Manual continuation of the 1 N npr=15 anchor PAST
% run_mintime_mee's automatic decadeMin stall guard (same escape used in
% Task 7b at N=1104, here at the much cheaper N=662). Resumes from the last
% cached automatic round (B_round02, defect 5.096e-01) and keeps chaining
% casadi_lt_mee calls with the SAME feasibility-selected barrier policy
% (warmTight iff benign status AND defect<1e-4) but with NO decadeMin floor
% enforcement -- any forward progress is retained, a true regression falls
% back to the best point so far and retries loose. Saves every round to its
% own file (task7c_manual_round%02d.mat) so a crash loses at most one round.
here = '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo';
cd(here);
addpath(here);

resDir = fullfile(here, 'results');
S25 = load(fullfile(resDir, 'MEE_ladder_T25.mat'));
prevAnchor = S25.rung.anchor;
thrustN = 1;
par = kepler_lt_params(thrustN, 1500, 2000);

% resume point: round02's outNew (best point on disk from the automatic driver)
S2 = load(fullfile(resDir, 'MEE_mintime_T10_npr15_B_round02.mat'));
out = S2.outNew;
x0state = out.X(:,1);
fprintf('RESUME from B_round02: defect=%.4e termErr=%.4e status=%s\n', ...
        out.maxDefect, out.termErr, out.ipoptStatus);

isGood = @(o) o.success && o.maxDefect < 1e-8 && o.termErr < 1e-8;

maxRoundsManual = 60;
maxWallS = 2.5*3600;   % leave margin inside the recipe's 3 h total budget
t0 = tic;
rnd = 0;
best = out;
lastDefects = [];
while ~isGood(out) && rnd < maxRoundsManual && toc(t0) < maxWallS
    rnd = rnd + 1;
    rf = fullfile(resDir, sprintf('task7c_manual_round%02d.mat', rnd));
    if isfile(rf)
        S = load(rf); outNew = S.outNew;
        fprintf('  [cached] manual round %d\n', rnd);
    else
        warmTight = any(strcmp(out.ipoptStatus, {'Solve_Succeeded', ...
            'Solved_To_Acceptable_Level', 'Maximum_Iterations_Exceeded'})) && out.maxDefect < 1e-4;
        maxIterR = 75;
        if out.maxDefect < 1e-3, maxIterR = 150; end   % cheap here; allow more once close
        rt0 = tic;
        outNew = casadi_lt_mee(linspace(0,1,size(out.X,2)).', out.X, out.U, out.dL, ...
            struct('par', par, 'mode', 'mintime', 'x0', x0state, 'maxIter', maxIterR, ...
            'warmTight', warmTight, 'printLevel', 3));
        rtElapsed = toc(rt0);
        save(rf, 'outNew', 'rtElapsed', 'warmTight', 'maxIterR');
    end
    fprintf(['MANUAL round %d (warmTight=%d, maxIter=%d): defect %.4e -> %.4e, termErr=%.4e, ' ...
             'status=%s (%.1f s elapsed total)\n'], rnd, warmTight, maxIterR, ...
            out.maxDefect, outNew.maxDefect, outNew.termErr, outNew.ipoptStatus, toc(t0));
    if outNew.maxDefect < best.maxDefect
        best = outNew;
    end
    out = outNew;   % always advance the warm-start chain (no floor gate here)
end

fprintf('\nMANUAL CONTINUATION STOPPED after %d round(s), %.1f s wall.\n', rnd, toc(t0));
fprintf('Final: defect=%.4e termErr=%.4e status=%s certified=%d\n', ...
        out.maxDefect, out.termErr, out.ipoptStatus, isGood(out));
fprintf('Best-ever: defect=%.4e termErr=%.4e status=%s\n', best.maxDefect, best.termErr, best.ipoptStatus);
save(fullfile(resDir, 'task7c_step1_manual_final.mat'), 'out', 'best', 'rnd');
