%% ARC_FAST2_dn  Pseudo-arclength continuation of the NEW fast family (anchor
% at sA 0.8671, 17.25 d) in the dn arrival direction, recording crossings of
% the 24-level ladder, so the sheet can be rebuilt with this family in it.
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
addpath(ind, '/Users/msc/Desktop/optimal_control/orbit_transfer/costate_common');
cd(ind);  capped_pool(2);
so = struct('thrustN', 0.070, 'ispS', 900, 'm0kg', 150, 'tauDRO', 1.0, 'NpTulip', 7, 'pmTulip', -1, 'sD', 0, ...
            'sA0', 0.8671, 'anchorMat', fullfile(ind, 'results', 'mintime_70mN_anchor_fast2.mat'));
[~, anc] = arclength_arrival('setup', so);
ao = so;  ao.direction = -1;  ao.sAStop = 0.8671 + (-1)*1.15;
ao.levels = 0.0754 + (-24:48)/24;  ao.nStep = 4000;  ao.deadlineSec = 6*3600;
ao.logFile = fullfile(ind, 'results', 'arrival_arc_fast2_dn_long.log');
t0 = tic;
A = arclength_arrival(anc, ao);
save(fullfile(ind, 'results', 'arrival_arc_fast2_dn_long.mat'), 'A', '-v7.3');
fprintf('ARC fast2 dn DONE: %d roots, sA %.4f -> %.4f, %d folds, %d crossings, stop = %s (%.0f s)\n', ...
    numel(A.q), A.q(1), A.q(end), numel(A.folds), numel(A.crossings), A.stop, toc(t0));
