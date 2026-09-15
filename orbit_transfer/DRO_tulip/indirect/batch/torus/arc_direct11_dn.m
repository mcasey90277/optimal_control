%% ARC_DIRECT18_dn  Pseudo-arclength continuation of the branch through the
% direct-found root at sA 0.4921 (17.96 d) in the dn arrival direction,
% recording crossings of the 24-level ladder. the direct probe at the fourth branch's fold (FINDINGS 68) found
% this root 2.3 d under the fast family; the arc says whether it joins the fourth branch's
% fast sheet past its 0.59 fold cluster or is a family of its own.
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
addpath(ind, '/Users/msc/Desktop/optimal_control/orbit_transfer/costate_common');
cd(ind);  capped_pool(2);
sA0 = 0.492066666666667;
so = struct('thrustN', 0.070, 'ispS', 900, 'm0kg', 150, 'tauDRO', 1.0, 'NpTulip', 7, 'pmTulip', -1, 'sD', 0, ...
            'sA0', sA0, 'anchorMat', fullfile(ind, 'results', 'mintime_70mN_anchor_direct11.mat'));
[~, anc] = arclength_arrival('setup', so);
ao = so;  ao.direction = -1;  ao.sAStop = sA0 + (-1)*1.15;
ao.levels = 0.0754 + (-24:48)/24;  ao.nStep = 4000;  ao.deadlineSec = 6*3600;
ao.logFile = fullfile(ind, 'results', 'arrival_arc_direct11_dn_long.log');
ao.partialFile = fullfile(ind, 'results', 'arrival_arc_direct11_dn_long.partial.mat');  ao.saveEvery = 50;
t0 = tic;
A = arclength_arrival(anc, ao);
save(fullfile(ind, 'results', 'arrival_arc_direct11_dn_long.mat'), 'A', '-v7.3');
if isfile(ao.partialFile), delete(ao.partialFile); end
fprintf('ARC direct11 dn DONE: %d roots, sA %.4f -> %.4f, %d folds, %d crossings, stop = %s (%.0f s)\n', ...
    numel(A.q), A.q(1), A.q(end), numel(A.folds), numel(A.crossings), A.stop, toc(t0));
