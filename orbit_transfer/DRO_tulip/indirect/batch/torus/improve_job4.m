%% FILL_HOLES_JOB (improve pass 4: the stuck cell (14,15) with a 900 s cap, then the chain)  Cell-by-cell direct solves at every hole of the round-6
% catalog (FINDINGS 70), each seeded from its nearest certified neighbour,
% harvested and certified; certified points land in
% results_fine_v6/fine_rib_direct_holes.mat for the next packaging.
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
addpath(ind, '/Users/msc/Desktop/optimal_control/orbit_transfer/costate_common');
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
cd(ind);
out = fill_holes_direct(fullfile(ind, 'results_fine_v6', 'costate_catalog_dro_tulip_70mN.mat'), ...
    struct('logFile', fullfile(ind, 'results_fine_v6', 'fill_holes_direct.log'), 'clearKm', 1900, 'improveDays', 2, 'maxCpuSec', 900, 'cells', [14 15]));
fprintf('FILL HOLES DONE: %d holes, %d tried, %d certified\n', out.nHoles, out.nTried, out.nCert);
