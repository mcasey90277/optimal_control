%% V2_CAMPAIGN_JOB  Round 2 of the 24x24 library under the corrected certifier:
% rebuild the arrival sheet (every crossing re-certified with the tight
% pointwise flight), reuse the 18 valid ribs copied in, and hand the rest
% (column 21 re-walk, plus any newly certified column) to the supervised
% campaign, which finalizes when done.
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
addpath(ind, '/Users/msc/Desktop/optimal_control/orbit_transfer/costate_common');
out = run_costate_library(struct('outDir', fullfile(ind, 'results_fine_v2'), 'nD', 24, 'nA', 24, 'sA0', 0.0754, ...
    'nWorkers', 2, 'launch', true, ...
    'run', struct('sheet', true, 'ribs', true, 'package', true, 'audit', true, 'sweep', true)));
fprintf('\nV2 ENTRY: state %s\n', out.state);  fprintf('  blocker: %s\n', out.blockers{:});
fprintf('V2 ENTRY DONE\n');
