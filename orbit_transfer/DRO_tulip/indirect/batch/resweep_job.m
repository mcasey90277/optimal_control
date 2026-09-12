% RESWEEP_JOB  70 mN chain, live path, BESIDE the shipped files (2026-09-11).
%   Sheet + package + audit + second-order SWEEP on; arcs and ribs reused
%   from results/. Outputs go to results_resweep/; the shipped results/ are
%   not touched. The sweep uses the v3 sidecar (resolved conj_spectrum:
%   located minima, shifted grids, UNRESOLVED class, endpoint refinement),
%   so every one of the 115 entries is re-measured. Pictures off (made
%   afterwards); deliverable off (Mike's decision).
%   Launch with batch/run_resweep.sh, never inline. Verdict = VERDICT.txt.
root = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
outDir = fullfile(root, 'results_resweep');
if ~isfolder(outDir), mkdir(outDir); end
diary(fullfile(outDir, 'resweep_diary.log'));  diary on
fprintf('JOB START %s\n', char(datetime('now')));
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
cd(root);
chainOverrides = struct('outDir', 'results_resweep', 'run', struct( ...
    'arcs', false, 'sheet', true, 'ribs', false, 'package', true, 'audit', true, ...
    'sweep', true, 'pictures', false, 'deliverable', false));
try
    build_70mN_library
    vtxt = 'CHAIN DONE';
catch ME
    vtxt = ['CHAIN FAILED: ' ME.identifier ' -- ' ME.message];
    fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
end
fid = fopen(fullfile(outDir, 'VERDICT.txt'), 'w');
fprintf(fid, '%s  %s\n', vtxt, char(datetime('now')));  fclose(fid);
fprintf('%s\n', vtxt);  diary off
