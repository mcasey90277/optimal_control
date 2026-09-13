%% FINE_LIBRARY_FINISH_JOB  Package, audit and sweep the 24 x 24 library, then
% render its full phase-sweep movie. Run by fine_library_autochain.sh once all
% 19 certified rib columns are on disk and no rib worker is alive.
%
% Everything goes through run_costate_library, so the same barrier applies:
% it refuses to package unless every column loads as a rib and no queue claim
% is live. The movie is rendered only from a catalog THIS call produced.
% A one-line-per-fact verdict goes to results_fine/FINISH_VERDICT.txt.
root = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
outDir = fullfile(root, 'results_fine');
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
addpath(root, fullfile(fileparts(fileparts(root)), 'costate_common'));
v = {sprintf('FINISH JOB START %s', char(datetime('now')))};
try
    out = run_costate_library(struct('outDir', outDir, 'nD', 24, 'nA', 24, 'launch', false, ...
        'run', struct('sheet', false, 'calibrate', false, 'ribs', true, ...
                      'package', true, 'audit', true, 'sweep', true)));
    v{end+1} = sprintf('LIBRARY state: %s  catalog: %s', out.state, out.catalog);
    for kb = 1:numel(out.blockers), v{end+1} = ['  blocker: ' out.blockers{kb}]; end
catch ME
    out = struct('state', 'threw', 'catalog', '');
    v{end+1} = sprintf('LIBRARY THREW: %s (%s)', ME.message, ME.identifier);
    fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
end
if strcmp(out.state, 'packaged')
    try
        tM = tic;
        M = movie_phase_sweep(struct('catMat', out.catalog, 'sD', 'all', 'sA', 'all', ...
            'order', 1, 'slow', 2, 'outStem', fullfile(outDir, 'sweep_full_library_24x24')));
        v{end+1} = sprintf('MOVIE: %d frames (%d certified, %d gaps) in %.0f s -> %s.{mp4,gif}', ...
                           numel(M.titles), M.nCertified, M.nGaps, toc(tM), ...
                           fullfile(outDir, 'sweep_full_library_24x24'));
    catch ME
        v{end+1} = sprintf('MOVIE THREW: %s (%s)', ME.message, ME.identifier);
    end
else
    v{end+1} = 'MOVIE not rendered: the library did not package in this call';
end
v{end+1} = sprintf('FINISH JOB END %s', char(datetime('now')));
fid = fopen(fullfile(outDir, 'FINISH_VERDICT.txt'), 'w');
fprintf(fid, '%s\n', v{:});  fclose(fid);
fprintf('%s\n', v{:});
