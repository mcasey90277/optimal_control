% FINE_SHEET_JOB  The arrival axis at double resolution (24 levels), built
% from the SAVED arcs by re-scanning their stored paths -- no arc is walked
% again. Writes to results_fine/; nothing existing is touched.
root = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
outDir = fullfile(root, 'results_fine');
if ~isfolder(outDir), mkdir(outDir); end
diary(fullfile(outDir, 'fine_sheet_diary.log'));  diary on
fprintf('JOB START %s\n', char(datetime('now')));
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
cd(root);  addpath(root, fullfile(fileparts(fileparts(root)), 'costate_common'));
capped_pool();
try
    S = build_arrival_sheet(struct('nA', 24, 'sA0', 0.0754, 'rescan', true, ...
                                   'out', fullfile(outDir, 'arrival_sheet_70mN_nA24.mat')));
    nc = nnz(isfinite(S.TF));
    fprintf('\nFINE SHEET: %d of %d arrival phases certified (coarse sheet had 11 of 12)\n', nc, numel(S.sA));
    for j = 1:numel(S.sA)
        c = S.cand{j};
        if isfinite(S.TF(j)), tfs = sprintf('%8.4f d', S.TF(j)); else, tfs = '       - '; end
        fprintf('  j %2d  sA %.4f  cand %2d  cert %2d  t_f %s\n', j, S.sA(j), numel(c), ...
                nnz([c.ok]), tfs);
    end
    vtxt = sprintf('FINE SHEET DONE: %d of %d columns certified', nc, numel(S.sA));
catch ME
    vtxt = ['FINE SHEET FAILED: ' ME.identifier ' -- ' ME.message];
    fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
end
fid = fopen(fullfile(outDir, 'VERDICT.txt'), 'w');
fprintf(fid, '%s  %s\n', vtxt, char(datetime('now')));  fclose(fid);
fprintf('%s\n', vtxt);  diary off
