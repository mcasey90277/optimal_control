% FINE_RIBS_JOB  The departure axis at double resolution (24 points) off
% every certified column of the 24-level arrival sheet.
%
% BATCHED BY COLUMN, because that is the unit that can be lost: one file per
% column, a column already on disk is skipped, so the job resumes for free
% and a column that stalls costs only its own row. A wall budget stops it
% cleanly between columns rather than mid-walk.
root = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
outDir = fullfile(root, 'results_fine');
diary(fullfile(outDir, 'fine_ribs_diary.log'));  diary on
fprintf('JOB START %s\n', char(datetime('now')));
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
cd(root);  addpath(root, fullfile(fileparts(fileparts(root)), 'costate_common'));
capped_pool();

sheetMat = fullfile(outDir, 'arrival_sheet_70mN_nA24.mat');
assert(isfile(sheetMat), 'the fine sheet is missing: %s', sheetMat);
L = load(sheetMat);  S = L.S;
cols = find(isfinite(S.TF));
nD = 24;                                  % <-- the doubled departure axis
budgetSec = 8*3600;
t0 = tic;  nDone = 0;  nSkip = 0;
fprintf('fine ribs: %d certified columns, nD = %d (%d points each)\n', numel(cols), nD, nD-1);
try
    for j = cols(:)'
        f = fullfile(outDir, sprintf('fine_rib_col%02d.mat', j));
        if isfile(f), nSkip = nSkip + 1; fprintf('col %2d: reused %s\n', j, f); continue, end
        if toc(t0) > budgetSec
            fprintf('budget reached with %d column(s) left; relaunch to continue\n', ...
                    numel(cols) - nDone - nSkip);
            break
        end
        fprintf('\n=== column %d of %d (sA %.4f, t_f %.4f d) ===\n', j, numel(cols), S.sA(j), S.TF(j));
        build_ribs(sheetMat, struct('only', j, 'direction', -1, 'nD', nD, ...
                                    'nPts', nD-1, 'wallSec', 900, 'out', f));
        nDone = nDone + 1;
    end
    left = numel(cols) - nDone - nSkip;
    vtxt = sprintf('FINE RIBS: %d walked, %d reused, %d left', nDone, nSkip, left);
    if left > 0, vtxt = [vtxt ' (RELAUNCH TO CONTINUE)']; end
catch ME
    vtxt = ['FINE RIBS FAILED: ' ME.identifier ' -- ' ME.message];
    fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
end
fid = fopen(fullfile(outDir, 'RIBS_VERDICT.txt'), 'w');
fprintf(fid, '%s  %s\n', vtxt, char(datetime('now')));  fclose(fid);
fprintf('%s\n', vtxt);  diary off
