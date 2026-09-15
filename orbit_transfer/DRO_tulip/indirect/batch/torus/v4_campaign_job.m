%% V4_CAMPAIGN_JOB  Round 4: the arrival sheet rebuilt from ALL arcs (both
% families) plus the exact-phase direct-found seeds, each column keeping its
% fastest certified solution; ribs reused from rounds 2/3 for every column
% whose spine is unchanged; the rest walked by the supervised campaign; the
% packager given both earlier rounds' ribs so per cell the fastest certified
% point of any family is kept. Log: results_fine_v4/round4.log
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
addpath(ind, '/Users/msc/Desktop/optimal_control/orbit_transfer/costate_common');
V4 = fullfile(ind, 'results_fine_v4');  V3 = fullfile(ind, 'results_fine_v3');  V2 = fullfile(ind, 'results_fine_v2');
if ~isfolder(V4), mkdir(V4); end
logF = fullfile(V4, 'round4.log');  lg = @(varargin) logline(logF, sprintf(varargin{:}));
tStar = 382981.289129055;
% stage 1: the sheet, from every arc in results/ (the fast2 arcs included) + the seeds
out1 = run_costate_library(struct('outDir', V4, 'nD', 24, 'nA', 24, 'sA0', 0.0754, 'launch', false, ...
    'run', struct('sheet', true, 'ribs', false, 'package', false, 'audit', false, 'sweep', false)));
lg('sheet stage: state %s', out1.state);
S4 = load(out1.sheet);  S4 = S4.S;  S3 = load(fullfile(V3, 'arrival_sheet_70mN_nA24.mat'));  S3 = S3.S;
changed = [];  reuse = [];
for j = 1:numel(S4.sA)
    a = S3.TF(j);  b = S4.TF(j);
    if isfinite(b) && (~isfinite(a) || abs(b - a) > 1e-6), changed(end+1) = j; tag = 'CHANGED'; %#ok<SAGROW>
    elseif isfinite(b), reuse(end+1) = j; tag = 'same'; %#ok<SAGROW>
    else, tag = 'none'; end
    lg('col %2d sA %.4f  v3 %s  v4 %s  %s', j, S4.sA(j), fmt_num(a, 8, 3), fmt_num(b, 8, 3), tag);
end
% reuse ribs for unchanged columns: prefer round 3's file, else round 2's
for j = reuse
    for src = {fullfile(V3, sprintf('fine_rib_col%02d.mat', j)), fullfile(V2, sprintf('fine_rib_col%02d.mat', j))}
        dst = fullfile(V4, sprintf('fine_rib_col%02d.mat', j));
        if isfile(src{1}) && ~isfile(dst), copyfile(src{1}, dst); lg('reused %s', src{1}); break, end
    end
end
lg('columns to walk: %s', mat2str(changed));
extra = {};
for V = {V2, V3}
    d_ = dir(fullfile(V{1}, 'fine_rib_col*.mat'));  extra = [extra, fullfile(V{1}, {d_.name})]; %#ok<AGROW>
end
out = run_costate_library(struct('outDir', V4, 'nD', 24, 'nA', 24, 'sA0', 0.0754, 'nWorkers', 4, 'launch', true, ...
    'extraRibFiles', {extra}, 'run', struct('sheet', false, 'ribs', true, 'package', true, 'audit', true, 'sweep', true)));
lg('ROUND 4 ENTRY: state %s; blockers: %s', out.state, strjoin(out.blockers, ' | '));
lg('ROUND 4 LAUNCH DONE');
function logline(f, s)
fid = fopen(f, 'a');  fprintf(fid, '%s %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), s);  fclose(fid);  fprintf('%s\n', s);
end
