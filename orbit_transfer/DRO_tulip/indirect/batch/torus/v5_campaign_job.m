%% V5_CAMPAIGN_JOB  Round 5: the arrival sheet rebuilt with the direct18 arcs
% (the fourth branch, FINDINGS 66) beside the other six, each column keeping
% its fastest certified solution; ribs reused from rounds 2-4 for every
% column whose spine is unchanged; the rest walked by the supervised
% campaign; the packager given all earlier rounds' ribs so per cell the
% fastest certified point of any family is kept. Log: results_fine_v5/round5.log
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
addpath(ind, '/Users/msc/Desktop/optimal_control/orbit_transfer/costate_common', '/Users/msc/Desktop/optimal_control/orbit_transfer/campaign_common');
V5 = fullfile(ind, 'results_fine_v5');  V4 = fullfile(ind, 'results_fine_v4');
V3 = fullfile(ind, 'results_fine_v3');  V2 = fullfile(ind, 'results_fine_v2');
if ~isfolder(V5), mkdir(V5); end
logF = fullfile(V5, 'round5.log');  lg = @(varargin) logline(logF, sprintf(varargin{:}));
assert(isfile(fullfile(ind, 'results', 'arrival_arc_direct18_up_long.mat')) && ...
       isfile(fullfile(ind, 'results', 'arrival_arc_direct18_dn_long.mat')), 'direct18 arcs not on disk yet');
% stage 1: the sheet, from every arc in results/ (direct18 included) + the seeds
out1 = run_costate_library(struct('outDir', V5, 'nD', 24, 'nA', 24, 'sA0', 0.0754, 'launch', false, ...
    'run', struct('sheet', true, 'ribs', false, 'package', false, 'audit', false, 'sweep', false)));
lg('sheet stage: state %s', out1.state);
S5 = load(out1.sheet);  S5 = S5.S;  S4 = load(fullfile(V4, 'arrival_sheet_70mN_nA24.mat'));  S4 = S4.S;
changed = [];  reuse = [];
for j = 1:numel(S5.sA)
    a = S4.TF(j);  b = S5.TF(j);
    if isfinite(b) && (~isfinite(a) || abs(b - a) > 1e-6), changed(end+1) = j; tag = 'CHANGED';
    elseif isfinite(b), reuse(end+1) = j; tag = 'same';
    else, tag = 'none'; end
    lg('col %2d sA %.4f  v4 %s  v5 %s  %s', j, S5.sA(j), fmt_num(a, 8, 3), fmt_num(b, 8, 3), tag);
end
% reuse ribs for unchanged columns: the newest round's file first
for j = reuse
    for src = {fullfile(V4, sprintf('fine_rib_col%02d.mat', j)), fullfile(V3, sprintf('fine_rib_col%02d.mat', j)), ...
               fullfile(V2, sprintf('fine_rib_col%02d.mat', j))}
        dst = fullfile(V5, sprintf('fine_rib_col%02d.mat', j));
        if isfile(src{1}) && ~isfile(dst), copyfile(src{1}, dst); lg('reused %s', src{1}); break, end
    end
end
lg('columns to walk: %s', mat2str(changed));
extra = {};
for V = {V2, V3, V4}
    d_ = dir(fullfile(V{1}, 'fine_rib_col*.mat'));  extra = [extra, fullfile(V{1}, {d_.name})];
end
out = run_costate_library(struct('outDir', V5, 'nD', 24, 'nA', 24, 'sA0', 0.0754, 'nWorkers', 4, 'launch', true, ...
    'extraRibFiles', {extra}, 'run', struct('sheet', false, 'ribs', true, 'package', true, 'audit', true, 'sweep', true)));
lg('ROUND 5 ENTRY: state %s; blockers: %s', out.state, strjoin(out.blockers, ' | '));
lg('ROUND 5 LAUNCH DONE');
function logline(f, s)
fid = fopen(f, 'a');  fprintf(fid, '%s %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), s);  fclose(fid);  fprintf('%s\n', s);
end
