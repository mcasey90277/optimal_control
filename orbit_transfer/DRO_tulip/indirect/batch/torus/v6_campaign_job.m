%% V6_CAMPAIGN_JOB  Round 6: the arrival sheet rebuilt with the direct11 arcs
% (the fifth branch, FINDINGS 68) beside the other twelve, each column keeping
% its fastest certified solution; ribs reused from rounds 2-5 for every
% column whose spine is unchanged; the rest walked by the supervised
% campaign; the packager given all earlier rounds' ribs so per cell the
% fastest certified point of any family is kept. Log: results_fine_v6/round6.log
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
addpath(ind, '/Users/msc/Desktop/optimal_control/orbit_transfer/costate_common', '/Users/msc/Desktop/optimal_control/orbit_transfer/campaign_common');
V6 = fullfile(ind, 'results_fine_v6');  V5 = fullfile(ind, 'results_fine_v5');  V4 = fullfile(ind, 'results_fine_v4');
V3 = fullfile(ind, 'results_fine_v3');  V2 = fullfile(ind, 'results_fine_v2');
if ~isfolder(V6), mkdir(V6); end
logF = fullfile(V6, 'round6.log');  lg = @(varargin) logline(logF, sprintf(varargin{:}));
assert(isfile(fullfile(ind, 'results', 'arrival_arc_direct11_up_long.mat')) && ...
       isfile(fullfile(ind, 'results', 'arrival_arc_direct11_dn_long.mat')), 'direct11 arcs not on disk yet');
% stage 1: the sheet, from every arc in results/ (direct11 included) + the seeds
out1 = run_costate_library(struct('outDir', V6, 'nD', 24, 'nA', 24, 'sA0', 0.0754, 'launch', false, ...
    'run', struct('sheet', true, 'ribs', false, 'package', false, 'audit', false, 'sweep', false)));
lg('sheet stage: state %s', out1.state);
S6 = load(out1.sheet);  S6 = S6.S;  S5 = load(fullfile(V5, 'arrival_sheet_70mN_nA24.mat'));  S5 = S5.S;
changed = [];  reuse = [];
for j = 1:numel(S6.sA)
    a = S5.TF(j);  b = S6.TF(j);
    if isfinite(b) && (~isfinite(a) || abs(b - a) > 1e-6), changed(end+1) = j; tag = 'CHANGED';
    elseif isfinite(b), reuse(end+1) = j; tag = 'same';
    else, tag = 'none'; end
    lg('col %2d sA %.4f  v5 %s  v6 %s  %s', j, S6.sA(j), fmt_num(a, 8, 3), fmt_num(b, 8, 3), tag);
end
% reuse ribs for unchanged columns: the newest round's file first
for j = reuse
    for src = {fullfile(V5, sprintf('fine_rib_col%02d.mat', j)), fullfile(V4, sprintf('fine_rib_col%02d.mat', j)), ...
               fullfile(V3, sprintf('fine_rib_col%02d.mat', j)), fullfile(V2, sprintf('fine_rib_col%02d.mat', j))}
        dst = fullfile(V6, sprintf('fine_rib_col%02d.mat', j));
        if isfile(src{1}) && ~isfile(dst), copyfile(src{1}, dst); lg('reused %s', src{1}); break, end
    end
end
lg('columns to walk: %s', mat2str(changed));
extra = {};
for V = {V2, V3, V4, V5}
    d_ = dir(fullfile(V{1}, 'fine_rib_col*.mat'));  extra = [extra, fullfile(V{1}, {d_.name})];
end
out = run_costate_library(struct('outDir', V6, 'nD', 24, 'nA', 24, 'sA0', 0.0754, 'nWorkers', 4, 'launch', true, ...
    'extraRibFiles', {extra}, 'run', struct('sheet', false, 'ribs', true, 'package', true, 'audit', true, 'sweep', true)));
lg('ROUND 6 ENTRY: state %s; blockers: %s', out.state, strjoin(out.blockers, ' | '));
lg('ROUND 6 LAUNCH DONE');
function logline(f, s)
fid = fopen(f, 'a');  fprintf(fid, '%s %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), s);  fclose(fid);  fprintf('%s\n', s);
end
