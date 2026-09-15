%% V3_SHEET_JOB  Round 3, stage 1 only: rebuild the arrival sheet with the seven
% direct-found certified solutions of the faster family seeded in beside the
% arcs' crossings. Each column keeps its fastest certified solution. Ribs are
% decided afterwards from which columns' spines changed.
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
addpath(ind, '/Users/msc/Desktop/optimal_control/orbit_transfer/costate_common');
out = run_costate_library(struct('outDir', fullfile(ind, 'results_fine_v3'), 'nD', 24, 'nA', 24, 'sA0', 0.0754, ...
    'launch', false, 'run', struct('sheet', true, 'ribs', false, 'package', false, 'audit', false, 'sweep', false)));
fprintf('\nV3 SHEET: state %s\n', out.state);  fprintf('  blocker: %s\n', out.blockers{:});
S = load(out.sheet);  S = S.S;  tStar = 382981.289129055;
S2 = load(fullfile(ind, 'results_fine_v2', 'arrival_sheet_70mN_nA24.mat'));  S2 = S2.S;
fprintf('\n col   sA     v2 t_f[d]   v3 t_f[d]   change\n');
for j = 1:numel(S.sA)
    a = S2.TF(j);  b = S.TF(j);  tag = '';
    if isfinite(b) && ~isfinite(a), tag = 'NEW'; elseif isfinite(b) && isfinite(a) && b < a - 1e-6, tag = sprintf('FASTER by %.2f d', a - b); elseif isfinite(a) && ~isfinite(b), tag = 'LOST'; end
    fprintf(' %2d  %.4f   %9s   %9s   %s\n', j, S.sA(j), fmt_num(a, 9, 3), fmt_num(b, 9, 3), tag);
end
fprintf('V3 SHEET DONE\n');
