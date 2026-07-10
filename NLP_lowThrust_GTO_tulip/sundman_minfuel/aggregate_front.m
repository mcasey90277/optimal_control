function [pts, pmp] = aggregate_front(makePlot)
% AGGREGATE_FRONT  Collect every stored min-fuel solution, PMP-verify each,
% and draw the HONEST Delta-V vs t_f front with three marker classes.
%
% Promotes the scratchpad combine_front.m pattern into the repo (cleanup plan
% Phase 0) and implements the honest-plot policy of
% HONEST_EVALUATION_DV_TF_FRONT.md:
%   class 1  feasible upper bound   (machine-tight dynamics+BCs, NOT certified)
%   class 2  direct-certified       (KKT-dual PMP check passes: local extremal)
%   class 3  direct+indirect        (independent adjoint/shooting match; none
%                                    yet -- populated by the ms_band certifier)
% The envelope is drawn ONLY through certified points (classes 2-3). A grey
% feasible point BELOW the envelope is information (a better local minimum
% exists there), never noise -- e.g. the 1.75x/2.523 point that dominates the
% certified 1.85x/2.667.
%
% Sources scanned (all optional): results/minfuel/minfuel_f*_*.mat (new
% layout), legacy ms_*.mat and tf_front_results.mat in the library root.
% Multiple solutions at the same factor are all kept and plotted.
%
% INPUTS:
%   makePlot - true (default): write results/plots/front_honest.png
% OUTPUTS:
%   pts - struct array per solution: .factor .tf_days .dV .switches .edge
%         .defect .class (1|2|3) .source
%   pmp - verify_tf_front output for the combined set
%
% REFERENCES:
%   [1] verify_tf_front.m (empirical-beta switching-law certification).
%   [2] HONEST_EVALUATION_DV_TF_FRONT.md ("Plan to the two goals").

here = fileparts(mfilename('fullpath'));  addpath(here);
cfg  = minfuel_config();
if nargin<1 || isempty(makePlot), makePlot = true; end

flds = {'factor','tf_days','dV','switches','edge','X','U','lamDef','primerAlignDeg'};
res  = struct([]);  src = {};

% --- new-layout results -----------------------------------------------------
if exist(cfg.dirs.minfuel,'dir')
    dd = dir(fullfile(cfg.dirs.minfuel, 'minfuel_f*.mat'));
    for k = 1:numel(dd)
        R = load(fullfile(dd(k).folder, dd(k).name));
        [res, src] = add_point(res, src, R.out, flds, dd(k).name);
    end
end
% --- legacy per-factor files (ms_<f>.mat carry `out`) ------------------------
dd = dir(fullfile(here, 'ms_*.mat'));
for k = 1:numel(dd)
    R = load(fullfile(dd(k).folder, dd(k).name));
    if isfield(R,'out'), [res, src] = add_point(res, src, R.out, flds, dd(k).name); end
end
% --- legacy up-pass front (struct array `results`) ---------------------------
fLeg = fullfile(here, 'tf_front_results.mat');
if isfile(fLeg)
    R = load(fLeg);
    for k = 1:numel(R.results)
        [res, src] = add_point(res, src, R.results(k), flds, 'tf_front_results.mat');
    end
end
if isempty(res), error('aggregate_front:noData','no stored solutions found'); end

% sort by factor NOW: verify_tf_front sorts internally (stable), so pre-sorting
% here keeps pmp(k) aligned with res(k)/src{k} below.
[~, ix] = sort([res.factor]);  res = res(ix);  src = src(ix);

% --- PMP-verify the combined set ---------------------------------------------
if ~exist(cfg.dirs.fronts,'dir'), mkdir(cfg.dirs.fronts); end
combF = fullfile(cfg.dirs.fronts, 'combined_front.mat');
results = res; %#ok<NASGU>
save(combF, 'results');
pmp = verify_tf_front(combF, false);

% --- classify ----------------------------------------------------------------
pts = struct('factor',{},'tf_days',{},'dV',{},'switches',{},'edge',{}, ...
             'defect',{},'class',{},'source',{});
for k = 1:numel(res)
    cls = 1;                                   % feasible upper bound
    if pmp(k).pmpPass, cls = 2; end            % direct-certified extremal
    % class 3 (direct+indirect) is set once the ms_band certifier stamps a
    % matching indirect solution -- field reserved, no data yet.
    pts(end+1) = struct('factor',res(k).factor,'tf_days',res(k).tf_days, ...
        'dV',res(k).dV,'switches',res(k).switches,'edge',res(k).edge, ...
        'defect',NaN,'class',cls,'source',src{k}); %#ok<AGROW>
end
fprintf('\naggregate_front: %d solutions | %d certified | %d feasible-only\n', ...
        numel(pts), sum([pts.class]>=2), sum([pts.class]==1));

% --- plot ---------------------------------------------------------------------
if makePlot
    if ~exist(cfg.dirs.plots,'dir'), mkdir(cfg.dirs.plots); end
    tStar = 382981.289129055;
    d=[pts.tf_days]; v=[pts.dV]; c=[pts.class];
    fig=figure('Color','w','Position',[100 100 880 520],'Visible','off');
    try, theme(fig,'light'); catch, end
    hold on; grid on; box on;
    plot(d(c==1), v(c==1), 'o', 'Color',[0.60 0.60 0.65], ...
         'MarkerFaceColor',[0.87 0.87 0.90],'MarkerSize',7);
    plot(d(c==2), v(c==2), 'o', 'Color',[0.10 0.45 0.15], ...
         'MarkerFaceColor',[0.20 0.65 0.25],'MarkerSize',9,'LineWidth',1.2);
    plot(d(c==3), v(c==3), 's', 'Color',[0.05 0.25 0.55], ...
         'MarkerFaceColor',[0.15 0.45 0.80],'MarkerSize',10,'LineWidth',1.2);
    % envelope through the best certified point per factor
    cf=[pts.class]>=2;
    if any(cf)
        fac=[pts.factor]; ufac=unique(fac(cf)); envD=zeros(size(ufac)); envV=zeros(size(ufac));
        for kk=1:numel(ufac)
            sel = cf & fac==ufac(kk);
            [envV(kk), pos] = min(v(sel)); dd2=d(sel); envD(kk)=dd2(pos);
        end
        plot(envD, envV, '-', 'Color',[0.10 0.45 0.15 0.5], 'LineWidth',1.5);
    end
    plot(cfg.tfMin*tStar/86400, 4.4665, 'ks', 'MarkerFaceColor','k','MarkerSize',9);
    text(cfg.tfMin*tStar/86400+0.4, 4.4665, 'min-time (4.4665, 0 sw)', ...
         'FontSize',9,'Color',[0.2 0.2 0.2]);
    xlabel('transfer time t_f (days)'); ylabel('\DeltaV (km/s)');
    title('Min-fuel \DeltaV vs t_f -- honest front (envelope through certified points only)');
    legend({'feasible upper bound','direct-certified extremal', ...
            'direct+indirect certified','certified envelope','min-time'}, ...
           'Location','northwest','Box','off');
    outP = fullfile(cfg.dirs.plots, 'front_honest.png');
    exportgraphics(fig, outP, 'Resolution',150); close(fig);
    fprintf('WROTE %s\n', outP);
end
end

% ---------------------------------------------------------------------------
function [res, src] = add_point(res, src, s, flds, sourceName)
% Append one solution record, keeping only the canonical fields.
q = struct();
for k = 1:numel(flds)
    if isfield(s, flds{k}), q.(flds{k}) = s.(flds{k}); else, q.(flds{k}) = []; end
end
if isempty(res), res = q; else, res(end+1) = q; end
src{end+1} = sourceName;
end
