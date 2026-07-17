function fig_basin_scatter()
% FIG_BASIN_SCATTER  Free-longitude local-minima scatter figure (Fig-18 analog).
%
% Scans results/ for every certified c_tf=1.2..3.0 sweep point at T_max=10 N,
% including the SUPERSEDED chain-seeded basins alongside the redo-set FRESH
% basins (run_ctf_sweep.m's basin-scatter finding: for c_tf in {1.2, 2.5,
% 3.0} a second cold full-schedule solve lands in a different, generally
% better, local optimum than the neighbor-chain seed -- this is the same
% phenomenon the paper's own Fig 18 documents). Each certified point is
% plotted individually (chain vs fresh marker styles), plus the single
% uncertified c_tf=1.2 chain attempt pulled from the run log (never cached,
% campaign rule "never cache uncertified" -- run_transfer.m), plus the
% best-of envelope line that run_ctf_sweep.m's front figure plots.
%
% INPUTS:  none
% OUTPUTS: none (writes results/fig_basin_scatter.png)
%
% REFERENCES:
%   [1] run_ctf_sweep.m header (basin-scatter finding, best-of envelope).
%   [2] results/M3_sweep.log line 59 (uncertified c_tf=1.2 chain point:
%       "DONE sweep_T100_c120: certified=0 revs=8.10 sw=10 edge=100.0%
%       mf=1322.31 kg dV=2.473 km/s apoBurn=0.92" -- uncertified results are
%       deliberately never saved to a .mat, so the log is the only source).
%   [3] paper Fig 18 (free-longitude local-minima scatter, the figure this
%       one is styled to match).
here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');

% --- certified points: {c_tf, chain tag, fresh tag} -------------------------
% Includes the 4 densify_front.m points (1.35, 1.75, 2.25, 2.75; Task 14
% cleanup pass) alongside the original 5. The densified points have no
% "_fresh" alternative on disk (only a single neighbor-seeded chain solve
% each, with the two non-monotone ones already re-seeded from the OTHER
% neighbor and kept as the best -- see densify_front.m / results/
% densify_run.log), so they show up here as chain-marker points.
ctfs = [1.20 1.35 1.50 1.75 2.00 2.25 2.50 2.75 3.00];
pts  = struct('cf', {}, 'src', {}, 'mf', {});
for cf = ctfs
    chainFn = fullfile(resDir, sprintf('sweep_T100_c%03d.mat', round(100*cf)));
    freshFn = fullfile(resDir, sprintf('sweep_T100_c%03d_fresh.mat', round(100*cf)));
    if isfile(chainFn)
        S = load(chainFn);
        if S.res.report.certified
            pts(end+1) = struct('cf', cf, 'src', 'chain', 'mf', S.res.report.m_f_kg); %#ok<AGROW>
        end
    end
    if isfile(freshFn)
        S = load(freshFn);
        if S.res.report.certified
            pts(end+1) = struct('cf', cf, 'src', 'fresh', 'mf', S.res.report.m_f_kg); %#ok<AGROW>
        end
    end
end

% --- best-of envelope (same 5 merged values as run_ctf_sweep's front) -------
envMf = nan(1, numel(ctfs));
for kc = 1:numel(ctfs)
    m = [pts([pts.cf] == ctfs(kc)).mf];
    envMf(kc) = max(m);
end

% --- the one uncertified point, hardcoded from the run log ------------------
% results/M3_sweep.log line 59: "DONE sweep_T100_c120: certified=0 revs=8.10
% sw=10 edge=100.0% mf=1322.31 kg dV=2.473 km/s apoBurn=0.92" -- run_transfer
% never saves an uncertified result (campaign rule), so this number exists
% ONLY in the log, not in any .mat file.
uncertCf = 1.20;
uncertMf = 1322.31;

% --- figure -------------------------------------------------------------------
fig = figure('Color','w','Visible','off');
try, theme(fig,'light'); catch, end
hold on; grid on; box on

hChain = []; hFresh = [];
chainMask = strcmp({pts.src}, 'chain');
freshMask = strcmp({pts.src}, 'fresh');
if any(chainMask)
    hChain = plot([pts(chainMask).cf], [pts(chainMask).mf], 's', ...
        'MarkerFaceColor',[0.10 0.35 0.85], 'MarkerEdgeColor','k', 'MarkerSize',9, 'LineStyle','none');
end
if any(freshMask)
    hFresh = plot([pts(freshMask).cf], [pts(freshMask).mf], '^', ...
        'MarkerFaceColor',[0.85 0.35 0.10], 'MarkerEdgeColor','k', 'MarkerSize',9, 'LineStyle','none');
end
hUncert = plot(uncertCf, uncertMf, 'o', 'MarkerFaceColor','none', ...
    'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize',10, 'LineWidth',1.6);
text(uncertCf + 0.04, uncertMf, 'uncertified (never cached -- from run log)', ...
    'FontSize',9, 'Color',[0.4 0.4 0.4]);

hEnv = plot(ctfs, envMf, '-', 'Color',[0.15 0.15 0.15], 'LineWidth',1.4);

legendH = [hEnv];  legendS = {'best-of envelope'};
if ~isempty(hChain), legendH(end+1) = hChain; legendS{end+1} = 'chain-seeded (certified)'; end
if ~isempty(hFresh), legendH(end+1) = hFresh; legendS{end+1} = 'fresh redo (certified)'; end
legendH(end+1) = hUncert; legendS{end+1} = 'uncertified (log only)';
legend(legendH, legendS, 'Location','southeast');

xlabel('c_{tf}');  ylabel('m_f [kg]');
title('Free-longitude local-minima scatter (T_{max}=10 N) -- cf. paper Fig 18');
xlim([1.05 3.10]);

exportgraphics(fig, fullfile(resDir, 'fig_basin_scatter.png'), 'Resolution', 150);
close(fig);
fprintf('WROTE %s\n', fullfile(resDir, 'fig_basin_scatter.png'));
end
