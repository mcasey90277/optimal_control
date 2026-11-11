function [T, info] = combined_timeline(leg1Mat, leg2Mat, opts)
%% Purpose:
%
%   THE COMBINED LEG-1 + LEG-2 TIMELINE of the lunar PNT mission, per
%   satellite: GTO wait, leg 1 (GTO -> DRO slot, 5 N), DRO coast to the
%   hand-off epoch T0, DRO wait, leg 2 (DRO -> tulip slot, 70 mN), on
%   station to the end of the leg-2 service window. Mission time is the
%   LEG-1 clock: leg-2 times (its own clock, t = 0 at the hand-off) are
%   shifted by T0. Satellite k of leg 1 is spacecraft k of leg 2.
%
%   COMPATIBILITY FIRST (Review Focus 4); nothing is drawn or returned
%   unless every check holds:
%     identity    a package file whose SHA-256 equals the leg-1 package
%                 identity stored in the leg-1 results (run_meta.identity.
%                 sha256) is found by leg1_package_for, HASH ONLY: the file
%                 at opts.catalogMat (default: the GTO_DRO package path),
%                 else run_meta.packageSnapshot, else a snapshot in
%                 inputs/ beside the leg-1 results or in this campaign's
%                 results/inputs; none -> refused unless
%                 opts.allowStalePackage (then info.packageMatch says so)
%     dro         the leg-1 DRO has period 1 (run_meta.tauD, 1e-9) and its
%                 MD5 at 64 phases (run_meta.droMd5) equals the leg-2
%                 fingerprint's droMd5 (the same recipe, catalog_identity)
%     droPhase0   the leg-2 fingerprint's droPhase0 equals (0:5)/6 and the
%                 leg-1 run_meta.droPhase0 (1e-12)
%     units       leg-1 run_meta.regime.tStar equals the leg-2 tStar
%     verdict     the leg-1 study's run_meta.verdict is exactly 'PASS'
%     incomplete  six leg-1 transfers, every one accepted with G.pass true;
%                 six leg-2 ones
%     leg2        the leg-2 results passed: its gate table (the leg-2
%                 study saves no verdict string; its VERDICT is PASS iff
%                 every gate row is PASS) is non-empty, every status is
%                 'PASS', and it holds the validation gate G0
%     epoch       T0 >= every leg-1 arrival (1e-12)
%   Each refusal throws combined_timeline:<check>.
%
%   Self-demo (nargin == 0): the saved results of both studies, the figure
%   written to results/mission_timeline.png, the table printed.
%
%% Inputs:
%
%  leg1Mat                  char                    lunar_pnt_leg1_study
%                                                   results (P1 [1 x 6] with
%                                                   .tDep .tArr .accepted
%                                                   .G.pass, T0, run_meta
%                                                   with .verdict)
%
%  leg2Mat                  char                    lunar_pnt_mission_study
%                                                   results (P [1 x 6] with
%                                                   .tDep .tArr,
%                                                   fingerprint, H.tEnd,
%                                                   gates {n x 4})
%
%  opts                     struct (optional)       .catalogMat package file
%                                                   [GTO_DRO/results/
%                                                   costate_catalog_gto_dro
%                                                   .mat]
%                                                   .allowStalePackage [false]
%                                                   .outPng [''] figure path
%                                                   ('' = no figure)
%                                                   .visible [false]
%                                                   .keepOpen [false] leave
%                                                   the figure open
%
%% Outputs:
%
%  T                        table [6 rows]          .sat, and [6 x 2] start /
%                                                   end pairs in DAYS of
%                                                   mission (leg-1) time:
%                                                   .gtoWait .leg1 .droCoast
%                                                   .droWait .leg2 .onStation
%
%  info                     struct                  .T0Days .endDays .tStar
%                                                   .packageMatch
%                                                   .packageFile ('' if none)
%                                                   .packageSource ('primary'
%                                                   | 'extra' | 'snapshot' |
%                                                   '') .checks
%                                                   (char, one line each)
%                                                   .fig (handle when
%                                                   keepOpen, else [])
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  M. Casey  fix round 1: verdict, G.pass and leg-2 gate checks (c) 11/11/2026
%  M. Casey  final fix wave: package found by SHA-256 with a
%            snapshot fallback (leg1_package_for)             (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(fileparts(mfilename('fullpath')));
if nargin == 0
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    leg1Mat = fullfile(here, 'results', 'lunar_pnt_leg1_study.mat');
    leg2Mat = fullfile(here, 'results', 'lunar_pnt_mission_study.mat');
    opts = struct('outPng', fullfile(here, 'results', 'mission_timeline.png'));
end
if nargin < 3, opts = struct(); end
known = {'catalogMat', 'allowStalePackage', 'outPng', 'visible', 'keepOpen'};
for f = fieldnames(opts).'
    if ~any(strcmp(f{1}, known))
        error('combined_timeline:opts', 'unknown option "%s"', f{1});
    end
end
catalogMat = fdef(opts, 'catalogMat', fullfile(fileparts(here), 'GTO_DRO', 'results', ...
                                               'costate_catalog_gto_dro.mat'));
allowStale = fdef(opts, 'allowStalePackage', false);
outPng     = fdef(opts, 'outPng', '');
visible    = fdef(opts, 'visible', false);
keepOpen   = fdef(opts, 'keepOpen', false);

%% Load:
A = load(leg1Mat, 'P1', 'T0', 'run_meta');
B = load(leg2Mat, 'P', 'fingerprint', 'H', 'gates');
rm = A.run_meta;  fp = B.fingerprint;
chk = {};

%% Compatibility:
snapX = {};
if isfield(rm, 'packageSnapshot') && ~isempty(rm.packageSnapshot), snapX = {char(rm.packageSnapshot)}; end
inDirs = unique({fullfile(fileparts(canon_path(leg1Mat)), 'inputs'), fullfile(here, 'results', 'inputs')}, 'stable');
[pkgFile, pkgWhy, pkgInfo] = leg1_package_for(rm.identity, struct('primary', catalogMat, 'extra', {snapX}, ...
                                                                  'inputsDirs', {inDirs}));
packageMatch = ~isempty(pkgFile);
if isfile(catalogMat), shaNow = file_sha256(catalogMat); else, shaNow = 'missing'; end
if ~packageMatch && ~allowStale
    error('combined_timeline:identity', ['the leg-1 results were seeded from package %s... but the file ' ...
          'now at %s is %s... and no snapshot matches (%s); refusing (opts.allowStalePackage to override)'], ...
          rm.identity.sha256(1:12), catalogMat, shaNow(1:min(12, end)), pkgWhy);
end
if packageMatch
    chk{end+1} = sprintf('identity: package SHA-256 %s... matches (%s: %s)', rm.identity.sha256(1:12), ...
                         pkgInfo.source, pkgFile);
else
    chk{end+1} = sprintf('identity: package SHA-256 %s... NOT FOUND (allowStalePackage; %s is %s...)', ...
                         rm.identity.sha256(1:12), catalogMat, shaNow(1:min(12, end)));
end
if ~(abs(rm.tauD - 1) < 1e-9)
    error('combined_timeline:dro', 'the leg-1 DRO period is %.12g ND, not 1', rm.tauD);
end
if ~strcmp(rm.droMd5, fp.droMd5)
    error('combined_timeline:dro', 'the leg-2 DRO (MD5 %s) is not the leg-1 DRO (MD5 %s)', fp.droMd5, rm.droMd5);
end
chk{end+1} = sprintf('dro: period 1, MD5 at 64 phases %s on both legs', rm.droMd5);
d0 = (0:5)/6;
if ~(isequal(size(fp.droPhase0), [1 6]) && max(abs(fp.droPhase0 - d0)) < 1e-12 && ...
     max(abs(rm.droPhase0(:).' - fp.droPhase0)) < 1e-12)
    error('combined_timeline:droPhase0', 'leg-2 droPhase0 [%s] is not (0:5)/6 / the leg-1 hand-off phases [%s]', ...
          num2str(fp.droPhase0, 6), num2str(rm.droPhase0(:).', 6));
end
chk{end+1} = 'droPhase0: (0:5)/6 on both legs';
tStar = rm.regime.tStar;
if ~(abs(tStar - fp.tStar) < 1e-6)
    error('combined_timeline:units', 'leg-1 tStar %.9g s differs from leg-2 tStar %.9g s', tStar, fp.tStar);
end
v1 = '';
if isfield(rm, 'verdict'), v1 = char(rm.verdict); end
if ~strcmp(v1, 'PASS')
    error('combined_timeline:verdict', 'the leg-1 study verdict is "%s", not PASS', v1);
end
gOk = arrayfun(@(q) isfield(q, 'G') && isstruct(q.G) && isfield(q.G, 'pass') && isequal(q.G.pass, true), A.P1);
if ~(numel(A.P1) == 6 && all([A.P1.accepted]) && all(gOk) && numel(B.P) == 6)
    error('combined_timeline:incomplete', ['need six accepted leg-1 transfers with G.pass and six leg-2 ones ' ...
          '(accepted %d of %d, G.pass %d, leg 2 %d)'], nnz([A.P1.accepted]), numel(A.P1), nnz(gOk), numel(B.P));
end
g2 = B.gates;
if ~(iscell(g2) && size(g2, 1) > 0 && all(strcmp(g2(:, 2), 'PASS')) && any(strcmp(g2(:, 1), 'G0')))
    if iscell(g2) && size(g2, 1) > 0
        badG = g2(~strcmp(g2(:, 2), 'PASS'), 1);
        why = sprintf('not PASS: [%s]; G0 present %d', strjoin(badG.', ' '), any(strcmp(g2(:, 1), 'G0')));
    else
        why = 'no gate table';
    end
    error('combined_timeline:leg2', 'the leg-2 results did not pass (%s)', why);
end
chk{end+1} = sprintf('verdicts: leg 1 PASS, six G.pass; leg 2 all %d gates PASS incl. G0', size(g2, 1));
tArr1 = [A.P1.tArr];
if ~(A.T0 >= max(tArr1) - 1e-12)
    error('combined_timeline:epoch', 'T0 %.12g ND precedes the last leg-1 arrival %.12g ND', A.T0, max(tArr1));
end
chk{end+1} = sprintf('units, six transfers per leg, T0 >= last leg-1 arrival (by %.4g ND)', A.T0 - max(tArr1));

%% Table (days, leg-1 clock):
dayNd = tStar/86400;
T0 = A.T0;  tEnd = T0 + B.H.tEnd;
td1 = [A.P1.tDep].';  ta1 = tArr1.';
td2 = T0 + [B.P.tDep].';  ta2 = T0 + [B.P.tArr].';
z6 = zeros(6, 1);
T = table((1:6).', [z6, td1]*dayNd, [td1, ta1]*dayNd, [ta1, T0 + z6]*dayNd, ...
          [T0 + z6, td2]*dayNd, [td2, ta2]*dayNd, [ta2, tEnd + z6]*dayNd, ...
          'VariableNames', {'sat', 'gtoWait', 'leg1', 'droCoast', 'droWait', 'leg2', 'onStation'});
info = struct('T0Days', T0*dayNd, 'endDays', tEnd*dayNd, 'tStar', tStar, ...
              'packageMatch', packageMatch, 'packageFile', pkgFile, 'packageSource', pkgInfo.source, ...
              'checks', {chk}, 'fig', []);

%% Figure:
if ~isempty(outPng) || keepOpen
    fig = figure('Theme', 'light', 'Color', 'w', 'Position', [60 60 1400 640], ...
                 'Visible', tern(visible, 'on', 'off'));
    closeFig = onCleanup(@() close_unless(fig, keepOpen));
    ax = axes(fig);  hold(ax, 'on');
    FS = 14;
    ph = {'gtoWait', 'GTO wait'; 'leg1', 'leg 1 (5 N)'; 'droCoast', 'DRO coast to T0'; ...
          'droWait', 'DRO wait'; 'leg2', 'leg 2 (70 mN)'; 'onStation', 'on station'};
    col = [0.70 0.70 0.70; 0.85 0.33 0.10; 0.20 0.40 0.80; 0.62 0.72 0.88; 0.93 0.60 0.20; 0.20 0.60 0.30];
    for kp = 1:size(ph, 1)
        seg = T.(ph{kp, 1});
        for k = 1:6
            hv = tern(k == 1, 'on', 'off');
            plot(ax, seg(k, :), [k k], '-', 'Color', col(kp, :), 'LineWidth', 16, ...
                 'DisplayName', ph{kp, 2}, 'HandleVisibility', hv);
        end
    end
    xline(ax, info.T0Days, '--k', 'LineWidth', 1.4, 'DisplayName', sprintf('hand-off T0 = %.2f d', info.T0Days));
    set(ax, 'YTick', 1:6, 'YTickLabel', compose('sat %d', 1:6), 'YDir', 'reverse', 'FontSize', FS, 'Color', 'w');
    xlim(ax, [0 info.endDays]);  ylim(ax, [0.4 6.6]);  grid(ax, 'on');
    xlabel(ax, 'mission time [d] (leg-1 clock)');
    legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', FS - 1, 'NumColumns', 4);
    title(ax, sprintf('Lunar PNT mission: GTO -> DRO (leg 1, done at %.2f d) -> tulip (leg 2, last arrival %.2f d)', ...
                      info.T0Days, max(T.leg2(:, 2))), 'FontSize', FS + 1);
    if ~isempty(outPng)
        exportgraphics(fig, outPng, 'Resolution', 150);
    end
    if keepOpen, info.fig = fig; end
end

if nargin == 0
    disp(T);
    fprintf('%s\n', chk{:});
end
end

% ------------------------------------------------------------------------
function close_unless(fig, keep)
% CLOSE_UNLESS  Close the figure this function opened unless kept.
%   INPUTS: fig handle; keep logical.  OUTPUTS: none.
if ~keep && isgraphics(fig), close(fig); end
end

% ------------------------------------------------------------------------
function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c logical; a, b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
