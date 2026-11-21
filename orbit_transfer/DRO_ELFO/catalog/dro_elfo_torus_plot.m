function H = dro_elfo_torus_plot(catalog, raanDeg, thrustN, opts)
%% Purpose:
%
%   THE PHASE TORUS of a DRO -> ELFO catalog at one thrust rung.
%
%   ONE SHEET (raanDeg a number): two panels over (arrival phase,
%   departure phase): the flight time in days of every cell holding an
%   entry, empty cells left blank (transparent over a hatched grey
%   background, so they read as "no entry", never as a colour of the
%   scale); and the verdict map, four classes: pass / refuted / not run /
%   no entry.
%
%   ALL SHEETS (raanDeg = 'all'): small multiples of the flight-time panel,
%   one tile per sheet of that DRO size (ordered by raanDeg), on ONE shared
%   colour scale.
%
%  ASSUMPTIONS / NOTES:
%
% - Sheets are matched exactly by dro_elfo_find_sheet (no fallback;
%   dro_elfo_torus_plot:noSheet). The DRO size is opts.droTau, or the
%   catalog's only one; a catalog with several DRO sizes needs opts.droTau
%   (dro_elfo_torus_plot:droTau).
% - thrustN must be a rung within 1e-9 N (dro_elfo_torus_plot:noRung).
% - Verdict codes in H.verdict: 1 pass, 0 refuted, -1 not run, NaN no
%   entry (a sheet without conj_pass gives -1 at every entry).
%
%% Inputs:
%
%  catalog                  struct | char           The catalog struct, or
%                                                   the path of its .mat
%
%  raanDeg                  double | char           ELFO plane orientation,
%                                                   deg [1 x 1], or 'all'
%
%  thrustN                  double                  Thrust rung, N [1 x 1]
%
%  opts                     struct (optional)       .droTau [the catalog's
%                                                   only DRO size] .out ['']
%                                                   PNG path .visible [true]
%
%% Outputs:
%
%  H                        struct                  [1 x 1] .fig .sD [1 x nD]
%                                                   .sA [1 x nA] .thrustN
%                                                   .raanDeg [1 x nS]
%                                                   .sheets [1 x nS]
%                                                   .tfDays [nD x nA x nS]
%                                                   days, NaN = no entry
%                                                   .verdict [nD x nA x nS]
%                                                   (1 / 0 / -1 / NaN)
%                                                   .clim [1 x 2] days
%
%% Revision History:
%  M. Casey                                                   (c) 11/21/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    catalog = fullfile(here, 'results', 'costate_catalog_dro_elfo.mat');
    if ~isfile(catalog), catalog = fullfile(here, 'results', 'pilot_chain', 'costate_catalog_dro_elfo.mat'); end
    if ~isfile(catalog), fprintf('dro_elfo_torus_plot: self-demo needs a catalog under results/\n');  H = [];  return, end
    H = dro_elfo_torus_plot(catalog, 0, 1);
    return
end
if nargin < 4 || isempty(opts), opts = struct(); end
if ischar(catalog) || isstring(catalog)
    L = load(char(catalog));  fn = fieldnames(L);
    assert(isscalar(fn), 'dro_elfo_torus_plot:catalog', 'expected one variable in %s', char(catalog));
    catalog = L.(fn{1});
end
if ~(isfield(catalog, 'sheets') && ~isempty(catalog.sheets) && isfield(catalog.sheets, 'sheet_key'))
    error('dro_elfo_torus_plot:catalog', 'the catalog has no sheets carrying sheet_key');
end
outPng  = getf(opts, 'out', '');
visible = getf(opts, 'visible', true);
keysAll = [catalog.sheets.sheet_key];
droTau  = getf(opts, 'droTau', []);
if isempty(droTau)
    uT = unique([keysAll.depSize]);
    if ~isscalar(uT)
        error('dro_elfo_torus_plot:droTau', 'the catalog holds %d DRO sizes (%s): pass opts.droTau', ...
              numel(uT), mat2str(uT, 6));
    end
    droTau = uT;
end
kr = find(abs(catalog.rungs_N - thrustN) < 1e-9, 1);
if isempty(kr)
    error('dro_elfo_torus_plot:noRung', 'thrust %g N is not a catalog rung -- rungs: %s N', ...
          thrustN, mat2str(catalog.rungs_N));
end

%% Sheets drawn:
allMode = (ischar(raanDeg) || isstring(raanDeg)) && strcmpi(raanDeg, 'all');
if allMode
    sel = find(arrayfun(@(k) strcmp(k.target, 'elfo') && abs(k.depSize - droTau) < 1e-9, keysAll));
    if isempty(sel)
        error('dro_elfo_torus_plot:noSheet', 'no ELFO sheet with depSize %.12g', droTau);
    end
    [~, ord] = sort(mod([keysAll(sel).raanDeg], 360));
    sel = sel(ord);
elseif isnumeric(raanDeg) && isscalar(raanDeg)
    sel = dro_elfo_find_sheet(catalog, droTau, raanDeg, 'dro_elfo_torus_plot');
else
    error('dro_elfo_torus_plot:raanDeg', 'raanDeg must be a number or ''all''');
end

%% Matrices:
sh1 = catalog.sheets(sel(1));
sD  = sh1.sD_frac(:).';  sA = sh1.sA_frac(:).';
nD  = numel(sD);  nA = numel(sA);  nS = numel(sel);
tf  = nan(nD, nA, nS);  vd = nan(nD, nA, nS);
d2day = catalog.constants.tStar_s/86400;
for kk = 1:nS
    sh  = catalog.sheets(sel(kk));
    if ~isequal(size(sh.has_solution, 1, 2), [nD nA])
        error('dro_elfo_torus_plot:grid', 'sheet %d has a %d x %d phase grid, sheet %d has %d x %d', ...
              sel(kk), size(sh.has_solution, 1), size(sh.has_solution, 2), sel(1), nD, nA);
    end
    has = logical(sh.has_solution(:, :, kr));
    ek  = sh.entry_index(:, :, kr);
    t1  = nan(nD, nA);  t1(has) = sh.z8(8, ek(has))*d2day;
    v1  = nan(nD, nA);
    if isfield(sh, 'conj_pass') && ~isempty(sh.conj_pass)
        cp = double(sh.conj_pass(:, :, kr));
        cp(cp ~= 0 & cp ~= 1) = -1;
        v1(has) = cp(has);
    else
        v1(has) = -1;
    end
    tf(:, :, kk) = t1;  vd(:, :, kk) = v1;
end
fin = tf(isfinite(tf));
if isempty(fin), cl = [0 1]; else, cl = [min(fin) max(fin)]; end
if cl(2) <= cl(1), cl = cl + [-0.5 0.5]; end
H = struct('fig', [], 'sD', sD, 'sA', sA, 'thrustN', catalog.rungs_N(kr), ...
           'raanDeg', [keysAll(sel).raanDeg], 'sheets', sel, 'tfDays', tf, 'verdict', vd, 'clim', cl);

%% Draw:
if visible, vis = 'on'; else, vis = 'off'; end
H.fig = figure('Position', [100 100 1280 720], 'Color', 'w', 'Visible', vis);
if ~allMode
    k  = sh1.sheet_key;
    ax = subplot(1, 2, 1, 'Parent', H.fig);
    tf_panel(ax, sA, sD, tf(:, :, 1), cl);
    title(ax, sprintf('t_f [d] at %g N: %d cells with an entry, %d empty (hatched)', ...
                      H.thrustN, nnz(isfinite(tf)), nnz(~isfinite(tf))));
    ax = subplot(1, 2, 2, 'Parent', H.fig);
    code = 4*ones(nD, nA);                            % 4 no entry
    v1 = vd(:, :, 1);
    code(v1 == -1) = 3;  code(v1 == 0) = 2;  code(v1 == 1) = 1;
    imagesc(ax, sA, sD, code);
    colormap(ax, [0.20 0.65 0.30; 0.85 0.15 0.15; 0.95 0.80 0.20; 0.92 0.92 0.92]);
    clim(ax, [0.5 4.5]);
    cb = colorbar(ax, 'Ticks', 1:4, 'TickLabels', {'pass', 'refuted', 'not run', 'no entry'});
    cb.Color = 'k';
    title(ax, sprintf('conjugate test: %d pass, %d refuted, %d not run', nnz(v1 == 1), nnz(v1 == 0), nnz(v1 == -1)));
    axes_style(ax, nD, nA);
    sgtitle(H.fig, sprintf('DRO %.4g ND -> ELFO raan %g deg, %g N', k.depSize, k.raanDeg, H.thrustN));
else
    nc = ceil(sqrt(nS*16/9));  nr = ceil(nS/nc);
    tl = tiledlayout(H.fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
    for kk = 1:nS
        ax = nexttile(tl);
        tf_panel(ax, sA, sD, tf(:, :, kk), cl, false);
        title(ax, sprintf('raan %g deg: %d cells', H.raanDeg(kk), nnz(isfinite(tf(:, :, kk)))), 'FontSize', 9);
    end
    cb = colorbar(ax);  cb.Layout.Tile = 'east';  cb.Label.String = 't_f [d]';  cb.Color = 'k';
    title(tl, sprintf('DRO %.4g ND -> ELFO, every plane orientation at %g N (one colour scale; hatched = no entry)', ...
                      droTau, H.thrustN));
end
drawnow;
if ~isempty(outPng)
    od = fileparts(outPng);
    if ~isempty(od) && ~isfolder(od), mkdir(od); end
    set(H.fig, 'PaperPositionMode', 'auto');
    print(H.fig, outPng, '-dpng', '-r0');
end
end

% ------------------------------------------------------------------------
function tf_panel(ax, sA, sD, t1, cl, withBar)
% TF_PANEL  One flight-time panel: entries coloured, empty cells
% transparent over a hatched grey background.  INPUTS: ax [axes]; sA
% [1 x nA]; sD [1 x nD]; t1 [nD x nA] days (NaN = empty); cl [1 x 2];
% withBar [logical, default true].  OUTPUTS: none.
if nargin < 6, withBar = true; end
nD = numel(sD);  nA = numel(sA);
hold(ax, 'on');
x0 = -0.5/nA;  x1 = 1 - 0.5/nA;  y0 = -0.5/nD;  y1 = 1 - 0.5/nD;
patch(ax, [x0 x1 x1 x0], [y0 y0 y1 y1], [0.80 0.80 0.80], 'EdgeColor', 'none');
for c = linspace(-1, 1, 25)                       % hatching under the empty cells
    plot(ax, [x0 x1], c + [x0 x1]*(y1 - y0)/(x1 - x0), '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.5);
end
hI = imagesc(ax, sA, sD, t1);
hI.AlphaData = double(isfinite(t1));
colormap(ax, parula);  clim(ax, cl);
if withBar, cb = colorbar(ax);  cb.Label.String = 't_f [d]';  cb.Color = 'k'; end
axes_style(ax, nD, nA);
end

% ------------------------------------------------------------------------
function axes_style(ax, nD, nA)
% AXES_STYLE  Limits, labels and colours of a torus panel.  INPUTS: ax
% [axes]; nD, nA [1 x 1] grid size.  OUTPUTS: none.
set(ax, 'YDir', 'normal', 'XLim', [-0.5/nA 1-0.5/nA], 'YLim', [-0.5/nD 1-0.5/nD], 'Box', 'on', ...
        'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'Layer', 'top');
xlabel(ax, 'ELFO arrival phase s_A [period]');  ylabel(ax, 'DRO departure phase s_D [period]');
set([ax.Title ax.XLabel ax.YLabel], 'Color', 'k');
end

% ------------------------------------------------------------------------
function v = getf(s, fld, d)
% GETF  Field with default.  INPUTS: s [struct]; fld [char]; d default.
% OUTPUTS: v the field value, or d when absent or empty.
if isfield(s, fld) && ~isempty(s.(fld)), v = s.(fld); else, v = d; end
end
