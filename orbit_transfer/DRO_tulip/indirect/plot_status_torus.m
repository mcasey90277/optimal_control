function H = plot_status_torus(cat, opts)
%% Purpose:
%
%   Draw the phase torus of the library: for every (departure phase,
%   arrival phase) cell, the flight time of the primary transfer, how many
%   ALTERNATIVE transfers the cell holds, and which cells hold a status-2
%   transfer (a conjugate point on the arc: not locally optimal). Rows of
%   the alternatives that sit off the grid (iD or iA = 0) are drawn as
%   markers at their phases and counted separately, never folded into a
%   cell. Three panels: 'tf', 'nalt', 'conj'.
%
%  ASSUMPTIONS / NOTES:
%
% • Cells are matched by the stored indices (iD, iA), not by comparing
%   phases. Phases are FRACTIONS of a period.
% • A catalog without .alternatives gives empty 'nalt'/'conj' panels
%   carrying a note.
% • The figure is created with 'Visible' set from opts.visible, so nothing
%   appears in batch.
%
%% Inputs:
%
%  cat                      struct                  phase catalog (.sheets(1)
%                                                   .sD_frac .sA_frac .has_solution
%                                                   .entry_index .z8, .constants.tStar_s,
%                                                   optional .alternatives)
%
%  opts                     struct                  optional: .out PNG path ('')
%                                                   .visible [true]
%                                                   .panels [{'tf','nalt','conj'}]
%
%% Outputs:
%
%  H                        struct                  .fig .nAltPerCell [nD x nA]
%                                                   .conjCell [nD x nA logical]
%                                                   .nOffGrid .nOffGridConj
%
%% Revision History:
%  M. Casey                                                   (c) 10/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'out'),     opts.out = '';                      end
if ~isfield(opts, 'visible'), opts.visible = true;                end
if ~isfield(opts, 'panels'),  opts.panels = {'tf', 'nalt', 'conj'}; end

S   = cat.sheets(1);
sD  = S.sD_frac(:).';   sA = S.sA_frac(:).';
nD  = numel(sD);        nA = numel(sA);

% primaries' flight time in days:
tf = nan(nD, nA);
has = logical(S.has_solution);
tf(has) = S.z8(8, S.entry_index(has)) * cat.constants.tStar_s / 86400;

% alternatives:
nAlt   = zeros(nD, nA);
conj   = false(nD, nA);
offD   = [];  offA = [];  offSt = [];
hasAlt = isfield(cat, 'alternatives') && ~isempty(cat.alternatives);
if hasAlt
    A  = cat.alternatives(:);
    aD = [A.iD].';  aA = [A.iA].';  aS = [A.status].';
    on = aD >= 1 & aA >= 1;
    for ka = find(on).'
        nAlt(aD(ka), aA(ka)) = nAlt(aD(ka), aA(ka)) + 1;
        if aS(ka) == 2, conj(aD(ka), aA(ka)) = true; end
    end
    offD = [A(~on).sD].';  offA = [A(~on).sA].';  offSt = aS(~on);
end
nOff     = numel(offD);
nOffConj = nnz(offSt == 2);

H = struct('fig', [], 'nAltPerCell', nAlt, 'conjCell', conj, 'nOffGrid', nOff, 'nOffGridConj', nOffConj);

if opts.visible, vis = 'on'; else, vis = 'off'; end
H.fig = figure('Position', [100 100 1280 720], 'Color', 'w', 'Visible', vis);
np = numel(opts.panels);
for kp = 1:np
    ax = subplot(1, np, kp, 'Parent', H.fig);
    switch opts.panels{kp}
        case 'tf'
            imagesc(ax, sA, sD, tf);  colormap(ax, parula);  cb = colorbar(ax);  cb.Label.String = 't_f [days]';  cb.Color = 'k';
            title(ax, sprintf('primary t_f [d], %d cells', nnz(has)));
        case 'nalt'
            imagesc(ax, sA, sD, nAlt);  colormap(ax, flipud(gray));  cb = colorbar(ax);  cb.Label.String = 'alternatives';  cb.Color = 'k';
            clim(ax, [0 max(1, max(nAlt(:)))]);  hold(ax, 'on');
            if nOff > 0, plot(ax, offA, offD, 'r.', 'MarkerSize', 10); end
            title(ax, sprintf('alternatives: %d on grid, %d off grid (red)', sum(nAlt(:)), nOff));
            if ~hasAlt, text(ax, 0.5, 0.5, 'no alternatives in this catalog', 'HorizontalAlignment', 'center'); end
        case 'conj'
            imagesc(ax, sA, sD, double(conj));  colormap(ax, [1 1 1; 0.85 0.15 0.15]);  clim(ax, [0 1]);  hold(ax, 'on');
            if nOffConj > 0, plot(ax, offA(offSt == 2), offD(offSt == 2), 'k.', 'MarkerSize', 10); end
            title(ax, sprintf('status 2: %d cells, %d off grid (black)', nnz(conj), nOffConj));
            if ~hasAlt, text(ax, 0.5, 0.5, 'no alternatives in this catalog', 'HorizontalAlignment', 'center'); end
        otherwise
            error('plot_status_torus:panel', 'unknown panel %s', opts.panels{kp});
    end
    set(ax, 'YDir', 'normal', 'XLim', [-0.5/nA 1-0.5/nA], 'YLim', [-0.5/nD 1-0.5/nD], 'Box', 'on', 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    set([ax.Title ax.XLabel ax.YLabel], 'Color', 'k');
    xlabel(ax, 'arrival phase s_A [period]');  ylabel(ax, 'departure phase s_D [period]');
end
drawnow;
if ~isempty(opts.out)
    exportgraphics(H.fig, opts.out, 'Resolution', 100);
end
end
