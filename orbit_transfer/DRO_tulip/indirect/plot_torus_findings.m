function I = plot_torus_findings(catMat, outPng)
%% Purpose:
%
%   The phase torus COLOUR-CODED BY WHAT WE FOUND: six panels over the same
%   (departure x arrival) grid, so the geometry, the cost and the strength of
%   the certificate can be read against each other.
%
%     1  transfer time        the design variable
%     2  Delta-V              what it costs
%     3  Legendre margin      min |lambda_v| -- how far from a degenerate control
%     4  switching margin     min Q_mt -- how far from wanting to coast
%     5  H6 margin            (c/T)/lambda_m(0) -- spurious-zero exclusion
%     6  rank margin          sigma_6 / measured error -- dim S = 1 certified
%
%   Panels 5 and 6 appear only once `second_order_pass` has written its
%   measurements; until then they say so rather than showing a blank grid
%   that might be mistaken for zeros.
%
%   Cells with nothing certified are GAPS in every panel, never zeros: on
%   this catalog 29 of 144 have no entry, and two whole arrival phases are
%   empty for reasons worth reading (one had both candidates refuted by the
%   conjugate test; the other stalls immediately in the departure walk).
%
%% Inputs:
%
%  catMat                   char                    catalog .mat
%  outPng                   char (optional)         save path
%
%% Outputs:
%
%  I                        struct                  .fig .nCert .nCells
%                                                   .panels (cellstr of what
%                                                   was drawable)
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, outPng = ''; end
L = load(catMat);  fn = fieldnames(L);  c = L.(fn{1});
s = c.sheets(1);
tStar = c.constants.tStar_s;  lStar = c.constants.lStar_km;
OK = s.has_solution(:,:,1);
nD = numel(s.sD_frac);  nA = numel(s.sA_frac);

tf = nan(nD, nA);  tf(OK) = s.tf_nd(OK)*tStar/86400;
dv = nan(nD, nA);
% THE shared propulsion conversion (costate_common/nd_propulsion)
ndp = nd_propulsion(c.rungs_N(1), c.thruster.isp_s, c.thruster.m0_kg, lStar, tStar);
cnd = ndp.cnd;   Tnd = ndp.Tnd;
mf = 1 - (Tnd/cnd)*s.tf_nd(:,:,1);
dv(OK) = cnd*log(1./mf(OK))*lStar/tStar;

panels = { tf, 'minimum transfer time  [days]', 'parula', false
           dv, '\DeltaV  [km/s]', 'parula', false
           grab(s,'gate_min_lamv',OK), 'Legendre margin   min |\lambda_v|', 'viridis', false
           grab(s,'gate_min_qmt',OK),  'switching margin   min Q_{mt}', 'viridis', false
           grab(s,'h6_margin',OK),     'H6 margin   (c/T) / \lambda_m(0)', 'copper', true
           grab(s,'lift_margin',OK),   'rank margin   \sigma_6 / measured error', 'copper', true };

fig = figure('Color', 'w', 'Position', [60 60 1500 880]);
I.panels = {};
for k = 1:6
    ax = subplot(2, 3, k, 'Parent', fig);
    Z = panels{k,1};
    if isempty(Z) || all(isnan(Z(:)))
        axis(ax, 'off');
        text(ax, 0.5, 0.5, {panels{k,2}, '', 'not yet measured', ...
             '(run second\_order\_pass with writeback)'}, ...
             'HorizontalAlignment', 'center', 'Color', [0.45 0.45 0.5], 'FontSize', 11);
        continue
    end
    I.panels{end+1} = panels{k,2};
    h = imagesc(ax, 1:nA, 1:nD, Z);
    set(h, 'AlphaData', ~isnan(Z));
    set(ax, 'Color', [0.92 0.92 0.94], 'YDir', 'normal', 'Layer', 'top');
    try, colormap(ax, panels{k,3}); catch, colormap(ax, parula); end
    if panels{k,4}, set(ax, 'ColorScale', 'log'); end
    cb = colorbar(ax);  cb.FontSize = 8;
    title(ax, panels{k,2}, 'FontWeight', 'normal', 'FontSize', 11);
    xticks(ax, 1:2:nA);  yticks(ax, 1:2:nD);
    xticklabels(ax, arrayfun(@(v) sprintf('%.2f', v), s.sA_frac(1:2:end), 'UniformOutput', false));
    yticklabels(ax, arrayfun(@(v) sprintf('%.2f', v), s.sD_frac(1:2:end), 'UniformOutput', false));
    if k > 3, xlabel(ax, 'arrival phase s_A'); end
    if mod(k, 3) == 1, ylabel(ax, 'departure phase s_D'); end
    grid(ax, 'on');  set(ax, 'GridColor', 'w', 'GridAlpha', 0.35);
    % mark the fastest cell on the time panel
    if k == 1
        [~, kb] = min(Z(:));  [ib, jb] = ind2sub(size(Z), kb);
        hold(ax, 'on');  plot(ax, jb, ib, 'o', 'MarkerSize', 13, 'LineWidth', 2, 'Color', [0 0.5 0]);
    end
end

I.nCert = nnz(OK);  I.nCells = numel(OK);  I.fig = fig;
sgtitle(fig, sprintf(['DRO \\rightarrow %d-petal tulip minimum time, %.0f mN / I_{sp} %g s / %g kg' ...
    '   --   %d of %d phase pairs certified,  t_f %.2f to %.2f d' ...
    '   (grey = no certified transfer)'], s.Np, c.rungs_N(1)*1000, c.thruster.isp_s, ...
    c.thruster.m0_kg, I.nCert, I.nCells, min(tf(:)), max(tf(:))), 'FontSize', 13);
if ~isempty(outPng), exportgraphics(fig, outPng, 'Resolution', 150); end
end

function Z = grab(s, f, OK)
% GRAB  A per-cell field masked to certified cells, or [] if absent.
% INPUTS: s; f; OK.  OUTPUTS: Z.
Z = [];
if ~isfield(s, f), return, end
Z = nan(size(OK));  V = s.(f);
Z(OK) = V(OK);
end
