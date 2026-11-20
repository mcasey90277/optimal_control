function out = plot_pair_transfer(T, opts)
%% Purpose:
%
%   ONE FIGURE FOR ONE SOLVED TRANSFER, for any orbit pair pair_transfer
%   handles. Six tiles in a 2 x 3 grid, white background:
%
%     1  rotating-frame 3-D view (Mm) sized to the departure orbit, the
%        arrival orbit and the transfer together (extent plus 10 percent,
%        axis equal); the transfer coloured by time; Earth / Moon drawn to
%        scale only when inside the limits; departure / arrival marked
%     2  zoom on the ARRIVAL orbit (limits from its extent and the last
%        part of the transfer within twice that extent), camera along the
%        orbit normal tilted 25 deg; no Moon assumption
%     3  thrust direction components (-lambda_v / |lambda_v|) against time
%     4  altitude above the Moon and the Earth (km, log axis), the gate
%        floors as dashed lines and the LOCATED minima (T.gates) marked
%     5  mass (kg) against time
%     6  text tile: the verdict table of the gates
%
%   The same function builds the verdict table for the user script:
%
%     L = plot_pair_transfer('verdict', T)
%
%   returns the table as a cellstr (one line per entry): one row per
%   blocking gate with its measured value, its threshold and PASS / FAIL,
%   then the reported (non-blocking) second-order lines, then one line
%   (two lines) with t_f, m_f and delta-V. The rows apply the same tests as
%   pair_transfer's pass logic against T.thresholds.
%
%   Refuses an unsolved transfer (T.ok false) with
%   plot_pair_transfer:notSolved.
%
%   Usage errors (no arguments, 'verdict' without T) raise
%   plot_pair_transfer:usage.
%
%   Self-demo (nargin == 0): none (a figure needs a solved transfer;
%   pair_transfer's own self-demo makes one).
%
%  ASSUMPTIONS / NOTES:
%
% - Needs get_family_orbit's path (pumpkyn, costate_common) for the two
%   orbits; nothing is solved here.
% - The figure is created with Visible = opts.visible and stays that way.
%
%% Inputs:
%
%  T                        struct                  pair_transfer result
%                                                   (.ok .flight .gates
%                                                   .thresholds .dep .arr
%                                                   .rv0 .rvf .m0kg ...),
%                                                   or the char 'verdict'
%
%  opts                     struct (optional)       .visible ['on']
%                                                   .file ['' = do not save]
%                                                   .title ['']
%                                                   (with 'verdict': the T
%                                                   struct)
%
%% Outputs:
%
%  out                      figure handle           the figure (default);
%                           or cellstr [L x 1]      the verdict lines in
%                                                   'verdict' mode
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    error('plot_pair_transfer:usage', 'plot_pair_transfer(T, opts) needs a solved transfer (see pair_transfer)');
end
if ischar(T) || isstring(T)
    if ~strcmp(char(T), 'verdict')
        error('plot_pair_transfer:mode', 'unknown mode "%s" (verdict)', char(T));
    end
    if nargin < 2
        error('plot_pair_transfer:usage', 'plot_pair_transfer(''verdict'', T) needs the transfer T as second argument');
    end
    out = verdict_lines(opts);
    return
end
if nargin < 2, opts = struct(); end
visible = fieldd(opts, 'visible', 'on');
file    = fieldd(opts, 'file', '');
ttl     = fieldd(opts, 'title', '');
if ~(isstruct(T) && isfield(T, 'ok') && T.ok)
    why = '';
    if isstruct(T) && isfield(T, 'reason'), why = T.reason; end
    error('plot_pair_transfer:notSolved', 'transfer was not solved (reason: "%s")', why);
end

muStar = 0.012150585609624;
lStar  = 389703.264829278;                    % km
tStar  = 382981.289129055;                    % s
rMoon  = 1737.4;                              % km
rEarth = 6378.137;                            % km
moonC  = [1 - muStar, 0, 0] * lStar / 1e3;    % Mm
earthC = [-muStar, 0, 0] * lStar / 1e3;       % Mm
d2nd   = tStar / 86400;                       % ND time -> days

t  = T.flight.t(:);
Y  = T.flight.Y;
al = T.flight.alpha;
days = t * d2nd;
rr   = Y(:, 1:3);
dMoon  = sqrt(sum(((rr - moonC / lStar * 1e3) * lStar).^2, 2)) - rMoon;     % km
dEarth = sqrt(sum(((rr - earthC / lStar * 1e3) * lStar).^2, 2)) - rEarth;   % km
mass = Y(:, 7) * T.m0kg;
rrM  = rr * lStar / 1e3;                      % Mm
depM = family_xyz(T.dep) * lStar / 1e3;
arrM = family_xyz(T.arr) * lStar / 1e3;
col  = [0.00 0.35 0.70; 0.85 0.33 0.10; 0.20 0.55 0.25];     % dep, arr, transfer accent
cmap = parula(256);
midc = cmap(128, :);
tEnd = days(end);

hf = figure('Visible', visible, 'Color', 'w', 'Position', [60 60 1600 950]);
tl = tiledlayout(hf, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
if ~isempty(ttl), title(tl, ttl, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k'); end

% ---- tile 1: both orbits and the transfer, limits from their combined extent
ax = nexttile(tl);  hold(ax, 'on');
allP = [depM; arrM; rrM];
ctr  = (min(allP, [], 1) + max(allP, [], 1)) / 2;
half = 0.5 * max(max(allP, [], 1) - min(allP, [], 1)) * 1.1;
draw_body(ax, earthC, rEarth / 1e3, [0.30 0.50 0.85], ctr, half);
draw_body(ax, moonC,  rMoon / 1e3,  [0.60 0.60 0.60], ctr, half);
hDep = plot3(ax, depM(:, 1), depM(:, 2), depM(:, 3), '-', 'Color', col(1, :), 'LineWidth', 1.5);
hArr = plot3(ax, arrM(:, 1), arrM(:, 2), arrM(:, 3), '-', 'Color', col(2, :), 'LineWidth', 1.5);
scatter3(ax, rrM(:, 1), rrM(:, 2), rrM(:, 3), 6, days, 'filled');
hTr  = plot3(ax, NaN, NaN, NaN, '-', 'Color', midc, 'LineWidth', 2);
hP0  = plot3(ax, T.rv0(1) * lStar / 1e3, T.rv0(2) * lStar / 1e3, T.rv0(3) * lStar / 1e3, 'o', 'MarkerSize', 9, ...
             'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k');
hP1  = plot3(ax, T.rvf(1) * lStar / 1e3, T.rvf(2) * lStar / 1e3, T.rvf(3) * lStar / 1e3, 's', 'MarkerSize', 9, ...
             'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k');
style3(ax, 'x (Mm)', 'y (Mm)', 'z (Mm)', 'Transfer (rotating frame)');
axis(ax, [ctr(1) + [-half half], ctr(2) + [-half half], ctr(3) + [-half half]]);
colormap(ax, cmap);  clim(ax, [0 tEnd]);
cb = colorbar(ax);  cb.Label.String = 'time (d)';  cb.Color = 'k';
legend(ax, [hDep hArr hTr hP0 hP1], {[T.dep.family ' orbit'], [T.arr.family ' orbit'], 'transfer', 'departure', 'arrival'}, ...
       'Location', 'northwest', 'FontSize', 8, 'Box', 'off', 'TextColor', 'k');

% ---- tile 2: zoom on the arrival orbit, seen face-on (tilted 25 deg)
ax = nexttile(tl);  hold(ax, 'on');
c2   = mean(arrM, 1);
ext  = max(sqrt(sum((arrM - c2).^2, 2)));
dist = sqrt(sum((rrM - c2).^2, 2));
k0   = find(dist > 2 * ext, 1, 'last');
if isempty(k0), k0 = 0; end
near = rrM(k0 + 1:end, :);
half2 = 1.1 * max([ext; sqrt(sum((near - c2).^2, 2))]);
draw_body(ax, moonC,  rMoon / 1e3,  [0.60 0.60 0.60], c2, half2);
draw_body(ax, earthC, rEarth / 1e3, [0.30 0.50 0.85], c2, half2);
plot3(ax, arrM(:, 1), arrM(:, 2), arrM(:, 3), '-', 'Color', col(2, :), 'LineWidth', 1.5);
scatter3(ax, rrM(:, 1), rrM(:, 2), rrM(:, 3), 6, days, 'filled');
hTr2 = plot3(ax, NaN, NaN, NaN, '-', 'Color', midc, 'LineWidth', 2);
hA2  = plot3(ax, T.rvf(1) * lStar / 1e3, T.rvf(2) * lStar / 1e3, T.rvf(3) * lStar / 1e3, 's', 'MarkerSize', 9, ...
             'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k');
hO2  = plot3(ax, NaN, NaN, NaN, '-', 'Color', col(2, :), 'LineWidth', 1.5);
style3(ax, 'x (Mm)', 'y (Mm)', 'z (Mm)', ['Arrival ' T.arr.family ' orbit, face-on']);
axis(ax, [c2(1) + [-half2 half2], c2(2) + [-half2 half2], c2(3) + [-half2 half2]]);
clim(ax, [0 tEnd]);  colormap(ax, cmap);
legend(ax, [hO2 hTr2 hA2], {[T.arr.family ' orbit'], 'transfer', 'arrival'}, 'Location', 'northwest', ...
       'FontSize', 8, 'Box', 'off', 'TextColor', 'k');
face_on_view(ax, arrM, c2, half2, 25);

% ---- tile 3: thrust direction
ax = nexttile(tl);  hold(ax, 'on');
plot(ax, days, al(:, 1), 'LineWidth', 1.3);
plot(ax, days, al(:, 2), 'LineWidth', 1.3);
plot(ax, days, al(:, 3), 'LineWidth', 1.3);
style2(ax, 'time (d)', 'thrust direction', 'Thrust direction (rotating frame)');
xlim(ax, [0 tEnd]);
legend(ax, {'\alpha_x', '\alpha_y', '\alpha_z'}, 'Location', 'best', 'FontSize', 9, 'Box', 'off', 'TextColor', 'k');

% ---- tile 4: altitudes, floors and the LOCATED minima
ax = nexttile(tl);  hold(ax, 'on');
g  = T.gates;
th = T.thresholds;
hM = plot(ax, days, max(dMoon, 1e-3), 'LineWidth', 1.5, 'Color', col(2, :));
hE = plot(ax, days, max(dEarth, 1e-3), 'LineWidth', 1.5, 'Color', col(1, :));
hs = [hM hE];  leg = {'above Moon', 'above Earth'};
vals = [dMoon(:); dEarth(:); th.floorKm];
hs(end + 1) = yline(ax, th.floorKm, '--', 'Color', col(2, :), 'LineWidth', 1.1);
leg{end + 1} = sprintf('Moon floor %g km', th.floorKm);
if ~isempty(th.floorEarthKm)
    hs(end + 1) = yline(ax, th.floorEarthKm, '--', 'Color', col(1, :), 'LineWidth', 1.1);
    leg{end + 1} = sprintf('Earth floor %g km', th.floorEarthKm);
    vals(end + 1) = th.floorEarthKm;
end
if isfinite(g.moonMinKm) && isfinite(g.moonMinT)
    hs(end + 1) = plot(ax, g.moonMinT * d2nd, g.moonMinKm, 'v', 'MarkerSize', 8, 'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k');
    leg{end + 1} = sprintf('located Moon min %.0f km', g.moonMinKm);
    vals(end + 1) = g.moonMinKm;
end
if isfinite(g.earthMinKm) && isfinite(g.earthMinT)
    hs(end + 1) = plot(ax, g.earthMinT * d2nd, g.earthMinKm, 'v', 'MarkerSize', 8, 'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k');
    leg{end + 1} = sprintf('located Earth min %.0f km', g.earthMinKm);
    vals(end + 1) = g.earthMinKm;
end
set(ax, 'YScale', 'log');
style2(ax, 'time (d)', 'altitude (km)', 'Distance to the bodies');
xlim(ax, [0 tEnd]);
vals = vals(vals > 0);
ylim(ax, [0.5 * min(vals), 2 * max([dMoon(:); dEarth(:)])]);
legend(ax, hs, leg, 'Location', 'best', 'FontSize', 8, 'Box', 'off', 'TextColor', 'k');

% ---- tile 5: mass
ax = nexttile(tl);  hold(ax, 'on');
plot(ax, days, mass, 'LineWidth', 1.5, 'Color', col(3, :));
style2(ax, 'time (d)', 'mass (kg)', 'Spacecraft mass');
xlim(ax, [0 tEnd]);

% ---- tile 6: verdict text
ax = nexttile(tl);
axis(ax, 'off');
L = verdict_lines(T);
text(ax, 0, 1, strjoin(L, newline), 'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'FontName', 'Menlo', 'FontSize', 7, 'Color', 'k', 'Interpreter', 'none', 'Tag', 'verdictText');

if ~isempty(file)
    exportgraphics(hf, file, 'Resolution', 150, 'BackgroundColor', 'white');
end
out = hf;
end

% ------------------------------------------------------------------------
function L = verdict_lines(T)
% VERDICT_LINES  the verdict table of a pair_transfer result as a cellstr.
%   One row per blocking gate (name, measured, threshold, PASS / FAIL),
%   the reported second-order lines, one summary line.  INPUTS: T struct
%   (.ok .reason .gates .thresholds .tf_nd .tf_d .mf_kg .dV_kms).
%   OUTPUTS: L cellstr [n x 1].
g  = T.gates;
th = T.thresholds;
L  = {sprintf('%-10s %-24s %-14s %s', 'gate', 'measured', 'threshold', 'verdict')};
L{end + 1} = sprintf('%-10s %-24s %-14s %s', 'ok', ternary(T.ok, 'engine accepted', T.reason), 'accepted', pf(T.ok));
L{end + 1} = sprintf('%-10s %-24s %-14s %s', 'flown', sprintf('%.3g km', g.flownKm), sprintf('< %g km', th.gateKm), ...
                     pf(g.flownKm < th.gateKm));
L{end + 1} = sprintf('%-10s %-24s %-14s %s', 'accDz', sprintf('%.3g', g.accDz), sprintf('< %g', th.accTol), ...
                     pf(g.accDz < th.accTol));
L{end + 1} = sprintf('%-10s %-24s %-14s %s', 'admissible', ternary(g.admissible, 'yes', ['no: ' g.admissibleReason]), ...
                     'yes', pf(g.admissible));
L{end + 1} = sprintf('%-10s %-24s %-14s %s', 'moon', sprintf('%.1f km (%d failed)', g.moonMinKm, g.moonNFailed), ...
                     sprintf('>= %g km', th.floorKm), pf(g.moonMinKm >= th.floorKm && g.moonNFailed == 0));
if isempty(th.floorEarthKm)
    L{end + 1} = sprintf('%-10s %-24s %-14s %s', 'earth', sprintf('%.1f km', g.earthMinKm), 'no floor set', '--');
else
    L{end + 1} = sprintf('%-10s %-24s %-14s %s', 'earth', sprintf('%.1f km (%d failed)', g.earthMinKm, g.earthNFailed), ...
                         sprintf('>= %g km', th.floorEarthKm), pf(g.earthMinKm >= th.floorEarthKm && g.earthNFailed == 0));
end
L{end + 1} = sprintf('%-10s %s', 'BLOCKING', ternary(g.pass, 'all gates PASS', ['FAIL: ' strjoin(g.failed, ', ')]));
L{end + 1} = 'reported (not blocking):';
L{end + 1} = sprintf('  second order: %s; conjugate test %s', g.secondStatus, conj_text(g.conjPass));
L{end + 1} = sprintf('  min |lambda_v| %.4g, Q_mt %.4g, dim S %g', g.minLamV, g.minQmt, g.dimS);
L{end + 1} = sprintf('t_f %.6f ND = %.4f d', T.tf_nd, T.tf_d);
L{end + 1} = sprintf('m_f %.3f kg, delta-V %.4f km/s', T.mf_kg, T.dV_kms);
L = L(:);
end

function s = conj_text(c)
% CONJ_TEXT  conjugate-test flag as text.  INPUTS: c (1 / 0 / NaN).
%   OUTPUTS: s char.
if isnan(c), s = 'not run'; elseif c, s = 'pass'; else, s = 'fail'; end
end

function s = pf(c)
% PF  PASS / FAIL text.  INPUTS: c logical.  OUTPUTS: s char.
s = ternary(logical(c), 'PASS', 'FAIL');
end

function s = ternary(c, a, b)
% TERNARY  a if c else b.  INPUTS: c logical, a, b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

function v = fieldd(s, name, d)
% FIELDD  field of s, or the default d when absent or empty.  INPUTS: s
%   struct, name char, d default.  OUTPUTS: v.
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = d; end
end

function xyz = family_xyz(spec)
% FAMILY_XYZ  rotating-frame states of one orbit of a pair_transfer
%   family spec (.family, .params).  INPUTS: spec struct.  OUTPUTS: xyz
%   [M x 3] (ND).
p = struct();
if isfield(spec, 'params') && ~isempty(spec.params), p = spec.params; end
[~, rv] = get_family_orbit(spec.family, p);
xyz = rv(:, 1:3);
end

function draw_body(ax, c, r, rgb, ctr, half)
% DRAW_BODY  sphere of radius r at c, only when its centre lies inside the cube of half-width half around ctr.
%   INPUTS: ax axes, c [1 x 3], r scalar, rgb [1 x 3], ctr [1 x 3], half scalar.  OUTPUTS: none.
if all(abs(c - ctr) <= half), sphere_at(ax, c, r, rgb); end
end

function face_on_view(ax, orb, c, half, tiltDeg)
% FACE_ON_VIEW  point the camera along the orbit's normal (sum of cross products of successive position
%   differences), tilted by tiltDeg about an in-plane axis so depth still reads.  INPUTS: ax axes, orb [M x 3],
%   c [1 x 3] target, half scalar, tiltDeg scalar.  OUTPUTS: none.
dd = diff(orb, 1, 1);
n  = sum(cross(dd(1:end - 1, :), dd(2:end, :), 2), 1);
n  = n / sqrt(sum(n.^2));
u  = dd(1, :) - (dd(1, :) * n.') * n;                  % in-plane direction
u  = u / sqrt(sum(u.^2));
v  = cosd(tiltDeg) * n + sind(tiltDeg) * u;
up = [0 0 1] - ([0 0 1] * v.') * v;
if sqrt(sum(up.^2)) < 1e-3, up = [0 1 0] - ([0 1 0] * v.') * v; end
set(ax, 'Projection', 'orthographic', 'CameraTarget', c, 'CameraPosition', c + 8 * half * v, 'CameraUpVector', up);
end

function sphere_at(ax, c, r, rgb)
% SPHERE_AT  draw a sphere of radius r at centre c.  INPUTS: ax axes, c
%   [1 x 3], r scalar, rgb [1 x 3].  OUTPUTS: none.
[sx, sy, sz] = sphere(24);
surf(ax, c(1) + r * sx, c(2) + r * sy, c(3) + r * sz, 'FaceColor', rgb, 'EdgeColor', 'none', ...
     'FaceLighting', 'gouraud', 'AmbientStrength', 0.6);
end

function style3(ax, xl, yl, zl, ttl)
% STYLE3  common 3-D axes dressing.  INPUTS: ax axes, labels and title
%   char.  OUTPUTS: none.
axis(ax, 'equal');  grid(ax, 'on');  box(ax, 'on');  view(ax, 3);
xlabel(ax, xl);  ylabel(ax, yl);  zlabel(ax, zl);
title(ax, ttl, 'FontWeight', 'bold', 'Color', 'k');
set(ax, 'FontSize', 10, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', 'GridColor', [0.3 0.3 0.3]);
camlight(ax, 'headlight');
end

function style2(ax, xl, yl, ttl)
% STYLE2  common 2-D axes dressing.  INPUTS: ax axes, labels and title
%   char.  OUTPUTS: none.
grid(ax, 'on');  box(ax, 'on');
xlabel(ax, xl);  ylabel(ax, yl);
title(ax, ttl, 'FontWeight', 'bold', 'Color', 'k');
set(ax, 'FontSize', 10, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', 'GridColor', [0.3 0.3 0.3]);
end
