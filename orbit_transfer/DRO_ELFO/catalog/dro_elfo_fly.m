function [F, hFig] = dro_elfo_fly(catalog, entry, opts)
%% Purpose:
%
%   ONE DRO -> ELFO catalog entry, flown and handed back: the time history
%   of state, mass and thrust direction in ND and in days / km / kg, the
%   arrival miss against the ELFO state at the entry's arrival phase, the
%   located minimum altitude above the Moon, the number of revolutions
%   about the Moon, and an optional 3D picture. Nothing is solved or
%   re-certified here.
%
%   The flight is ONE call of pumpkyn's own propagator
%   (pumpkyn.cr3bp.tfMinProp, all-burn minimum time) from z8, started at
%   the DRO state at the entry's departure phase with mass fraction 1 --
%   the catalog stores no junction states. Endpoints follow the catalog's
%   orbit_reconstruction rule: get_family_orbit on the sheet's dep/arr
%   family and parameters at the catalog's muStar, then phase_state.
%
%  ASSUMPTIONS / NOTES:
%
% - The entry may be a dro_elfo_pick result, a dro_elfo_transfers_at row
%   or a dro_elfo_extremes row: each carries z8, sD, sA, thrustN and its
%   sheet_key, and the sheet is resolved from the key by
%   dro_elfo_find_sheet (exact match, no fallback), never from an index.
% - Minimum altitude: costate_common/located_min_altitude on the true
%   min-time flow (the range-rate root between samples), Moon radius
%   1737.4 km -- the same measure the chain's clearance sweep applies.
% - Revolutions: the unwrapped angle of the Moon-centred position in the
%   rotating x-y plane, start to end, over 2*pi (signed; positive is
%   counter-clockwise seen from +z).
% - The ELFO is drawn over one Kepler period from its locus; it is a fixed
%   curve of osculating states in the rotating frame, as in the catalog.
%
%% Inputs:
%
%  catalog                  struct | char           The catalog struct, or
%                                                   the path of its .mat
%
%  entry                    struct                  [1 x 1] a pick /
%                                                   transfers_at / extremes
%                                                   row: .z8 [8 x 1] .sD
%                                                   .sA .thrustN .sheet_key
%
%  opts                     struct (optional)       .plot [false] 3D picture
%                                                   .visible [true] figure
%                                                   visibility
%                                                   .out [''] PNG path
%                                                   (implies .plot)
%                                                   .title [''] extra line
%
%% Outputs:
%
%  F                        struct                  [1 x 1] .t [n x 1] ND
%                                                   .tDays [n x 1]
%                                                   .Y [n x 14] ND
%                                                   .rKm [n x 3] barycentric
%                                                   .rMoonKm [n x 3]
%                                                   Moon-centred .vKms
%                                                   [n x 3] .massKg [n x 1]
%                                                   .u [n x 3] unit thrust
%                                                   direction .tfNd .tfDays
%                                                   .dvKms .propellantKg
%                                                   .finalMassKg .flyKm
%                                                   (arrival position miss)
%                                                   .flyVms (velocity miss,
%                                                   m/s) .minMoonAltKm
%                                                   .tMinMoonDays .moonRevs
%                                                   .sD .sA .thrustN .sheet
%                                                   .sheet_key .z8 [8 x 1]
%                                                   .rv0 .rvf [6 x 1] ND
%                                                   .orbits (.tD [mD x 1]
%                                                   .rvD [mD x 6] .tA
%                                                   [mA x 1] .rvA [mA x 6])
%
%  hFig                     handle                  The figure when drawn,
%                                                   else []
%
%% Revision History:
%  M. Casey                                                   (c) 11/21/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    catalog = fullfile(here, 'results', 'costate_catalog_dro_elfo.mat');
    if ~isfile(catalog), catalog = fullfile(here, 'results', 'pilot_chain', 'costate_catalog_dro_elfo.mat'); end
    if ~isfile(catalog), fprintf('dro_elfo_fly: self-demo needs a catalog under results/\n');  F = [];  hFig = [];  return, end
    L = load(catalog);  fn = fieldnames(L);  catalog = L.(fn{1});
    k = catalog.sheets(1).sheet_key;
    entry = dro_elfo_pick(catalog, k.depSize, k.raanDeg, 0.25, 0.5, 1);
    [F, hFig] = dro_elfo_fly(catalog, entry, struct('plot', true));
    fprintf('dro_elfo_fly: t_f %.3f d, arrival miss %.2e km, min Moon altitude %.1f km, %.2f revs\n', ...
            F.tfDays, F.flyKm, F.minMoonAltKm, F.moonRevs);
    return
end
if nargin < 3 || isempty(opts), opts = struct(); end
outPng  = getf(opts, 'out', '');
doPlot  = getf(opts, 'plot', false) || ~isempty(outPng);
visible = getf(opts, 'visible', true);
if ischar(catalog) || isstring(catalog)
    L = load(char(catalog));  fn = fieldnames(L);
    assert(isscalar(fn), 'dro_elfo_fly:catalog', 'expected one variable in %s', char(catalog));
    catalog = L.(fn{1});
end
if ~(isstruct(entry) && isscalar(entry))
    error('dro_elfo_fly:entry', 'the entry must be ONE struct (index a transfers_at result: T(k))');
end
need = {'z8', 'sD', 'sA', 'thrustN', 'sheet_key'};
for kk = 1:numel(need)
    if ~isfield(entry, need{kk}) || isempty(entry.(need{kk}))
        error('dro_elfo_fly:entry', 'the entry needs .%s (a dro_elfo_pick, transfers_at or extremes row)', need{kk});
    end
end
if numel(entry.z8) ~= 8
    error('dro_elfo_fly:entry', 'the entry''s z8 must have 8 elements, got %d', numel(entry.z8));
end
ks = dro_elfo_find_sheet(catalog, entry.sheet_key.depSize, entry.sheet_key.raanDeg, 'dro_elfo_fly');
sh = catalog.sheets(ks);

%% Orbits and endpoints (the catalog's orbit_reconstruction rule):
mu    = catalog.constants.muStar;
tStar = catalog.constants.tStar_s;
lStar = catalog.constants.lStar_km;
m0kg  = catalog.thruster.m0_kg;
cnd   = catalog.thruster.c_nd;
dp = sh.dep_params;  dp.muStar = mu;
ap = sh.arr_params;  ap.muStar = mu;
[tD, rvD] = get_family_orbit(sh.dep_family, dp);
[tA, rvA] = get_family_orbit(sh.arr_family, ap);
xD  = phase_state(tD, rvD);
xA  = phase_state(tA, rvA);
sD  = mod(entry.sD, 1);  sA = mod(entry.sA, 1);
rv0 = xD(sD);  rv0 = reshape(rv0(1:6), 1, 6);      % [1 x 6] ND
rvf = xA(sA);  rvf = reshape(rvf(1:6), 1, 6);
Tnd = (entry.thrustN/m0kg)*tStar^2/(lStar*1000);

%% Fly:
z8 = entry.z8(:);
[t, Y] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(:); 1; z8(1:7)], Tnd, cnd, mu);
lv  = Y(:, 11:13);
mf  = Y(end, 7);
moonC = [1 - mu, 0, 0];
Cm  = located_min_altitude(t, Y, moonC, 1737.4, -inf, ...
                           struct('lStar', lStar, 'phys', struct('Tnd', Tnd, 'cnd', cnd, 'muStar', mu)));
ang = unwrap(atan2(Y(:, 2), Y(:, 1) - moonC(1)));

F = struct('t', t, 'tDays', t*tStar/86400, 'Y', Y, 'rKm', Y(:, 1:3)*lStar, ...
           'rMoonKm', (Y(:, 1:3) - moonC)*lStar, 'vKms', Y(:, 4:6)*lStar/tStar, ...
           'massKg', Y(:, 7)*m0kg, 'u', -lv ./ sqrt(sum(lv.^2, 2)), ...
           'tfNd', z8(8), 'tfDays', z8(8)*tStar/86400, 'dvKms', cnd*log(1/mf)*lStar/tStar, ...
           'propellantKg', (1 - mf)*m0kg, 'finalMassKg', mf*m0kg, ...
           'flyKm', sqrt(sum((Y(end, 1:3) - rvf(1:3)).^2))*lStar, ...
           'flyVms', sqrt(sum((Y(end, 4:6) - rvf(4:6)).^2))*lStar/tStar*1000, ...
           'minMoonAltKm', Cm.hMinKm, 'tMinMoonDays', Cm.tAtMin*tStar/86400, ...
           'moonRevs', (ang(end) - ang(1))/(2*pi), 'sD', sD, 'sA', sA, ...
           'thrustN', entry.thrustN, 'sheet', ks, 'sheet_key', sh.sheet_key, 'z8', z8, ...
           'rv0', rv0(:), 'rvf', rvf(:), 'orbits', struct('tD', tD, 'rvD', rvD, 'tA', tA, 'rvA', rvA));

%% Plot (rotating frame, Moon-centred, km):
hFig = [];
if ~doPlot, return, end
if visible, vis = 'on'; else, vis = 'off'; end
hFig = figure('Color', 'w', 'Visible', vis, 'Position', [80 80 960 720]);
ax = axes(hFig);  hold(ax, 'on');
dro  = (rvD(:, 1:3) - moonC)*lStar;
elfo = (rvA(:, 1:3) - moonC)*lStar;
r0   = (rv0(1:3) - moonC)*lStar;
rf   = (rvf(1:3) - moonC)*lStar;
P    = F.rMoonKm;
[xs, ys, zs] = sphere(48);
surf(ax, 1737.4*xs, 1737.4*ys, 1737.4*zs, 'FaceColor', [0.62 0.62 0.62], 'EdgeColor', 'none', ...
     'FaceLighting', 'gouraud', 'DisplayName', 'Moon (to scale)');
plot3(ax, dro(:, 1), dro(:, 2), dro(:, 3), 'k--', 'LineWidth', 0.9, 'DisplayName', 'DRO');
plot3(ax, elfo(:, 1), elfo(:, 2), elfo(:, 3), 'r--', 'LineWidth', 0.9, 'DisplayName', 'ELFO');
patch(ax, [P(:, 1); NaN], [P(:, 2); NaN], [P(:, 3); NaN], [F.tDays; NaN], 'FaceColor', 'none', ...
      'EdgeColor', 'interp', 'LineWidth', 1.6, 'DisplayName', 'transfer');
plot3(ax, r0(1), r0(2), r0(3), 'o', 'MarkerFaceColor', [0.1 0.7 0.2], 'MarkerEdgeColor', 'k', ...
      'MarkerSize', 8, 'DisplayName', 'departure');
plot3(ax, rf(1), rf(2), rf(3), 's', 'MarkerFaceColor', [0.85 0.1 0.1], 'MarkerEdgeColor', 'k', ...
      'MarkerSize', 8, 'DisplayName', 'arrival');
colormap(ax, parula);  clim(ax, [0 max(F.tDays(end), eps)]);
cb = colorbar(ax);  cb.Label.String = 'time since departure [d]';
light(ax, 'Position', [1 -1 1]);
axis(ax, 'equal');  grid(ax, 'on');  view(ax, -35, 25);
xlabel(ax, 'x from Moon [km]');  ylabel(ax, 'y from Moon [km]');  zlabel(ax, 'z from Moon [km]');
ttl = {sprintf('DRO %.4g ND -> ELFO raan %g deg, %.3g N, s_D %.4f, s_A %.4f (rotating frame)', ...
               sh.sheet_key.depSize, sh.sheet_key.raanDeg, entry.thrustN, sD, sA), ...
       sprintf('t_f %.3f d, \\DeltaV %.4f km/s, %.3f kg, min Moon alt %.0f km, miss %.1e km', ...
               F.tfDays, F.dvKms, F.propellantKg, F.minMoonAltKm, F.flyKm)};
extra = getf(opts, 'title', '');
if ~isempty(extra), ttl{end+1} = extra; end
title(ax, ttl, 'FontSize', 10);
legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal');
drawnow;
if ~isempty(outPng)
    od = fileparts(outPng);
    if ~isempty(od) && ~isfolder(od), mkdir(od); end
    print(hFig, outPng, '-dpng', '-r150');
end
end

% ------------------------------------------------------------------------
function v = getf(s, fld, d)
% GETF  Field with default.  INPUTS: s [struct]; fld [char]; d default.
% OUTPUTS: v the field value, or d when absent or empty.
if isfield(s, fld) && ~isempty(s.(fld)), v = s.(fld); else, v = d; end
end
