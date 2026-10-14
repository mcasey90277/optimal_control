function [F, hFig] = catalog_transfer(cat_, entry, opts)
%% Purpose:
%
%   ONE library transfer, flown and handed back: the time history of the
%   state, thrust direction and mass of a chosen entry, in days, km and kg
%   as well as ND, with an optional 3D plot in the rotating Earth-Moon
%   frame. The flight uses pumpkyn's own propagator (pumpkyn.cr3bp.tfMinProp)
%   segment by segment from the entry's stored multiple-shooting junction
%   states (the first junction IS the departure state), or, when it has
%   none, in one flight from z8 started at the departure phase of the
%   departure orbit. Nothing is solved or re-checked here.
%
%  ASSUMPTIONS / NOTES:
%
% • Junctions sit on the uniform grid linspace(0, t_f, K+1): column k is the
%   augmented state [r; v; m; lam_r; lam_v; lam_m] at the start of segment k.
% • Minimum-time, all-burn: thrust direction u = -lam_v/|lam_v|, throttle 1.
% • Orbits come from the sheet's own recipe (dep_family / arr_family) at the
%   catalog's mass ratio; a sheet without those keys uses the DRO / tulip
%   rule (deliverable-3 catalogs). The arrival miss F.flyKm is measured
%   against a spline of the arrival orbit at the entry's arrival phase.
% • Phases are FRACTIONS of a period, wrapping.
%
%% Inputs:
%
%  cat_                     struct                  the catalog
%
%  entry                    struct                  a transfers_at row (z8,
%                                                   sD, sA, junctions,
%                                                   status, ...); OR a
%                                                   costate_catalog_extremes
%                                                   row (dep_frac / arr_frac);
%                                                   OR struct('iD', , 'iA', )
%                                                   = the primary of that
%                                                   grid cell; OR
%                                                   struct('z8', , 'sD', ,
%                                                   'sA', ) with optional
%                                                   .junctions
%
%  opts                     struct (optional)       .method ['auto'] |
%                                                   'junctions' | 'z8';
%                                                   .sheet [1]; .thrustN
%                                                   [rungs_N(1)]; .plot
%                                                   [false]; .visible [true]
%                                                   figure visibility;
%                                                   .title [''] extra title
%                                                   line
%
%% Outputs:
%
%  F                        struct                  .t [n x 1] ND, .tDays,
%                                                   .Y [n x 14] ND, .rKm
%                                                   [n x 3], .vKms [n x 3],
%                                                   .u [n x 3] unit thrust
%                                                   direction, .massKg,
%                                                   .tfNd .tfDays .dvKms
%                                                   .propellantKg
%                                                   .finalMassKg .flyKm
%                                                   .method .status
%                                                   .statusName .sD .sA .z8
%                                                   .rv0 .rvf (ND) .orbits
%                                                   (.tD .rvD .tT .rvT)
%
%  hFig                     handle                  the figure when
%                                                   opts.plot, else []
%
%% Revision History:
%  M. Casey                                                   (c) 10/14/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the fastest primary of whichever catalog sits in the current
   %folder, flown from its junctions and plotted:
     D = dir('costate_catalog_*.mat');
     assert(~isempty(D), 'demo: no costate_catalog_*.mat in the current folder');
       L = load(D(1).name);  fn = fieldnames(L);  cat_ = L.(fn{1});
   fprintf('demo catalog: %s\n', D(1).name);
      sh = cat_.sheets(1);
     tfg = sh.tf_nd(:, :, 1);  tfg(~sh.has_solution(:, :, 1)) = inf;
   [~, kBest] = min(tfg(:));  [iD, iA] = ind2sub(size(tfg), kBest);
       F = catalog_transfer(cat_, struct('iD', iD, 'iA', iA), struct('plot', true));
   fprintf('fastest primary: cell (%d, %d), t_f %.3f d, Delta-V %.4f km/s, arrival miss %.3f km (%s flight)\n', ...
           iD, iA, F.tfDays, F.dvKms, F.flyKm, F.method);
     return
end
if nargin < 3 || isempty(opts), opts = struct(); end
method  = getf(opts, 'method', 'auto');
ks      = getf(opts, 'sheet', 1);
thrustN = getf(opts, 'thrustN', cat_.rungs_N(1));
doPlot  = getf(opts, 'plot', false);
visible = getf(opts, 'visible', true);
assert(any(strcmp(method, {'auto', 'junctions', 'z8'})), 'catalog_transfer:method', ...
       'opts.method must be auto, junctions or z8');
sh = cat_.sheets(ks);

%% Resolve the entry form:
if isfield(entry, 'iD') && isfield(entry, 'iA') && ~isfield(entry, 'z8')
    T  = transfers_at(cat_, sh.sD_frac(entry.iD), sh.sA_frac(entry.iA));
    ip = find(strcmp({T.kind}, 'primary'), 1);
    if isempty(ip), error('catalog_transfer:noPrimary', 'cell (%d, %d) has no primary', entry.iD, entry.iA); end
    entry = T(ip);
end
if isfield(entry, 'dep_frac') && ~isfield(entry, 'sD'), entry.sD = entry.dep_frac; end
if isfield(entry, 'arr_frac') && ~isfield(entry, 'sA'), entry.sA = entry.arr_frac; end
if ~isfield(entry, 'z8') || numel(entry.z8) ~= 8
    error('catalog_transfer:entry', 'the entry needs a z8 of 8 elements');
end
if ~isfield(entry, 'sD') || ~isfield(entry, 'sA')
    error('catalog_transfer:entry', 'the entry needs phase fractions sD and sA');
end
z8 = entry.z8(:);
J  = [];  if isfield(entry, 'junctions'), J = entry.junctions; end
hasJ = ~isempty(J) && size(J, 1) == 14 && all(isfinite(J(:)));
if strcmp(method, 'junctions') && ~hasJ
    error('catalog_transfer:noJunctions', 'the entry carries no 14 x K junction states');
end
useJ = hasJ && ~strcmp(method, 'z8');

%% Orbits (the sheet's own recipe at the catalog's mass ratio) and physics:
mu = cat_.constants.muStar;
if isfield(sh, 'dep_family'), df = sh.dep_family;  dp = sh.dep_params;
else,                         df = 'dro';          dp = struct('tau', sh.tauDRO); end
if isfield(sh, 'arr_family'), af = sh.arr_family;  ap = sh.arr_params;
else,                         af = 'tulip';        ap = struct('Np', sh.Np, 'pm', -1); end
dp.muStar = mu;  ap.muStar = mu;
[tD, rvD] = get_family_orbit(df, dp);
[tT, rvT] = get_family_orbit(af, ap);
tStar = cat_.constants.tStar_s;  lStar = cat_.constants.lStar_km;
m0kg  = cat_.thruster.m0_kg;     cnd   = cat_.thruster.c_nd;
Tnd   = (thrustN/m0kg)*tStar^2/(lStar*1000);
sD = mod(entry.sD, 1);  sA = mod(entry.sA, 1);
rv0 = interp1(tD/tD(end), rvD, sD, 'spline');  rv0 = rv0(1:6);
rvf = interp1(tT/tT(end), rvT, sA, 'spline');  rvf = rvf(1:6);

%% Fly:
tf = z8(8);
if useJ
    K = size(J, 2);  h = tf/K;  t = [];  Y = [];
    for k = 1:K
        [tk, Yk] = pumpkyn.cr3bp.tfMinProp(h, J(:, k), Tnd, cnd, mu);
        if k > 1, tk = tk(2:end);  Yk = Yk(2:end, :); end
        t = [t; (k - 1)*h + tk];  Y = [Y; Yk];
    end
    method = 'junctions';
else
    [t, Y] = pumpkyn.cr3bp.tfMinProp(tf, [rv0(:); 1; z8(1:7)], Tnd, cnd, mu);
    method = 'z8';
end
lv = Y(:, 11:13);
mf = Y(end, 7);

%% Assemble:
F = struct('t', t, 'tDays', t*tStar/86400, 'Y', Y, 'rKm', Y(:, 1:3)*lStar, ...
           'vKms', Y(:, 4:6)*lStar/tStar, 'u', -lv ./ sqrt(sum(lv.^2, 2)), 'massKg', Y(:, 7)*m0kg, ...
           'tfNd', tf, 'tfDays', tf*tStar/86400, 'dvKms', cnd*log(1/mf)*lStar/tStar, ...
           'propellantKg', (1 - mf)*m0kg, 'finalMassKg', mf*m0kg, ...
           'flyKm', sqrt(sum((Y(end, 1:3) - rvf(1:3)).^2))*lStar, ...
           'method', method, 'status', NaN, 'statusName', 'not annotated', ...
           'sD', sD, 'sA', sA, 'z8', z8, 'rv0', rv0(:), 'rvf', rvf(:), ...
           'orbits', struct('tD', tD, 'rvD', rvD, 'tT', tT, 'rvT', rvT));
if isfield(entry, 'status') && ~isempty(entry.status) && isfinite(entry.status)
    F.status = double(entry.status);
    Kst = status_key();  kk = find(Kst.codes == F.status, 1);
    if ~isempty(kk), F.statusName = Kst.names{kk}; end
end
if isfield(entry, 'statusName') && ~isempty(entry.statusName), F.statusName = entry.statusName; end

%% Plot:
hFig = [];
if ~doPlot, return, end
if visible, vis = 'on'; else, vis = 'off'; end
hFig = figure('Color', 'w', 'Visible', vis, 'Position', [80 80 900 700]);
plot3(rvD(:, 1), rvD(:, 2), rvD(:, 3), 'k--', 'LineWidth', 0.8); hold on
plot3(rvT(:, 1), rvT(:, 2), rvT(:, 3), 'r--', 'LineWidth', 0.8);
plot3(Y(:, 1), Y(:, 2), Y(:, 3), 'b-', 'LineWidth', 1.4);
plot3(rv0(1), rv0(2), rv0(3), 'g.', 'MarkerSize', 20);
plot3(rvf(1), rvf(2), rvf(3), 'r.', 'MarkerSize', 20);
plot3(1 - mu, 0, 0, 'ko', 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerSize', 7);
plot3(-mu, 0, 0, 'bo', 'MarkerFaceColor', [0.3 0.5 1.0], 'MarkerSize', 9);
axis equal; grid on; view(-35, 25);
xlabel('x [ND]'); ylabel('y [ND]'); zlabel('z [ND]');
ttl = {sprintf('%s -> %s transfer, sD = %.4f, sA = %.4f, %.3g N', upper(df), af, sD, sA, thrustN), ...
       sprintf('t_f %.3f d, \\DeltaV %.4f km/s, %.3f kg propellant  [%s, %s flight]', ...
               F.tfDays, F.dvKms, F.propellantKg, F.statusName, method)};
extra = getf(opts, 'title', '');
if ~isempty(extra), ttl{end+1} = extra; end
title(ttl);
legend('departure orbit', 'arrival orbit', 'transfer', 'depart', 'arrive', 'Moon', 'Earth', 'Location', 'best');
rotate3d on
end

% ------------------------------------------------------------------------
function v = getf(s, fld, d)
% GETF  Field with default.  INPUTS: s; fld; d.  OUTPUTS: v.
if isfield(s, fld) && ~isempty(s.(fld)), v = s.(fld); else, v = d; end
end
