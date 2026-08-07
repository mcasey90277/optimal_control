function [eMin, eMax] = costate_catalog_extremes(cat_, filt, doPlot, metric)
%% Purpose:
%
%   Finds the EXTREME transfers in the multi-orbit costate CATALOG -- by
%   Delta-V (default) or by flight time -- and reports each extreme's full
%   coordinates: which sheet (DRO period, petal count), which phase pair
%   (nondimensional time AND days), which thrust rung, plus flight time,
%   Delta-V and propellant. Optionally flies both and plots them side by
%   side, each above its own sheet's orbits.
%
%   The catalog is COMPACT, so Delta-V is derived on the fly (exact for
%   these all-burn minimum-time solutions):
%       mf = 1 - Tmax_nd*tf_nd/c_nd,   dV = c_nd*ln(1/mf)*lStar/tStar.
%
%   RESTRICT THE SEARCH with the filter -- comparisons are most meaningful
%   with thrust held fixed (at one thrust, min time and min Delta-V are the
%   SAME transfer; they only separate across thrust levels), and often with
%   one sheet held fixed (isolating pure phasing effects).
%
%% Inputs:
%
%  cat_                     struct                  costate_catalog_dro_tulip
%
%  filt                     struct                  Optional restrictions,
%                                                   any subset of:
%                                                   .tauDRO  DRO period (ND)
%                                                   .Np      petal count
%                                                   .thrustN rung (N)
%                                                   [] or omitted = search
%                                                   the whole catalog
%
%  doPlot                   logical                 true: fly both extremes
%                                                   and draw side by side
%                                                   (default false)
%
%  metric                   char                    'deltaV' (default) or
%                                                   'time'
%
%% Outputs:
%
%  eMin                     struct                  The extreme entries:
%  eMax                                             .tauDRO .Np .thrustN
%                                                   .dep_frac/_nd/_days
%                                                   .arr_frac/_nd/_days
%                                                   .tf_nd/.tf_days
%                                                   .deltaV_kms
%                                                   .propellant_kg
%                                                   .z8 [8x1] (tfMin-ready)
%
%% Revision History:
%  M. Casey                                                   (c) 08/06/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: pure phasing effect -- one sheet, one thrust, with the plot:
       L = load('costate_catalog_dro_tulip.mat');
    cat_ = L.costate_catalog_dro_tulip;
[eMin,eMax] = costate_catalog_extremes(cat_, ...
                  struct('tauDRO',2.0,'Np',7,'thrustN',5), true);
     return;
end
if ~exist('filt','var') || isempty(filt), filt = struct(); end
if ~exist('doPlot','var'), doPlot = false; end
if ~exist('metric','var'), metric = 'deltaV'; end
assert(any(strncmpi(metric, {'deltaV','time'}, 1)), ...
       'metric must be ''deltaV'' or ''time'', got ''%s''', metric);

   tStar = cat_.constants.tStar_s;
   lStar = cat_.constants.lStar_km;
     cnd = cat_.thruster.c_nd;
    m0kg = cat_.thruster.m0_kg;
   d2day = tStar/86400;
     ndT = @(TN) (TN/m0kg)*tStar^2/(lStar*1000);

%% Walk every stored entry, apply the filter, rank by the metric:
    best = [];  worst = [];
   vBest = inf;  vWorst = -inf;
for ks = 1:numel(cat_.sheets)
    sh = cat_.sheets(ks);
    if isfield(filt,'tauDRO') && abs(sh.tauDRO - filt.tauDRO) > 1e-9, continue, end
    if isfield(filt,'Np')     && sh.Np ~= filt.Np,                    continue, end
    for kr = 1:numel(cat_.rungs_N)
        TN = cat_.rungs_N(kr);
        if isfield(filt,'thrustN') && abs(TN - filt.thrustN) > 1e-9, continue, end
        for iD = 1:numel(sh.sD_frac)
            for iA = 1:numel(sh.sA_frac)
                if ~sh.has_solution(iD,iA,kr), continue, end
                tf = sh.tf_nd(iD,iA,kr);
                mf = 1 - ndT(TN)*tf/cnd;
                if ~(mf > 0 && mf <= 1), continue, end   % corrupt entry
                dV = cnd*log(1/mf)*lStar/tStar;
                if strncmpi(metric,'t',1), v = tf; else, v = dV; end
                if v < vBest
                    vBest = v;  best  = pack(sh, iD, iA, TN, tf, dV, mf, m0kg, d2day, ...
                                             sh.z8(:,sh.entry_index(iD,iA,kr)));
                end
                if v > vWorst
                    vWorst = v;  worst = pack(sh, iD, iA, TN, tf, dV, mf, m0kg, d2day, ...
                                              sh.z8(:,sh.entry_index(iD,iA,kr)));
                end
            end
        end
    end
end
assert(~isempty(best), 'no catalog entries match the filter');
eMin = best;  eMax = worst;

%% Report:
if strncmpi(metric,'t',1), mName = 'flight time'; else, mName = 'Delta-V'; end
fprintf('\n%s extremes%s:\n', mName, filtstr(filt));
report('MINIMUM', eMin);
report('MAXIMUM', eMax);

if ~doPlot, return, end

%% Side-by-side plot, each extreme above its own sheet's orbits:
      mu = cat_.constants.muStar;
figure('Color','w','Position',[60 60 1200 560]);
ttl = {sprintf('MINIMUM %s', mName), sprintf('MAXIMUM %s', mName)};
for kp = 1:2
    if kp == 1, e = eMin; else, e = eMax; end
    [tD, rvD] = get_family_orbit(e.dep_family, e.dep_params);
    tauT = 2*pi*(e.Np-2)/(e.Np-1);
    [~, rvT0] = pumpkyn.cr3bp.getTulip(tauT, e.Np, -1);
    rvT0 = pumpkyn.cr3bp.cont_np(rvT0, tauT, mu, 1e-12);
    [tT, rvT] = pumpkyn.cr3bp.prop(tauT, rvT0, mu);
    rv0 = interp1(tD, rvD, e.dep_frac*tD(end), 'spline');
    rvf = interp1(tT, rvT, e.arr_frac*tT(end), 'spline');
    [~, y] = pumpkyn.cr3bp.tfMinProp(e.z8(8), [rv0(1:6)'; 1; e.z8(1:7)], ...
                 ndT(e.thrustN), cnd, mu);
    subplot(1,2,kp);
    plot3(y(:,1), y(:,2), y(:,3), 'b', 'LineWidth', 1.2); hold on;
    plot3(rvD(:,1), rvD(:,2), rvD(:,3), 'k--');
    plot3(rvT(:,1), rvT(:,2), rvT(:,3), 'r--');
    plot3(rv0(1), rv0(2), rv0(3), '.g', 'MarkerSize', 18);
    plot3(rvf(1), rvf(2), rvf(3), '.r', 'MarkerSize', 18);
    axis equal; grid on; view(-35, 25);
    xlabel('x [ND]'); ylabel('y [ND]'); zlabel('z [ND]');
    title({ttl{kp}, ...
        sprintf('%s \\tau=%.2f, N_p=%d, %.2f N', upper(e.dep_family), ...
                e.tauDRO, e.Np, e.thrustN), ...
        sprintf('\\DeltaV %.4f km/s,  t_f %.3f d,  %.2f kg prop', ...
                e.deltaV_kms, e.tf_days, e.propellant_kg), ...
        sprintf('depart %.3f d, arrive %.3f d', e.dep_days, e.arr_days)});
    if kp == 1
        legend('transfer','departure','tulip','depart','arrive','Location','best');
    end
end
rotate3d on
end

% ------------------------------------------------------------------------
function e = pack(sh, iD, iA, TN, tf, dV, mf, m0kg, d2day, z8)
% PACK  Assemble one extreme's report struct.
% INPUTS: sheet struct; indices; thrust; tf (ND); dV (km/s); final mass
%         fraction; m0 (kg); ND->days factor; z8.  OUTPUTS: e struct.
Pd = sh.tauDRO;  Pa = sh.period_tulip_nd;
% Departure reconstruction recipe: multi-family catalogs carry it per sheet;
% the original DRO catalog predates the fields, so fall back to its rule
if isfield(sh, 'dep_family')
    df = sh.dep_family;  dp = sh.dep_params;
else
    df = 'dro';  dp = struct('tau', sh.tauDRO);
end
e = struct('tauDRO',sh.tauDRO, 'dep_family',df, 'dep_params',dp, ...
    'Np',sh.Np, 'thrustN',TN, ...
    'dep_frac',sh.sD_frac(iD), 'dep_nd',sh.sD_frac(iD)*Pd, ...
    'dep_days',sh.sD_frac(iD)*Pd*d2day, ...
    'arr_frac',sh.sA_frac(iA), 'arr_nd',sh.sA_frac(iA)*Pa, ...
    'arr_days',sh.sA_frac(iA)*Pa*d2day, ...
    'tf_nd',tf, 'tf_days',tf*d2day, 'deltaV_kms',dV, ...
    'propellant_kg',(1-mf)*m0kg, 'z8',z8);
end

% ------------------------------------------------------------------------
function report(tag, e)
% REPORT  Print one extreme in both unit systems.
% INPUTS: tag [char]; e extreme struct.  OUTPUTS: none.
fprintf('  %s  dV = %.4f km/s,  t_f = %.5f ND = %.3f d   [%s tau=%.2f, Np=%d, %.2f N, %.2f kg]\n', ...
        tag, e.deltaV_kms, e.tf_nd, e.tf_days, upper(e.dep_family), ...
        e.tauDRO, e.Np, e.thrustN, ...
        e.propellant_kg);
fprintf('      departure phase: %.6f ND = %.4f days  (fraction %.4f)\n', ...
        e.dep_nd, e.dep_days, e.dep_frac);
fprintf('      arrival   phase: %.6f ND = %.4f days  (fraction %.4f)\n', ...
        e.arr_nd, e.arr_days, e.arr_frac);
end

% ------------------------------------------------------------------------
function s = filtstr(filt)
% FILTSTR  Human-readable filter description.  INPUTS: filt.  OUTPUTS: s.
s = '';
if isfield(filt,'tauDRO'), s = [s sprintf(' tau=%.2f', filt.tauDRO)]; end
if isfield(filt,'Np'),     s = [s sprintf(' Np=%d', filt.Np)]; end
if isfield(filt,'thrustN'),s = [s sprintf(' %.2g N', filt.thrustN)]; end
if isempty(s), s = ' over the WHOLE catalog';
else, s = [' (restricted to' s ')'];
end
end
