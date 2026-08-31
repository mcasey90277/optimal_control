function [eMin, eMax] = costate_catalog_extremes(cat_, filt, doPlot, metric)
%% Purpose:
%
%   Finds the EXTREME transfers in the multi-orbit costate CATALOG -- by
%   Delta-V (default) or by flight time -- and reports each extreme's full
%   coordinates: which sheet (departure key, petal count), which phase pair
%   (nondimensional time AND days), which thrust rung, plus flight time,
%   Delta-V and propellant. Optionally flies both and plots them side by
%   side, each above its own sheet's orbits.
%
%   GTO-CATALOG COPY (2026-08-31): this is costate_catalog_gto_tulip's OWN
%   copy of the shared extremes tool (DRO_tulip/indirect). Two GTO-specific
%   fixes, both required because sheets(k).tauDRO holds orientDeg (DEGREES)
%   for this catalog, not a period (the Stage B task-6 sheet-key override):
%     (1) the DEPARTURE PERIOD used for phase-day conversion is the sheet's
%         TRUE GTO Kepler period, reconstructed from dep_params.sma_km
%         (formula inlined verbatim from get_family_orbit.m's 'gto' case)
%         -- the shared copy's Pd = sh.tauDRO would report a departure
%         phase of 0 days at every orientation 0 and nonsense at the rest;
%     (2) reports and plot titles label the sheet as 'GTO orient=270 deg',
%         not 'GTO tau=270.00'. The filter accepts .orientDeg as a synonym
%         for .tauDRO.
%   Period-keyed catalogs fall through to the shared behavior unchanged.
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
%  cat_                     struct                  A compact costate catalog
%                                                   (any family; this copy's
%                                                   home is gto_tulip)
%
%  filt                     struct                  Optional restrictions,
%                                                   any subset of:
%                                                   .orientDeg / .tauDRO
%                                                     departure key (deg for
%                                                     this catalog)
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
%                                                   .depLab .arrLab
%                                                   .dep_frac/_nd/_days
%                                                   .arr_frac/_nd/_days
%                                                   .tf_nd/.tf_days
%                                                   .deltaV_kms
%                                                   .propellant_kg
%                                                   .z8 [8x1] (tfMin-ready)
%
%% Revision History:
%  M. Casey                                                   (c) 08/31/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: pure phasing effect -- one sheet, one thrust, with the plot.
   %Runs on whichever compact catalog sits here (or in results/):
     F = dir('costate_catalog_*.mat');
     if isempty(F), F = dir(fullfile('results', 'costate_catalog_*.mat')); end
     assert(~isempty(F), ...
            'demo: no costate_catalog_*.mat here or in results/');
       L = load(fullfile(F(1).folder, F(1).name));
      fn = fieldnames(L);
    cat_ = L.(fn{1});
   fprintf('demo catalog: %s\n', F(1).name);
      sh = cat_.sheets(ceil(numel(cat_.sheets)/2));
   if isnan(sh.Np)
       fl = struct('tauDRO',sh.tauDRO,'tauArr',sh.tau_arr,'thrustN',5);
   else
       fl = struct('tauDRO',sh.tauDRO,'Np',sh.Np,'thrustN',5);
   end
[eMin,eMax] = costate_catalog_extremes(cat_, fl, true);
     return;
end
if ~exist('filt','var') || isempty(filt), filt = struct(); end
if ~exist('doPlot','var'), doPlot = false; end
if ~exist('metric','var'), metric = 'deltaV'; end
assert(any(strncmpi(metric, {'deltaV','time'}, 1)), ...
       'metric must be ''deltaV'' or ''time'', got ''%s''', metric);
if isfield(filt, 'orientDeg'), filt.tauDRO = filt.orientDeg; end

   tStar = cat_.constants.tStar_s;
   lStar = cat_.constants.lStar_km;
     cnd = cat_.thruster.c_nd;
    m0kg = cat_.thruster.m0_kg;
   d2day = tStar/86400;
     ndT = @(TN) (TN/m0kg)*tStar^2/(lStar*1000);
     muE = 6.67384e-20*(1 - cat_.constants.muStar)*(5.9736E24 + 7.35E22);

%% Walk every stored entry, apply the filter, rank by the metric:
    best = [];  worst = [];
   vBest = inf;  vWorst = -inf;
for ks = 1:numel(cat_.sheets)
    sh = cat_.sheets(ks);
    if isfield(filt,'tauDRO') && abs(sh.tauDRO - filt.tauDRO) > 1e-9, continue, end
    if isfield(filt,'Np')
        assert(~isnan(sh.Np), 'costate_catalog_extremes:NpFilter', ...
            'this catalog''s arrival family is not tulip: filter by .tauArr, not .Np');
        if sh.Np ~= filt.Np, continue, end
    end
    if isfield(filt,'tauArr') && isfield(sh,'tau_arr') && ...
       abs(sh.tau_arr - filt.tauArr) > 1e-9,                          continue, end
    Pd = dep_period_nd(sh, muE, tStar);      % TRUE departure period (ND)
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
                    vBest = v;  best  = pack(sh, iD, iA, TN, tf, dV, mf, m0kg, ...
                                             d2day, Pd, ...
                                             sh.z8(:,sh.entry_index(iD,iA,kr)));
                end
                if v > vWorst
                    vWorst = v;  worst = pack(sh, iD, iA, TN, tf, dV, mf, m0kg, ...
                                              d2day, Pd, ...
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
    [tT, rvT] = get_family_orbit(e.arr_family, e.arr_params);
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
        sprintf('%s, %s, %.2f N', e.depLab, e.arrLab, e.thrustN), ...
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
function Pd = dep_period_nd(sh, muE, tStar)
% DEP_PERIOD_ND  The sheet's departure-orbit period in ND time. For the
% 'gto' pseudo-family the sheet key tauDRO is orientDeg (NOT a period), so
% the true Kepler period is rebuilt from dep_params.sma_km (formula from
% get_family_orbit.m's 'gto' case); every other family keys sheets on the
% period itself.
% INPUTS: sh sheet struct; muE km^3/s^2; tStar s.  OUTPUTS: Pd [ND].
if isfield(sh, 'dep_family') && strcmpi(sh.dep_family, 'gto')
    Pd = 2*pi*sqrt((sh.dep_params.sma_km)^3/muE)/tStar;
else
    Pd = sh.tauDRO;
end
end

% ------------------------------------------------------------------------
function e = pack(sh, iD, iA, TN, tf, dV, mf, m0kg, d2day, Pd, z8)
% PACK  Assemble one extreme's report struct.
% INPUTS: sheet struct; indices; thrust; tf (ND); dV (km/s); final mass
%         fraction; m0 (kg); ND->days factor; Pd departure period (ND,
%         family-correct); z8.  OUTPUTS: e struct.
Pa = sh.period_tulip_nd;
% Departure reconstruction recipe: multi-family catalogs carry it per sheet;
% the original DRO catalog predates the fields, so fall back to its rule
if isfield(sh, 'dep_family')
    df = sh.dep_family;  dp = sh.dep_params;
else
    df = 'dro';  dp = struct('tau', sh.tauDRO);
end
if strcmpi(df, 'gto')
    dLab = sprintf('GTO orient=%g deg', sh.tauDRO);
else
    dLab = sprintf('%s \\tau=%.2f', upper(df), sh.tauDRO);
end
% Arrival recipe: v2 catalogs carry it; v1 catalogs are tulip-arrival
if isfield(sh, 'arr_family')
    af = sh.arr_family;  ap = sh.arr_params;
else
    af = 'tulip';  ap = struct('Np', sh.Np, 'pm', -1);
end
if strcmpi(af, 'tulip')
    aLab = sprintf('Np=%d', sh.Np);
else
    aLab = sprintf('arr \\tau=%.2f', ap.tau);
end
e = struct('tauDRO',sh.tauDRO, 'dep_family',df, 'dep_params',dp, ...
    'arr_family',af, 'arr_params',ap, 'depLab',dLab, 'arrLab',aLab, ...
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
fprintf('  %s  dV = %.4f km/s,  t_f = %.5f ND = %.3f d   [%s, %s, %.2f N, %.2f kg]\n', ...
        tag, e.deltaV_kms, e.tf_nd, e.tf_days, ...
        strrep(e.depLab,'\tau','tau'), strrep(e.arrLab,'\tau','tau'), ...
        e.thrustN, e.propellant_kg);
fprintf('      departure phase: %.6f ND = %.4f days  (fraction %.4f)\n', ...
        e.dep_nd, e.dep_days, e.dep_frac);
fprintf('      arrival   phase: %.6f ND = %.4f days  (fraction %.4f)\n', ...
        e.arr_nd, e.arr_days, e.arr_frac);
end

% ------------------------------------------------------------------------
function s = filtstr(filt)
% FILTSTR  Human-readable filter description.  INPUTS: filt.  OUTPUTS: s.
s = '';
if isfield(filt,'tauDRO'), s = [s sprintf(' dep-key=%g', filt.tauDRO)]; end
if isfield(filt,'Np'),     s = [s sprintf(' Np=%d', filt.Np)]; end
if isfield(filt,'tauArr'), s = [s sprintf(' arr tau=%.2f', filt.tauArr)]; end
if isfield(filt,'thrustN'),s = [s sprintf(' %.2g N', filt.thrustN)]; end
if isempty(s), s = ' over the WHOLE catalog';
else, s = [' (restricted to' s ')'];
end
end
