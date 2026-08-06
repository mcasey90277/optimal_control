function [eMin, eMax] = costate_lib_extremes(lib, thrustN, doPlot, metric)
%% Purpose:
%
%   Finds the EXTREME transfers in a costate library -- either by Delta-V or
%   by flight time -- and reports the phasing of each in both nondimensional
%   time and days. Optionally flies both and plots them side by side, so the
%   geometric difference between a well-phased departure and a badly-phased
%   one is visible rather than merely numerical.
%
%   AT A FIXED THRUST THE TWO METRICS AGREE EXACTLY. These transfers burn
%   continuously, so the final mass is m_f = 1 - T*t_f/c and
%   dV = c*ln(1/m_f): a strictly increasing function of t_f when thrust is
%   held fixed. The cheapest phasing at a rung is therefore also the fastest,
%   and both metric options return the same entries. They diverge only in a
%   whole-library search across thrust levels, where the fastest transfer
%   (highest thrust) is nearly the most expensive one.
%
%   Restrict the search to one thrust rung by passing thrustN; omit it (or
%   pass []) to search every entry in the library. Restricting is usually
%   what you want, because Delta-V depends strongly on thrust -- a 15 N
%   transfer costs several times a 1 N transfer regardless of phasing -- so
%   an unrestricted search mostly reports the highest and lowest rungs.
%
%  ASSUMPTIONS / NOTES:
%
% • Delta-V comes from each entry's stored deltaV_kms, which is exact for
%   these all-burn minimum-time solutions: dV = c*ln(m0/mf).
% • Plotting needs pumpkyn/pumpkynPie on the path; the report alone does not.
%
%% Inputs:
%
%  lib                      struct                  costate_lib_dro_tulip_v2
%                                                   (or v1) loaded from .mat
%
%  thrustN                  double                  Restrict to this thrust
%                                                   rung (N). [] or omitted
%                                                   searches all entries
%
%  doPlot                   logical                 true: fly both and draw
%                                                   them side by side
%                                                   (default false)
%
%  metric                   char                    'deltaV' (default) ranks
%                                                   by Delta-V; 'time' ranks
%                                                   by flight time. Same
%                                                   answer at fixed thrust
%                                                   (see above)
%
%% Outputs:
%
%  eMin                     struct                  Library entry with the
%                                                   LOWEST value of the
%                                                   chosen metric
%
%  eMax                     struct                  Library entry with the
%                                                   HIGHEST value
%
%                                                   Each carries phases
%                                                   (frac/nd/days), thrust,
%                                                   lambda0, tf, deltaV_kms,
%                                                   propellant_kg, z8
%
%% Revision History:
%  M. Casey                                                   (c) 08/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: extremes at the 5 N rung, with the side-by-side plot:
       L = load('costate_lib_dro_tulip_v2.mat');
     lib = L.costate_lib_dro_tulip_v2;
[eMin,eMax] = costate_lib_extremes(lib, 5, true, 'deltaV');
     return;
end
if ~exist('thrustN','var'), thrustN = []; end
if ~exist('doPlot','var'),  doPlot  = false; end
if ~exist('metric','var'),  metric  = 'deltaV'; end

%% Select the population and the quantity being ranked:
      dv = [lib.entries.deltaV_kms];
      tN = [lib.entries.thrust_N];
if strncmpi(metric, 't', 1)
     rank = [lib.entries.tf_days];   mName = 'flight time';  mUnit = 'days';
else
     rank = dv;                      mName = 'Delta-V';      mUnit = 'km/s';
end
     sel = true(size(dv));
if ~isempty(thrustN)
     sel = abs(tN - thrustN) < 1e-9;
    if ~any(sel)
        error('costate_lib_extremes:noRung', ...
              'no entries at %.3f N; rungs are %s', thrustN, ...
              num2str(unique(tN)));
    end
end
     idx = find(sel);
[~, kLo] = min(rank(idx));
[~, kHi] = max(rank(idx));
    eMin = lib.entries(idx(kLo));
    eMax = lib.entries(idx(kHi));

%% Report:
if isempty(thrustN)
    fprintf('\n%s extremes over the whole library (%d entries):\n', ...
            mName, numel(idx));
else
    fprintf('\n%s extremes at %.2f N (%d entries):\n', mName, thrustN, numel(idx));
end
report('MINIMUM', eMin, lib);
report('MAXIMUM', eMax, lib);
fprintf('  ratio max/min = %.2f (%s, %s)\n\n', ...
        rank(idx(kHi))/rank(idx(kLo)), mName, mUnit);

if ~doPlot, return, end

%% Rebuild the orbits once, then fly both transfers:
      mu = lib.constants.muStar;
    tauD = lib.departure_params.tau;
[~,rvD0] = pumpkynPie.cr3bp.getDRO(tauD);
    rvD0 = pumpkyn.cr3bp.cont_np(rvD0, tauD, mu, 1e-12);
[tD,rvD] = pumpkyn.cr3bp.prop(tauD, rvD0, mu);
    tauT = lib.arrival_params.tau;
[~,rvT0] = pumpkyn.cr3bp.getTulip(tauT, lib.arrival_params.Np, lib.arrival_params.pm);
    rvT0 = pumpkyn.cr3bp.cont_np(rvT0, tauT, mu, 1e-12);
[tT,rvT] = pumpkyn.cr3bp.prop(tauT, rvT0, mu);

figure('Color','w','Position',[60 60 1200 560]);
if strncmpi(metric,'t',1)
    ttl = {'MINIMUM flight time', 'MAXIMUM flight time'};
else
    ttl = {'MINIMUM \DeltaV', 'MAXIMUM \DeltaV'};
end
for kp = 1:2
    if kp == 1, e = eMin; else, e = eMax; end
    rv0 = interp1(tD, rvD, e.departure_phase_frac*tD(end), 'spline');
    rvf = interp1(tT, rvT, e.arrival_phase_frac*tT(end),  'spline');
    Tnd = (e.thrust_N/lib.thruster.m0_kg)*lib.constants.tStar_s^2 / ...
          (lib.constants.lStar_km*1000);
    [~,y] = pumpkyn.cr3bp.tfMinProp(e.tf_nd, [rv0(1:6)'; 1; e.lambda0], ...
                Tnd, lib.thruster.c_nd, mu);
    subplot(1,2,kp);
    plot3(y(:,1), y(:,2), y(:,3), 'b', 'LineWidth', 1.2); hold on;
    plot3(rvD(:,1), rvD(:,2), rvD(:,3), 'k--');
    plot3(rvT(:,1), rvT(:,2), rvT(:,3), 'r--');
    plot3(rv0(1), rv0(2), rv0(3), '.g', 'MarkerSize', 18);
    plot3(rvf(1), rvf(2), rvf(3), '.r', 'MarkerSize', 18);
    axis equal; grid on; view(-35, 25);
    xlabel('x [ND]'); ylabel('y [ND]'); zlabel('z [ND]');
    title({sprintf('%s', ttl{kp}), ...
           sprintf('\\DeltaV = %.4f km/s,  t_f = %.3f d,  %.2f N,  %.2f kg prop', ...
                   e.deltaV_kms, e.tf_days, e.thrust_N, e.propellant_kg), ...
           sprintf('depart %.3f d,  arrive %.3f d', ...
                   e.departure_phase_days, e.arrival_phase_days)});
    if kp == 1
        legend('transfer','DRO','tulip','depart','arrive','Location','best');
    end
end
rotate3d on

end

% ------------------------------------------------------------------------
function report(tag, e, lib)
% REPORT  Print one extreme entry in both unit systems.
% INPUTS: tag [char]; e entry struct; lib library struct.  OUTPUTS: none.
fprintf(['  %s  dV = %.4f km/s  (t_f = %.5f ND = %.3f d, %.2f N, ', ...
         '%.2f kg propellant)\n'], tag, e.deltaV_kms, e.tf_nd, e.tf_days, ...
        e.thrust_N, e.propellant_kg);
fprintf('      departure phase: %.6f ND = %.4f days  (fraction %.4f of the %.3f-d DRO)\n', ...
        e.departure_phase_nd, e.departure_phase_days, ...
        e.departure_phase_frac, lib.departure_params.period_days);
fprintf('      arrival   phase: %.6f ND = %.4f days  (fraction %.4f of the %.3f-d tulip)\n', ...
        e.arrival_phase_nd, e.arrival_phase_days, ...
        e.arrival_phase_frac, lib.arrival_params.period_days);
end
