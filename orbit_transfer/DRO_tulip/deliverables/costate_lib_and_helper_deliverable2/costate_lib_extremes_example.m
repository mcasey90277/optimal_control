%% costate_lib_extremes_example
%
%   Worked example for costate_lib_extremes: find the cheapest and the most
%   expensive transfers a costate library holds, in Delta-V, and see what
%   distinguishes them.
%
%   Four questions are asked in turn:
%     1. At ONE thrust level, which phasing is cheapest and which is dearest?
%        (This is the useful comparison: phasing alone, thrust held fixed.)
%     2. Across the WHOLE library, what is the cheapest transfer available
%        at any thrust?
%     3. How does the spread between best and worst phasing vary with thrust?
%     4. The same search ranked by FLIGHT TIME instead of Delta-V.
%
%   A caution worth internalizing: at a FIXED thrust the two rankings are
%   identical. These transfers burn continuously, so dV = c*ln(1/(1-T*t_f/c))
%   rises monotonically with t_f -- the cheapest phasing at a rung is also
%   the fastest. The metrics separate only across thrust levels, where the
%   quickest transfer of all is nearly the most expensive.
%
%   Files needed on the path:
%     costate_lib_dro_tulip_v2.mat   the library
%     costate_lib_extremes.m         the search + side-by-side plot
%     pumpkyn / pumpkynPie           orbit families and propagation
%
%   Delta-V is exact for these solutions rather than estimated: minimum-time
%   transfers burn continuously, so dV = c*ln(m0/mf) with mf fixed by thrust
%   and flight time.
%
%  M. Casey                                                     08/05/2026

%% 1. Load the library (EDIT libFile for your machine):
 libFile = 'costate_lib_dro_tulip_v2.mat';   %sits beside this script
       L = load(libFile);
     lib = L.costate_lib_dro_tulip_v2;

%% 2. Extremes at ONE thrust level, with the side-by-side trajectory plot:
 thrustN = 5.0;                       % any rung in lib.thruster.thrust_rungs_N
[eMin, eMax] = costate_lib_extremes(lib, thrustN, true);

%% 3. Extremes across the entire library (all thrusts):
%   Expect this to pick out the lowest and highest RUNGS rather than a
%   phasing effect -- thrust dominates Delta-V.
[gMin, gMax] = costate_lib_extremes(lib, [], false);
fprintf('cheapest anywhere: %.4f km/s at %.2f N;  dearest: %.4f km/s at %.2f N\n\n', ...
        gMin.deltaV_kms, gMin.thrust_N, gMax.deltaV_kms, gMax.thrust_N);

%% 4. Ranked by FLIGHT TIME instead -- at one rung (identical answer, as
%   explained above) and then across the whole library (a different answer):
[tMinRung, tMaxRung] = costate_lib_extremes(lib, thrustN, false, 'time');
fprintf('at %.1f N, fastest phasing matches the cheapest: %d\n\n', ...
        thrustN, isequal(tMinRung.z8, eMin.z8));

[tMin, tMax] = costate_lib_extremes(lib, [], true, 'time');
fprintf(['fastest transfer anywhere: %.3f d at %.2f N (costs %.4f km/s)\n', ...
         'slowest transfer anywhere: %.3f d at %.2f N (costs %.4f km/s)\n\n'], ...
        tMin.tf_days, tMin.thrust_N, tMin.deltaV_kms, ...
        tMax.tf_days, tMax.thrust_N, tMax.deltaV_kms);

%% 5. How much does phasing alone matter, rung by rung?
    rungs = lib.thruster.thrust_rungs_N;
   spread = nan(size(rungs));
    dvLo  = nan(size(rungs));
    dvHi  = nan(size(rungs));
fprintf('thrust    dV_min     dV_max    spread   (km/s)\n');
for k = 1:numel(rungs)
    try
        [a, b] = evalc_quiet(lib, rungs(k));
        dvLo(k) = a;  dvHi(k) = b;  spread(k) = b - a;
        fprintf('%6.2f N  %8.4f  %8.4f  %8.4f\n', rungs(k), a, b, b-a);
    catch
        fprintf('%6.2f N  (no entries)\n', rungs(k));
    end
end

figure('Color','w');
plot(rungs, dvLo, '-o', 'LineWidth',1.3); hold on
plot(rungs, dvHi, '-s', 'LineWidth',1.3);
fill([rungs, fliplr(rungs)], [dvLo, fliplr(dvHi)], [0.2 0.4 0.8], ...
     'FaceAlpha',0.12, 'EdgeColor','none');
grid on; set(gca,'XScale','log');
xlabel('thrust [N]  (log scale)'); ylabel('\DeltaV [km/s]');
title({'DRO \rightarrow tulip minimum-time transfers', ...
       'best and worst phasing at each thrust; band = the cost of phasing'});
legend('best phasing','worst phasing','spread','Location','northwest');

%% ------------------------------------------------------------------------
function [dvLo, dvHi] = evalc_quiet(lib, thrustN)
% EVALC_QUIET  Delta-V extremes at one rung, without the printed report.
% INPUTS: lib structure; thrustN rung (N).  OUTPUTS: dvLo, dvHi [km/s].
      dv = [lib.entries.deltaV_kms];
      tN = [lib.entries.thrust_N];
     sel = abs(tN - thrustN) < 1e-9;
if ~any(sel), error('no entries at this rung'); end
    dvLo = min(dv(sel));
    dvHi = max(dv(sel));
end
