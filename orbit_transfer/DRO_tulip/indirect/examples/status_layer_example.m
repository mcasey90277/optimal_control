%% status_layer_example
%
%   Worked example for the OPTIMALITY-STATUS layer of
%   costate_catalog_dro_tulip_70mN.mat. Every transfer the library records
%   carries a status (4 sufficient, 3 necessary only, 2 conjugate point
%   found, 1 neither) and a one-line reason; a phase cell can hold its
%   primary transfer plus ALTERNATIVES. This script takes the cell with the
%   most status-2 alternatives, ranks everything recorded there by flight
%   time, flies the primary and its fastest status-2 alternative from their
%   multiple-shooting junctions, shows what status 2 means, confirms that
%   pumpkyn tfMin still accepts the alternative (a conjugate-point transfer
%   is a converged extremal, just not a locally optimal one), and draws the
%   status map of the whole phase torus.
%
%   An orbit transfer does not have to be optimal to be useful: the status
%   says what was ESTABLISHED about each one, and a mission designer can
%   take a faster status-2 alternative knowing exactly what that means.
%
%   FILES: the catalog is loaded BY BARE FILENAME, so run this from the
%   folder that holds it; every function called is in that same folder (or
%   is pumpkyn / pumpkynPie). No absolute paths, no addpath.
%
%   Run time: about 10 seconds.
%
%  M. Casey                                                     10/13/2026

if ~usejava('desktop'), figVis = 'off'; else, figVis = 'on'; end

%% 1. Load the catalog:
 catFile = 'costate_catalog_dro_tulip_70mN.mat';   % sits beside this script
       L = load(catFile);
    cat_ = L.costate_catalog_dro_tulip_70mN;
      sh = cat_.sheets(1);
       D = costate_lib_describe(cat_);
fprintf('[1] verdict: %d primaries + %d alternatives recorded (status counts primaries [%s], alternatives [%s] for codes 4 3 2 1)\n', ...
        D.nEntries, D.nAlternatives, num2str(D.statusCounts(:,1).'), num2str(D.statusCounts(:,2).'));

%% 2. The sample cell: the one that holds the most status-2 alternatives:
[prim, alt2] = sample_status_pair(cat_);
 sD = prim.sD;                      % departure phase [fraction of the DRO period]
 sA = prim.sA;                      % arrival phase   [fraction of the tulip period]
fprintf('[2] verdict: cell sD = %.4f (%.2f d), sA = %.4f (%.2f d); primary %.2f d, fastest status-2 alternative %.2f d\n', ...
        sD, sD*sh.tauDRO*cat_.constants.tStar_s/86400, sA, ...
        sA*sh.period_tulip_nd*cat_.constants.tStar_s/86400, prim.tfDays, alt2.tfDays);

%% 3. Everything recorded at the cell, ranked by flight time:
 T = transfers_at(cat_, sD, sA);
fprintf('  %-12s %-22s %8s %10s   %s\n', 'kind', 'status', 't_f [d]', 'dV [km/s]', 'reason (first 60 chars)');
for k = 1:numel(T)
    fprintf('  %-12s %d %-20s %8.2f %10.4f   %s\n', T(k).kind, T(k).status, T(k).statusName, ...
            T(k).tfDays, T(k).dvKms, T(k).reason(1:min(60, numel(T(k).reason))));
end
fprintf('[3] verdict: %d transfers at this cell, t_f from %.2f d to %.2f d\n', ...
        numel(T), T(1).tfDays, T(end).tfDays);

%% 4. Fly the primary and the status-2 alternative from their junctions:
 B = catalog_endpoints(cat_, 1, cat_.rungs_N(1));
Fp = fly_catalog_entry(cat_, prim, struct('B', B));
Fa = fly_catalog_entry(cat_, alt2, struct('B', B));
fprintf(' primary (status %d, %s):\n', Fp.status, Fp.statusName);
print_transfer_summary(Fp);
fprintf(' alternative (status %d, %s):\n', Fa.status, Fa.statusName);
print_transfer_summary(Fa);
fprintf('[4] verdict: primary arrives within %.3f km / %.4f m/s (junction defect %.1e km); alternative within %.3f km / %.4f m/s (defect %.1e km)\n', ...
        Fp.flyKm, Fp.flyVms, max([Fp.defectKm 0]), Fa.flyKm, Fa.flyVms, max([Fa.defectKm 0]));

%% 5. The two thrust programs: |u| components and mass against days:
 figure('Color', 'w', 'Visible', figVis, 'Position', [100 100 900 600]);
subplot(2, 1, 1); hold on
plot(Fp.tDays, Fp.u(:, 1), 'b-',  Fp.tDays, Fp.u(:, 2), 'b--',  Fp.tDays, Fp.u(:, 3), 'b:');
plot(Fa.tDays, Fa.u(:, 1), 'r-',  Fa.tDays, Fa.u(:, 2), 'r--',  Fa.tDays, Fa.u(:, 3), 'r:');
grid on; xlabel('time [days]'); ylabel('thrust direction u_x, u_y, u_z [-]');
title('primary (blue) and status-2 alternative (red): unit thrust direction components');
subplot(2, 1, 2); hold on
plot(Fp.tDays, Fp.massKg, 'b-', Fa.tDays, Fa.massKg, 'r-');
grid on; xlabel('time [days]'); ylabel('mass [kg]');
legend(sprintf('primary, %.2f d', Fp.tfDays), sprintf('alternative, %.2f d', Fa.tfDays), 'Location', 'southwest');
fprintf('[5] verdict: primary burns %.3f kg in %.2f d, alternative %.3f kg in %.2f d\n', ...
        Fp.propellantKg, Fp.tfDays, Fa.propellantKg, Fa.tfDays);

%% 6. What status 2 means, and why this alternative has it:
  K = status_key();
 k2 = find(K.codes == alt2.status, 1);
fprintf(' status %d = %s: %s\n', K.codes(k2), K.names{k2}, K.meaning{k2});
fprintf(' this alternative''s reason, in full: %s\n', alt2.reason);
fprintf('[6] verdict: %.2f d %s than the primary, and NOT locally optimal (a conjugate point lies on the arc)\n', ...
        abs(alt2.tfDays - prim.tfDays), ternary(alt2.tfDays < prim.tfDays, 'faster', 'slower'));

%% 7. tfMin accepts the alternative unchanged (it is a converged extremal):
[acc, dz, ~, tinfo] = tfmin_accepts(cat_, alt2, struct('B', B));
fprintf('[7] verdict: |dz| = %.2e (tolerance %.0e) -> %s (tfMin took %.1f s)\n', ...
        dz, tinfo.tolDz, ternary(acc, 'ACCEPTED unchanged', 'NOT accepted'), tinfo.wallSec);

%% 8. The status map of the whole phase torus:
 H = plot_status_torus(cat_, struct('visible', strcmp(figVis, 'on')));
fprintf('[8] verdict: torus drawn, %d cells with a status-2 transfer, %d off-grid alternatives (%d status 2)\n', ...
        sum(H.conjCell(:)), H.nOffGrid, H.nOffGridConj);

fprintf('ALL SECTIONS DONE\n');
if ~usejava('desktop'), close all; end

function s = ternary(c, a, b)
% TERNARY  a if c is true, b otherwise (local helper for the verdict line).
if c, s = a; else, s = b; end
end
