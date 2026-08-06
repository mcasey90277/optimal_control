%% costate_catalog_example
%
%   Worked example for costate_catalog_dro_tulip.mat: the MULTI-ORBIT
%   minimum-time costate catalog. Pick five coordinates -- DRO period, tulip
%   petal count, departure phase, arrival phase, thrust -- and the catalog
%   returns the flight time (interpolated between thrust rungs) plus a
%   converged z8 seed. Then rebuild both orbits from the catalog's recipes,
%   fly the seed, hand it to pumpkyn.cr3bp.tfMin, and plot.
%
%   The catalog is COMPACT (data-minimized): it stores only canonical
%   nondimensional quantities; days, Delta-V and masses are derived on the
%   fly using the formulas in cat.derive. This example shows each derivation
%   as it uses it.
%
%   HONESTY CONTRACT: whenever the returned solution differs from the
%   request (nearest sheet, snapped phases, nearest solved pair, different
%   rung), costate_catalog_pick PRINTS A WARNING saying exactly what you
%   are getting. Try changing the request below to an exact stored point
%   and watch the warnings disappear.
%
%   Files needed on the path:
%     costate_catalog_dro_tulip.mat   the catalog
%     costate_catalog_pick.m          the five-coordinate lookup
%     pumpkyn / pumpkynPie            orbit families, propagation, tfMin
%
%  M. Casey                                                     08/06/2026

%% 1. Load the catalog (EDIT catFile for your machine):
 catFile = 'costate_catalog_dro_tulip.mat';   %sits beside this script
       L = load(catFile);
    cat_ = L.costate_catalog_dro_tulip;
      mu = cat_.constants.muStar;
   tStar = cat_.constants.tStar_s;

fprintf('catalog: %d sheets, %d entries; DRO periods [%s] ND; petals [%s]\n', ...
        numel(cat_.sheets), cat_.n_entries, ...
        num2str(unique([cat_.sheets.tauDRO])), ...
        num2str(unique([cat_.sheets.Np])));

%% 2. Choose the five coordinates of a transfer:
  tauDRO = 2.0;       % DRO period (ND; 1 ND = 4.43 d) -- sheets at 0.5,1,2,3
      Np = 7;         % tulip petal count               -- sheets at 5,7,9,12
 depDays = 3.0;       % departure phase (days past the DRO reference)
 arrDays = 11.6;      % arrival phase (days past the tulip reference)
 thrustN = 5.0;       % thrust (N), anywhere in the rung range

[tf_nd, z8, info] = costate_catalog_pick(cat_, tauDRO, Np, depDays, arrDays, thrustN);

% derived quantities, straight from cat_.derive:
 tf_days = tf_nd * tStar/86400;
     Tnd = (info.delivered.thrustN/cat_.thruster.m0_kg) * tStar^2 / ...
           (cat_.constants.lStar_km*1000);
     cnd = cat_.thruster.c_nd;
      mf = 1 - Tnd*z8(8)/cnd;
  dV_kms = cnd*log(1/mf) * cat_.constants.lStar_km/tStar;
fprintf('flight time %.4f ND = %.3f days;  seed rung %.1f N;  dV %.4f km/s\n', ...
        tf_nd, tf_days, info.delivered.thrustN, dV_kms);

%% 3. Rebuild the DELIVERED sheet's orbits from the catalog's recipes:
      shD = info.delivered.tauDRO;
[~, rvD0] = pumpkynPie.cr3bp.getDRO(shD);
     rvD0 = pumpkyn.cr3bp.cont_np(rvD0, shD, mu, 1e-12);
 [tD,rvD] = pumpkyn.cr3bp.prop(shD, rvD0, mu);

     tauT = 2*pi*(info.delivered.Np-2)/(info.delivered.Np-1);
[~, rvT0] = pumpkyn.cr3bp.getTulip(tauT, info.delivered.Np, -1);
     rvT0 = pumpkyn.cr3bp.cont_np(rvT0, tauT, mu, 1e-12);
 [tT,rvT] = pumpkyn.cr3bp.prop(tauT, rvT0, mu);

%% 4. Endpoint states at the DELIVERED phases (info tells the truth):
      Pd = shD * tStar/86400;
      Pa = tT(end) * tStar/86400;
     rv0 = interp1(tD, rvD, mod(info.delivered.depDays/Pd,1)*tD(end), 'spline');
     rvf = interp1(tT, rvT, mod(info.delivered.arrDays/Pa,1)*tT(end), 'spline');

%% 5a. Fly the transfer straight from the catalog costates (no solver):
 [tau,y] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(1:6)'; 1; z8(1:7)], Tnd, cnd, mu);
  missKm = norm(y(end,1:3) - rvf(1:3)) * cat_.constants.lStar_km;
fprintf('flown arrival miss: %.3f km\n', missKm);

%% 5b. ...or hand it to tfMin (accepts the seed unchanged, ~1 s):
       z = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), z8, Tnd, cnd, mu);
fprintf('tfMin seed change |z - z8| = %.2e\n', norm(z - z8));

%% 6. Plot:
figure('Color','w');
plot3(y(:,1), y(:,2), y(:,3), 'b', 'LineWidth', 1.2); hold on;
plot3(rvD(:,1), rvD(:,2), rvD(:,3), 'k--');
plot3(rvT(:,1), rvT(:,2), rvT(:,3), 'r--');
plot3(rv0(1), rv0(2), rv0(3), '.g', 'MarkerSize', 18);
plot3(rvf(1), rvf(2), rvf(3), '.r', 'MarkerSize', 18);
axis equal; grid on;
xlabel('x [ND]'); ylabel('y [ND]'); zlabel('z [ND]');
title(sprintf('DRO(\\tau=%.1f) \\rightarrow tulip(N_p=%d), %.1f N, t_f = %.2f d', ...
      info.delivered.tauDRO, info.delivered.Np, info.delivered.thrustN, tf_days));
legend('transfer','DRO','tulip','depart','arrive');

%% 7. The period axis at a glance: best flight time vs DRO period at 5 N:
figure('Color','w'); hold on
kr5 = find(cat_.rungs_N == 5, 1);
for k = 1:numel(cat_.sheets)
    sh = cat_.sheets(k);
    if sh.Np ~= info.delivered.Np, continue, end
    tf5 = sh.tf_nd(:,:,kr5);
    tf5 = tf5(sh.has_solution(:,:,kr5));
    if isempty(tf5), continue, end
    plot(sh.tauDRO*tStar/86400, min(tf5)*tStar/86400, 'ob', 'MarkerFaceColor','b');
    plot(sh.tauDRO*tStar/86400, max(tf5)*tStar/86400, 'sr', 'MarkerFaceColor','r');
end
grid on
xlabel('DRO period [days]'); ylabel('t_f at 5 N [days]');
title(sprintf('best (blue) and worst (red) phasing vs DRO period, N_p = %d', ...
      info.delivered.Np));
