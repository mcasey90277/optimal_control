%% costate_catalog_dpo_example
%
%   Worked example for costate_catalog_dpo_tulip.mat: the DPO -> TULIP
%   minimum-time costate catalog. Pick five coordinates -- DPO period,
%   tulip petal count, departure phase, arrival phase, thrust -- and the
%   catalog returns the flight time (interpolated between thrust rungs)
%   plus a converged z8 seed. Then rebuild both orbits from the catalog's
%   recipes, fly the seed, hand it to pumpkyn.cr3bp.tfMin, and plot.
%
%   The catalog is COMPACT (data-minimized): it stores only canonical
%   nondimensional quantities; days, Delta-V and masses are derived on the
%   fly using the formulas in cat.derive. This example shows each
%   derivation as it uses it.
%
%   SCHEMA NOTE: sheets(k).tauDRO is a LEGACY field name shared with the
%   DRO catalog so the same pickers work on both -- here it holds the DPO
%   period (ND). The full departure reconstruction recipe rides in
%   sheets(k).dep_family / dep_params (the DPO recipe is its period).
%
%   HONESTY CONTRACT: whenever the returned solution differs from the
%   request (nearest sheet, snapped phases, nearest solved pair, different
%   rung), costate_catalog_pick PRINTS A WARNING saying exactly what you
%   are getting. Try changing the request below to an exact stored point
%   and watch the warnings disappear.
%
%   Files needed on the path:
%     costate_catalog_dpo_tulip.mat  the catalog
%     costate_catalog_pick.m          the five-coordinate lookup
%     get_family_orbit.m              family name + params -> orbit
%     pumpkyn / pumpkynPie            orbit families, propagation, tfMin
%
%  M. Casey                                                     08/07/2026

%% 1. Load the catalog (sits beside this script):
 catFile = 'costate_catalog_dpo_tulip.mat';
       L = load(catFile);
    cat_ = L.costate_catalog_dpo_tulip;
      mu = cat_.constants.muStar;
   tStar = cat_.constants.tStar_s;

fprintf('catalog: %d sheets, %d entries; DPO periods [%s] ND; petals [%s]\n', ...
        numel(cat_.sheets), cat_.n_entries, ...
        num2str(unique([cat_.sheets.tauDRO])), ...
        num2str(unique([cat_.sheets.Np])));

%% 2. Choose the five coordinates of a transfer:
  tauDPO = 3.0;       % DPO period (ND; 1 ND = 4.43 d) -- sheets at 1,2,3,4
      Np = 9;         % tulip petal count               -- sheets at 5,7,9,12
 depDays = 4.0;       % departure phase (days past the DPO reference)
 arrDays = 10.0;      % arrival phase (days past the tulip reference)
 thrustN = 5.0;       % thrust (N), anywhere in the rung range

[tf_nd, z8, info] = costate_catalog_pick(cat_, tauDPO, Np, depDays, arrDays, thrustN);

% derived quantities, straight from cat_.derive. NOTE the honest split:
% tf_nd is interpolated at tfThrustN; the z8 seed is a STORED converged
% vector at seedThrustN (z8 is never interpolated) -- fly it at that rung.
 tf_days = tf_nd * tStar/86400;
   seedN = info.delivered.seedThrustN;
     Tnd = (seedN/cat_.thruster.m0_kg) * tStar^2 / ...
           (cat_.constants.lStar_km*1000);
     cnd = cat_.thruster.c_nd;
      mf = 1 - Tnd*z8(8)/cnd;
  dV_kms = cnd*log(1/mf) * cat_.constants.lStar_km/tStar;
fprintf(['flight time %.4f ND = %.3f days (at %.2f N);  seed rung %.1f N;  ', ...
         'seed dV %.4f km/s\n'], tf_nd, tf_days, info.delivered.tfThrustN, ...
        seedN, dV_kms);

%% 3. Rebuild the DELIVERED sheet's orbits from the catalog's recipes:
% the sheet index the picker chose:
      ks = find(abs([cat_.sheets.tauDRO] - info.delivered.tauDRO) < 1e-9 & ...
                [cat_.sheets.Np] == info.delivered.Np, 1);
      sh = cat_.sheets(ks);
[tD,rvD] = get_family_orbit(sh.dep_family, sh.dep_params);   % the DPO

     tauT = 2*pi*(sh.Np-2)/(sh.Np-1);
[~, rvT0] = pumpkyn.cr3bp.getTulip(tauT, sh.Np, sh.pm);
     rvT0 = pumpkyn.cr3bp.cont_np(rvT0, tauT, mu, 1e-12);
 [tT,rvT] = pumpkyn.cr3bp.prop(tauT, rvT0, mu);

%% 4. Endpoint states at the DELIVERED phases (info tells the truth):
      Pd = tD(end) * tStar/86400;
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
title(sprintf('dpo(\\tau=%.2f) \\rightarrow tulip(N_p=%d), %.1f N, t_f = %.2f d', ...
      info.delivered.tauDRO, info.delivered.Np, seedN, z8(8)*tStar/86400));
legend('transfer','DPO','tulip','depart','arrive');

%% 7. The period axis at a glance: best flight time vs DPO period at 5 N:
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
xlabel('DPO period [days]'); ylabel('t_f at 5 N [days]');
title(sprintf('best (blue) and worst (red) phasing vs DPO period, N_p = %d', ...
      info.delivered.Np));
