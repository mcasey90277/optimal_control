%% costate_catalog_gto_example
%
%   Worked example for costate_catalog_gto_tulip.mat: the GTO -> tulip
%   minimum-time costate catalog. Pick five coordinates -- GTO departure
%   ORIENTATION (degrees), tulip petal count, departure phase, arrival
%   phase, thrust -- and the catalog returns the flight time (interpolated
%   between thrust rungs) plus a converged z8 seed. Then rebuild both
%   orbits from the catalog's recipes, fly the seed, hand it to
%   pumpkyn.cr3bp.tfMin, and plot.
%
%   TWO THINGS ARE DIFFERENT FROM THE PERIOD-KEYED CATALOGS (DRO/HALO/DPO):
%     - The departure key is an ORIENTATION ANGLE, not a period: orientDeg
%       is the angle from the Earth->Moon line to the GTO's perigee
%       direction (rotating-frame rotation sense), sheets at {0,90,180,270}.
%       The physical GTO never changes (350 x 35,786 km); only its
%       orientation relative to the Moon does.
%     - Departure phase is a TIME fraction since perigee over one Kepler
%       period (~0.44 d), reconstructed algebraically -- no propagation.
%
%   The catalog is COMPACT (data-minimized): it stores only canonical
%   nondimensional quantities; days, Delta-V and masses are derived on the
%   fly using the formulas in cat_.derive. This example shows each
%   derivation as it uses it.
%
%   HONESTY CONTRACT: whenever the returned solution differs from the
%   request (nearest sheet, snapped phases, nearest solved pair, different
%   rung), costate_catalog_pick PRINTS A WARNING saying exactly what you
%   are getting. Try changing the request below to an exact stored point
%   and watch the warnings disappear.
%
%   Files needed on the path:
%     costate_catalog_gto_tulip.mat   the catalog (here or in results/)
%     costate_catalog_pick.m          THIS FOLDER'S five-coordinate lookup
%                                     (orientDeg-keyed; NOT the generic
%                                     period-keyed copy -- see README)
%     get_family_orbit.m              the 'gto' pseudo-family + tulip recipes
%     pumpkyn / pumpkynPie            propagation, tfMin
%
%  M. Casey                                                     08/31/2026

%% 1. Load the catalog (beside this script, or in results/):
 catFile = 'costate_catalog_gto_tulip.mat';
if ~exist(catFile, 'file')
 catFile = fullfile('results', 'costate_catalog_gto_tulip.mat');
end
       L = load(catFile);
    cat_ = L.costate_catalog_gto_tulip;
      mu = cat_.constants.muStar;
   tStar = cat_.constants.tStar_s;

fprintf('catalog: %d sheets, %d entries; orientations [%s] deg; petals [%s]\n', ...
        numel(cat_.sheets), cat_.n_entries, ...
        num2str(unique([cat_.sheets.tauDRO])), ...
        num2str(unique([cat_.sheets.Np])));

%% 2. Choose the five coordinates of a transfer:
orientDeg = 90;       % departure orientation (deg) -- sheets at 0,90,180,270
       Np = 7;        % tulip petal count           -- sheets at 5,7,9,12
  depDays = 0.15;     % departure phase (days past GTO perigee, period ~0.44 d)
  arrDays = 3.0;      % arrival phase (days past the tulip reference)
  thrustN = 10.0;     % thrust (N), anywhere in the rung range [15..5]

[tf_nd, z8, info] = costate_catalog_pick(cat_, orientDeg, Np, ...
                                         depDays, arrDays, thrustN);

% derived quantities, straight from cat_.derive (the z8 seed lives at the
% nearest STORED rung -- info.delivered.seedThrustN -- so mass and Delta-V
% below are derived AT THAT RUNG, consistent with z8(8)):
 tf_days = tf_nd * tStar/86400;
     Tnd = (info.delivered.seedThrustN/cat_.thruster.m0_kg) * tStar^2 / ...
           (cat_.constants.lStar_km*1000);
     cnd = cat_.thruster.c_nd;
      mf = 1 - Tnd*z8(8)/cnd;
  dV_kms = cnd*log(1/mf) * cat_.constants.lStar_km/tStar;
fprintf('flight time %.4f ND = %.3f days;  seed rung %.1f N;  dV %.4f km/s\n', ...
        tf_nd, tf_days, info.delivered.seedThrustN, dV_kms);

%% 3. Rebuild the DELIVERED sheet's orbits from the catalog's recipes:
      ks = find(abs([cat_.sheets.tauDRO] - info.delivered.tauDRO) < 1e-9 & ...
                [cat_.sheets.Np] == info.delivered.Np, 1);
      sh = cat_.sheets(ks);
[tG, rvG] = get_family_orbit('gto', sh.dep_params);    % algebraic Kepler locus
[tT, rvT] = get_family_orbit('tulip', sh.arr_params);  % getTulip->cont_np->prop

%% 4. Endpoint states at the DELIVERED phases (info tells the truth).
%  The GTO phase-day conversion uses the TRUE Kepler period (from
%  dep_params.sma_km), never sh.tauDRO (which is orientDeg here):
     muE = 6.67384e-20*(1 - mu)*(5.9736E24 + 7.35E22);
      Pd = 2*pi*sqrt((sh.dep_params.sma_km)^3/muE)/86400;   % days
      Pa = tT(end) * tStar/86400;
%  Interpolate with 'spline' (as the campaign's audit does): the min-time
%  flow is SENSITIVE -- a 0.03 km linear-interp departure error grows to
%  ~500 km at arrival; spline endpoints reproduce the audit's <0.01 km:
     rv0 = interp1(tG, rvG, mod(info.delivered.depDays/Pd,1)*tG(end), 'spline');
     rvf = interp1(tT, rvT, mod(info.delivered.arrDays/Pa,1)*tT(end), 'spline');

%% 5a. Fly the transfer straight from the catalog costates (no solver):
 [tau,y] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(1:6)'; 1; z8(1:7)], Tnd, cnd, mu);
  missKm = sqrt(sum((y(end,1:3) - rvf(1:3)).^2)) * cat_.constants.lStar_km;
fprintf('flown arrival miss: %.3f km\n', missKm);

%% 5b. ...or hand it to tfMin (accepts the seed unchanged, ~1 s):
       z = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), z8, Tnd, cnd, mu);
fprintf('tfMin seed change |z - z8| = %.2e\n', sqrt(sum((z - z8).^2)));

%% 6. Plot:
figure('Color','w');
plot3(y(:,1), y(:,2), y(:,3), 'b', 'LineWidth', 1.2); hold on;
plot3(rvG(:,1), rvG(:,2), rvG(:,3), 'k--');
plot3(rvT(:,1), rvT(:,2), rvT(:,3), 'r--');
plot3(rv0(1), rv0(2), rv0(3), '.g', 'MarkerSize', 18);
plot3(rvf(1), rvf(2), rvf(3), '.r', 'MarkerSize', 18);
axis equal; grid on;
xlabel('x [ND]'); ylabel('y [ND]'); zlabel('z [ND]');
title(sprintf('GTO(orient=%g%c) \\rightarrow tulip(N_p=%d), %.1f N, t_f = %.2f d', ...
      info.delivered.tauDRO, char(176), info.delivered.Np, ...
      info.delivered.seedThrustN, tf_days));
legend('transfer','GTO','tulip','depart','arrive');

%% 7. The orientation axis at a glance: best/worst flight time vs
%  orientation at 10 N, for the delivered petal count -- the pi-dip
%  geometry (180 deg = apogee toward the Moon) is the sparse quadrant, but
%  where it DOES solve, the flight times are ordinary:
figure('Color','w'); hold on
kr10 = find(cat_.rungs_N == 10, 1);
for k = 1:numel(cat_.sheets)
    sh_k = cat_.sheets(k);
    if sh_k.Np ~= info.delivered.Np, continue, end
    tf10 = sh_k.tf_nd(:,:,kr10);
    tf10 = tf10(sh_k.has_solution(:,:,kr10));
    if isempty(tf10), continue, end
    plot(sh_k.tauDRO, min(tf10)*tStar/86400, 'ob', 'MarkerFaceColor','b');
    plot(sh_k.tauDRO, max(tf10)*tStar/86400, 'sr', 'MarkerFaceColor','r');
end
grid on;  xticks(unique([cat_.sheets.tauDRO]));
xlabel('departure orientation [deg]'); ylabel('t_f at 10 N [days]');
title(sprintf('best (blue) and worst (red) phasing vs orientation, N_p = %d', ...
      info.delivered.Np));
