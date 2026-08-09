%% costate_catalog_hh_example
%
%   Worked example for the L1 <-> L2 HALO-TO-HALO minimum-time costate
%   catalogs -- the first catalogs on the SCHEMA-V2 arrival-period axis.
%   Because BOTH ends are halos, there is no petal count: the picker's
%   third argument is the ARRIVAL PERIOD (ND), and sheets are selected by
%   nearest (departure period, arrival period). Everything else works
%   exactly as in the tulip catalogs.
%
%   TWO catalogs ship in this folder (no symmetry maps one direction onto
%   the other):
%     costate_catalog_halo_halo.mat     L1 -> L2   (1,952 entries)
%     costate_catalog_halo_halo_B.mat   L2 -> L1   (2,096 entries)
%
%   HONESTY CONTRACT: whenever the returned solution differs from the
%   request (nearest sheet, snapped phases, nearest solved pair, different
%   rung), costate_catalog_pick PRINTS A WARNING saying exactly what you
%   are getting.
%
%   Files needed on the path:
%     the two catalog .mats            (sit beside this script)
%     costate_catalog_pick.m           five-coordinate lookup
%     get_family_orbit.m               family name + params -> orbit
%     catalog_schema.m                 validator + named derive registry
%     pumpkyn / pumpkynPie             orbit families, propagation, tfMin
%
%  M. Casey                                                     08/09/2026
%
%% 1. Load the L1 -> L2 catalog and validate it against the schema:
 catFile = 'costate_catalog_halo_halo.mat';
       L = load(catFile);
    cat_ = L.costate_catalog_halo_halo;
      mu = cat_.constants.muStar;
   tStar = cat_.constants.tStar_s;
problems = catalog_schema('validate', cat_);
assert(isempty(problems), 'catalog failed schema validation');
fprintf('catalog: schema v%d, %d sheets, %d entries\n', ...
        cat_.schema, numel(cat_.sheets), cat_.n_entries);
fprintf('L1 departures [%s] ND x L2 arrivals [%s] ND\n', ...
        num2str(unique([cat_.sheets.tau_dep])), ...
        num2str(unique([cat_.sheets.tau_arr])));

%% 2. Choose the five coordinates (NOTE: arg 3 is the ARRIVAL PERIOD):
  tauL1 = 2.7433;     % L1 departure period (ND) -- sheets at 1.8037, 2.7433
  tauL2 = 2.8;        % L2 arrival period (ND)   -- sheets at 1.75...3.4
 depDays = 6.0;       % departure phase (days past the L1 halo reference)
 arrDays = 3.0;       % arrival phase (days past the L2 halo reference)
 thrustN = 2.0;       % thrust (N), anywhere in the rung range

[tf_nd, z8, info] = costate_catalog_pick(cat_, tauL1, tauL2, depDays, arrDays, thrustN);

% derived quantities through the NAMED registry (one home for formulas):
   seedN = info.delivered.seedThrustN;
 tf_days = catalog_schema('derive', cat_, 'days', struct('t_nd', tf_nd));
  dV_kms = catalog_schema('derive', cat_, 'deltaV_kms', ...
               struct('tf_nd', z8(8), 'thrustN', seedN));
fprintf(['flight time %.4f ND = %.3f days (at %.2f N);  seed rung %.1f N;  ', ...
         'seed dV %.4f km/s\n'], tf_nd, tf_days, info.delivered.tfThrustN, ...
        seedN, dV_kms);

%% 3. Rebuild BOTH halos from the sheet's recipes:
      ks = find(abs([cat_.sheets.tau_dep] - info.delivered.tauDRO) < 1e-9 & ...
                abs([cat_.sheets.tau_arr] - info.delivered.arrKey) < 1e-9, 1);
      sh = cat_.sheets(ks);
[tD,rvD] = get_family_orbit(sh.dep_family, sh.dep_params);   % the L1 halo
[tA,rvA] = get_family_orbit(sh.arr_family, sh.arr_params);   % the L2 halo

%% 4. Endpoint states at the DELIVERED phases:
      Pd = sh.tau_dep * tStar/86400;
      Pa = sh.tau_arr * tStar/86400;
     rv0 = interp1(tD/tD(end), rvD, mod(info.delivered.depDays/Pd,1), 'spline');
     rvf = interp1(tA/tA(end), rvA, mod(info.delivered.arrDays/Pa,1), 'spline');

%% 5a. Fly the transfer straight from the catalog costates (no solver):
     Tnd = catalog_schema('derive', cat_, 'thrust_nd', struct('thrustN', seedN));
     cnd = cat_.thruster.c_nd;
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
plot3(rvA(:,1), rvA(:,2), rvA(:,3), 'r--');
plot3(rv0(1), rv0(2), rv0(3), '.g', 'MarkerSize', 18);
plot3(rvf(1), rvf(2), rvf(3), '.r', 'MarkerSize', 18);
axis equal; grid on;
xlabel('x [ND]'); ylabel('y [ND]'); zlabel('z [ND]');
title(sprintf('L1(\\tau=%.2f) \\rightarrow L2(\\tau=%.2f), %.1f N, t_f = %.2f d', ...
      sh.tau_dep, sh.tau_arr, seedN, z8(8)*tStar/86400));
legend('transfer','L1 halo','L2 halo','depart','arrive');

%% 7. The reverse direction lives in its own catalog (no symmetry):
LB = load('costate_catalog_halo_halo_B.mat');
cB = LB.costate_catalog_halo_halo_B;
[tfB, ~, infoB] = costate_catalog_pick(cB, tauL2, tauL1, arrDays, depDays, thrustN, false);
fprintf(['reverse (L2 %.2f -> L1 %.4f): tf = %.4f ND = %.3f d ', ...
         '(forward was %.3f d)\n'], infoB.delivered.tauDRO, ...
        infoB.delivered.arrKey, tfB, ...
        catalog_schema('derive', cB, 'days', struct('t_nd', tfB)), tf_days);
