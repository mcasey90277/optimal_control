%% costate_catalog_example
%
%   Worked example for costate_catalog_dro_tulip_70mN.mat: the 70 mN
%   DRO -> tulip minimum-time costate library. Pick five coordinates (DRO
%   period, tulip petal count, departure phase, arrival phase, thrust), get
%   the converged costate seed z8 the library holds there, fly it from the
%   departure orbit to the tulip, and hand it to pumpkyn's tfMin to confirm
%   that tfMin accepts it unchanged.
%
%   The catalog is COMPACT: days, Delta-V and propellant are derived from
%   the canonical nondimensional numbers with the formulas in cat_.derive,
%   and section 3 shows each derivation as it is used.
%
%   HONESTY CONTRACT: whenever the returned solution differs from the
%   request (a snapped phase pair, another sheet, another rung),
%   costate_catalog_pick PRINTS A WARNING saying exactly what you are
%   getting. The request below is deliberately off-grid, so you will see
%   one; ask for an exact stored phase and the warning disappears.
%
%   FILES: this script loads the catalog BY BARE FILENAME, so run it from
%   the folder that holds the catalog; every function it calls is in that
%   same folder (or is pumpkyn / pumpkynPie, which must be on the path).
%   Nothing here needs an absolute path or an addpath.
%
%   Run time: about 10 seconds.
%
%  M. Casey                                                     10/13/2026

if ~usejava('desktop'), figVis = 'off'; else, figVis = 'on'; end

%% 1. Load the catalog and describe it (counts, grids in days, coverage):
 catFile = 'costate_catalog_dro_tulip_70mN.mat';   % sits beside this script
       L = load(catFile);
    cat_ = L.costate_catalog_dro_tulip_70mN;
   tStar = cat_.constants.tStar_s;
 lStarKm = cat_.constants.lStar_km;
       D = costate_lib_describe(cat_);
fprintf('[1] verdict: %d entries, %.1f %% of the phase grid solved, thrust %.3f N\n', ...
        D.nEntries, D.coveragePct, cat_.rungs_N(1));

%% 2. Choose the five coordinates and look the transfer up:
  tauDep = 1;         % DRO period [ND]  (this library has one sheet: 1 ND = 4.43 d)
      Np = 7;         % tulip petal count [-]
 depDays = 0.5;       % departure phase [days past the DRO reference]
 arrDays = 3.0;       % arrival phase [days past the tulip reference]
 thrustN = 0.07;      % thrust [N]

[tfNd, z8, info] = costate_catalog_pick(cat_, tauDep, Np, depDays, arrDays, thrustN);
fprintf('[2] verdict: asked dep %.3f d / arr %.3f d, delivered dep %.3f d / arr %.3f d (warned = %d)\n', ...
        depDays, arrDays, info.delivered.depDays, info.delivered.arrDays, info.warned);

%% 3. Derive days, Delta-V and propellant from the canonical numbers:
     sh = cat_.sheets(1);
  tfDays = tfNd*tStar/86400;                                       % [days]
     Tnd = (info.delivered.tfThrustN/cat_.thruster.m0_kg)*tStar^2/(lStarKm*1000);
     cnd = cat_.thruster.c_nd;
      mf = 1 - Tnd*tfNd/cnd;                                       % final mass fraction [-]
  dvKms  = cnd*log(1/mf)*lStarKm/tStar;                            % [km/s]
 propKg  = (1 - mf)*cat_.thruster.m0_kg;                           % [kg]
fprintf('[3] verdict: t_f = %.3f d, Delta-V = %.4f km/s, propellant = %.3f kg of %.1f kg (rung %.3f N)\n', ...
        tfDays, dvKms, propKg, cat_.thruster.m0_kg, info.delivered.seedThrustN);

%% 4. The endpoints: orbits rebuilt from the catalog's own keys:
 [B, orbits] = catalog_endpoints(cat_, 1, thrustN);
         sD = info.delivered.depDays*86400/tStar/sh.tauDRO;        % phase FRACTION [-]
         sA = info.delivered.arrDays*86400/tStar/sh.period_tulip_nd;
        rv0 = B.stateD(sD);
        rvf = B.stateA(sA);
fprintf('[4] verdict: departure orbit %d samples, arrival orbit %d samples; start r = [%.4f %.4f %.4f] ND, target r = [%.4f %.4f %.4f] ND\n', ...
        numel(orbits.tD), numel(orbits.tT), rv0(1:3), rvf(1:3));

%% 5. Fly the transfer from the catalog costates (no solver) and summarize:
 F = fly_catalog_entry(cat_, struct('z8', z8, 'sD', sD, 'sA', sA), struct('B', B));
print_transfer_summary(F);
fprintf('[5] verdict: flown arrival miss %.3f km and %.4f m/s (method %s)\n', ...
        F.flyKm, F.flyVms, F.method);

%% 6. Does pumpkyn tfMin accept the seed unchanged? (|dz| = |z_tfMin - z8|):
% |dz| = 0 exactly means tfMin's first residual was already at machine level,
% so the solver did not move: the seed is accepted unchanged.
[acc, dz, ~, tinfo] = tfmin_accepts(cat_, struct('z8', z8, 'sD', sD, 'sA', sA), struct('B', B));
fprintf('    (|dz| = 0 exactly = tfMin''s first residual was already at machine level; the solver did not move)\n');
fprintf('[6] verdict: |dz| = %.2e (tolerance %.0e) -> %s (tfMin took %.1f s)\n', ...
        dz, tinfo.tolDz, ternary(acc, 'ACCEPTED unchanged', 'NOT accepted'), tinfo.wallSec);

%% 7. Plot the transfer in 3D (rotating Earth-Moon frame):
 P = plot_transfer_3d(F, F.B, struct('visible', strcmp(figVis, 'on')));
fprintf('[7] verdict: figure drawn, start miss %.3f km, end miss %.3f km, t_f %.2f d\n', ...
        P.startMissKm, P.endMissKm, P.tfDays);

fprintf('ALL SECTIONS DONE\n');
if ~usejava('desktop'), close all; end

function s = ternary(c, a, b)
% TERNARY  a if c is true, b otherwise (local helper for the verdict line).
if c, s = a; else, s = b; end
end
