%% costate_library_example
%
%   Worked example for the USER KIT of costate_catalog_dro_tulip_70mN.mat,
%   the 70 mN DRO -> tulip minimum-time costate library: what it holds, how
%   to pick a transfer, how to read what the library established about it,
%   how to hand it to pumpkyn's tfMin, and how to watch it.
%
%     1. describe the library (grids, coverage, status counts)
%     2. pick by phases in DAYS (the honesty contract) and derive the units
%     3. hand the costates to pumpkyn.cr3bp.tfMin at the pick's endpoints
%     4. every transfer recorded at that cell, with its status
%     5. the status legend, and the fastest conjugate-point transfer
%     6. the extremes of the whole library, flight time, with alternatives
%     7. the status map of the phase torus
%     8. a side-by-side movie: a cell's primary and its fastest alternative
%
%   Every number is derived from the catalog's canonical nondimensional
%   values with the formulas in cat_.derive; nothing here solves, flies with
%   our own engine, or re-certifies an entry.
%
%   FILES: the catalog is loaded BY BARE FILENAME, so run this from the
%   folder that holds it; every function called is in that same folder (or
%   is pumpkyn / pumpkynPie, which must be on the path). No absolute paths,
%   no addpath.
%
%   Run time: about 15 seconds (the movie takes most of it).
%
%  M. Casey                                                     10/14/2026

if ~usejava('desktop'), figVis = 'off'; else, figVis = 'on'; end

%% 1. Load the catalog and describe it:
 catFile = 'costate_catalog_dro_tulip_70mN.mat';   % sits beside this script
       L = load(catFile);
    cat_ = L.costate_catalog_dro_tulip_70mN;
      sh = cat_.sheets(1);
   tStar = cat_.constants.tStar_s;
 lStarKm = cat_.constants.lStar_km;
       D = costate_lib_describe(cat_);
fprintf('[1] verdict: %d primaries + %d alternatives, %.1f %% of the %d x %d phase grid solved, thrust %.3f N\n', ...
        D.nEntries, D.nAlternatives, D.coveragePct, numel(sh.sD_frac), numel(sh.sA_frac), cat_.rungs_N(1));

%% 2. Pick by phases in DAYS; the picker says exactly what it delivers:
  tauDep = 1;         % DRO period [ND]  (this library has one sheet: 1 ND = 4.43 d)
      Np = 7;         % tulip petal count [-]
 depDays = 0.5;       % departure phase [days past the DRO reference]
 arrDays = 3.0;       % arrival phase [days past the tulip reference]
 thrustN = 0.07;      % thrust [N]
[tfNd, z8, info] = costate_catalog_pick(cat_, tauDep, Np, depDays, arrDays, thrustN);
%  The units are derived at the thrust the picker DELIVERED (info.delivered),
%  which is the library's rung whatever thrust was asked for:
   ndOf  = @(TN) (TN/cat_.thruster.m0_kg)*tStar^2/(lStarKm*1000);   % thrust acceleration [ND]
     Tnd = ndOf(info.delivered.tfThrustN);                         % the rung of the delivered t_f
     cnd = cat_.thruster.c_nd;                                     % exhaust velocity [ND]
      mf = 1 - Tnd*tfNd/cnd;                                       % final mass fraction [-]
   dvKms = cnd*log(1/mf)*lStarKm/tStar;                            % [km/s]
  propKg = (1 - mf)*cat_.thruster.m0_kg;                           % [kg]
fprintf('[2] verdict: asked dep %.3f d / arr %.3f d at %.3f N, delivered dep %.3f d / arr %.3f d at %.3f N (warned = %d): t_f %.3f d, Delta-V %.4f km/s, %.3f kg propellant\n', ...
        depDays, arrDays, thrustN, info.delivered.depDays, info.delivered.arrDays, info.delivered.tfThrustN, ...
        info.warned, tfNd*tStar/86400, dvKms, propKg);

%% 3. Hand the costates to pumpkyn tfMin at the pick's endpoints:
%  The endpoints are states on the two periodic orbits at the delivered
%  phases (fractions of each period); the orbits come from the catalog's
%  own recipe at its mass ratio. tfMin starts from z8 at the thrust the
%  seed belongs to and returns the root it converges to; |dz| = 0 means the
%  first residual was already at machine level. (Needs the Optimization
%  Toolbox through pumpkyn's fsolve; nothing else in this example does.)
      sD = info.delivered.depDays*86400/tStar/sh.tauDRO;           % phase FRACTION [-]
      sA = info.delivered.arrDays*86400/tStar/sh.period_tulip_nd;
      pD = sh.dep_params;  pD.muStar = cat_.constants.muStar;
      pA = sh.arr_params;  pA.muStar = cat_.constants.muStar;
[tD, rvD] = get_family_orbit(sh.dep_family, pD);
[tT, rvT] = get_family_orbit(sh.arr_family, pA);
     rv0 = interp1(tD/tD(end), rvD, sD, 'spline');                 % [1 x 6] ND
     rvf = interp1(tT/tT(end), rvT, sA, 'spline');
 TndSeed = ndOf(info.delivered.seedThrustN);                       % the rung the z8 seed belongs to
  tic;  zTf = pumpkyn.cr3bp.tfMin(rv0, rvf, z8, TndSeed, cnd, cat_.constants.muStar);  wall = toc;
      dz = sqrt(sum((zTf(:) - z8(:)).^2));
fprintf('[3] verdict: tfMin returned in %.1f s with |dz| = %.2e (t_f %.4f d -> %.4f d)\n', ...
        wall, dz, z8(8)*tStar/86400, zTf(8)*tStar/86400);

%% 4. Every transfer recorded at the delivered cell, ranked by flight time:
 T = transfers_at(cat_, sD, sA);
fprintf('  %-12s %-22s %8s %10s   %s\n', 'kind', 'status', 't_f [d]', 'dV [km/s]', 'reason (first 60 chars)');
for k = 1:numel(T)
    fprintf('  %-12s %d %-20s %8.2f %10.4f   %s\n', T(k).kind, T(k).status, T(k).statusName, ...
            T(k).tfDays, T(k).dvKms, T(k).reason(1:min(60, numel(T(k).reason))));
end
fprintf('[4] verdict: %d transfer(s) at sD = %.4f, sA = %.4f\n', numel(T), sD, sA);

%% 5. The status legend, and the fastest conjugate-point transfer in the library:
 K = status_key();
for k = find(K.codes >= 1)                % the four statuses a shipped entry carries (0 and -1 are builder-internal)
    fprintf('  status %d = %-22s %s\n', K.codes(k), K.names{k}, K.meaning{k});
end
 f2 = struct('includeAlternatives', true, 'status', 2);
[e2, ~] = costate_catalog_extremes(cat_, f2, false, 'time');
fprintf('[5] verdict: fastest status-2 transfer: %.2f d at sD = %.4f, sA = %.4f (%s); its reason: %s\n', ...
        e2.tf_days, e2.dep_frac, e2.arr_frac, e2.kind, e2.reason(1:min(80, numel(e2.reason))));

%% 6. The extremes of the whole library (flight time), alternatives included:
 fAll = struct('includeAlternatives', true);
[eMin, eMax] = costate_catalog_extremes(cat_, fAll, strcmp(figVis, 'on'), 'time');
fprintf('[6] verdict: fastest %.2f d (%s, status %d), slowest %.2f d (%s, status %d)\n', ...
        eMin.tf_days, eMin.kind, eMin.status, eMax.tf_days, eMax.kind, eMax.status);

%% 7. The status map of the phase torus:
 H = plot_status_torus(cat_, struct('visible', strcmp(figVis, 'on')));
fprintf('[7] verdict: torus drawn, %d cells with a status-2 transfer, %d off-grid alternatives (%d status 2)\n', ...
        sum(H.conjCell(:)), H.nOffGrid, H.nOffGridConj);

%% 8. A movie: the cell with the most alternatives, its primary beside its fastest alternative
%  (60 frames here; the shipped status_pair_movie has 240 -- pass
%  struct('pair', ..., 'nFrames', 240) for the full length):
 [nMost, kCell] = max(H.nAltPerCell(:));
 assert(nMost > 0, 'costate_library_example:noAlternatives', 'no cell of this library holds an alternative');
 [iD, iA]   = ind2sub(size(H.nAltPerCell), kCell);
 Tc = transfers_at(cat_, sh.sD_frac(iD), sh.sA_frac(iA));
 ip = find(strcmp({Tc.kind}, 'primary'), 1);
 ia = find(strcmp({Tc.kind}, 'alternative'), 1);          % rows are ranked by t_f: the first is the fastest
 assert(~isempty(ip) && ~isempty(ia), 'costate_library_example:cell', ...
        'cell (%d, %d) holds %d alternatives by index but transfers_at found no primary/alternative pair at its phases', iD, iA, nMost);
 mo = struct('pair', {{Tc(ip), Tc(ia)}}, 'labels', {{'primary', 'fastest alternative'}}, 'nFrames', 60);
 stem = costate_catalog_extremes_movies(cat_, [], 'time', 'example_pair_movie', mo);
fprintf('[8] verdict: movie %s.mp4 / .gif of cell (%d, %d): primary %.2f d (status %d) vs alternative %.2f d (status %d)\n', ...
        stem, iD, iA, Tc(ip).tfDays, Tc(ip).status, Tc(ia).tfDays, Tc(ia).status);

fprintf('ALL SECTIONS DONE\n');
if ~usejava('desktop'), close all; end
