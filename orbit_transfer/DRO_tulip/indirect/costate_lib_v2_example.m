%% costate_lib_v2_example
%
%   Worked example for costate_lib_dro_tulip_v2.mat: the DRO -> tulip
%   minimum-time costate library WITH A THRUST AXIS. Load the library, choose
%   a departure phase, an arrival phase, and a thrust level anywhere in the
%   library's range; the flight time is interpolated between the bracketing
%   thrust rungs and a converged costate seed comes back from the nearest
%   rung. Then rebuild both orbits from the library's own parameters, fly the
%   transfer, optionally hand the seed to pumpkyn.cr3bp.tfMin, and plot.
%
%   Files needed on the path:
%     costate_lib_dro_tulip_v2.mat   the library (990 entries)
%     costate_lib_pick.m             phase + thrust lookup helper
%     pumpkyn / pumpkynPie           orbit families, propagation, tfMin
%
%   Every library entry is a converged minimum-time PMP solution verified
%   three ways: multiple-shooting residual ~1e-12, end-to-end PMP flight
%   arriving at the target, and acceptance by tfMin (which returns the seed
%   unchanged in about a second). So step 5a needs no solver at all.
%
%   Thrust levels available: see lib.thruster.thrust_rungs_N. Isp and initial
%   mass are fixed for this library (lib.thruster.isp_s / .m0_kg); the
%   dimensionless thrust for any entry follows from
%   Tmax_nd = (thrust_N/m0_kg)*tStar^2/(lStar_km*1000).
%
%  M. Casey                                                     08/05/2026

%% 1. Load the library (EDIT libFile for your machine):
 libFile = ['/Users/msc/Desktop/optimal_control/orbit_transfer/', ...
            'DRO_tulip/direct/results/costate_lib_dro_tulip_v2.mat'];
       L = load(libFile);
     lib = L.costate_lib_dro_tulip_v2;
      mu = lib.constants.muStar;

fprintf('library: %d entries, thrust rungs [%s] N, Isp %d s, m0 %d kg\n', ...
        lib.n_entries, num2str(lib.thruster.thrust_rungs_N), ...
        lib.thruster.isp_s, lib.thruster.m0_kg);

%% 2. Choose phasing (DAYS) and thrust (N):
 depDays = 20.0;      %departure: days past the DRO reference point
 arrDays = 4.0;     %arrival:   days past the tulip reference point
 thrustN = 1.0;      %anywhere within the library's rung range

[tf_days, e, bracket] = costate_lib_pick(lib, depDays, arrDays, thrustN);

fprintf('%.1f N: min flight time %.3f days (interpolated between %.1f N -> %.3f d and %.1f N -> %.3f d)\n', ...
        thrustN, tf_days, bracket.thrust_N(1), bracket.tf_days(1), ...
        bracket.thrust_N(2), bracket.tf_days(2));
fprintf('seed taken from the %.1f N rung at dep %.3f d, arr %.3f d\n', ...
        e.thrust_N, e.departure_phase_days, e.arrival_phase_days);

%% 3. Rebuild the two orbits from the library's recorded parameters:
    tauD = lib.departure_params.tau;             %DRO selected by PERIOD (ND)
[~,rvD0] = pumpkynPie.cr3bp.getDRO(tauD);
    rvD0 = pumpkyn.cr3bp.cont_np(rvD0, tauD, mu, 1e-12);
[tD,rvD] = pumpkyn.cr3bp.prop(tauD, rvD0, mu);

    tauT = lib.arrival_params.tau;
[~,rvT0] = pumpkyn.cr3bp.getTulip(tauT, lib.arrival_params.Np, ...
                                  lib.arrival_params.pm);
    rvT0 = pumpkyn.cr3bp.cont_np(rvT0, tauT, mu, 1e-12);
[tT,rvT] = pumpkyn.cr3bp.prop(tauT, rvT0, mu);

%% 4. Endpoint states at the entry's phases, and its dimensionless thrust:
     rv0 = interp1(tD, rvD, e.departure_phase_frac*tD(end), 'spline');
     rvf = interp1(tT, rvT, e.arrival_phase_frac*tT(end),  'spline');
     Tnd = (e.thrust_N/lib.thruster.m0_kg) * lib.constants.tStar_s^2 / ...
           (lib.constants.lStar_km*1000);
     cnd = lib.thruster.c_nd;

%% 5a. Fly the transfer directly from the library costates (no solver):
 [tau,y] = pumpkyn.cr3bp.tfMinProp(e.tf_nd, [rv0(1:6)'; 1; e.lambda0], ...
                Tnd, cnd, mu);
  missKm = norm(y(end,1:3) - rvf(1:3)) * lib.constants.lStar_km;
  mFinal = y(end,7) * lib.thruster.m0_kg;
fprintf('flown arrival miss: %.3f km;  propellant used: %.2f kg\n', ...
        missKm, lib.thruster.m0_kg - mFinal);

%% 5b. ...or hand the seed to tfMin (accepts it unchanged, ~1 s):
       z = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), e.z8, Tnd, cnd, mu);
fprintf('tfMin seed change |z - z8| = %.2e\n', norm(z - e.z8));

%% 6. Plot the transfer:
figure('Color','w');
plot3(y(:,1), y(:,2), y(:,3), 'b'); hold on;
plot3(rvD(:,1), rvD(:,2), rvD(:,3), 'k--');
plot3(rvT(:,1), rvT(:,2), rvT(:,3), 'r--');
plot3(rv0(1), rv0(2), rv0(3), '.g', 'MarkerSize', 15);
plot3(rvf(1), rvf(2), rvf(3), '.r', 'MarkerSize', 15);
axis equal; grid on;
xlabel('x [ND]'); ylabel('y [ND]'); zlabel('z [ND]');
title(sprintf('DRO \\rightarrow tulip minimum time, %.1f N, t_f = %.3f days', ...
      e.thrust_N, e.tf_days));
legend('transfer', 'DRO', 'tulip', 'depart', 'arrive');

%% 7. The phasing torus at this thrust, with the selected point starred:
%
%   Both phases wrap, so the space of (departure, arrival) choices is a
%   torus. Two views of the same map are drawn -- the flat unwrap and the
%   genuine 3-D torus, which is rotatable (drag it with the mouse):
%
%     GREEN  the library holds a verified entry at this thrust
%     RED    the pair solves at some other thrust, but not this one
%     GREY   the pair has no entry at any thrust
%     YELLOW STAR  the phase pair selected above
%
      kr = find(lib.grid.thrust_N == e.thrust_N, 1);
fprintf('at %.2f N, %d of %d phase pairs have solutions\n', ...
        e.thrust_N, nnz(lib.grid.has_solution(:,:,kr)), ...
        numel(lib.grid.departure_phase_frac)*numel(lib.grid.arrival_phase_frac));

plot_costate_torus(lib, e.thrust_N, depDays, arrDays);

%% 8. ...and the flight-time surface at this thrust, for contrast:
figure('Color','w');
imagesc(lib.grid.arrival_phase_days, lib.grid.departure_phase_days, ...
        lib.grid.tf_days(:,:,kr));
set(gca,'YDir','normal');  colorbar;  axis square
xlabel('arrival phase along the tulip [days]');
ylabel('departure phase along the DRO [days]');
title(sprintf('minimum flight time [days] at %.2f N (white = no solution)', ...
      e.thrust_N));
