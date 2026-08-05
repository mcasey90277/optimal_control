%% costate_lib_example
%
%   Worked example for costate_lib_dro_tulip.mat: load the library, pick
%   the entry nearest a requested (departure, arrival) phasing, rebuild
%   both orbits from the library's recorded parameters, fly the transfer
%   straight from the stored costates, optionally hand the seed to
%   pumpkyn.cr3bp.tfMin, and plot the result.
%
%   Files needed on the path:
%     costate_lib_dro_tulip.mat     the library
%     costate_lib_nearest.m         phase-pair lookup helper
%     pumpkyn / pumpkynPie          orbit families, propagation, tfMin
%
%   Every library entry is a converged minimum-time PMP solution
%   (multiple-shooting residual ~1e-11; flown arrival miss ~meters), so
%   step 5a needs no solver at all, and in step 5b tfMin accepts the seed
%   essentially unchanged in about a second.
%
%  M. Casey                                                    08/04/2026

%% 1. Load the library:
       L = load('costate_lib_dro_tulip.mat');
     lib = L.costate_lib_dro_tulip;
      mu = lib.constants.muStar;

%% 2. Pick the entry nearest the phasing you want (phases in DAYS):
 depDays = 1.0;      %departure: days past the DRO reference point
 arrDays = 10.0;     %arrival:   days past the tulip reference point
       e = costate_lib_nearest(lib, depDays, arrDays);

fprintf('using entry at dep %.2f d, arr %.2f d;  transfer time %.2f days\n', ...
        e.departure_phase_days, e.arrival_phase_days, e.tf_days);

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

%% 4. Endpoint states at the entry's phases:
     rv0 = interp1(tD, rvD, e.departure_phase_frac*tD(end), 'spline');
     rvf = interp1(tT, rvT, e.arrival_phase_frac*tT(end),  'spline');

%% 5a. Fly the transfer directly from the library costates (no solver):
 [tau,y] = pumpkyn.cr3bp.tfMinProp(e.tf_nd, [rv0(1:6), 1, e.lambda0'], ...
                lib.thruster.Tmax_nd, lib.thruster.c_nd, mu);
  missKm = norm(y(end,1:3) - rvf(1:3)) * lib.constants.lStar_km;
fprintf('flown arrival miss: %.3f km\n', missKm);

%% 5b. ...or hand it to tfMin (accepts the seed unchanged, ~1 s):
       z = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), e.z8, ...
               lib.thruster.Tmax_nd, lib.thruster.c_nd, mu);
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
title(sprintf('DRO \\rightarrow tulip minimum-time transfer, t_f = %.2f days', ...
      e.tf_days));
legend('transfer', 'DRO', 'tulip', 'depart', 'arrive');
