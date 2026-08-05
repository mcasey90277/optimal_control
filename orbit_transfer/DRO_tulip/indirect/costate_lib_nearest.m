function e = costate_lib_nearest(lib, depPhaseDays, arrPhaseDays)
%% Purpose:
%
%   Selects the costate-library entry nearest a requested pair of orbit
%   phases (departure phase along the DRO, arrival phase along the tulip).
%   Both phases are periodic, so distance is measured on the phasing torus
%   (with wraparound), each axis in units of its own orbit period. The
%   returned entry holds a converged minimum-time PMP solution: the field
%   z8 = [lambda_r(3); lambda_v(3); lambda_m; tf] can be handed directly
%   to pumpkyn.cr3bp.tfMin (it is accepted unchanged) or flown with
%   pumpkyn.cr3bp.tfMinProp.
%
%  ASSUMPTIONS / NOTES:
%
% • Phases are measured in days past each orbit's reference point (the
%   state returned by getDRO / getTulip), modulo the orbit period.
% • Nearest-grid-entry lookup; for off-grid phases, feed the returned z8
%   to tfMin as the initial guess and let it converge the small correction.
% • lib.grid.has_solution is the at-a-glance availability map of the
%   torus; lib.grid.entry_index maps grid cells to lib.entries rows.
%
%% Inputs:
%
%  lib                      struct                  costate_lib_dro_tulip
%                                                   structure loaded from
%                                                   costate_lib_dro_tulip
%                                                   .mat
%
%  depPhaseDays             double                  Departure phase along
%                                                   the DRO, days past
%                                                   the reference point
%
%  arrPhaseDays             double                  Arrival phase along
%                                                   the tulip, days past
%                                                   the reference point
%
%% Outputs:
%
%  e                        struct                  Nearest library entry:
%                                                   -------
%                                                   lambda0        [7 x 1]
%                                                   tf_nd, tf_days
%                                                   z8             [8 x 1]
%                                                   departure_phase_*
%                                                   arrival_phase_*
%                                                   ms_residual
%                                                   flight_miss_km
%
%% Revision History:
%  M. Casey                                                   (c) 08/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: load the library and pull the entry nearest a phase pair:
       L = load('costate_lib_dro_tulip.mat');
     lib = L.costate_lib_dro_tulip;
       e = costate_lib_nearest(lib, 1.0, 10.0);
     fprintf('nearest entry: dep %.3f d, arr %.3f d, tf = %.3f d\n', ...
             e.departure_phase_days, e.arrival_phase_days, e.tf_days);
     fprintf('z8 = [%s]\n', sprintf('%.6f ', e.z8));
     return;
end

      Pd = lib.departure_params.period_days;
      Pa = lib.arrival_params.period_days;
      fd = mod(depPhaseDays/Pd, 1);
      fa = mod(arrPhaseDays/Pa, 1);

    best = inf;
       e = [];
for k = 1:lib.n_entries
      dd = abs(lib.entries(k).departure_phase_frac - fd);
      dd = min(dd, 1-dd);                          %torus wraparound
      da = abs(lib.entries(k).arrival_phase_frac - fa);
      da = min(da, 1-da);
    if dd^2 + da^2 < best
        best = dd^2 + da^2;
           e = lib.entries(k);
    end
end

end
