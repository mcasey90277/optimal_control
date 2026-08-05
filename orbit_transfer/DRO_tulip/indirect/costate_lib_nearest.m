function e = costate_lib_nearest(lib, depPhaseDays, arrPhaseDays)
% COSTATE_LIB_NEAREST  Pick the library entry nearest a requested phase pair.
%
% Phases are periodic: distance is computed on the torus (wraparound), in
% units of each orbit's period. Returns the nearest entry with a solution.
%
% INPUTS:
%   lib          - costate_lib_dro_tulip structure (loaded from the .mat)
%   depPhaseDays - requested departure phase along the DRO [days, scalar]
%   arrPhaseDays - requested arrival phase along the tulip [days, scalar]
%
% OUTPUTS:
%   e - the nearest lib.entries element (.z8 is the tfMin-ready seed)
%
% REFERENCES:
%   [1] build_costate_lib.m -- the library format.

Pd = lib.departure_params.period_days;
Pa = lib.arrival_params.period_days;
fd = mod(depPhaseDays/Pd, 1);
fa = mod(arrPhaseDays/Pa, 1);
best = inf;  e = [];
for k = 1:lib.n_entries
    dd = abs(lib.entries(k).departure_phase_frac - fd);  dd = min(dd, 1-dd);
    da = abs(lib.entries(k).arrival_phase_frac   - fa);  da = min(da, 1-da);
    if dd^2 + da^2 < best
        best = dd^2 + da^2;
        e = lib.entries(k);
    end
end
end
