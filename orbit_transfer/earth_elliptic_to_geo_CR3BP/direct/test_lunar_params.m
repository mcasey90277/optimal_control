% TEST_LUNAR_PARAMS  Unit conversions + physical sanity of the Moon spec.
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths;
par  = kepler_lt_params(10, 1500, 2000);
pert = lunar_params(par, 0, 1);
% (a) canonical-unit roundtrips against physical values
assert(abs(pert.DM * par.LU_km - 384400) < 1e-6, 'DM roundtrip [km]');
assert(abs(pert.muM/par.mu - 4902.800/398600.47) < 1e-9, 'mass-ratio 0.0123');
periodDays = (2*pi/pert.nM) * par.TU_s / 86400;
assert(abs(periodDays - 27.32) < 0.05, 'sidereal month ~27.32 d');
% (b) sidereal rate consistent with two-body circular orbit about barycenter
assert(abs(pert.nM - sqrt((par.mu + pert.muM)/pert.DM^3)) < 1e-14, 'nM formula');
% (c) tidal acceleration at GEO radius ~ 7.3e-6 m/s^2 (spec sec 7)
aTide = 2*pert.muM*1/pert.DM^3 * par.AU_ms2;       % r = 1 LU = GEO radius
assert(abs(aTide - 7.3e-6) < 0.3e-6, 'lunar tide at GEO ~7.3e-6 m/s^2');
% (d) phi0/gain pass-throughs
p2 = lunar_params(par, 1.25, 0.5);
assert(p2.phi0 == 1.25 && p2.gain == 0.5, 'phi0/gain stored');
fprintf('test_lunar_params: ALL PASS\n');
