% TEST_PARAMS  Unit checks for kepler_lt_params (canonical units + paper BCs).
p = kepler_lt_params(10, 1500, 2000);
assert(abs(p.TU_s - 13713.8) < 1.0,      'TU wrong');
assert(abs(p.VU_kms - 3.0747) < 1e-3,    'VU wrong');
assert(abs(p.Tmax - 0.029735) < 1e-5,    'nondim thrust wrong');
assert(abs(p.c - 6.3790) < 1e-3,         'nondim exhaust velocity wrong');
assert(p.mu == 1 && p.pSund == 1.5);
% paper initial-orbit geometry in these units
P0 = 11625/p.LU_km;  e0 = 0.75;  a0 = P0/(1-e0^2);
assert(abs(P0 - 0.275703) < 1e-5);
assert(abs(a0*(1+e0) - 1.102810) < 1e-5);   % apogee ~46,500 km
assert(abs(a0*(1-e0) - 0.157544) < 1e-5);   % perigee ~6,643 km
assert(abs(2*atand(0.0612) - 7.0052) < 1e-3); % inclination ~7 deg
fprintf('test_params: ALL PASS\n');
