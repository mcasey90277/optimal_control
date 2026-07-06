r1n = 1; r2n = 1.2; mu = 1; dth = pi/2;
A = sin(dth)*sqrt(r1n*r2n/(1-cos(dth)))   % expect sqrt(1.2)
[t0, y0] = lambert_tof(0, r1n, r2n, A, mu)
lambert_tof(-5, r1n, r2n, A, mu)          % expect NaN
lambert_tof(20, r1n, r2n, A, mu)
zz = linspace(-20, (2*pi)^2 - 1e-6, 2000);
tt = lambert_tof(zz, r1n, r2n, A, mu);
sum(diff(tt(~isnan(tt))) <= 0)            % monotonicity violations
