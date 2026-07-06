r1n = 1; r2n = 1.2; mu = 1; dth = pi/2;
A = sum(sin(dth)*sqrt(r1n*r2n/(1-cos(dth))));
eps = 1e-12;
[t_inner, y_inner] = lambert_tof((2*pi)^2 + eps, r1n, r2n, A, mu);
disp(['y near multirev left edge: ', num2str(y_inner)]);
disp(['t near multirev left edge: ', num2str(t_inner)]);
