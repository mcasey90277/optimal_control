r1n = 1; r2n = 1.2; mu = 1; dth = pi/2;
A = sin(dth)*sqrt(r1n*r2n/(1-cos(dth)));
zhi = -3.2; % Infeasible!
zlo = -3.1; % Feasible!
[~, y1] = lambert_tof(zhi, r1n, r2n, A, mu);
[~, y2] = lambert_tof(zlo, r1n, r2n, A, mu);
disp(['y at zhi: ', num2str(y1)]);
disp(['y at zlo: ', num2str(y2)]);
