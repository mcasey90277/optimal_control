r1n = 1; r2n = 1.2; mu = 1; dth = pi/2;
A = sum(sin(dth)*sqrt(r1n*r2n/(1-cos(dth))));
eps = logspace(-1,-12, 10);
zL = (2*pi - eps).^2;
zR = (2*pi + eps).^2;
[tL, yL] = lambert_tof(zL, r1n, r2n, A, mu);
[tR, yR] = lambert_tof(zR, r1n, r2n, A, mu);
disp('Left side limit for y (eps -> 0):');
disp(yL);
disp('Right side limit for y:');
disp(yR);
disp('Left side t (eps -> 0):')
disp(tL);
disp('Right side t:')
disp(tR);
