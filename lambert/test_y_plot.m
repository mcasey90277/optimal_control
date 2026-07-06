r1n = 1; r2n = 1.2; mu = 1; dth = pi/2;
A = sin(dth)*sqrt(r1n*r2n/(1-cos(dth)));
z = linspace(-20, 20, 1000);
[t, y] = lambert_tof(z, r1n, r2n, A, mu);

A2 = sum(sin(3*pi/2)*sqrt(r1n*r2n/(1-cos(3*pi/2))));  % A < 0 case
[t2, y2] = lambert_tof(z, r1n, r2n, A2, mu);

disp('A > 0 case: Is y increasing?');
disp(all(diff(y) >= 0));
disp('Where is y < 0?');
disp(z(find(y < 0, 1, 'last')));

disp('A < 0 case: Is y increasing?');
disp(all(diff(y2) >= 0));
disp('Where is y2 < 0?');
disp(z(find(y2 < 0, 1, 'last')));
