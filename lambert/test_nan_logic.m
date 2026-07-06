r1n = 1; r2n = 1.2; mu = 1; dth = pi/2;
A = sum(sin(dth)*sqrt(r1n*r2n/(1-cos(dth))));
% Feasible region is z > -3.18
% Let's pretend our bracket was zlo = -5, zhi = -2.
% Is it possible to have t(-2) < dt?
% t(-2) ~ 0.2. What if dt = 0.5?
% If zlo = -5 (NaN), zhi = -2 (t=0.2 < 0.5). That bracket has no root.

% Let's test the text's logic. If zm is NaN, we treat as "too long" (zhi = zm).
disp('If we do zhi = zm when zm is NaN:');
disp('We are moving the UPPER limit to the LEFT.');
disp('Since infeasible is on the LEFT, we are moving entirely into the inefasible region!');
