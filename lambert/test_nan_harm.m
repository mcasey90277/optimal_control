r1n = 1; r2n = 1.2; mu = 1; dth = pi/2;
A = sum(sin(dth)*sqrt(r1n*r2n/(1-cos(dth))));
% Feasible region is z > -3.18.
% t goes to 0 as z -> -3.18 from bounded right side
zfeas = -3.1832;

% Imagine a bracket [zlo, zhi] = [-4, 0]
% -4 is INFEASIBLE. So t(-4) = NaN.
% 0 is FEASIBLE. t(0) = 1.13.
% If dt = 1.0, and our bracket is [-4, 0].
% zm = -2. t(-2) = 0.2 < 1.0. zlo = -2. Next zm = -1. t(-1) = 0.5 < 1.0. 
% Wait, what if bracket is [-4, -3.1]?
% zm = -3.55. t(-3.55) = NaN. 
disp('If zm is NaN (-3.55), does the author logic say zhi = zm?');
disp('If zhi = zm = -3.55, new bracket is [-4, -3.55], both NaN!');
disp('The root is at z > -3.18. We LOST the root! The logic is wrong!');
