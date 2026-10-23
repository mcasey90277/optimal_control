function ok = test_layouts()
%% Purpose:
%
%   The layouts reproduce the paper's lower bounds (its Table 1 J_LB
%   column: 0.200, 0.400, 0.800, 1.156, 2.297, 4.704 at T = 10), targets
%   are antipodal, the 2D layouts have z = 0, pair_list is complete and
%   ordered, and malformed requests are refused by identifier.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
T = 10;

Ms   = [4 8 16 32 64 128];
JLBp = [0.200 0.400 0.800 1.156 2.297 4.704];
rs   = [0.02 0.02 0.02 0.02 0.01 0.01];
for k = 1:numel(Ms)
    [x0, xT, geom] = drone_layouts('paper', Ms(k));
    J = lower_bound(x0, xT, T);
    ok = chk(ok, abs(J - JLBp(k)) < 6e-4, sprintf('M = %3d: J_LB = %.4f (paper %.3f)', Ms(k), J, JLBp(k)));
    ok = chk(ok, max(max(abs(xT + x0))) < 1e-15, sprintf('M = %3d: targets antipodal', Ms(k)));
    ok = chk(ok, all(x0(3,:) == 0) && all(xT(3,:) == 0) && geom.d == 2, sprintf('M = %3d: planar, d = 2', Ms(k)));
    ok = chk(ok, geom.r == rs(k) && geom.M == Ms(k) && geom.eps == 1e-5 && geom.l == 1e-4 && isempty(geom.ball), ...
             sprintf('M = %3d: geom r = %.2f, eps, l, no ball', Ms(k), rs(k)));
end
[x0, ~, geom] = drone_layouts('paper', 8, 0.05);
ok = chk(ok, geom.r == 0.05 && abs(sqrt(sum(x0(:,1).^2)) - 0.5) < 1e-15, 'explicit r is honoured; circle radius 0.5');

[x0, xT, geom] = drone_layouts('headon');
ok = chk(ok, isequal(x0, [-0.4 0.4; 0 0; 0 0]) && isequal(xT, [0.4 -0.4; 0 0; 0 0]) && geom.M == 2 && geom.r == 0.02, ...
         'headon: the two-agent oracle layout');

[x0, xT, geom] = drone_layouts('sphere', 16);
rad = sqrt(sum(x0.^2, 1));
ok = chk(ok, max(abs(rad - 0.4)) < 1e-12 && max(max(abs(xT + x0))) < 1e-15 && geom.d == 3, 'sphere: radius 0.4, antipodal, d = 3');
ok = chk(ok, isstruct(geom.ball) && isequal(geom.ball.c, [0; 0; 0]) && geom.ball.R == 0.1, 'sphere: ball obstacle c = 0, R = 0.1');
dmin = inf;
for a = 1:15
    for b = a+1:16
        dmin = min(dmin, sqrt(sum((x0(:,a) - x0(:,b)).^2)));
    end
end
ok = chk(ok, dmin > 2*geom.r, sprintf('sphere: initial points clear of each other (min %.3f > 2r)', dmin));

P = pair_list(5);
ok = chk(ok, isequal(P, [1 1 1 1 2 2 2 3 3 4; 2 3 4 5 3 4 5 4 5 5]), 'pair_list(5): the 10 pairs, lexicographic');
ok = chk(ok, isequal(size(pair_list(1)), [2 0]), 'pair_list(1) is [2 x 0]');

try
    drone_layouts('rings', 31);
    ok = chk(ok, false, 'rings with M not divisible refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'drone_layouts:M'), 'rings with M not divisible refused by drone_layouts:M');
end
try
    drone_layouts('hexagon', 6);
    ok = chk(ok, false, 'unknown layout refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'drone_layouts:which'), 'unknown layout refused by drone_layouts:which');
end

if ok, fprintf('TEST_LAYOUTS: ALL PASS\n'); else, fprintf('TEST_LAYOUTS: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
