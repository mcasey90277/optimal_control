function ok = test_direct_oracle()
%% Purpose:
%
%   THE closed-form oracle for the direct solver: two agents head-on,
%   (-0.4,0) <-> (0.4,0), r = 0.02, T = 10. In centre-of-mass and relative
%   coordinates the cost is integral |cdot|^2 + 1/4 |wdot|^2, c rests, and
%   w = x1 - x2 follows the shortest path from (-0.8,0) to (0.8,0) outside
%   the disk |w| < 0.04 at constant speed:
%       ell = 2 sqrt(D^2 - R^2) + R (pi - 2 acos(R/D)),  D = 0.8, R = 0.04
%       J*  = ell^2 / (4 T) = 0.0641602   (J_LB = 0.064)
%   The NLP knows none of this. Its optimum sits ABOVE J* by O(h^2): with
%   clearance enforced at the Hermite midpoints as well as the nodes, the
%   polygonal path must circumscribe the 2r circle (a chord's midpoint lies
%   inside it). Bars: 5e-4 at N = 100, 5e-5 at N = 400.
%
%   Also: the z rows are exactly zero (2D), the defects are re-evaluated,
%   and a bump of zero is refused.
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
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));

[x0, xT, geom] = drone_layouts('headon');
T = 10;  D = 0.8;  R = 2*geom.r;
ell   = 2*sqrt(D^2 - R^2) + R*(pi - 2*acos(R/D));
Jstar = ell^2/(4*T);
ok = chk(ok, abs(Jstar - 0.06416013) < 5e-8, sprintf('closed form J* = %.7f', Jstar));

out100 = solve_drones_direct(x0, xT, geom, struct('N', 100, 'T', T, 'mode', 'hard'));
ok = chk(ok, out100.success, sprintf('N = 100 solve succeeded (%s, %d iterations, %.1f s)', out100.ipoptStatus, out100.iterations, out100.wall));
rel100 = (out100.J - Jstar)/Jstar;
ok = chk(ok, abs(rel100) < 5e-4, sprintf('N = 100: J = %.7f, (J - J*)/J* = %+.2e (bar 5e-4)', out100.J, rel100));
ok = chk(ok, rel100 >= -1e-6, 'N = 100: the discrete optimum does not fall below J* (midpoint clearance circumscribes the arc)');
ok = chk(ok, out100.J >= out100.JLB, sprintf('J >= J_LB = %.4f', out100.JLB));
ok = chk(ok, max(max(abs(out100.X(3:3:end,:)))) == 0 && max(max(abs(out100.V(3:3:end,:)))) == 0, '2D: z rows exactly zero');
ok = chk(ok, out100.maxDefect < 1e-9, sprintf('defects re-evaluated: %.1e', out100.maxDefect));
ok = chk(ok, out100.minHnodes >= -1e-9, sprintf('nodal clearance %.2e >= 0', out100.minHnodes));

out400 = solve_drones_direct(x0, xT, geom, struct('N', 400, 'T', T, 'mode', 'hard'));
rel400 = (out400.J - Jstar)/Jstar;
ok = chk(ok, out400.success && abs(rel400) < 5e-5, sprintf('N = 400: (J - J*)/J* = %+.2e (bar 5e-5), %.1f s', rel400, out400.wall));
ok = chk(ok, abs(rel100)/abs(rel400) > 3, sprintf('N 100 -> 400: gap shrinks %.1fx (mixed-order terms; > 3 required)', abs(rel100)/abs(rel400)));

try
    solve_drones_direct(x0, xT, geom, struct('bump', 0));
    ok = chk(ok, false, 'bump = 0 refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'solve_drones_direct:bump'), 'bump = 0 refused by solve_drones_direct:bump');
end

if ok, fprintf('TEST_DIRECT_ORACLE: ALL PASS\n'); else, fprintf('TEST_DIRECT_ORACLE: FAILURE\n'); end
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
