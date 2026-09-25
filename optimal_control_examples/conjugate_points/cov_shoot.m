function S = cov_shoot(prob, p, tEnd)
%% Purpose:
%
%   Fly ONE Euler-Lagrange extremal from the left end, y(a) = ya, with
%   initial slope y'(a) = p, together with its Jacobi field h:
%
%       y'' = g(t, y, y'),             y(a) = ya,  y'(a) = p,
%       h'' = g_y h + g_yp h',         h(a) = 0,   h'(a) = 1.
%
%   h is the derivative of the trajectory with respect to the initial
%   slope, h(t) = dy(t)/dp: a neighbour started with slope p + delta is
%   y(t) + delta h(t) to first order. A zero of h after a is therefore a
%   point where the neighbours refocus -- a CONJUGATE POINT of a.
%
%   The flight stops early (S.tStop < tEnd, S.ok = false) if the state
%   leaves the domain: non-finite, or |y|, |y'| above 1e6.
%
%% Inputs:
%
%  prob                     struct                  cov_problem output
%  p                        scalar                  initial slope y'(a)
%  tEnd                     scalar                  final time, > a
%
%% Outputs:
%
%  S                        struct                  .sol (ode45 solution,
%                                                   states [y y' h h'], use
%                                                   deval) .p .tStop .ok
%                                                   .tConj (zeros of h in
%                                                   (a, tStop], sorted)
%                                                   .yEnd (y at tStop)
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    %Demo: the oscillator's Jacobi field is sin(t); first zero at pi
    prob = cov_problem('yp^2 - y^2', 0, 4, 0, 1);
    S = cov_shoot(prob, 1/sin(4), 4);
    fprintf('conjugate point(s): %s   (pi = %.10f)\n', mat2str(S.tConj, 12), pi);
    return
end
a = prob.a;
assert(tEnd > a, 'cov_shoot:tEnd', 'tEnd must exceed a');
big = 1e6;
rhs  = @(t, z) [z(2); prob.g(t, z(1), z(2)); z(4); ...
                prob.gy(t, z(1), z(2))*z(3) + prob.gyp(t, z(1), z(2))*z(4)];
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, 'Events', @events);
z0 = [prob.ya; p; 0; 1];

S = struct('sol', [], 'p', p, 'tStop', a, 'ok', false, 'tConj', [], 'yEnd', NaN);
warnState = warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
restore = onCleanup(@() warning(warnState));
try
    sol = ode45(@(t, z) guard(rhs(t, z)), [a tEnd], z0, opts);
catch
    return                                   % the field left its domain
end
S.sol = sol;
S.tStop = sol.x(end);
S.yEnd = sol.y(1, end);
S.ok = abs(S.tStop - tEnd) <= 1e-12*max(1, abs(tEnd));
if isfield(sol, 'xe') && ~isempty(sol.xe)
    keep = sol.ie == 2 & sol.xe > a + 1e-9*(tEnd - a);
    S.tConj = sort(sol.xe(keep));
end

    function [value, isterminal, direction] = events(~, z)
        % 1: leaving the domain (terminal); 2: a zero of h (recorded only)
        value      = [big - max(abs(z(1:2))); z(3)];
        isterminal = [1; 0];
        direction  = [0; 0];
    end
end

% ---------------------------------------------------------------------------
function dz = guard(dz)
% GUARD  Refuse a non-finite field value so the solver stops instead of
% integrating NaN. INPUTS/OUTPUTS: dz [4x1].
if ~all(isfinite(dz))
    error('cov_shoot:domain', 'the Euler-Lagrange field is not finite here');
end
end
