function [x, seam, dx] = phase_state(t, y, opts)
%% Purpose:
%
%   One periodic orbit in, two closures out: the state at a PHASE FRACTION
%   of the period, and its derivative with respect to that fraction.
%
%     [xA, seam, dxA] = phase_state(tTulip, rvTulip);
%     rvf = xA(0.0754);           % the arrival endpoint
%     dxA(0.0754)                 % what a continuation in arrival phase needs
%
%   This is the endpoint rule every campaign writes by hand -- a dozen
%   interp1(..., 'spline') sites across the catalogs, plus the private
%   copies in transfer_study and arclength_arrival. Built on periodic_pp,
%   so the interpolant is C1 across the seam at s = 0 and the derivative is
%   the derivative of the function the boundary condition actually matches.
%
%  ASSUMPTIONS / NOTES:
%
% • The table covers ONE period starting at t = 0; t(end) IS the period.
% • The phase wraps: s, s+1 and s-1 give the same state, bitwise.
% • dx is per unit PHASE (it carries the period), not per unit time.
%
%% Inputs:
%
%  t                        [m x 1] or [1 x m]      sample times of one
%                                                   period, t(1) = 0
%  y                        [m x n]                 samples, one row per time
%  opts                     struct (optional)       passed to periodic_pp
%                                                   (.scheme, .onMissing)
%
%% Outputs:
%
%  x                        fhandle                 x(s) -> [n x 1] state at
%                                                   phase fraction s
%  seam                     struct                  periodic_pp's seam
%                                                   (.value .deriv .scheme)
%  dx                       fhandle                 dx(s) -> [n x 1], d x / d s
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
            thD = linspace(0, 2*pi, 61).';
    [xD, sD, dxD] = phase_state(thD, [cos(thD), sin(thD)]);
    fprintf('circle at s = 0.25: x = [%.3f %.3f], dx/ds = [%.3f %.3f] (analytic [%.3f %.3f])\n', ...
            xD(0.25), dxD(0.25), -2*pi*sin(2*pi*0.25), 2*pi*cos(2*pi*0.25));
    [x, seam, dx] = deal(xD, sD, dxD);
    return;
end

if nargin < 3, opts = struct(); end
[pp, seam, dpp] = periodic_pp(t, y, opts);
                t = t(:).';
if abs(t(1)) > eps(t(end))
    error('phase_state:tStart', ...
          'the table must cover one period starting at t = 0 (t(1) = %g)', t(1));
end
           period = t(end);

                x = @(s) ppval(pp,  mod(s, 1)*period);
               dx = @(s) period*ppval(dpp, mod(s, 1)*period);
end
