function ok = assert_periodic_orbit(tau, rv, tol, throwOnFail)
%% Purpose:
%
%   PERIODICITY GUARD for a propagated orbit, made a single-home helper
%   (migration #4): an interpolated family seed that cont_np could not
%   truly converge produces a non-closing "orbit" whose downstream metrics
%   are garbage (measured in the survey: periselene below the lunar
%   surface, distances of 1e10 km). Every consumer of get_family_orbit
%   that cannot tolerate a junk orbit should call this.
%
%  ASSUMPTIONS / NOTES:
%
% • NEGATED comparison, deliberately: a NaN closure (integrator blow-up
%   mid-arc) must fail this test too, and NaN passes any '>' comparison.
%   ~(err < tol) is NaN-safe; (err >= tol) is not.
%
%% Inputs:
%
%  tau                      [M x 1]                 Time along the orbit
%                                                   (unused; kept so the
%                                                   call site reads
%                                                   assert_periodic_orbit(
%                                                   tau, rv, ...) straight
%                                                   from get_family_orbit)
%
%  rv                       [M x 6]                 States over one period
%
%  tol                      double                  Closure tolerance (ND)
%                                                   [default 1e-6]
%
%  throwOnFail              logical                 error() on failure
%                                                   [default true]; false
%                                                   returns ok = false
%
%% Outputs:
%
%  ok                       logical                 Orbit closes and is
%                                                   finite
%
%% Revision History:
%  M. Casey                                                   (c) 08/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3 || isempty(tol), tol = 1e-6; end
if nargin < 4, throwOnFail = true; end

perErr = norm(rv(end,1:6) - rv(1,1:6));
ok = (perErr < tol) && ~any(~isfinite(rv(:)));
if ~ok && throwOnFail
    error('assert_periodic_orbit:open', ...
        'orbit does not close: |closure| = %.3e (tol %.1e), finite = %d', ...
        perErr, tol, ~any(~isfinite(rv(:))));
end
end
