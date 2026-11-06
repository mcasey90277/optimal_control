function [res, tStar] = symmetry_residual(x1, T, mu, opts)
%% Purpose:
%
%   How far the mirror image of a state is from the orbit through it:
%   min over t in [0, T] of max|phi(x1, t) - mirror_state(x1)|. Zero (to
%   integrator accuracy) when the periodic orbit through x1 is x-axis
%   symmetric. The check is meaningful only for an OFF-PLANE x1: a state
%   with y = vx = vz = 0 is its own image and the residual is zero at
%   t = 0 whatever the orbit does. A coarse grid over one period locates
%   the minimum, fminbnd refines it on the bracketing interval. This does
%   not test periodicity; the caller checks closure separately.
%
%% Inputs:
%
%  x1                       [6 x 1]                 An off-plane state on
%                                                   the orbit
%
%  T                        double                  Period (ND)
%
%  mu                       double                  Mass ratio
%
%  opts                     struct (optional)       .nGrid [400] coarse
%                                                   grid points over [0, T]
%
%% Outputs:
%
%  res                      double                  The minimum residual
%                                                   (max-abs over the six
%                                                   components, ND)
%
%  tStar                    double                  The minimizing time
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: pumpkyn's own L1 halo, from a point 0.37 T along it
     x0 = [9.2294247990951594E-1; 0; 2.1592634764719765E-1; 4.2453615136710976E-13; 1.2786666755130519E-1; -2.3343675990673871E-12];
     T0 = 1.8036821222727220;  mu0 = 1.215058560962404E-2;
     [~, X] = pumpkyn.cr3bp.prop([0 0.37*T0], x0.', mu0);
     [res, tStar] = symmetry_residual(X(end, 1:6).', T0, mu0);
     fprintf('residual %.2e at t = %.4f\n', res, tStar);
     return
end

if nargin < 4 || isempty(opts), opts = struct(); end
if isfield(opts, 'nGrid') && ~isempty(opts.nGrid), nGrid = opts.nGrid; else, nGrid = 400; end
x1  = x1(:);
Rx  = mirror_state(x1);
tg  = linspace(0, T, nGrid);
[~, X] = pumpkyn.cr3bp.prop(tg, x1.', mu);                 % [nGrid x 6]
d   = max(abs(X(:, 1:6).' - Rx), [], 1);                  % [1 x nGrid]
[~, k] = min(d);
lo  = tg(max(k - 1, 1));  hi = tg(min(k + 1, nGrid));
f   = @(t) resid_at(t, x1, Rx, mu);
[tStar, res] = fminbnd(f, lo, hi, optimset('TolX', 1e-12));
if d(k) < res, res = d(k); tStar = tg(k); end               % never worse than the grid
end

% ------------------------------------------------------------------------
function r = resid_at(t, x1, Rx, mu)
%% Purpose:
%
%   The residual at one time (t = 0 needs no propagation).
%
if t <= 0
    r = max(abs(x1 - Rx));
else
    [~, X] = pumpkyn.cr3bp.prop([0 t], x1.', mu);
    r = max(abs(X(end, 1:6).' - Rx));
end
end
