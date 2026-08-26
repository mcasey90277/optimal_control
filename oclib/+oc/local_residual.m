function [dX, out] = local_residual(X, tGrid, rhs, opts)
%% Purpose:
%
%   THE local-residual engine (oclib move 3): the per-interval TRUE
%   continuous-time error of a direct-transcription solution -- "defect
%   is not accuracy" made executable at interval granularity. For every
%   interval, restart the integrator FROM THE TRANSCRIPTION'S OWN LEFT
%   NODE, integrate the true dynamics (with the caller's reconstructed
%   control inside the rhs closure), and report the miss at the right
%   node:
%
%       dX(:,k) = Phi(t_{k+1}; t_k, X(:,k)) - X(:,k+1)
%
%   This is the number that decides whether a solution is PHYSICAL: on
%   DRO->tulip a defect of 1.4e-14 coexisted with a local error of 1.5
%   (ratio 1e7). Sibling of oc.fly_control, which CARRIES the state
%   (global G1b/G2 flight); this engine RESETS it (local G1 residual).
%   Both leave dynamics and control reconstruction to the caller: those
%   are domain policy.
%
%  ASSUMPTIONS / NOTES:
%
% • dX is returned RAW, per state component. Splitting into position /
%   velocity / mass errors and converting to km, m/s is the CALLER's job
%   -- state layouts differ across campaigns (CR3BP [rv; m] vs MEE
%   [p f g h k L; m]), and norming mixed components is the exact bug an
%   external review caught in the original gate.
% • The independent variable is whatever tGrid is in (time, longitude,
%   Sundman s) -- the engine never interprets it.
% • Integrator tolerances must sit well below the residuals being
%   measured, or the measurement reports its own integrator error
%   (defaults 1e-12 / 1e-14, the G1 convention).
%
%% Inputs:
%
%  X                        [nx x N+1]              Transcription states at
%                                                   the nodes
%
%  tGrid                    [1 x N+1]               Node values of the
%                                                   independent variable
%
%  rhs                      fhandle                 dz = rhs(t, z): true
%                                                   dynamics with the
%                                                   reconstructed control
%                                                   inside
%
%  opts                     struct (optional)
%   .solver                 fhandle                 @ode113 (default)
%   .RelTol                 double                  [1e-12]
%   .AbsTol                 double                  [1e-14]
%
%% Outputs:
%
%  dX                       [nx x N]                Per-interval miss at
%                                                   the right node (flown
%                                                   minus transcription)
%
%  out                      struct                  .tMid [1 x N] interval
%                                                   midpoints, .kWorst
%                                                   (index of max column
%                                                   norm -- convenience
%                                                   only, dimensionless)
%
%% Revision History:
%  M. Casey                                                   (c) 08/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 4, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
solver = d('solver', @ode113);
oo = odeset('RelTol', d('RelTol', 1e-12), 'AbsTol', d('AbsTol', 1e-14));

t = tGrid(:).';
N = numel(t) - 1;
assert(size(X,2) == N+1, 'oc:local_residual:shape', ...
       'X has %d columns for %d node times', size(X,2), N+1);

dX = zeros(size(X,1), N);
for k = 1:N
    [~, Z] = solver(rhs, [t(k) t(k+1)], X(:,k), oo);
    dX(:,k) = Z(end,:).' - X(:,k+1);
end

if nargout > 1
    [~, kW] = max(sqrt(sum(dX.^2, 1)));
    out = struct('tMid', 0.5*(t(1:N)+t(2:N+1)), 'kWorst', kW);
end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
