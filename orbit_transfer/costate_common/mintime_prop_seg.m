function [yh, PHI] = mintime_prop_seg(dt, y0, needSTM, Tmax, c, muStar)
%% Purpose:
%
%   One multiple-shooting SEGMENT of the minimum-time PMP flow via pumpkyn
%   tfMinProp, with the 14x14 state transition matrix when asked. The
%   shared propagator closure of the min-time bindings (ms_tfmin and
%   ms_tfmin_hom): extracted 2026-09-08 on its second consumer, per the
%   library rule that code reused twice moves to costate_common.
%
%   MASS-DEPLETION GUARD (2026-11-10). The min-time mass flow is
%   mDot = -u*Tmax/c with u in {0, 1}, so
%
%       m(dt) >= y0(7) - (Tmax/c)*dt
%
%   is an exact lower bound on the segment's end mass. When that bound is
%   <= 0 the segment MAY reach m = 0, where the thrust acceleration T/m
%   blows up, ode45's step collapses, and pumpkyn.cr3bp.tfMinProp restarts
%   the failed integration at the same instant forever (its while loop has
%   no progress check; GTO -> DRO pilot 2026-11-10: one residual evaluation
%   ran 55 min). Such a propagation is REFUSED with
%   mintime_prop_seg:massDepleted before tfMinProp is called; oc.ms_bvp's
%   residual catches the throw and rejects the iterate. At a root the bound
%   is at most the next junction's mass (or the final mass), which is
%   positive in every catalog regime, so a converged solution never trips it.
%
%% Inputs:
%
%  dt                       double                  segment duration, ND
%  y0                       [14 x 1]                [r; v; m; lam]
%  needSTM                  logical                 also return PHI
%  Tmax, c, muStar          double                  as tfMinProp
%
%% Outputs:
%
%  yh                       [14 x 1]                state at dt
%  PHI                      [14 x 14] or []         STM d y(dt) / d y0
%
%% Revision History:
%  M. Casey                                                   (c) 09/08/2026
%  M. Casey  mass-depletion guard (massDepleted)               (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~(isscalar(dt) && isreal(dt) && isfinite(dt) && dt > 0)
    error('mintime_prop_seg:duration', 'segment duration must be a positive finite real scalar');
end
mFloor = y0(7) - Tmax/c*dt;                     % exact end-mass lower bound (header)
if ~(mFloor > 0)
    error('mintime_prop_seg:massDepleted', ...
          'segment may deplete the mass: m0 %.6g - (Tmax/c) dt = %.6g <= 0', y0(7), mFloor);
end
if needSTM, y0 = [y0; reshape(eye(14), [], 1)]; end
[tau, Yout] = pumpkyn.cr3bp.tfMinProp(dt, y0, Tmax, c, muStar);
% an integrator that stops early without throwing hands back a short,
% finite trajectory; a consumer that takes the last row then reports a
% spectrum at the wrong time (Astra review #3, 2026-09-11)
if isempty(tau) || abs(tau(end) - dt) > 1e-9*max(dt, 1) || size(Yout, 2) < 14*(1 + 13*needSTM)
    error('mintime_prop_seg:incomplete', ...
          'propagation stopped at %.9g of the requested %.9g (%d columns)', ...
          tern(isempty(tau), NaN, tau(end)), dt, size(Yout, 2));
end
yh = Yout(end, 1:14)';
if needSTM, PHI = reshape(Yout(end, 15:210), 14, 14); else, PHI = []; end
if ~all(isfinite(yh)) || (needSTM && ~all(isfinite(PHI(:))))
    error('mintime_prop_seg:nonfinite', 'propagation returned non-finite state or STM');
end
end

function v = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end
