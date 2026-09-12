function [yh, PHI] = mintime_prop_seg(dt, y0, needSTM, Tmax, c, muStar)
%% Purpose:
%
%   One multiple-shooting SEGMENT of the minimum-time PMP flow via pumpkyn
%   tfMinProp, with the 14x14 state transition matrix when asked. The
%   shared propagator closure of the min-time bindings (ms_tfmin and
%   ms_tfmin_hom): extracted 2026-09-08 on its second consumer, per the
%   library rule that code reused twice moves to costate_common.
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
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~(isscalar(dt) && isreal(dt) && isfinite(dt) && dt > 0)
    error('mintime_prop_seg:duration', 'segment duration must be a positive finite real scalar');
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
