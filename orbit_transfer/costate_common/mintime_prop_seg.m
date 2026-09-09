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

if needSTM, y0 = [y0; reshape(eye(14), [], 1)]; end
[~, Yout] = pumpkyn.cr3bp.tfMinProp(dt, y0, Tmax, c, muStar);
yh = Yout(end, 1:14)';
if needSTM, PHI = reshape(Yout(end, 15:210), 14, 14); else, PHI = []; end
end
