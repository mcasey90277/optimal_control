function F = mintime_rhs_point(y, Tmax, c, muStar)
%% Purpose:
%
%   The 14-state minimum-time PMP vector field at a point (pumpkyn
%   tfMinEoM, state + costate rows). Shared by ms_tfmin and ms_tfmin_hom;
%   extracted 2026-09-08 on its second consumer.
%
%% Inputs:
%
%  y                        [14 x 1]                [r; v; m; lam]
%  Tmax, c, muStar          double                  as tfMinEoM
%
%% Outputs:
%
%  F                        [14 x 1]                dy/dt
%
%% Revision History:
%  M. Casey                                                   (c) 09/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

F = pumpkyn.cr3bp.tfMinEoM(0, [y(1:14); reshape(eye(14), [], 1)], Tmax, c, muStar);
F = F(1:14);
end
