function u = ctrl_quad(ua, um, ub, w)
%% Purpose:
%
%   Lagrange quadratic control reconstruction through the node/midpoint/
%   node samples of a Hermite-Simpson interval, at normalized time w.
%   Extracted verbatim from certify_dro_mintime/local_q (migration #4);
%   shared by flown_control_error and true_min_altitude.
%
%% Inputs:
%
%  ua                       [4 x 1]                 Control at w = 0
%
%  um                       [4 x 1]                 Control at w = 1/2
%
%  ub                       [4 x 1]                 Control at w = 1
%
%  w                        double                  Normalized time in the
%                                                   interval [0, 1]
%
%% Outputs:
%
%  u                        [4 x 1]                 Reconstructed control
%
%% Revision History:
%  M. Casey                                                   (c) 08/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

u = (2*w-1).*(w-1).*ua + 4*w.*(1-w).*um + w.*(2*w-1).*ub;
end
