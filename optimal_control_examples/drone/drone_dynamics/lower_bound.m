function JLB = lower_bound(x0, xT, T)
%% Purpose:
%
%   The unconstrained lower bound on the multi-agent minimum-energy cost
%   J = integral sum_i 1/2 |v_i|^2 dt over fixed T: straight lines at
%   constant speed, J_LB = sum_i |xT_i - x0_i|^2 / (2 T). Every clearance
%   constraint can only raise the cost above this.
%
%% Inputs:
%
%  x0                       [3 x M]                 Initial positions
%
%  xT                       [3 x M]                 Final positions
%
%  T                        double                  Horizon
%
%% Outputs:
%
%  JLB                      double                  The bound
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the paper's M = 4 circle, J_LB = 0.200:
     [a, b] = drone_layouts('paper', 4);
     fprintf('J_LB(M = 4) = %.4f\n', lower_bound(a, b, 10));
     return
end

JLB = sum(sum((xT - x0).^2))/(2*T);
end
