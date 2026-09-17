function [F, G] = cartpole_field(x, p)
%% Purpose:
%
%   The cart-pole dynamics, split so the control enters affinely:
%   xdot = F(x) + G(x) u. The split is what makes the Pontryagin control
%   explicit -- with a running cost u^2, stationarity gives u* = -lam'G/2 --
%   so this file, not the PMP field, is the one home for the physics.
%
%   Same equations and constants as the direct example
%   (ex2_cart_pole_swing_up/try2/cart_accel, pendulum_accel), which
%   tests/test_cartpole_field pins.
%
%  ASSUMPTIONS / NOTES:
%
% • COMPLEX-STEP SAFE, and must stay so: no abs, no norm, no max/min, no
%   branch on a state value. The STM's Jacobian is taken by complex step
%   through this function.
%
%% Inputs:
%
%  x                        [4 x 1]                 [q1; q2; q1dot; q2dot]:
%                                                   cart position (m),
%                                                   pendulum angle (rad, 0 =
%                                                   hanging), and their rates
%
%  p                        struct                  .m1 cart mass (kg), .m2
%                                                   bob mass (kg), .L
%                                                   pendulum length (m), .g
%                                                   gravity (m/s^2)
%
%% Outputs:
%
%  F                        [4 x 1]                 Drift: xdot at u = 0
%
%  G                        [4 x 1]                 Control column: d(xdot)/du
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the drift and control column at the hanging equilibrium:
     pd = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
     [F, G] = cartpole_field([0; 0; 0; 0], pd);
     fprintf('at rest, hanging:  F = [%s],  G = [%s]\n', ...
             sprintf('%.4f ', F), sprintf('%.4f ', G));
     return
end

q1d = x(3);  q2d = x(4);
s = sin(x(2));  c = cos(x(2));
D1 = p.m1 + p.m2*(1 - c^2);
D2 = p.L*(p.m1 + p.m2)*(1 - (p.m2/(p.m1 + p.m2))*c^2);

F = [q1d;
     q2d;
     (p.L*p.m2*s*q2d^2 + p.m2*p.g*c*s)/D1;
     (p.L*p.m2*c*s*q2d^2 + (p.m1 + p.m2)*p.g*s)/D2];
G = [0; 0; 1/D1; c/D2];
end
