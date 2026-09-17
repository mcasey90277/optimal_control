function [dy, u] = cartpole_pmp_rhs(y, p)
%% Purpose:
%
%   The Pontryagin field of the cart-pole minimum-effort swing-up at fixed
%   final time: the 8-state flow y = [x; lam] whose boundary-value problem
%   IS the optimality condition.
%
%   With J = integral u^2 dt and xdot = F(x) + G(x) u, the Hamiltonian is
%       H = u^2 + lam'(F + G u),
%   so dH/du = 2u + lam'G = 0 gives the control in closed form,
%       u* = -lam'G / 2,
%   and d2H/du2 = 2 > 0 makes it the minimiser everywhere -- there is no
%   switching structure, which is exactly why this problem is the clean
%   demonstration that the shooting engine is generic. The costate obeys
%   lamdot = -dH/dx at u*; the running cost has no state dependence, so that
%   is -A(x,u*)' lam with A the state Jacobian.
%
%  ASSUMPTIONS / NOTES:
%
% • Complex-step safe: A comes from the GENERATED cartpole_state_jac, never
%   from a complex step, so a caller may complex-step through this function
%   (that is how the STM's Jacobian is built).
% • u* is substituted before the costate equation is evaluated. Doing so is
%   legitimate precisely because dH/du = 0 there: the total and partial
%   x-derivatives of H agree at the stationary control.
%
%% Inputs:
%
%  y                        [8 x 1]                 [x; lam]: the four
%                                                   states then the four
%                                                   costates
%
%  p                        struct                  .m1 .m2 .L .g
%
%% Outputs:
%
%  dy                       [8 x 1]                 [xdot; lamdot]
%
%  u                        double                  The PMP control at y
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the field and its control at a representative point:
     here = fileparts(mfilename('fullpath'));
     addpath(here, fullfile(fileparts(here), 'cartpole_common'));
     pd = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
     [dy, u] = cartpole_pmp_rhs([0; 0.5; 0; 0; 1; -2; 0.5; 0.3], pd);
     fprintf('u* = %.6f,  dy = [%s]\n', u, sprintf('%.4f ', dy));
     return
end

x = y(1:4);  lam = y(5:8);
[F, G] = cartpole_field(x, p);
u = -(lam.'*G)/2;
A = cartpole_state_jac(x, u, p);
dy = [F + G*u; -A.'*lam];
end
