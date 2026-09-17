function p = cartpole_params()
%% Purpose:
%
%   THE cart-pole's physical constants, in one place. Every example that
%   solves this plant -- minimum energy, minimum time, minimum fuel -- and
%   every test of it reads them from here rather than repeating literals.
%
%   The same rule as the dynamics themselves: one home. A sign error survived
%   in this plant until 2026-09-17 because the equations existed in several
%   copies and every test compared one copy against another.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  p                        struct                  .m1 cart mass (kg), .m2
%                                                   bob mass (kg), .L
%                                                   pendulum length (m), .g
%                                                   gravity (m/s^2)
%
%% Revision History:
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0 && nargout == 0
   %Demo: the constants, and the hanging small-oscillation period they imply:
     pd = cartpole_params();
     w  = sqrt((pd.m1 + pd.m2)*pd.g/(pd.L*pd.m1));
     fprintf('m1 %g kg, m2 %g kg, L %g m, g %g m/s^2  ->  hanging period %.4f s\n', ...
             pd.m1, pd.m2, pd.L, pd.g, 2*pi/w);
     return
end

p = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
end
