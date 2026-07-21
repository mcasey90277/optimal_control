function [psi, sol] = shoot_residual(lam0, x0, xf, tf, mu)
% SHOOT_RESIDUAL  Terminal-state miss for the indirect (shooting) method.
%
%   Integrates the coupled state-costate system from [x0; lam0] over
%   [0, tf] and returns the rendezvous residual psi = x(tf) - xf. Solving
%   psi(lam0) = 0 for lam0 (4 equations, 4 unknowns) solves the min-energy
%   BVP: x(0) fixed (4 conditions) + x(tf) fixed (4 conditions) = 2n = 8.
%
% INPUTS:
%   lam0 - Initial costate guess [4x1]
%   x0   - Initial state (on departure orbit) [4x1]
%   xf   - Target state (rendezvous point on arrival orbit) [4x1]
%   tf   - Fixed final time [scalar]
%   mu   - Gravitational parameter [scalar]
%
% OUTPUTS:
%   psi  - Terminal miss x(tf) - xf [4x1]
%   sol  - ode45 solution structure (deval-able), for post-processing
%
% REFERENCES:
%   [1] Bryson & Ho, "Applied Optimal Control", Ch. 7 (indirect methods).

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
sol  = ode45(@(t, z) ocp_dynamics(t, z, mu), [0 tf], [x0; lam0], opts);
zf   = deval(sol, tf);
psi  = zf(1:4) - xf;
end
