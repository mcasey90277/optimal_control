function [lam0, sol, J, psi_norm] = solve_indirect(x0, xf, tf, mu, lam0_guess)
% SOLVE_INDIRECT  Min-energy transfer by indirect shooting on the costate BVP.
%
%   Finds the initial costate lam0 such that the coupled state-costate flow
%   hits the rendezvous target: psi(lam0) = x(tf) - xf = 0 (fsolve). Then
%   evaluates the realized cost J = int 0.5*||u||^2 dt with u = -lam_v on a
%   dense time grid.
%
% INPUTS:
%   x0         - Initial state [4x1]
%   xf         - Target state [4x1]
%   tf         - Fixed final time [scalar]
%   mu         - Gravitational parameter [scalar]
%   lam0_guess - Initial costate guess [4x1] (zeros(4,1) works for the
%                near-Hohmann tutorial instance)
%
% OUTPUTS:
%   lam0     - Converged initial costate [4x1]
%   sol      - ode45 solution structure for the converged trajectory
%   J        - Min-energy cost 0.5 * int ||u||^2 dt [scalar]
%   psi_norm - Terminal miss norm ||x(tf) - xf|| at the solution [scalar]
%
% REFERENCES:
%   [1] Bryson & Ho, "Applied Optimal Control", Ch. 7.

opts = optimoptions('fsolve', 'Display', 'off', ...
                    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, ...
                    'OptimalityTolerance', 1e-10);
lam0 = fsolve(@(l) shoot_residual(l, x0, xf, tf, mu), lam0_guess, opts);

[psi, sol] = shoot_residual(lam0, x0, xf, tf, mu);
psi_norm   = norm(psi);

% Cost on a dense grid: u(t) = -lam_v(t)
tq = linspace(0, tf, 2001);
zq = deval(sol, tq);
uq = -zq(7:8, :);
J  = trapz(tq, 0.5 * sum(uq.^2, 1));
end
