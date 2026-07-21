function G = gravity_gradient(r_vec, mu)
% GRAVITY_GRADIENT  Jacobian of the planar two-body acceleration.
%
%   G = d g / d r_vec = (mu / r^5) * (3 * r_vec * r_vec' - r^2 * I),  2x2,
%   symmetric. This is the matrix that drives both the state variational
%   equation (delta-r'' = G delta-r) and Lawden's primer equation (p'' = G p).
%
% INPUTS:
%   r_vec - Position vector [2x1]
%   mu    - Gravitational parameter [scalar]
%
% OUTPUTS:
%   G     - Gravity-gradient matrix [2x2], symmetric
%
% REFERENCES:
%   [1] Lawden, "Optimal Trajectories for Space Navigation", 1963.
%   [2] myLatex/notes/cov_pmp_orbit_transfer.tex (companion theory note).

r = norm(r_vec);
G = (mu / r^5) * (3 * (r_vec * r_vec') - r^2 * eye(2));
end
