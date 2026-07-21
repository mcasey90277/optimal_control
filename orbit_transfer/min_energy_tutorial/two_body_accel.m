function g = two_body_accel(r_vec, mu)
% TWO_BODY_ACCEL  Planar two-body gravitational acceleration.
%
%   g(r) = -mu * r_vec / ||r_vec||^3
%
% INPUTS:
%   r_vec - Position vector [2x1]
%   mu    - Gravitational parameter [scalar]
%
% OUTPUTS:
%   g     - Gravitational acceleration [2x1]
%
% REFERENCES:
%   [1] Bate, Mueller, White, "Fundamentals of Astrodynamics", Ch. 1.

r = norm(r_vec);
g = -mu * r_vec / r^3;
end
