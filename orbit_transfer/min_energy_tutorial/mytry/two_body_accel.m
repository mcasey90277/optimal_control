% function g = two_body_accel(r_vec, mu)
function g = two_body_accel(r_vec, mu)
% PURPOSE: calculate gravitational acceleration from a central body on a satellite
% INPUTS:
%   r_vec: 2x1 position vector
%   mu: central bodies gravitational parameter = GM
% OUTPUTS:
% g: 2x1 gravitational acceleration vector components

% radius 
r = norm(r_vec,2);

% gravitaional acceleration
g = -mu .* r_vec/r^3;

