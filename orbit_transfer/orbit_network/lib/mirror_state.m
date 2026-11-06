function R = mirror_state(x)
%% Purpose:
%
%   The x-axis (x-z plane) symmetry of the CR3BP applied to a state:
%   (x, y, z, vx, vy, vz) -> (x, -y, z, -vx, vy, -vz). With time reversed
%   the image of a ballistic trajectory is a ballistic trajectory, and the
%   image of a symmetric periodic orbit is the same orbit. For a THRUSTING
%   arc the thrust acceleration is reflected by diag(1, -1, 1) and the
%   mass runs backward (a mass-gain model), so a mirrored transfer is a
%   geometric object only -- see the spec's "hypothetical reverse graph".
%
%% Inputs:
%
%  x                        [6 x n]                 States as columns
%
%% Outputs:
%
%  R                        [6 x n]                 Mirrored states
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a state on the symmetry plane is its own image; an off-plane one is not
     x0 = [0.85; 0; 0; 0; 0.5; 0];  x1 = [0.8; 0.1; 0.02; -0.1; 0.4; 0.01];
     fprintf('on-plane : %s -> %s\n', mat2str(x0.'), mat2str(mirror_state(x0).'));
     fprintf('off-plane: %s -> %s\n', mat2str(x1.'), mat2str(mirror_state(x1).'));
     R = mirror_state(x1);
     return
end

if size(x, 1) ~= 6
    error('orbit_network:mirror_state:shape', 'x must be [6 x n] (got %d x %d)', size(x, 1), size(x, 2));
end
R = [1; -1; 1; -1; 1; -1].*x;
end
