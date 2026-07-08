function [F, A, B] = lt_dynamics(X, W, Tmax, c, muStar)
% LT_DYNAMICS  Max-thrust low-thrust CR3BP dynamics with analytic Jacobians.
%
% State: x = [r(3); v(3); m] in the rotating nondimensional barycentric
% frame with mass fraction m = m/m0. Control: w = thrust direction vector,
% constrained to the unit sphere (w'*w = 1) by the transcription. The
% throttle is FIXED at 1 (always burn): for this min-time transfer the
% indirect solution's switching function S = -||lambda_v||c/m - lambda_m
% stays strictly negative (max -2.5), so full thrust is optimal throughout
% and carrying a throttle variable only creates active-bound trouble for
% interior-point methods. Thrust acceleration = Tmax*w/m; mass flow =
% -Tmax/c (constant).
%
% INPUTS:
%   X      - states [7xM] (M nodes evaluated at once)
%   W      - thrust direction controls [3xM]
%   Tmax   - max thrust acceleration at m = 1 (ND) [scalar]
%   c      - exhaust velocity (ND) [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   F      - state derivatives [7xM]
%   A      - (optional) df/dx per node [7x7xM]
%   B      - (optional) df/dw per node [7x3xM]
%
% REFERENCES:
%   [1] Zhang et al., "Low-Thrust Minimum-Fuel Optimization in the
%       Circular Restricted Three-Body Problem," JGCD 38(8), 2015.

M  = size(X, 2);
r  = X(1:3, :);
v  = X(4:6, :);
m  = X(7, :);

% Distances to Earth (at [-muStar,0,0]) and Moon (at [1-muStar,0,0])
dd = [r(1,:) + muStar;     r(2,:); r(3,:)];
rr = [r(1,:) - 1 + muStar; r(2,:); r(3,:)];
d1 = sqrt(sum(dd.^2, 1));
r1 = sqrt(sum(rr.^2, 1));
d3 = d1.^3;  r3 = r1.^3;

gr = [r(1,:); r(2,:); zeros(1,M)] ...
     - (1 - muStar).*dd./d3 - muStar.*rr./r3;
hv = [2*v(2,:); -2*v(1,:); zeros(1,M)];

F = [v;
     gr + hv + Tmax.*W./m;
     (-Tmax/c)*ones(1, M)];

if nargout > 1
    % Vectorized Jacobian assembly (implicit expansion over the node
    % dimension; per-node loops are measurably slower at M ~ 1e4).
    d5 = d1.^5;  r5 = r1.^5;

    % Gravity-gradient + centrifugal block G, all nodes: 3x3xM
    ddOuter = reshape(dd, 3, 1, M).*reshape(dd, 1, 3, M);
    rrOuter = reshape(rr, 3, 1, M).*reshape(rr, 1, 3, M);
    I3  = eye(3);
    G   = diag([1, 1, 0]) ...
          - (1 - muStar)*(I3./reshape(d3, 1, 1, M) - 3*ddOuter./reshape(d5, 1, 1, M)) ...
          -      muStar *(I3./reshape(r3, 1, 1, M) - 3*rrOuter./reshape(r5, 1, 1, M));

    A = zeros(7, 7, M);
    A(1:3, 4:6, :) = repmat(I3, 1, 1, M);
    A(4:6, 1:3, :) = G;
    A(4:6, 4:6, :) = repmat([0 2 0; -2 0 0; 0 0 0], 1, 1, M);
    A(4:6, 7,   :) = reshape(-Tmax.*W./m.^2, 3, 1, M);

    B = zeros(7, 3, M);
    B(4:6, 1:3, :) = I3.*reshape(Tmax./m, 1, 1, M);
end
end
