function [F, A, B] = lt_dynamics_throttle(X, U, Tmax, c, muStar)
% LT_DYNAMICS_THROTTLE  Low-thrust CR3BP dynamics with variable throttle.
%
% The general 4-control form used by the MIN-FUEL transcription: control
% u = [w(3); s] with thrust-direction vector w and throttle s coupled by
% the transcription's cone constraint w'*w = s^2 (so ||w|| = s at any
% feasible point). Thrust acceleration = Tmax*w/m; mass flow = -Tmax*s/c.
% The cone must be an EQUALITY: with an inequality the optimizer can set
% s > ||w|| and burn propellant without thrusting (the "ballast exploit" --
% profitable whenever lighter-is-better, which min-fuel makes explicit).
%
% INPUTS:
%   X      - states [7xM]: [r(3); v(3); m] per column
%   U      - controls [4xM]: [w(3); s] per column
%   Tmax   - max thrust acceleration at m = 1 (ND) [scalar]
%   c      - exhaust velocity (ND) [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   F      - state derivatives [7xM]
%   A      - (optional) df/dx per node [7x7xM]
%   B      - (optional) df/du per node [7x4xM]
%
% REFERENCES:
%   [1] Zhang et al., JGCD 38(8), 2015.
%   [2] Betts, SIAM 2010, Ch. 4.

M  = size(X, 2);
r  = X(1:3, :);
v  = X(4:6, :);
m  = X(7, :);
w  = U(1:3, :);
s  = U(4, :);

dd = [r(1,:) + muStar;     r(2,:); r(3,:)];
rr = [r(1,:) - 1 + muStar; r(2,:); r(3,:)];
d1 = sqrt(sum(dd.^2, 1));
r1 = sqrt(sum(rr.^2, 1));
d3 = d1.^3;  r3 = r1.^3;

gr = [r(1,:); r(2,:); zeros(1,M)] ...
     - (1 - muStar).*dd./d3 - muStar.*rr./r3;
hv = [2*v(2,:); -2*v(1,:); zeros(1,M)];

F = [v;
     gr + hv + Tmax.*w./m;
     -(Tmax/c).*s];

if nargout > 1
    d5 = d1.^5;  r5 = r1.^5;
    ddOuter = reshape(dd, 3, 1, M).*reshape(dd, 1, 3, M);
    rrOuter = reshape(rr, 3, 1, M).*reshape(rr, 1, 3, M);
    I3 = eye(3);
    G  = diag([1, 1, 0]) ...
         - (1 - muStar)*(I3./reshape(d3, 1, 1, M) - 3*ddOuter./reshape(d5, 1, 1, M)) ...
         -      muStar *(I3./reshape(r3, 1, 1, M) - 3*rrOuter./reshape(r5, 1, 1, M));

    A = zeros(7, 7, M);
    A(1:3, 4:6, :) = repmat(I3, 1, 1, M);
    A(4:6, 1:3, :) = G;
    A(4:6, 4:6, :) = repmat([0 2 0; -2 0 0; 0 0 0], 1, 1, M);
    A(4:6, 7,   :) = reshape(-Tmax.*w./m.^2, 3, 1, M);

    B = zeros(7, 4, M);
    B(4:6, 1:3, :) = I3.*reshape(Tmax./m, 1, 1, M);
    B(7,   4,   :) = -Tmax/c;
end
end
