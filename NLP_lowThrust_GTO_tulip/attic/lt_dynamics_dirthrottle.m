function [F, A, B] = lt_dynamics_dirthrottle(X, U, Tmax, c, muStar)
% LT_DYNAMICS_DIRTHROTTLE  Low-thrust CR3BP dynamics, direction+throttle form.
%
% Cone-ELIMINATED control parameterization for min-fuel: control
% u = [alpha(3); s] with alpha a UNIT thrust direction (enforced by a
% separate ||alpha|| = 1 constraint in the NLP) and s the throttle in [0,1].
% Thrust acceleration = s*Tmax*alpha/m; mass flow = -s*Tmax/c.
%
% Why this instead of the (w, s) cone form w'*w = s^2: there, a coast
% (s -> 0) forces w -> 0, which makes the thrust DIRECTION degenerate and
% wedges the solver at bang-bang. Here alpha stays a well-defined unit
% vector at a coast while s -> 0 -- the direction and throttle are decoupled.
%
% INPUTS:
%   X      - states [7xM]: [r(3); v(3); m] per column
%   U      - controls [4xM]: [alpha(3); s] per column
%   Tmax   - max thrust acceleration at m = 1 (ND) [scalar]
%   c      - exhaust velocity (ND) [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   F - state derivatives [7xM]
%   A - (optional) df/dx per node [7x7xM]
%   B - (optional) df/du per node [7x4xM]
%
% REFERENCES:
%   [1] Betts, SIAM 2010, Ch. 4.

M = size(X, 2);
r = X(1:3, :);  v = X(4:6, :);  m = X(7, :);
alpha = U(1:3, :);  s = U(4, :);

dd = [r(1,:) + muStar;     r(2,:); r(3,:)];
rr = [r(1,:) - 1 + muStar; r(2,:); r(3,:)];
d1 = sqrt(sum(dd.^2, 1));  r1 = sqrt(sum(rr.^2, 1));
d3 = d1.^3;  r3 = r1.^3;

gr = [r(1,:); r(2,:); zeros(1,M)] - (1 - muStar).*dd./d3 - muStar.*rr./r3;
hv = [2*v(2,:); -2*v(1,:); zeros(1,M)];

sm = s.*Tmax./m;                       % 1xM throttle*Tmax/m
F = [v;
     gr + hv + sm.*alpha;
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
    A(4:6, 7,   :) = reshape(-s.*Tmax./m.^2.*alpha, 3, 1, M);   % d(sm*alpha)/dm

    B = zeros(7, 4, M);
    B(4:6, 1:3, :) = I3.*reshape(sm, 1, 1, M);                  % d(sm*alpha)/dalpha
    B(4:6, 4,   :) = reshape((Tmax./m).*alpha, 3, 1, M);        % d(sm*alpha)/ds
    B(7,   4,   :) = -Tmax/c;
end
end
