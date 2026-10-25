function out = sleigh_quantities(t, Y, p)
%% Purpose:
%
%   Everything derived from a reduced solution: the blade and center-of-
%   mass paths, the energy E = m/2 v^2 + Jb/2 w^2, the blade's constraint
%   force lambda = m J v w / Jb, the invariant C = v^2 + (Jb/m) w^2, the
%   exact solution in its well-conditioned form (b = sqrt(Jb/m), S =
%   sqrt(C), alpha = asinh(v0/(b |w0|)), z = alpha + k S t: v = S tanh z,
%   w = sign(w0) (S/b) sech z; constants when a = 0 or w0 = 0), the
%   asymptote sign(a) S (v0 when a = 0 or w0 = 0), the exact total heading
%   change theta_inf = theta_0 + sign(w0) (b/|a|) acos(sign(a) v0/S), and
%   the full-state arrays [q; qd] and [qd; qdd] of the center-of-mass
%   description that the symbolic builder is checked against.
%
%% Inputs:
%
%  t                        [1 x N]                 Times (s)
%
%  Y                        [5 x N]                 Reduced solution
%                                                   columns [x_b; y_b;
%                                                   th; v; w]
%
%  p                        struct                  From sleigh_params
%
%% Outputs:
%
%  out                      struct                  .rb .rc [2 x N]; .E
%                                                   .lambda .vExact
%                                                   .wExact [1 x N]; .C;
%                                                   .asymptote; .thetaInf
%                                                   (NaN when a = 0 or
%                                                   w0 = 0); .Yfull [6
%                                                   x N] = [xc; yc; th;
%                                                   xcd; ycd; thd];
%                                                   .Yfulld [6 x N] =
%                                                   [qd; qdd] from
%                                                   sleigh_rhs
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the defaults, five seconds
     p = sleigh_params(1, 0.1, 0.2);
     t = 0:0.01:5;
     [~, yy] = ode113(@(a, b) sleigh_rhs(a, b, p), t, [0; 0; 0; 1; 2], odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
     out = sleigh_quantities(t, yy.', p);
     fprintf('C = %.6f, asymptote %.6f, v(end) %.6f, max |v - vExact| %.1e, energy spread %.1e\n', ...
             out.C, out.asymptote, yy(end, 4), max(abs(yy(:, 4).' - out.vExact)), max(out.E) - min(out.E));
     return
end
t  = t(:).';
th = Y(3, :);  v = Y(4, :);  w = Y(5, :);
e1 = [cos(th); sin(th)];
e2 = [-sin(th); cos(th)];

out.rb     = Y(1:2, :);
out.rc     = out.rb + p.a*e1;
vc         = v.*e1 + p.a*w.*e2;                      % center-of-mass velocity
out.vc     = vc;                                     % [2 x N] center-of-mass velocity
out.E      = p.m/2*v.^2 + p.Jb/2*w.^2;
out.lambda = p.m*p.J*v.*w/p.Jb;
out.C      = v(1)^2 + (p.Jb/p.m)*w(1)^2;

yd = sleigh_rhs(t, Y, p);
vd = yd(4, :);  wd = yd(5, :);
ac = (vd - p.a*w.^2).*e1 + (v.*w + p.a*wd).*e2;     % d/dt (v e1 + a w e2)
out.Yfull  = [out.rc; th; vc; w];
out.Yfulld = [vc; w; ac; wd];

s = t - t(1);
if p.a ~= 0 && w(1) ~= 0
    % the well-conditioned form: alpha = asinh(v0/(b |w0|)) never sees the
    % rounded ratio v0/sqrt(C) = +-1 that the atanh form produces for tiny w0
    b  = sqrt(p.Jb/p.m);
    S  = sqrt(out.C);
    al = asinh(v(1)/(b*abs(w(1))));
    z  = al + p.k*S*s;
    out.vExact    = S*tanh(z);
    out.wExact    = sign(w(1))*(S/b)*sech(z);
    out.asymptote = sign(p.a)*S;
    out.thetaInf  = th(1) + sign(w(1))*(b/abs(p.a))*acos(sign(p.a)*v(1)/S);   % total heading change, exact
else
    out.vExact    = v(1)*ones(size(s));              % knife edge or straight ride: constants
    out.wExact    = w(1)*ones(size(s));
    out.asymptote = v(1);
    out.thetaInf  = NaN;                             % the heading keeps turning (a = 0) or never changes (w0 = 0)
end
end
