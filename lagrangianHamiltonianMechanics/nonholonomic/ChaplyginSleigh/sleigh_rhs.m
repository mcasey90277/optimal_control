function yd = sleigh_rhs(t, y, p)
%% Purpose:
%
%   The hand-derived reduced equations of the Chaplygin sleigh on the
%   BLADE-point state y = [x_b; y_b; th; v; w]: v is the forward speed of
%   the blade (equal to the forward speed of the center of mass), w the
%   yaw rate. From m rdd_c = lambda e2 and J wd = -a lambda with rd_c = v
%   e1 + a w e2:
%
%       xd_b = v cos th,  yd_b = v sin th,  thd = w,
%       vd = a w^2,       wd = -k v w,      k = m a / (J + m a^2).
%
%   No symbolic call: mech.lagrange_dalembert is this file's independent
%   check (test_sleigh_oracle D1). Vectorized over columns of y.
%
%% Inputs:
%
%  t                        [scalar or 1 x N]       Time (unused;
%                                                   autonomous)
%
%  y                        [5 x N]                 [x_b; y_b; th; v; w]
%
%  p                        struct                  From sleigh_params
%
%% Outputs:
%
%  yd                       [5 x N]                 Time derivative
%
%% References:
%   [1] Bloch, "Nonholonomic Mechanics and Control," 2nd ed., Springer,
%       2015, Sec. 1.7 (the Chaplygin sleigh).
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: five seconds from the defaults
     p = sleigh_params(1, 0.1, 0.2);
     [tt, yy] = ode113(@(a, b) sleigh_rhs(a, b, p), [0 5], [0; 0; 0; 1; 2], odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
     fprintf('t = %.1f s: v = %.6f, w = %.6f\n', tt(end), yy(end, 4), yy(end, 5));
     yd = sleigh_rhs(0, yy(end, :).', p);
     return
end
th = y(3, :);  v = y(4, :);  w = y(5, :);
yd = [v.*cos(th); v.*sin(th); w; p.a*w.^2; -p.k*v.*w];
end
