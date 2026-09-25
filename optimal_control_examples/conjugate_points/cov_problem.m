function prob = cov_problem(Fstr, a, b, ya, yb)
%% Purpose:
%
%   Turn a calculus-of-variations problem typed as text into the functions
%   every other file here needs:
%
%       minimise  J[y] = int_a^b F(t, y, y') dt,   y(a) = ya,  y(b) = yb.
%
%   The Lagrangian is written in the variables t, y and yp (yp = y'). The
%   Euler-Lagrange equation d/dt F_yp = F_y is expanded and solved for y'':
%
%       y'' = g(t,y,y') = (F_y - F_ypt - F_ypy y') / F_ypyp,
%
%   which needs F_ypyp ~= 0 (the Legendre coefficient). Its linearisation
%   about an extremal is the Jacobi equation, written here in first-order
%   form:  h'' = g_y h + g_yp h'.  The second variation uses the three
%   coefficients P = F_ypyp, R = F_ypy, Q0 = F_yy:
%
%       d^2 J[eta] = int_a^b ( P eta'^2 + 2 R eta eta' + Q0 eta^2 ) dt.
%
%   All derivatives are symbolic (Symbolic Math Toolbox), then turned into
%   vectorised function handles; nothing downstream differentiates.
%
%% Inputs:
%
%  Fstr                     char/string             Lagrangian in t, y, yp,
%                                                   e.g. 'yp^2 - y^2'
%  a, b                     scalar                  interval, a < b
%  ya, yb                   scalar                  boundary values
%
%% Outputs:
%
%  prob                     struct                  .F .g .gy .gyp .P .R .Q0
%                                                   (handles of (t,y,yp),
%                                                   elementwise) .a .b .ya
%                                                   .yb .Fstr .gStr (the
%                                                   Euler-Lagrange right-hand
%                                                   side as text)
%
%% References:
%
%   [1] I. M. Gelfand and S. V. Fomin, "Calculus of Variations,"
%       Prentice-Hall, 1963, Ch. 5 (the second variation, Jacobi's
%       equation, conjugate points).
%   [2] D. Liberzon, "Calculus of Variations and Optimal Control Theory,"
%       Princeton University Press, 2012, Sec. 2.6.
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    %Demo: the harmonic oscillator, whose conjugate point sits at pi
    prob = cov_problem('yp^2 - y^2', 0, 4, 0, 1);
    fprintf('y'''' = %s\n', prob.gStr);
    return
end
assert(isscalar(a) && isscalar(b) && isfinite(a) && isfinite(b) && b > a, ...
       'cov_problem:interval', 'need finite a < b');
assert(isscalar(ya) && isscalar(yb) && isfinite(ya) && isfinite(yb), ...
       'cov_problem:bc', 'need finite boundary values');

syms t y yp real
try
    F = str2sym(char(Fstr));
catch err
    error('cov_problem:parse', 'could not parse F = "%s": %s', char(Fstr), err.message);
end
extra = setdiff(symvar(F), [t y yp]);
if ~isempty(extra)
    error('cov_problem:vars', ['F may use only t, y and yp; found: %s'], ...
          strjoin(string(extra), ', '));
end

Fy    = diff(F, y);
Fyp   = diff(F, yp);
Fypyp = diff(Fyp, yp);
Fypy  = diff(Fyp, y);
Fypt  = diff(Fyp, t);
Fyy   = diff(Fy, y);
if isAlways(Fypyp == 0, 'Unknown', 'false')
    error('cov_problem:degenerate', ...
          'F_ypyp is identically zero: the Euler-Lagrange equation is not second order');
end

g   = (Fy - Fypt - Fypy*yp)/Fypyp;
gy  = diff(g, y);
gyp = diff(g, yp);

vars = {t, y, yp};
prob = struct();
prob.F   = vec(matlabFunction(F,     'Vars', vars));
prob.g   = vec(matlabFunction(g,     'Vars', vars));
prob.gy  = vec(matlabFunction(gy,    'Vars', vars));
prob.gyp = vec(matlabFunction(gyp,   'Vars', vars));
prob.P   = vec(matlabFunction(Fypyp, 'Vars', vars));
prob.R   = vec(matlabFunction(Fypy,  'Vars', vars));
prob.Q0  = vec(matlabFunction(Fyy,   'Vars', vars));
prob.a = a;  prob.b = b;  prob.ya = ya;  prob.yb = yb;
prob.Fstr = char(Fstr);
prob.gStr = char(simplify(g, 'Steps', 20));
end

% ---------------------------------------------------------------------------
function h = vec(f)
% VEC  Make a generated handle return an array the size of its inputs even
% when the expression is constant (matlabFunction returns a scalar then).
% INPUTS: f handle of (t,y,yp).  OUTPUTS: h handle of (t,y,yp).
h = @(t, y, yp) f(t, y, yp) + zeros(size(t + y + yp));
end
