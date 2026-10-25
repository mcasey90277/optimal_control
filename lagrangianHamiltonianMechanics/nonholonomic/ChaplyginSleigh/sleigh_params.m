function p = sleigh_params(m, J, a)
%% Purpose:
%
%   The Chaplygin sleigh's constants and the two derived ones every other
%   file uses: Jb = J + m a^2 (inertia about the blade point) and k = m a
%   / Jb (the decay rate in wd = -k v w). a is SIGNED: a > 0 puts the
%   center of mass ahead of the blade in the forward (+e1) direction.
%
%% Inputs:
%
%  m                        [scalar]                Mass (kg), > 0 [1]
%
%  J                        [scalar]                Inertia about the
%                                                   center of mass (kg
%                                                   m^2), > 0 [0.1]
%
%  a                        [scalar]                Blade to center of
%                                                   mass (m), signed, finite
%                                                   [0.2]
%
%% Outputs:
%
%  p                        struct                  .m .J .a .Jb .k
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the defaults
     p = sleigh_params(1, 0.1, 0.2);
     disp(p);
     return
end
if ~(isscalar(m) && isfinite(m) && m > 0 && isscalar(J) && isfinite(J) && J > 0 && isscalar(a) && isfinite(a))
    error('sleigh_params:value', 'need scalar m > 0, J > 0 and finite a (got m = %g, J = %g, a = %g)', m, J, a);
end
p.m  = m;
p.J  = J;
p.a  = a;
p.Jb = J + m*a^2;
p.k  = m*a/p.Jb;
end
