function [U, U1, U2] = barrier_shape(h, epsB, l)
%% Purpose:
%
%   The C1 log barrier of Zhang, Liu, Darbon & Karniadakis (arXiv:
%   2610.30388, eq. 3b), elementwise in h, with its first two derivatives:
%
%       U  = -eps log h                                    h > l
%          = -eps log l + (eps/2) (((h - 2l)/l)^2 - 1)     h <= l
%       U' = -eps/h                | (eps/l^2)(h - 2l)
%       U''=  eps/h^2              |  eps/l^2
%
%   Both branches meet at h = l in U, U' AND U''. The quadratic branch keeps
%   everything finite through h = 0 (interpenetration), so an optimizer or
%   an integrator never evaluates log of a non-positive gap.
%
%  ASSUMPTIONS / NOTES:
%
% • Complex-step safe: the branch is on real(h). Route 1's CasADi barrier
%   expression is a SECOND COPY of these two formulas (CasADi needs
%   if_else); its test pins that copy to barrier_terms.
%
%% Inputs:
%
%  h                        [any]                   Gaps (may be complex)
%
%  epsB                     double                  Barrier weight eps
%
%  l                        double                  Switch gap
%
%% Outputs:
%
%  U, U1, U2                [size(h)]               U, dU/dh, d2U/dh2
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the shape across the switch:
     hh = linspace(-1e-4, 1e-3, 400);
     figure('color', [1 1 1]);  plot(hh, barrier_shape(hh, 1e-5, 1e-4));  grid on
     xlabel('h');  ylabel('U_{\epsilon,l}(h)');
     return
end

logBranch = real(h) > l;
U  = zeros(size(h), 'like', h);
U1 = zeros(size(h), 'like', h);
U2 = zeros(size(h), 'like', h);
hl = h(logBranch);
U(logBranch)  = -epsB*log(hl);
U1(logBranch) = -epsB./hl;
U2(logBranch) =  epsB./hl.^2;
hq = h(~logBranch);
U(~logBranch)  = -epsB*log(l) + (epsB/2)*(((hq - 2*l)/l).^2 - 1);
U1(~logBranch) = (epsB/l^2)*(hq - 2*l);
U2(~logBranch) = (epsB/l^2)*ones(size(hq), 'like', hq);
end
