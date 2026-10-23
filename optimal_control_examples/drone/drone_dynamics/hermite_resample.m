function [Xf, Vf, tf] = hermite_resample(X, V, tN, nSub)
%% Purpose:
%
%   Cubic-Hermite reconstruction of a trajectory from node values X and
%   node derivatives V, sampled nSub times per interval. Under the
%   trapezoid rule with xdot = v (the KINEMATIC model) v is linear and x
%   is quadratic, and the cubic Hermite reproduces that solution exactly
%   because the trapezoid defect kills the cubic coefficient
%   (2(x_k - x_{k+1}) + h(v_k + v_{k+1}) = 0). For pointmass/rigid6dof
%   position rows it is a plain cubic Hermite on node values and
%   derivatives, NOT the transcription's own interpolant. Used by the clearance gate to look BETWEEN the nodes, where
%   a nodal constraint says nothing.
%
%% Inputs:
%
%  X                        [nx x N+1]              Node values
%
%  V                        [nx x N+1]              Node derivatives
%
%  tN                       [1 x N+1]               Node times
%
%  nSub                     double                  Sub-intervals per
%                                                   interval
%
%% Outputs:
%
%  Xf                       [nx x N nSub + 1]       Values on the fine grid
%
%  Vf                       [nx x N nSub + 1]       Derivatives there
%
%  tf                       [1 x N nSub + 1]        The fine grid
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a sine from 5 nodes:
     t = linspace(0, 2*pi, 6);
     [xf, ~, tt] = hermite_resample(sin(t), cos(t), t, 20);
     fprintf('max error vs sin on the fine grid: %.2e\n', max(abs(xf - sin(tt))));
     return
end

N  = numel(tN) - 1;
nx = size(X, 1);
nF = N*nSub + 1;
Xf = zeros(nx, nF);  Vf = zeros(nx, nF);  tf = zeros(1, nF);
s  = (0:nSub-1)/nSub;
for k = 1:N
    hk = tN(k+1) - tN(k);
    h00 = 2*s.^3 - 3*s.^2 + 1;   h10 = s.^3 - 2*s.^2 + s;
    h01 = -2*s.^3 + 3*s.^2;      h11 = s.^3 - s.^2;
    d00 = 6*s.^2 - 6*s;          d10 = 3*s.^2 - 4*s + 1;
    d01 = -6*s.^2 + 6*s;         d11 = 3*s.^2 - 2*s;
    cols = (k-1)*nSub + (1:nSub);
    Xf(:, cols) = X(:,k)*h00 + hk*V(:,k)*h10 + X(:,k+1)*h01 + hk*V(:,k+1)*h11;
    Vf(:, cols) = (X(:,k)*d00 + hk*V(:,k)*d10 + X(:,k+1)*d01 + hk*V(:,k+1)*d11)/hk;
    tf(cols) = tN(k) + hk*s;
end
Xf(:, end) = X(:, end);  Vf(:, end) = V(:, end);  tf(end) = tN(end);
end
