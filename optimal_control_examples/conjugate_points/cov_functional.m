function J = cov_functional(prob, sol, eta, nEl)
%% Purpose:
%
%   Evaluate J[y] = int_a^b F(t, y, y') dt for an extremal, or for the
%   extremal plus a piecewise-linear variation:
%
%       y(t) = y0(t) + eta(t),   eta given at nEl+1 equally spaced nodes,
%                                eta(a) = eta(b) = 0.
%
%   Quadrature: 4-point Gauss-Legendre on each of nEl equal elements, with
%   y0 and y0' taken from the ODE solution by deval (exact to the ODE's
%   tolerance) and eta, eta' exact on each element (eta' is constant there).
%   The same rule is used for the extremal alone and for every comparison
%   curve, so their DIFFERENCE is not polluted by quadrature error.
%
%% Inputs:
%
%  prob                     struct                  cov_problem output
%  sol                      struct                  ode45 solution of the
%                                                   extremal (cov_shoot .sol)
%  eta                      [(nEl+1)x1]             nodal variation [default
%                                                   zeros]
%  nEl                      scalar                  elements [default 400, or
%                                                   numel(eta)-1]
%
%% Outputs:
%
%  J                        scalar                  the functional
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3 || isempty(eta)
    if nargin < 4 || isempty(nEl), nEl = 400; end
    eta = zeros(nEl + 1, 1);
else
    nEl = numel(eta) - 1;
end
eta = eta(:);
[tg, wg, e1, xi] = gauss_grid(prob.a, prob.b, nEl);
z = deval(sol, tg);
he = (prob.b - prob.a)/nEl;
etaG  = (1 - xi).*eta(e1) + xi.*eta(e1 + 1);
detaG = (eta(e1 + 1) - eta(e1))/he;
J = sum(wg.*prob.F(tg, z(1,:).' + etaG, z(2,:).' + detaG));
end

% ---------------------------------------------------------------------------
function [tg, wg, e1, xi] = gauss_grid(a, b, nEl)
% GAUSS_GRID  4-point Gauss-Legendre points on nEl equal elements.
% INPUTS: a, b, nEl. OUTPUTS: tg, wg [4nEl x 1] points and weights; e1 the
% left node index of each point's element; xi its local coordinate in [0,1].
x = [-0.861136311594053; -0.339981043584856; 0.339981043584856; 0.861136311594053];
w = [ 0.347854845137454;  0.652145154862546; 0.652145154862546; 0.347854845137454];
he = (b - a)/nEl;
left = a + he*(0:nEl-1);
xi = repmat((x + 1)/2, nEl, 1);
e1 = kron((1:nEl).', ones(4, 1));
tg = left(e1).' + he*xi;
wg = repmat(w*he/2, nEl, 1);
end
