function V = cov_second_variation(prob, ext, nEl, epsList)
%% Purpose:
%
%   The second variation of J along an extremal, as a quadratic form on
%   variations eta with eta(a) = eta(b) = 0:
%
%       d^2 J[eta] = int_a^b ( P eta'^2 + 2 R eta eta' + Q0 eta^2 ) dt,
%       P = F_ypyp,  R = F_ypy,  Q0 = F_yy   (evaluated on the extremal),
%
%   discretised with piecewise-linear finite elements into K (stiffness)
%   and M (the L2 mass matrix). The generalised eigenvalues K v = lambda M v
%   are an INDEPENDENT instrument for the conjugate point: with P > 0,
%
%       number of negative eigenvalues = number of conjugate points in (a,b)
%
%   (the Morse index theorem). A negative eigenvalue also HANDS YOU a
%   direction that lowers J: the eigenvector eta*. This file evaluates
%   Delta J(eps) = J[y0 + eps eta*] - J[y0] exactly (same quadrature as
%   cov_functional) for each eps, next to the quadratic prediction
%   (eps^2/2) eta*' K eta*.
%
%% Inputs:
%
%  prob                     struct                  cov_problem output
%  ext                      struct                  one cov_extremals entry
%  nEl                      scalar                  elements [default 400]
%  epsList                  [1xk]                   amplitudes for Delta J
%                                                   [default 0.1*[-1 -0.5
%                                                   0.5 1] * max(1,|y0|max)]
%
%% Outputs:
%
%  V                        struct                  .lambda (sorted, [nEl-1])
%                                                   .nNeg (count below
%                                                   -tolNeg) .tolNeg .t
%                                                   (nodes) .eta (lowest
%                                                   mode, max|eta| = 1)
%                                                   .eps .dJ (exact) .dJquad
%                                                   (eps^2/2 eta'K eta) .Pmin
%                                                   (min of P on the extremal:
%                                                   the Legendre check)
%
%% References:
%
%   [1] I. M. Gelfand and S. V. Fomin, "Calculus of Variations,"
%       Prentice-Hall, 1963, Ch. 5.
%   [2] J. Milnor, "Morse Theory," Princeton University Press, 1963,
%       Part III (the index theorem).
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3 || isempty(nEl), nEl = 400; end
a = prob.a;  b = prob.b;
sol = ext.S.sol;

% 4-point Gauss on each element: points, weights, local coordinate
x = [-0.861136311594053; -0.339981043584856; 0.339981043584856; 0.861136311594053];
w = [ 0.347854845137454;  0.652145154862546; 0.652145154862546; 0.347854845137454];
he = (b - a)/nEl;
tn = linspace(a, b, nEl + 1).';                % exact end nodes
xi = (x + 1)/2;
nG = numel(x);
tg = reshape(tn(1:nEl).' + he*xi, [], 1);       % element-major, nG per element
z = deval(sol, tg);
Pg = prob.P(tg, z(1,:).', z(2,:).');
Rg = prob.R(tg, z(1,:).', z(2,:).');
Qg = prob.Q0(tg, z(1,:).', z(2,:).');

% element matrices, assembled (basis 1-xi and xi; derivatives -1/he, 1/he)
phi  = [1 - xi, xi];                           % nG x 2
dphi = [-1, 1]/he;                             % 1 x 2
I = zeros(4*nEl, 1);  Jc = I;  Kv = I;  Mv = I;
n = 0;
for e = 1:nEl
    gi = (e-1)*nG + (1:nG);
    wp = w*he/2;
    Ke = zeros(2);  Me = zeros(2);
    for i = 1:2
        for j = 1:2
            Ke(i,j) = sum(wp.*( Pg(gi)*dphi(i)*dphi(j) ...
                              + Rg(gi).*(phi(:,i)*dphi(j) + phi(:,j)*dphi(i)) ...
                              + Qg(gi).*phi(:,i).*phi(:,j) ));
            Me(i,j) = sum(wp.*phi(:,i).*phi(:,j));
        end
    end
    idx = [e, e+1];
    for i = 1:2
        for j = 1:2
            n = n + 1;
            I(n) = idx(i);  Jc(n) = idx(j);  Kv(n) = Ke(i,j);  Mv(n) = Me(i,j);
        end
    end
end
K = sparse(I, Jc, Kv, nEl+1, nEl+1);
M = sparse(I, Jc, Mv, nEl+1, nEl+1);
in = 2:nEl;                                    % Dirichlet: drop the end nodes
K = full(K(in, in));  M = full(M(in, in));
K = (K + K.')/2;  M = (M + M.')/2;
[Vec, D] = eig(K, M);
[lambda, ord] = sort(real(diag(D)));
Vec = Vec(:, ord);
tolNeg = 1e-8*max(1, max(abs(lambda)));

eta = [0; Vec(:,1); 0];
[~, im] = max(abs(eta));
eta = eta/eta(im);                             % max |eta| = 1, positive there

if nargin < 4 || isempty(epsList)
    ys = deval(sol, tn);
    epsList = 0.1*[-1 -0.5 0.5 1]*max(1, max(abs(ys(1,:))));
end
J0 = cov_functional(prob, sol, zeros(nEl+1, 1));
q = eta(in).'*K*eta(in);
dJ = zeros(size(epsList));
for k = 1:numel(epsList)
    dJ(k) = cov_functional(prob, sol, epsList(k)*eta) - J0;
end

V = struct('lambda', lambda, 'nNeg', nnz(lambda < -tolNeg), 'tolNeg', tolNeg, ...
           't', tn, 'eta', eta, 'eps', epsList, 'dJ', dJ, ...
           'dJquad', 0.5*epsList.^2*q, 'Pmin', min(Pg));
end
