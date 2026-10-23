function [U, gU, HU] = barrier_terms(x, pairs, geom)
%% Purpose:
%
%   The summed barrier of one agent configuration, with its gradient and
%   sparse Hessian: U(x) = sum_pairs U(h_ab) + sum_i U(h_ball,i). For one
%   pair with rho = |x_a - x_b|, n = (x_a - x_b)/rho, h = rho - 2r:
%
%       grad_a =  U'(h) n,   grad_b = -U'(h) n
%       H_aa   =  U''(h) n n' + U'(h) (I - n n')/rho = H_bb,   H_ab = H_ba = -H_aa
%
%   and the ball term has the a-block only (c fixed). This is the
%   objective-independent penalty SHAPE on geometry; the objective that
%   integrates it, and the PMP field A = [0 -I; -HU 0] that reads it, live
%   in the routes. Complex-step safe (branch on real(h) inside
%   barrier_shape; no magnitudes).
%
%% Inputs:
%
%  x                        [3M x 1]                One configuration,
%                                                   agent-major
%
%  pairs                    [2 x P]                 From pair_list
%
%  geom                     struct                  .r .eps .l .ball
%
%% Outputs:
%
%  U                        double                  The summed barrier
%
%  gU                       [3M x 1]                dU/dx
%
%  HU                       [3M x 3M] sparse        d2U/dx2
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the M = 4 circle at t = 0 (gaps 0.67 and 0.96):
     [a, ~, g] = drone_layouts('paper', 4);
     [U, gU] = barrier_terms(a(:), pair_list(4), g);
     fprintf('U = %.3e, |gU| = %.3e\n', U, sqrt(sum(gU.^2)));
     return
end

M = numel(x)/3;
P = size(pairs, 2);
[h, rho, n, hB] = clearance_terms(x(:), pairs, geom);
[Uq, U1, U2] = barrier_shape(h, geom.eps, geom.l);
U  = sum(Uq);
gU = zeros(3*M, 1, 'like', x);
needH = nargout > 2;
if needH
    rows = zeros(36*P + 9*M, 1);  cols = rows;  vals = zeros(36*P + 9*M, 1, 'like', x);
    nz = 0;
end
E3 = eye(3);
for q = 1:P
    a = pairs(1, q);  b = pairs(2, q);
    ia = 3*(a-1) + (1:3);  ib = 3*(b-1) + (1:3);
    nq = n(:, q);
    gq = U1(q)*nq;
    gU(ia) = gU(ia) + gq;
    gU(ib) = gU(ib) - gq;
    if needH
        Hq = U2(q)*(nq*nq.') + (U1(q)/rho(q))*(E3 - nq*nq.');
        [rows, cols, vals, nz] = addBlock(rows, cols, vals, nz, ia, ia,  Hq);
        [rows, cols, vals, nz] = addBlock(rows, cols, vals, nz, ib, ib,  Hq);
        [rows, cols, vals, nz] = addBlock(rows, cols, vals, nz, ia, ib, -Hq);
        [rows, cols, vals, nz] = addBlock(rows, cols, vals, nz, ib, ia, -Hq);
    end
end
if ~isempty(hB)
    [Ub, U1b, U2b] = barrier_shape(hB, geom.eps, geom.l);
    U = U + sum(Ub);
    for a = 1:M
        ia = 3*(a-1) + (1:3);
        dv = x(ia) - geom.ball.c;
        rb = hB(a) + geom.ball.R;
        nb = dv/rb;
        gU(ia) = gU(ia) + U1b(a)*nb;
        if needH
            Hb = U2b(a)*(nb*nb.') + (U1b(a)/rb)*(E3 - nb*nb.');
            [rows, cols, vals, nz] = addBlock(rows, cols, vals, nz, ia, ia, Hb);
        end
    end
end
if needH
    HU = sparse(rows(1:nz), cols(1:nz), vals(1:nz), 3*M, 3*M);
end
end

% ------------------------------------------------------------------------
function [rows, cols, vals, nz] = addBlock(rows, cols, vals, nz, ia, ib, Hq)
%% Purpose:
%
%   Append one 3x3 block to the sparse triplet list.
%
[C, R] = meshgrid(ib, ia);
idx = nz + (1:9);
rows(idx) = R(:);  cols(idx) = C(:);  vals(idx) = Hq(:);
nz = nz + 9;
end
