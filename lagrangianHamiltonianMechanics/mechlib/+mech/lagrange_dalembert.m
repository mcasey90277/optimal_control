function sys = lagrange_dalembert(q, qd, L, A, Q)
%% Purpose:
%
%   Lagrange-d'Alembert equations for a mechanical system whose
%   constraints are linear in velocity, A(q) qd = 0 (nonholonomic, or
%   holonomic written in velocity form), from a SYMBOLIC Lagrangian. With
%   p = dL/dqd, M = dp/dqd and h = (dp/dq) qd - dL/dq, the equations
%
%       M qdd + h = Q + A' lambda,     A qdd + Ad qd = 0
%
%   are solved as the block system [M A'; A 0] [qdd; -lambda] = [Q - h;
%   -Ad qd], where Ad = sum_j (dA/dq_j) qd_j. The result is returned as
%   numeric function handles (matlabFunction) and as the symbolic reduced
%   expressions, so a document can print them. Physics never enters this
%   file: the caller supplies L, A and Q.
%
%% Inputs:
%
%  q                        [n x 1] sym             Generalized coordinates
%
%  qd                       [n x 1] sym             Generalized velocities
%
%  L                        sym scalar              Lagrangian L(q, qd)
%
%  A                        [k x n] sym             Constraint matrix,
%                                                   A(q) qd = 0; k = 0 is
%                                                   sym(zeros(0, n))
%
%  Q                        [n x 1] sym (optional)  Generalized forces
%                                                   Q(q, qd); zeros
%
%% Outputs:
%
%  sys                      struct                  .f(t, y): first-order
%                                                   RHS, y = [q; qd],
%                                                   returns [qd; qdd]
%                                                   [2n x 1]; .lambda(q,
%                                                   qd) [k x 1]; .energy(q,
%                                                   qd) = p' qd - L;
%                                                   .constraint(q, qd) = A
%                                                   qd [k x 1]; .qddSym [n
%                                                   x 1]; .lambdaSym [k x
%                                                   1]; .n; .k
%
%% References:
%   [1] Bloch, "Nonholonomic Mechanics and Control," 2nd ed., Springer,
%       2015, Sec. 1.4 (Lagrange-d'Alembert principle).
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the planar pendulum as a constrained particle; lambda = -T/l
     syms x y xd yd real
     m = 1;  g = 9.81;  l = 1;
     sys = mech.lagrange_dalembert([x; y], [xd; yd], m/2*(xd^2 + yd^2) - m*g*y, [x, y]);
     v = 0.7;
     lam = sys.lambda([0; -l], [v; 0]);
     fprintf('pendulum at the bottom, v = %.2f: lambda = %.6f, -T/l = %.6f\n', v, lam, -(m*g + m*v^2/l)/l);
     disp('reduced accelerations:');  disp(sys.qddSym);
     return
end

q = q(:);  qd = qd(:);
n = numel(q);
if nargin < 5 || isempty(Q), Q = sym(zeros(n, 1)); end
Q = Q(:);
if isempty(A), A = sym(zeros(0, n)); end
A = sym(A);
if ~isa(L, 'sym')
    error('mech:lagrange_dalembert:shape', 'L must be symbolic (a sym expression in q and qd)');
end
if numel(qd) ~= n || numel(Q) ~= n || ~isscalar(L) || size(A, 2) ~= n
    error('mech:lagrange_dalembert:shape', ...
          'q and qd must be [n x 1], Q [n x 1], L scalar and A [k x n] (n = %d; got qd %d, Q %d, A columns %d)', ...
          n, numel(qd), numel(Q), size(A, 2));
end
k = size(A, 1);

p = jacobian(L, qd).';                      % generalized momenta [n x 1]
M = jacobian(p, qd);                        % mass matrix [n x n]
if isAlways(simplify(det(M)) == 0, 'Unknown', 'false')
    error('mech:lagrange_dalembert:M', 'the mass matrix d2L/dqd2 is singular');
end
if k > 0 && rank(A) < k
    error('mech:lagrange_dalembert:A', 'the constraint matrix A has rank %d < k = %d', rank(A), k);
end
h  = jacobian(p, q)*qd - jacobian(L, q).';  % Euler-Lagrange terms without M qdd
Ad = sym(zeros(k, n));
for jj = 1:n
    Ad = Ad + diff(A, q(jj))*qd(jj);
end

K   = [M, A.'; A, sym(zeros(k, k))];
rhs = [Q - h; -Ad*qd];
sol = K \ rhs;
qdd = simplify(sol(1:n));
lam = simplify(-sol(n+1:end));
E   = simplify(p.'*qd - L);

fq = matlabFunction(qdd, 'Vars', {q, qd});
fE = matlabFunction(E,   'Vars', {q, qd});
if k > 0
    fl = matlabFunction(lam, 'Vars', {q, qd});
    fA = matlabFunction(A,   'Vars', {q});
    sys.lambda     = @(qq, qqd) reshape(fl(qq(:), qqd(:)), k, 1);
    sys.constraint = @(qq, qqd) reshape(fA(qq(:))*qqd(:), k, 1);
else
    sys.lambda     = @(qq, qqd) zeros(0, 1);
    sys.constraint = @(qq, qqd) zeros(0, 1);
end
sys.f         = @(t, y) [y(n+1:2*n); reshape(fq(y(1:n), y(n+1:2*n)), n, 1)];
sys.energy    = @(qq, qqd) fE(qq(:), qqd(:));
sys.qddSym    = qdd;
sys.lambdaSym = lam;
sys.n         = n;
sys.k         = k;
end
