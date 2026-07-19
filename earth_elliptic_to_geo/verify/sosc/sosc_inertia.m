function IN = sosc_inertia(H, A, tol)
% SOSC_INERTIA  Inertia of the KKT matrix [H A'; A 0] via sparse LDL^T; decide
% the subspace second-order condition (reduced Hessian PD on null(A)).
%
% INPUTS:
%   H   - Lagrangian Hessian [n x n sparse]
%   A   - active-constraint Jacobian [m_a x n sparse]
%   tol - tolerance struct from sosc_defaults() (uses tol.inertiaZero)
%
% OUTPUTS:
%   IN - struct:
%        .npos        - count of positive eigenvalues of KKT matrix
%        .nneg        - count of negative eigenvalues
%        .nzero       - count of (near-)zero eigenvalues
%        .expected    - [n, m_a, 0], the target inertia signature
%        .subspaceOK  - true iff [npos nneg nzero] == expected (bool)
%        .redMinEig   - NaN placeholder; non-gating curvature margin,
%                       not computed here (reserved for later enhancement)
%
% REFERENCES:
%   [1] process/DESIGN_sosc.md sec 4.5.
%   [2] Nocedal & Wright, "Numerical Optimization," 2nd ed., Thm 16.3:
%       inertia(KKT) = (n, m_a, 0) <=> reduced Hessian PD on null(A),
%       when A has full row rank.

n  = size(H,1);
ma = size(A,1);

K = [H, A.'; A, sparse(ma,ma)];
K = (K + K.') / 2;                              % symmetrize numerically

[~, D, ~] = ldl(K, 'vector');                   % D block-diagonal (1x1 & 2x2)

scale = max(1, normest(K));
zt = tol.inertiaZero * scale;

[npos, nneg, nzero] = count_inertia(D, zt);

IN.npos       = npos;
IN.nneg       = nneg;
IN.nzero      = nzero;
IN.expected   = [n, ma, 0];
IN.subspaceOK = isequal([npos nneg nzero], [n ma 0]);
IN.redMinEig  = NaN;   % non-gating placeholder; not computed here
end
