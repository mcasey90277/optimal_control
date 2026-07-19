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

npos = 0; nneg = 0; nzero = 0;
i = 1; nD = size(D,1);
while i <= nD
    if i < nD && D(i+1,i) ~= 0
        % 2x2 block: classify both eigenvalues of the symmetrized block
        b = full(D(i:i+1, i:i+1));
        ev = eig((b + b.') / 2);
        for e = ev.'
            if e > zt
                npos = npos + 1;
            elseif e < -zt
                nneg = nneg + 1;
            else
                nzero = nzero + 1;
            end
        end
        i = i + 2;
    else
        % 1x1 block
        e = D(i,i);
        if e > zt
            npos = npos + 1;
        elseif e < -zt
            nneg = nneg + 1;
        else
            nzero = nzero + 1;
        end
        i = i + 1;
    end
end

IN.npos       = npos;
IN.nneg       = nneg;
IN.nzero      = nzero;
IN.expected   = [n, ma, 0];
IN.subspaceOK = isequal([npos nneg nzero], [n ma 0]);
IN.redMinEig  = NaN;   % non-gating placeholder; not computed here
end
