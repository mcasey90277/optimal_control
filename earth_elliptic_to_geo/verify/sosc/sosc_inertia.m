function IN = sosc_inertia(H, A, tol)
% SOSC_INERTIA  Inertia of the KKT matrix [H A'; A 0]; decide the subspace
% second-order condition (reduced Hessian PD on null(A)). Gold-standard
% dense `eig` is the primary method when the KKT is small enough
% (size-guarded by tol.maxEigDim); `ldl`-pivot-sign inertia is a size-guard
% fallback and is flagged non-robust (sec 11.4).
%
% INPUTS:
%   H   - Lagrangian Hessian [n x n sparse]
%   A   - active-constraint Jacobian [m_a x n sparse]
%   tol - tolerance struct from sosc_defaults() (uses tol.inertiaZero,
%         tol.maxEigDim)
%
% OUTPUTS:
%   IN - struct:
%        .npos        - count of positive eigenvalues of KKT matrix
%        .nneg        - count of negative eigenvalues
%        .nzero       - count of (near-)zero eigenvalues
%        .expected    - [n, m_a, 0], the target inertia signature
%        .subspaceOK  - true iff [npos nneg nzero] == expected (bool)
%                       (reported only; no longer drives the verdict, sec 11.4)
%        .method      - 'eig' (gold-standard, size <= tol.maxEigDim) or
%                        'ldl' (size guard fallback, sec 11.4)
%        .robust      - true iff .method=='eig' (a validated inertia);
%                        false for the 'ldl' fallback (sec 11.4, 11.5 rule 2)
%        .rankA       - structural rank of A (sprank), r
%        .red         - reduced-Hessian inertia struct .npos .nneg .nzero
%                       via the Gould decomposition (sec 11.4)
%        .redConsistent - true iff red inertia is self-consistent (bool)
%        .redMinEig   - NaN placeholder; non-gating curvature margin,
%                       not computed here (reserved for later enhancement)
%
% REFERENCES:
%   [1] process/DESIGN_sosc.md sec 4.5, 11.4.
%   [2] Nocedal & Wright, "Numerical Optimization," 2nd ed., Thm 16.3:
%       inertia(KKT) = (n, m_a, 0) <=> reduced Hessian PD on null(A),
%       when A has full row rank.
%   [3] Gould, N.I.M., "On practical conditions for the existence and
%       uniqueness of solutions to the general equality quadratic
%       programming problem," Math. Prog. 32 (1985): the reduced-Hessian
%       inertia from the full KKT inertia and rank(A):
%       inertia(Z'HZ) = (npos-r, nneg-r, nzero-(m_a-r)), r = rank(A).

n  = size(H,1);
ma = size(A,1);

K = [H, A.'; A, sparse(ma,ma)];
K = (K + K.') / 2;                              % symmetrize numerically
nk = size(K,1);

scale = max(1, normest(K));
zt = tol.inertiaZero * scale;

if nk <= tol.maxEigDim
    % Gold standard: dense eig, Sylvester-exact, robust across a wide zt
    % window on the near-singular bang-bang KKT (sec 11.4).
    ev = eig(full(K));
    npos  = sum(ev >  zt);
    nneg  = sum(ev < -zt);
    nzero = sum(abs(ev) <= zt);
    IN.method = 'eig';
    IN.robust = true;
else
    % Size guard: too large for dense eig. Fall back to ldl-pivot-sign
    % inertia, which is UNRELIABLE on this problem class (sec 11.4) --
    % flag as non-robust so sosc_decide defers to INCONCLUSIVE rather than
    % trusting an unvalidated inertia.
    [~, D, ~] = ldl(K, 'vector');                % D block-diagonal (1x1 & 2x2)
    [npos, nneg, nzero] = count_inertia(D, zt);
    IN.method = 'ldl';
    IN.robust = false;
end

IN.npos       = npos;
IN.nneg       = nneg;
IN.nzero      = nzero;
IN.expected   = [n, ma, 0];
IN.subspaceOK = isequal([npos nneg nzero], [n ma 0]);  % reported only (sec 11.4)

% Reduced-Hessian inertia via the Gould decomposition (sec 11.4). r is the
% structural rank of the active Jacobian; the reduced Hessian Z'HZ (H on
% null(A)) has inertia (npos-r, nneg-r, nzero-(m_a-r)).
r = sprank(A);
red.npos  = npos  - r;
red.nneg  = nneg  - r;
red.nzero = nzero - (ma - r);
IN.rankA  = r;
IN.red    = struct('npos', red.npos, 'nneg', red.nneg, 'nzero', red.nzero);
IN.redConsistent = (red.npos + red.nneg + red.nzero == n - r) ...
                   && red.nneg >= 0 && red.nzero >= 0;

IN.redMinEig  = NaN;   % non-gating placeholder; not computed here
end
