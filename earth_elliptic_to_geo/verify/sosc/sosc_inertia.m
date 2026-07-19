function IN = sosc_inertia(H, A, tol)
% SOSC_INERTIA  Reduced-Hessian inertia via a DIRECT dense null-space basis
% (FINAL method, process/DESIGN_sosc.md sec 12.1). Forms Z = null(full(A)) (an
% orthonormal basis of null(A), where A is the strongly-active Jacobian), the
% reduced Hessian RH = Z'HZ, and counts eigenvalue signs of RH over a zt-band
% to expose whether the near-zero-eigenvalue classification is threshold-stable.
% This SUPERSEDES the Gould-decomposition + sprank + single-zt method of sec
% 11.4 (which gave spurious FAILs when RH had near-flat directions of
% unresolvable sign). No sprank, no rank subtraction: Z gives the rank exactly.
%
% INPUTS:
%   H   - Lagrangian Hessian [n x n sparse]
%   A   - strongly-active-constraint Jacobian [m_a x n sparse]
%   tol - tolerance struct from sosc_defaults() (uses tol.inertiaZero,
%         tol.maxNullDim)
%
% OUTPUTS:
%   IN - struct:
%        .red        - reduced-Hessian inertia at the tightest zt=1e-9*s:
%                      struct .npos .nneg .nzero (counts of eig(RH))
%        .nnegBand   - [1x4] #{ev < -zt} over ztr = [1e-9 1e-8 1e-7 1e-6]*s
%        .sensStable - true iff nnegBand is constant across the band (the
%                      negative-count is threshold-insensitive => trustworthy)
%        .robust     - true iff the reduced inertia was computed (n<=maxNullDim)
%        .method     - 'reduced-eig' (direct) or 'scale-skip' (n>maxNullDim)
%        .rankA      - exact rank r = n - size(Z,2) from the null-space basis
%        .redMinEig  - min(eig(RH)); a REAL reported curvature margin (NaN skip)
%
% REFERENCES:
%   [1] process/DESIGN_sosc.md sec 12.1 (FINAL inertia method), sec 12.2.
%   [2] Nocedal & Wright, "Numerical Optimization," 2nd ed., Thm 16.3:
%       reduced Hessian Z'HZ PD on null(A) <=> second-order sufficiency on the
%       (subspace) critical cone.

n  = size(H,1);
ma = size(A,1);

if n > tol.maxNullDim
    % Size guard (sec 12.1): dense null-space intractable -> do not form it.
    IN.robust     = false;
    IN.method     = 'scale-skip';
    IN.red        = struct('npos', NaN, 'nneg', NaN, 'nzero', NaN);
    IN.nnegBand   = nan(1,4);
    IN.sensStable = false;
    IN.rankA      = NaN;
    IN.redMinEig  = NaN;
    return;
end

Z = null(full(A));                 % orthonormal basis of null(A), n x (n-r)
r = n - size(Z,2);                 % exact rank of A from the basis

RH = Z.' * full(H) * Z;
RH = (RH + RH.') / 2;              % symmetrize numerically
ev = sort(eig(RH));

s   = max(1, normest(full(H)));    % Hessian scale
ztr = [1e-9 1e-8 1e-7 1e-6];       % relative zt-band
nnegBand = zeros(1,numel(ztr));
for k = 1:numel(ztr)
    nnegBand(k) = sum(ev < -ztr(k)*s);
end

zt0 = tol.inertiaZero * s;         % tightest band point
IN.red = struct( ...
    'npos',  sum(ev >  zt0), ...
    'nneg',  sum(ev < -zt0), ...
    'nzero', sum(abs(ev) <= zt0));

IN.nnegBand   = nnegBand;
IN.sensStable = all(nnegBand == nnegBand(1));
IN.robust     = true;
IN.method     = 'reduced-eig';
IN.rankA      = r;
IN.redMinEig  = min(ev);
end
