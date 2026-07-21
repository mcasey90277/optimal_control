function [cIneq, ceq, gIneq, gCeq] = nlp_constraints_minfuel(Z, sigma, tf, Tmax, c, muStar)
% NLP_CONSTRAINTS_MINFUEL  Defects + throttle-cone equalities, fixed tf.
%
% Constraint function for the MIN-FUEL direct transcription: the general
% 4-control form (thrust-direction vector w + throttle s) on a fixed,
% generally nonuniform mesh. tf is a PARAMETER here (fixed-time problem),
% so the decision vector is Z = [X(:); U(:)], X [7x(N+1)], U [4x(N+1)],
% nZ = 11*(N+1). Equalities:
%   defects  d_k = x_{k+1} - x_k - h_k/2*(f_k + f_{k+1}) = 0   (7N eqs),
%            h_k = tf*(sigma_{k+1}-sigma_k)  (constants),
%   cone     g_k = w_k'*w_k - s_k^2 = 0                        (N+1 eqs).
% The cone equality prevents the ballast exploit (see LT_DYNAMICS_THROTTLE).
%
% INPUTS:
%   Z      - decision vector [11*(N+1) x 1]
%   sigma  - normalized node times, increasing, sigma(1)=0, sigma(end)=1
%            [(N+1)x1]
%   tf     - fixed transfer time (ND) [scalar]
%   Tmax   - max thrust acceleration at m = 1 (ND) [scalar]
%   c      - exhaust velocity (ND) [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   cIneq  - [] (no inequalities)
%   ceq    - equality constraints [7N + (N+1) x 1]
%   gIneq  - []
%   gCeq   - gradient of ceq, SPARSE [nZ x (7N + N+1)]
%
% REFERENCES:
%   [1] Betts, SIAM 2010, Ch. 4.

sigma  = sigma(:);
N      = numel(sigma) - 1;
nNodes = N + 1;

X = reshape(Z(1:7*nNodes), 7, nNodes);
U = reshape(Z(7*nNodes + (1:4*nNodes)), 4, nNodes);

dSig = diff(sigma).';
h    = tf.*dSig;

if nargout > 2
    [F, A, B] = lt_dynamics_throttle(X, U, Tmax, c, muStar);
else
    F = lt_dynamics_throttle(X, U, Tmax, c, muStar);
end

defects = X(:, 2:end) - X(:, 1:end-1) ...
          - (F(:, 1:end-1) + F(:, 2:end)).*(h./2);

w    = U(1:3, :);
s    = U(4, :);
cone = sum(w.^2, 1) - s.^2;

cIneq = [];
ceq   = [defects(:); cone(:)];

if nargout > 2
    gIneq = [];

    nZ   = 11*nNodes;
    nDef = 7*N;
    uOff = 7*nNodes;

    % per segment: 49 + 49 + 28 + 28 = 154 triplets; cone 4/node
    nTrip = 154*N + 4*nNodes;
    rowT  = zeros(nTrip, 1);
    colT  = zeros(nTrip, 1);
    valT  = zeros(nTrip, 1);

    I7    = eye(7);
    rows7 = (1:7).';
    rXk   = repmat(rows7, 7, 1);   cXk = repelem((1:7).', 7);
    rUk   = repmat(rows7, 4, 1);   cUk = repelem((1:4).', 7);
    segRow = [rXk; rXk; rUk; rUk];

    for k = 1:N
        hk   = h(k);
        base = (k-1)*154;
        rOff = (k-1)*7;

        M1 = -I7 - (hk/2)*A(:,:,k);
        M2 =  I7 - (hk/2)*A(:,:,k+1);
        B1 = -(hk/2)*B(:,:,k);
        B2 = -(hk/2)*B(:,:,k+1);

        rowT(base+1:base+154)    = rOff + segRow;
        colT(base+1:base+49)     = (k-1)*7 + cXk;
        colT(base+50:base+98)    =  k   *7 + cXk;
        colT(base+99:base+126)   = uOff + (k-1)*4 + cUk;
        colT(base+127:base+154)  = uOff +  k   *4 + cUk;
        valT(base+1:base+154)    = [M1(:); M2(:); B1(:); B2(:)];
    end

    base = 154*N;
    coneRows = nDef + repelem((1:nNodes).', 4);
    coneCols = uOff + reshape((0:nNodes-1)*4 + (1:4).', [], 1);
    coneVals = reshape([2*w; -2*s], [], 1);
    rowT(base+1:end) = coneRows;
    colT(base+1:end) = coneCols;
    valT(base+1:end) = coneVals;

    Jceq = sparse(rowT, colT, valT, nDef + nNodes, nZ);
    gCeq = Jceq.';
end
end
