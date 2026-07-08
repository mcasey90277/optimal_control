function [cIneq, ceq, gIneq, gCeq] = nlp_constraints(Z, sigma, Tmax, c, muStar)
% NLP_CONSTRAINTS  Trapezoidal defects + unit-sphere control equalities.
%
% Nonlinear constraint function for the direct transcription of the
% min-time low-thrust transfer on a (generally NONUNIFORM) mesh. The node
% times are t_k = tf*sigma_k with FIXED normalized fractions sigma
% (sigma_1 = 0, sigma_{N+1} = 1), so the free final time scales every
% segment width h_k = tf*(sigma_{k+1}-sigma_k) smoothly. A nonuniform,
% dynamics-density-matched mesh is essential here: the GTO spiral's perigee
% passes are ~1000x faster than its apogee coasts, and a uniform mesh
% either wastes nodes or (fatally) under-resolves perigee.
%
% Decision vector:
%   Z = [X(:); W(:); tf],  X [7x(N+1)], W [3x(N+1)],
% nZ = 10*(N+1) + 1. Equality constraints:
%   defects  d_k = x_{k+1} - x_k - h_k/2*(f_k + f_{k+1}) = 0   (7N eqs)
%   sphere   g_k = w_k'*w_k - 1 = 0                            (N+1 eqs)
% The throttle is fixed at 1 (always burn; justified post-hoc by the
% indirect switching function staying strictly negative), so the only
% control is the unit thrust direction. This keeps every bound inactive
% at the warm start -- essential for interior-point warm starting.
%
% INPUTS:
%   Z      - decision vector [10*(N+1)+1 x 1]
%   sigma  - normalized node times, increasing, sigma(1)=0, sigma(end)=1
%            [(N+1)x1]
%   Tmax   - max thrust acceleration at m = 1 (ND) [scalar]
%   c      - exhaust velocity (ND) [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   cIneq  - inequality constraints (none) []
%   ceq    - equality constraints [7N + (N+1) x 1]
%   gIneq  - gradient of cIneq (none) []
%   gCeq   - gradient of ceq, SPARSE [nZ x (7N + N+1)] (fmincon layout:
%            one column per constraint)
%
% REFERENCES:
%   [1] Betts, "Practical Methods for Optimal Control and Estimation Using
%       Nonlinear Programming," 2nd ed., SIAM, 2010 (Ch. 4).

sigma  = sigma(:);
N      = numel(sigma) - 1;
nNodes = N + 1;
[X, W, tf] = unpack_z(Z, N);

dSig = diff(sigma).';                    % [1xN] normalized widths
h    = tf.*dSig;                         % [1xN] segment widths

if nargout > 2
    [F, A, B] = lt_dynamics(X, W, Tmax, c, muStar);
else
    F = lt_dynamics(X, W, Tmax, c, muStar);   % value-only: skip Jacobians
end

% --- defects (vectorized) -------------------------------------------------
defects = X(:, 2:end) - X(:, 1:end-1) ...
          - (F(:, 1:end-1) + F(:, 2:end)).*(h./2);

% --- unit-sphere control constraint ---------------------------------------
cone = sum(W.^2, 1) - 1;

cIneq = [];
ceq   = [defects(:); cone(:)];

if nargout > 2
    gIneq = [];

    nZ    = 10*nNodes + 1;
    nDef  = 7*N;
    uOff  = 7*nNodes;                    % offset of U block within Z
    tfCol = nZ;                          % index of tf within Z

    % Preallocate triplets: per segment 49+49+21+21+7 = 147; sphere 3/node.
    % Index patterns are identical for every segment, so build the
    % per-segment row/column templates ONCE and add offsets in the loop
    % (a per-block meshgrid helper costs real time at N ~ 1e4).
    nTrip = 147*N + 3*nNodes;
    rowT  = zeros(nTrip, 1);
    colT  = zeros(nTrip, 1);
    valT  = zeros(nTrip, 1);

    I7    = eye(7);
    rows7 = (1:7).';
    % templates for one segment, k = 1 (columns relative to that segment)
    rXk   = repmat(rows7, 7, 1);          cXk  = repelem((1:7).', 7);
    rUk   = repmat(rows7, 3, 1);          cUk  = repelem((1:3).', 7);
    segRow = [rXk; rXk; rUk; rUk; rows7];
    segCol = [cXk; 7 + cXk; zeros(21,1); zeros(21,1); zeros(7,1)];  % X parts

    for k = 1:N
        hk    = h(k);
        base  = (k-1)*147;
        rOff  = (k-1)*7;

        M1 = -I7 - (hk/2)*A(:,:,k);      % d d_k / d x_k
        M2 =  I7 - (hk/2)*A(:,:,k+1);    % d d_k / d x_{k+1}
        B1 = -(hk/2)*B(:,:,k);           % d d_k / d w_k
        B2 = -(hk/2)*B(:,:,k+1);         % d d_k / d w_{k+1}
        dTf = -(dSig(k)/2)*(F(:,k) + F(:,k+1));   % d d_k / d tf

        rowT(base+1:base+147) = rOff + segRow;
        colT(base+1:base+ 98) = (k-1)*7 + segCol(1:98);
        colT(base+ 99:base+119) = uOff + (k-1)*3 + cUk;
        colT(base+120:base+140) = uOff +  k   *3 + cUk;
        colT(base+141:base+147) = tfCol;
        valT(base+1:base+147) = [M1(:); M2(:); B1(:); B2(:); dTf];
    end

    base = 147*N;
    coneRows = nDef + repelem((1:nNodes).', 3);
    coneCols = uOff + reshape((0:nNodes-1)*3 + (1:3).', [], 1);
    rowT(base+1:end) = coneRows;
    colT(base+1:end) = coneCols;
    valT(base+1:end) = 2*W(:);

    Jceq = sparse(rowT, colT, valT, nDef + nNodes, nZ);
    gCeq = Jceq.';                       % fmincon wants nZ x nConstraints
end
end
