function [cIneq, ceq, gIneq, gCeq] = nlp_constraints_minfuel_hs(Z, sigma, tf, Tmax, c, muStar)
% NLP_CONSTRAINTS_MINFUEL_HS  Separated Hermite-Simpson min-fuel constraints.
%
% 4th-order (Hermite-Simpson) collocation of the cone-eliminated min-fuel
% problem (direction+throttle control via LT_DYNAMICS_DIRTHROTTLE). Each
% segment carries a MIDPOINT with its own state x_c and control u_c as
% decision variables (the SEPARATED form -- keeps the Jacobian free of
% chain-rule terms through the interpolated midpoint). Two vector defects
% per segment:
%   interpolation (Hermite):  x_c - [1/2(x_k+x_{k+1}) + h/8(f_k - f_{k+1})] = 0
%   collocation (Simpson):    x_{k+1} - x_k - h/6(f_k + 4 f_c + f_{k+1})    = 0
% plus unit-direction constraints ||alpha|| = 1 at every node AND midpoint.
%
% Decision vector layout (N segments, nN = N+1 nodes):
%   Z = [ X(:);  U(:);  Xc(:);  Uc(:) ]
%       X  [7 x nN]   endpoint states     (offset oX  = 0)
%       U  [4 x nN]   endpoint controls   (offset oU  = 7 nN)
%       Xc [7 x N]    midpoint states     (offset oXc = 11 nN)
%       Uc [4 x N]    midpoint controls   (offset oUc = 11 nN + 7 N)
%   nZ = 11 nN + 11 N. Controls are [alpha(3); s].
%
% Constraint ordering: [D1(:) (7N); D2(:) (7N); unitNode (nN); unitMid (N)],
% total 16N + 1.
%
% REFERENCES:
%   [1] Betts, SIAM 2010, Ch. 4 (separated Hermite-Simpson).
%   [2] Kelly, SIAM Review 59(4), 2017.

sigma = sigma(:);
N  = numel(sigma) - 1;
nN = N + 1;
oX = 0;  oU = 7*nN;  oXc = 7*nN + 4*nN;  oUc = oXc + 7*N;
nZ = oUc + 4*N;

X  = reshape(Z(oX  + (1:7*nN)), 7, nN);
U  = reshape(Z(oU  + (1:4*nN)), 4, nN);
Xc = reshape(Z(oXc + (1:7*N )), 7, N );
Uc = reshape(Z(oUc + (1:4*N )), 4, N );

h = (tf.*diff(sigma)).';                       % 1 x N segment widths

if nargout > 2
    [Fk, Ak, Bk] = lt_dynamics_dirthrottle(X,  U,  Tmax, c, muStar);
    [Fc, Ac, Bc] = lt_dynamics_dirthrottle(Xc, Uc, Tmax, c, muStar);
else
    Fk = lt_dynamics_dirthrottle(X,  U,  Tmax, c, muStar);
    Fc = lt_dynamics_dirthrottle(Xc, Uc, Tmax, c, muStar);
end

D1 = Xc - 0.5*(X(:,1:end-1) + X(:,2:end)) - (h/8).*(Fk(:,1:end-1) - Fk(:,2:end));
D2 = X(:,2:end) - X(:,1:end-1) - (h/6).*(Fk(:,1:end-1) + 4*Fc + Fk(:,2:end));
unitNode = sum(U(1:3,:).^2, 1) - 1;
unitMid  = sum(Uc(1:3,:).^2, 1) - 1;

cIneq = [];  gIneq = [];
ceq   = [D1(:); D2(:); unitNode(:); unitMid(:)];

if nargout > 2
    nTrip = 434*N + 3*(2*N + 1);
    rowT = zeros(nTrip, 1);  colT = zeros(nTrip, 1);  valT = zeros(nTrip, 1);
    ptr  = 0;
    I7   = eye(7);

    for k = 1:N
        hk = h(k);
        colXk  = oX  + (k-1)*7 + (1:7);   colXk1 = oX  + k*7 + (1:7);
        colUk  = oU  + (k-1)*4 + (1:4);   colUk1 = oU  + k*4 + (1:4);
        colXc  = oXc + (k-1)*7 + (1:7);   colUc  = oUc + (k-1)*4 + (1:4);
        rD1 = (k-1)*7 + (1:7);
        rD2 = 7*N + (k-1)*7 + (1:7);
        Akk = Ak(:,:,k);  Ak1 = Ak(:,:,k+1);  Bkk = Bk(:,:,k);  Bk1 = Bk(:,:,k+1);
        Ack = Ac(:,:,k);  Bck = Bc(:,:,k);

        addB(rD1, colXk , -0.5*I7 - (hk/8)*Akk);
        addB(rD1, colUk , -(hk/8)*Bkk);
        addB(rD1, colXk1, -0.5*I7 + (hk/8)*Ak1);
        addB(rD1, colUk1,  (hk/8)*Bk1);
        addB(rD1, colXc ,  I7);

        addB(rD2, colXk , -I7 - (hk/6)*Akk);
        addB(rD2, colUk , -(hk/6)*Bkk);
        addB(rD2, colXk1,  I7 - (hk/6)*Ak1);
        addB(rD2, colUk1, -(hk/6)*Bk1);
        addB(rD2, colXc , -(2*hk/3)*Ack);
        addB(rD2, colUc , -(2*hk/3)*Bck);
    end

    % unit-direction rows: d(alpha'alpha - 1)/d(alpha) = 2 alpha
    for kk = 1:nN
        addB(14*N + kk, oU + (kk-1)*4 + (1:3), (2*U(1:3,kk)).');
    end
    for kk = 1:N
        addB(14*N + nN + kk, oUc + (kk-1)*4 + (1:3), (2*Uc(1:3,kk)).');
    end

    Jceq = sparse(rowT(1:ptr), colT(1:ptr), valT(1:ptr), 16*N + 1, nZ);
    gCeq = Jceq.';
end

    function addB(rows, cols, Mblk)
        rr = numel(rows);  cc = numel(cols);  nn = rr*cc;
        rowT(ptr+1:ptr+nn) = repmat(rows(:), cc, 1);
        colT(ptr+1:ptr+nn) = repelem(cols(:), rr);
        valT(ptr+1:ptr+nn) = Mblk(:);
        ptr = ptr + nn;
    end
end
