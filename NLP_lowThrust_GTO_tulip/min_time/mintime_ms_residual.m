function [R, J, info] = mintime_ms_residual(z, prob, wantJac)
% MINTIME_MS_RESIDUAL  Multiple-shooting residual + block Jacobian for the
% min-TIME CR3BP TPBVP (free final time, always-burn PMP), built on pumpkyn's
% analytic-STM min-time propagator.
%
% Single shooting through ~13 revs floors at ||R||~1e-3 (the ~1e6 STM-product
% sensitivity). Multiple shooting splits the transfer into M arcs whose
% per-arc STMs are well-conditioned; the tiny continuity gaps absorb the
% amplification, so the root-find reaches machine precision (the same reason
% MS beat the wall for the energy problem).
%
% Unknowns  z = [lambda0(7); Y_2..Y_M (14 each); tf]  (dim 14(M-1)+8).
%   Y_k = [r(3); v(3); m; lambda_r(3); lambda_v(3); lambda_m]  (14), arc-k start
%   (arc 1 starts at [rv0; 1; lambda0]).
% Residual  R = [ continuity Yend_k - Y_{k+1} (14 each, k=1..M-1);
%                 rendezvous Yend_M(1:6) - rvf (6);
%                 lambda_m(tf) = Yend_M(14)  (1);
%                 H(tf) = 0                   (1) ]           (dim 14(M-1)+8).
% Square system: 8 free scalars (lambda0,tf) <-> 8 terminal conditions.
%
% INPUTS:
%   z    - decision vector [ (14(M-1)+8) x 1 ]
%   prob - .rv0 .rvf [1x6] .sig [(M+1)x1] .M .Tmax .c .muStar
%   wantJac - build J (and STMs) if true [logical]
%
% OUTPUTS:
%   R    - residual [ (14(M-1)+8) x 1 ]
%   J    - sparse Jacobian dR/dz (empty if ~wantJac)
%   info - .maxCont .termErr .grazed .nSwitch (switch count across arcs; the
%          min-time arc should be all-burn, nSwitch=0, or J is approximate)
%
% REFERENCES:
%   [1] pumpkyn.cr3bp.tfMinProp/tfMinEoM (Koblick) - min-time field + 14x14 STM.
%   [2] ztl_ms_residual (energy MS analog); Zhang JGCD 38(8) 2015 (STM shooting).

nY = 14;  M = prob.M;  sig = prob.sig(:);
Tmax = prob.Tmax;  c = prob.c;  mu = prob.muStar;
tf = z(end);
lam0 = z(1:7);

% arc start states: column k is the start of arc k (k=1..M)
Ystart = zeros(nY, M);
Ystart(:,1) = [prob.rv0(:); 1; lam0(:)];
for k = 2:M
    Ystart(:,k) = z(7 + nY*(k-2) + (1:nY));
end

Yend = zeros(nY, M);
if wantJac, PHI = cell(1,M);  Wtf = zeros(nY, M); end
nSwitch = 0;
for k = 1:M
    dt = (sig(k+1) - sig(k)) * tf;
    if wantJac
        y0 = [Ystart(:,k); reshape(eye(nY), [], 1)];
        [~, y] = pumpkyn.cr3bp.tfMinProp(dt, y0, Tmax, c, mu);
        Yend(:,k) = y(end, 1:nY).';
        PHI{k}    = reshape(y(end, nY+1:nY+nY^2), nY, nY);
        yDot = pumpkyn.cr3bp.tfMinEoM(0, Yend(:,k), Tmax, c, mu);
        Wtf(:,k)  = yDot(:) * (sig(k+1) - sig(k));   % dYend/dtf (arcs scale with tf)
    else
        [~, y] = pumpkyn.cr3bp.tfMinProp(dt, Ystart(:,k), Tmax, c, mu);
        Yend(:,k) = y(end, 1:nY).';
    end
end

% --- residual ---------------------------------------------------------------
Rc = Yend(:,1:M-1) - Ystart(:,2:M);              % continuity (14 x (M-1))
[~, Ht, dHdy] = pumpkyn.cr3bp.tfMinEoM(0, Yend(:,M), Tmax, c, mu);
Rt = [Yend(1:6,M) - prob.rvf(:); Yend(14,M); Ht];  % terminal (8 x 1)
R  = [Rc(:); Rt];

info = struct('maxCont', max([abs(Rc(:)); 0]), 'termErr', norm(Rt), ...
              'grazed', false, 'nSwitch', nSwitch);

% --- Jacobian ---------------------------------------------------------------
J = [];
if wantJac
    nZ = 14*(M-1) + 8;   nR = nZ;
    J = sparse(nR, nZ);
    colY = @(k) 7 + nY*(k-2) + (1:nY);   % z-columns of node Y_k (k=2..M)
    Inn = speye(nY);
    for k = 1:M-1
        rows = nY*(k-1) + (1:nY);
        if k == 1
            J(rows, 1:7) = PHI{1}(:, 8:14);          % dYend_1/dlambda0
        else
            J(rows, colY(k)) = PHI{k};               % dYend_k/dY_k
        end
        J(rows, colY(k+1)) = -Inn;                   % d(-Y_{k+1})
        J(rows, nZ) = Wtf(:,k);                      % dYend_k/dtf
    end
    rowsT = 14*(M-1) + (1:8);
    selM  = PHI{M}([1:6, 14], :);                    % rendezvous + lambda_m rows
    J(rowsT, colY(M)) = [selM; dHdy * PHI{M}];       % + H row
    J(rowsT, nZ) = [Wtf([1:6,14], M); dHdy * Wtf(:,M)];
end
end
