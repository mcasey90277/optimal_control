function [R, J, info] = ztl_ms_residual_sun(z, prob, wantJ)
% ZTL_MS_RESIDUAL_SUN  Sundman-regularized multiple-shooting residual + exact
% block Jacobian (per SUN_BUILD.md).
%
% Unknowns z = [lam0(7); Y_2(15); ...; Y_M(15); tauF(1)], dim 15M-7.
% Nodes at fixed sigma_k=(k-1)/M; node 1 = [rv0;1;lam0;0]. Arc k integrates
% Y_k over [sigma_k, sigma_{k+1}] (Sundman, scaled by tauF) -> F_k, STM
% Phi_k (15x15), tauF-sensitivity w_k (15x1).
%
% R = [ F_k - Y_{k+1}   (k=1..M-1, continuity, 15) ;
%       F_M(1:6)-rvf ; F_M(14) ; F_M(15)-tf   (terminal, 8) ]
% Jacobian blocks: continuity k -> d/dlam0=Phi_1(:,8:14) or d/dY_k=Phi_k,
%   d/dY_{k+1}=-I(15), d/dtauF=w_k; terminal -> d/dY_M=[Phi_M(rows 1:6,14,15)],
%   d/dtauF=w_M(rows 1:6,14,15).
%
% INPUTS:
%   z, prob (.rv0 .rvf .tf .sig [1x(M+1)] .M .P), wantJ [true]
% OUTPUTS:
%   R, J [15M-7 sq], info (.maxCont .termErr .nSegMax .grazed)

if nargin < 3, wantJ = true; end
M = prob.M;  P = prob.P;  sN = prob.sig;
n = 15*M - 7;
tauF = z(end);

Y = cell(1,M);
Y{1} = [prob.rv0(:); 1; z(1:7); 0];
for k = 2:M
    Y{k} = z(7 + 15*(k-2) + (1:15));
end

F = cell(1,M);  Phi = cell(1,M);  W = cell(1,M);
nSegMax = 0;  grazed = false;
for k = 1:M
    o = ztl_flow_sun(Y{k}, tauF, [sN(k) sN(k+1)], P, wantJ);
    F{k} = o.Yf;  Phi{k} = o.PHI;  W{k} = o.w;
    nSegMax = max(nSegMax, o.nSegs);  grazed = grazed || (o.flag==1);
end

R = zeros(n,1);  maxCont = 0;
for k = 1:M-1
    rk = F{k} - Y{k+1};
    R(15*(k-1) + (1:15)) = rk;
    maxCont = max(maxCont, norm(rk));
end
termR = [F{M}(1:6) - prob.rvf(:); F{M}(14); F{M}(15) - prob.tf];
R(15*(M-1) + (1:8)) = termR;

J = [];
if wantJ
    J = zeros(n,n);
    colLam0 = 1:7;
    colY = @(k) 7 + 15*(k-2) + (1:15);      % valid k=2..M
    colTau = n;
    for k = 1:M-1
        rows = 15*(k-1) + (1:15);
        if k == 1
            J(rows, colLam0) = Phi{1}(:, 8:14);
        else
            J(rows, colY(k)) = Phi{k};
        end
        J(rows, colY(k+1)) = -eye(15);
        J(rows, colTau) = W{k};
    end
    rowsT = 15*(M-1) + (1:8);
    sel = [1 2 3 4 5 6 14 15];
    J(rowsT, colY(M)) = Phi{M}(sel, :);
    J(rowsT, colTau)  = W{M}(sel);
end

info = struct('maxCont',maxCont,'termErr',norm(termR),'nSegMax',nSegMax,'grazed',grazed);
end
