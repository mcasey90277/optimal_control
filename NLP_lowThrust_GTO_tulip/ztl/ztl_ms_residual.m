function [R, J, info] = ztl_ms_residual(z, prob, wantJ)
% ZTL_MS_RESIDUAL  Multiple-shooting residual + exact block Jacobian for the
% ramp-family PMP BVP, using Z0's per-arc variational STMs.
%
% FORMULATION (Z1). Split [0, tf] into M arcs at node times
% tNodes = [t_0=0, t_1, ..., t_M=tf]. Shooting nodes carry the augmented
% state y = [r(3); v(3); m; lam_r(3); lam_v(3); lam_m] in R^14:
%   node 1 at t_0: Y_1 = [rv0; 1; lam0]   (r,v,m FIXED; lam0 the 7 unknowns)
%   nodes 2..M   : Y_k full 14 unknowns each
% Arc k integrates Y_k from tNodes(k) to tNodes(k+1) via ztl_flow -> F_k,
% with STM Phi_k = d F_k / d Y_k.
%
% UNKNOWNS  z = [lam0(7); Y_2(14); ...; Y_M(14)]           dim 14M-7
% RESIDUAL  R = [ F_k - Y_{k+1}     (k=1..M-1, continuity, 14 each);
%                 F_M(1:6) - rvf ;  F_M(14)  (terminal BC, 7) ]  dim 14M-7
% Square system. The single-shooting ill-conditioning (full-trajectory STM)
% is UNROLLED into a block-bidiagonal Jacobian whose blocks are the SHORT
% per-arc STMs -- the reason multiple shooting beats single shooting on
% sensitive BVPs (Stoer & Bulirsch; Betts).
%
% JACOBIAN (block-bidiagonal, dense assembly -- 14M-7 is small):
%   continuity block k rows:
%     k=1     : d/d lam0 = Phi_1(:,8:14) [14x7] ; d/d Y_2 = -I
%     k=2..M-1: d/d Y_k  = Phi_k [14x14]        ; d/d Y_{k+1} = -I
%   terminal rows: d/d Y_M = [Phi_M(1:6,:); Phi_M(14,:)] [7x14]
%
% INPUTS:
%   z     - unknown vector [14M-7 x 1]
%   prob  - struct: .rv0 [1x6] .rvf [1x6] .tNodes [1x(M+1)] .M .P (ztl P)
%   wantJ - (optional) assemble J [default true]
%
% OUTPUTS:
%   R    - residual [14M-7 x 1]
%   J    - Jacobian [14M-7 x 14M-7] (empty if ~wantJ)
%   info - struct: .maxCont (max continuity residual) .termErr (norm of the
%          7 terminal rows) .nSegMax (max regime segments over arcs)
%          .grazed (any arc flagged a graze)
%
% REFERENCES:
%   [1] Stoer & Bulirsch, Intro. to Numerical Analysis, sec. 7.3 (MS BVP).
%   [2] Betts, Practical Methods for Optimal Control, 2010 (MS transcription).
%   [3] ztl_flow.m (per-arc flow + variational STM); Z0_BUILD.md.

if nargin < 3, wantJ = true; end
M = prob.M;  P = prob.P;  tN = prob.tNodes;
n = 14*M - 7;

% --- unpack nodes -----------------------------------------------------------
Y = cell(1, M);
Y{1} = [prob.rv0(:); 1; z(1:7)];
for k = 2:M
    Y{k} = z(7 + 14*(k-2) + (1:14));
end

% --- integrate arcs ---------------------------------------------------------
F = cell(1, M);  Phi = cell(1, M);
nSegMax = 0;  grazed = false;
for k = 1:M
    o = ztl_flow(Y{k}, [tN(k) tN(k+1)], P, wantJ);
    F{k} = o.yf;  Phi{k} = o.PHI;
    nSegMax = max(nSegMax, o.nSegs);
    grazed = grazed || (o.flag == 1);
end

% --- residual ---------------------------------------------------------------
R = zeros(n, 1);
maxCont = 0;
for k = 1:M-1
    rk = F{k} - Y{k+1};
    R(14*(k-1) + (1:14)) = rk;
    maxCont = max(maxCont, norm(rk));
end
termR = [F{M}(1:6) - prob.rvf(:); F{M}(14)];
R(14*(M-1) + (1:7)) = termR;

% --- Jacobian ---------------------------------------------------------------
J = [];
if wantJ
    J = zeros(n, n);
    colLam0 = 1:7;
    colY = @(k) 7 + 14*(k-2) + (1:14);      % valid for k = 2..M
    for k = 1:M-1
        rows = 14*(k-1) + (1:14);
        if k == 1
            J(rows, colLam0) = Phi{1}(:, 8:14);
        else
            J(rows, colY(k)) = Phi{k};
        end
        J(rows, colY(k+1)) = -eye(14);
    end
    rowsT = 14*(M-1) + (1:7);
    J(rowsT, colY(M)) = [Phi{M}(1:6, :); Phi{M}(14, :)];
end

info = struct('maxCont', maxCont, 'termErr', norm(termR), ...
              'nSegMax', nSegMax, 'grazed', grazed);
end
