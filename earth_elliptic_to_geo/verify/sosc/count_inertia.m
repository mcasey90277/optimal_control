function [npos, nneg, nzero] = count_inertia(D, zt)
% COUNT_INERTIA  Count positive/negative/zero eigenvalues of a block-diagonal
% LDL^T factor (1x1 and 2x2 blocks), applying a zero-pivot threshold.
%
% INPUTS:
%   D  - block-diagonal factor from ldl(K,'vector'); a 2x2 block at rows
%        i,i+1 is indicated by D(i+1,i) ~= 0 [n x n, may be sparse]
%   zt - zero-pivot threshold (already scaled by the caller) [scalar]
%
% OUTPUTS:
%   npos  - count of eigenvalues classified positive (> zt) [scalar]
%   nneg  - count of eigenvalues classified negative (< -zt) [scalar]
%   nzero - count of eigenvalues classified (near-)zero [scalar]
%
% REFERENCES:
%   [1] process/DESIGN_sosc.md sec 4.5.
%   [2] Nocedal & Wright, "Numerical Optimization," 2nd ed., Thm 16.3.

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
end
