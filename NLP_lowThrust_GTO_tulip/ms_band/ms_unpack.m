function [lam0, yJ] = ms_unpack(Z, prob)
% MS_UNPACK  Recover per-arc initial states from the unknown vector.
%
% INPUTS:
%   Z    - unknown vector [(14M-7)x1] (see MS_PACK)
%   prob - problem struct from MS_PROBLEM; prob.tJ must be set [1x(M+1)]
%
% OUTPUTS:
%   lam0 - initial costates [7x1]
%   yJ   - initial augmented state of each arc [14xM]; column 1 is
%          [prob.rv0; prob.m0; lam0], columns 2..M come from Z

M    = numel(prob.tJ) - 1;
Z    = Z(:);
lam0 = Z(1:7);
yJ   = zeros(14, M, 'like', Z);          % 'like' keeps complex-step alive
yJ(:, 1) = [prob.rv0; prob.m0; lam0];
if M > 1
    yJ(:, 2:M) = reshape(Z(8:end), 14, M-1);
end
end
