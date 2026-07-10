function [lam0, yJ] = sms_unpack(Z, prob)
% SMS_UNPACK  Recover per-arc initial 16-dim states from the unknown vector.
%
% Node 1 has known state r,v = rv0, m = 1, t = 0; its 8 costates are the
% first unknowns. Interior joints carry the full 16-dim augmented state.
%
% INPUTS:
%   Z    - unknown vector [(16M-8)x1] (see SMS_PACK)
%   prob - problem struct from SMS_PROBLEM; prob.sJ must be set [1x(M+1)]
%
% OUTPUTS:
%   lam0 - initial costates [8x1]
%   yJ   - initial augmented state of each arc [16xM]; column 1 is
%          [prob.rv0; prob.m0; 0; lam0], columns 2..M come from Z

M    = numel(prob.sJ) - 1;
Z    = Z(:);
lam0 = Z(1:8);
yJ   = zeros(16, M, 'like', Z);          % 'like' keeps complex-step alive
yJ(:, 1) = [prob.rv0; prob.m0; 0; lam0];
if M > 1
    yJ(:, 2:M) = reshape(Z(9:end), 16, M-1);
end
end
