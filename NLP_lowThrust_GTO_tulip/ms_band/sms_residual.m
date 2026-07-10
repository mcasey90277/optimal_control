function [R, J] = sms_residual(Z, prob)
% SMS_RESIDUAL  Sundman-domain MS residual (+ Jacobian), 16-dim system.
%
% Unknowns Z = [lam0(8); Y_2(16); ...; Y_M(16)] with Y_k the augmented
% state [r;v;m;t;lamR;lamV;lamM;lamT] at interior joint k. Arc k integrates
% SMS_EOM over sigma in [sJ(k), sJ(k+1)] (sigf FIXED). Residual stacks the
% 16(M-1) joint continuity defects and the 8 terminal conditions
% [rv(sigf) - rvf (6); lamM(sigf) (1); t(sigf) - tf (1)] (fixed tf enforced
% on the carried time state; free final mass). Square system, 16M-8.
%
% INPUTS:
%   Z    - unknowns [(16M-8)x1]
%   prob - struct from SMS_PROBLEM with sJ set [1x(M+1)]
%
% OUTPUTS:
%   R - residual [(16M-8)x1]
%   J - (optional) sparse d R / d Z [(16M-8)x(16M-8)], block bidiagonal,
%       per-arc complex-step STMs (SMS_JACOBIAN_CS)
%
% REFERENCES:
%   [1] Zhang, Topputo, Bernelli-Zazzera, Zhao, JGCD 38(8), 2015.
%   [2] Martins, Sturdza, Alonso, ACM TOMS 29(3), 2003 (complex step).
%   [3] .superpowers/sdd/task-S1-brief.md (MS structure).

M = numel(prob.sJ) - 1;
if M < 2, error('sms_residual:M', 'need M >= 2 arcs'); end
[~, yJ] = sms_unpack(Z, prob);

yEnd = zeros(16, M, 'like', Z);
for k = 1:M
    yEnd(:, k) = propagate_arc(yJ(:, k), prob.sJ(k), prob.sJ(k+1), prob);
end

R = zeros(16*M - 8, 1, 'like', Z);
for k = 1:M-1
    R(16*(k-1)+(1:16)) = yEnd(:, k) - yJ(:, k+1);
end
yf = yEnd(:, M);
R(16*(M-1)+(1:8)) = [yf(1:6) - prob.rvf; yf(15); yf(8) - prob.tf];

if nargout > 1
    J = sms_jacobian_cs(yJ, prob, M);
end
end

% -------------------------------------------------------------------------
function ye = propagate_arc(y0, s0, s1, prob)
% One arc of the Sundman-domain PMP dynamics (complex-step safe end state).
[~, Y] = ode113(@(s, y) sms_eom(s, y, prob.Tmax, prob.c, prob.muStar, ...
                prob.epsSmooth, prob.pSund), [s0 s1], y0, prob.odeOpts);
ye = Y(end, :).';
end
