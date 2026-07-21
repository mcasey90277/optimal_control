function [R, J] = ms_residual(Z, prob)
% MS_RESIDUAL  Multiple-shooting residual (+ Jacobian) for smoothed min-fuel PMP.
%
% Unknowns Z = [lambda0(7); y_2(14); ...; y_M(14)] with y_k the augmented
% state [r;v;m;lambda_r;lambda_v;lambda_m] at interior arc joint k. Arc k
% integrates lt_pmp_eom_minfuel over [tJ(k), tJ(k+1)]. Residual stacks the
% 14(M-1) joint continuity defects and the 7 terminal conditions
% [rv(tf)-rvf; lambda_m(tf)] (fixed tf, free final mass). Square system.
%
% INPUTS:
%   Z    - unknowns [(14M-7)x1]
%   prob - struct from MS_PROBLEM with tJ set [1x(M+1)]
%
% OUTPUTS:
%   R - residual [(14M-7)x1]
%   J - (optional) sparse d R / d Z [(14M-7)x(14M-7)], block bidiagonal,
%       per-arc complex-step STMs (Task 4)
%
% REFERENCES:
%   [1] Zhang, Topputo, Bernelli-Zazzera, Zhao, JGCD 38(8), 2015.
%   [2] Martins, Sturdza, Alonso, ACM TOMS 29(3), 2003 (complex step).

M = numel(prob.tJ) - 1;
if M < 2, error('ms_residual:M', 'need M >= 2 arcs'); end
[~, yJ] = ms_unpack(Z, prob);

yEnd = zeros(14, M, 'like', Z);
for k = 1:M
    yEnd(:, k) = propagate_arc(yJ(:, k), prob.tJ(k), prob.tJ(k+1), prob);
end

R = zeros(14*M - 7, 1, 'like', Z);
for k = 1:M-1
    R(14*(k-1)+(1:14)) = yEnd(:, k) - yJ(:, k+1);
end
yf = yEnd(:, M);
R(14*(M-1)+(1:7)) = [yf(1:6) - prob.rvf; yf(14)];

if nargout > 1
    J = ms_jacobian_cs(yJ, prob, M);    % Task 4
end
end

% -------------------------------------------------------------------------
function ye = propagate_arc(y0, t0, t1, prob)
% One arc of the smoothed PMP dynamics (complex-step safe end state).
[~, Y] = ode113(@(t, y) lt_pmp_eom_minfuel(t, y, prob.Tmax, prob.c, ...
                prob.muStar, prob.epsSmooth), [t0 t1], y0, prob.odeOpts);
ye = Y(end, :).';
end
