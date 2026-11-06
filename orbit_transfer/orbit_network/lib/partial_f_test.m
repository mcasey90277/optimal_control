function out = partial_f_test(y, X0, X1)
%% Purpose:
%
%   Nested ordinary-least-squares models and the partial F test: does the
%   larger model X1 explain y better than the smaller model X0 by more than
%   chance? F = ((RSS0 - RSS1)/(p1 - p0)) / (RSS1/(n - p1)), and the
%   upper-tail p-value through the regularized incomplete beta function,
%   p = I_{d2/(d2 + d1 F)}(d2/2, d1/2), so no toolbox is needed. For the
%   larger model it also returns each coefficient's standard error, t
%   statistic, two-sided p-value (I_{d2/(d2 + t^2)}(d2/2, 1/2)) and 95%
%   confidence interval (t quantile through the inverse incomplete beta).
%   Column spaces are compared, not columns: X0 nested in X1 means every
%   column of X0 lies in the span of X1's columns. The design must be of
%   full column rank (the degrees of freedom are column counts only then).
%   Policies: RSS0 = RSS1 = 0 (both exact) gives F = 0, p = 1; a negative
%   RSS0 - RSS1 within 1e-12 RSS0 of zero is roundoff and is clamped to 0,
%   a larger one is refused (:nested) because a nested model cannot fit
%   worse than its superset.
%
%% Inputs:
%
%  y                        [n x 1]                 Response (finite)
%
%  X0                       [n x p0]                Smaller model's design
%                                                   matrix (include the
%                                                   intercept column)
%
%  X1                       [n x p1]                Larger model's design
%                                                   matrix, p1 > p0, full
%                                                   column rank
%
%% Outputs:
%
%  out                      struct                  .b0 [p0 x 1], .b1 [p1
%                                                   x 1] least-squares
%                                                   coefficients; .RSS0,
%                                                   .RSS1; .df1 = p1 - p0,
%                                                   .df2 = n - p1; .F; .p;
%                                                   .R2_0, .R2_1; and for
%                                                   the larger model the
%                                                   per-coefficient .se1,
%                                                   .t1, two-sided .pt1
%                                                   [p1 x 1] and .ci1 [p1
%                                                   x 2] (95%)
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  M. Casey  rank check, exact-fit policy, roundoff clamp, QR covariance,
%            confidence intervals (Astra review)          (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a 10-point line plus noise; the quadratic term should not help
     x  = (1:10).';
     y  = 2 + 0.5*x + [0.1; -0.2; 0.15; -0.05; 0.1; -0.1; 0.05; -0.15; 0.2; -0.1];
     out = partial_f_test(y, [ones(10, 1), x], [ones(10, 1), x, x.^2]);
     fprintf('F = %.4f  p = %.4f  (the quadratic term is not supported: p > 0.05)\n', out.F, out.p);
     fprintf('quadratic coefficient %.4g, 95%% CI [%.4g, %.4g]\n', out.b1(3), out.ci1(3, 1), out.ci1(3, 2));
     return
end

y  = y(:);
n  = numel(y);
p0 = size(X0, 2);  p1 = size(X1, 2);
if size(X0, 1) ~= n || size(X1, 1) ~= n
    error('orbit_network:partial_f_test:shape', 'X0 and X1 must have numel(y) = %d rows (got %d, %d)', n, size(X0, 1), size(X1, 1));
end
if ~all(isfinite(y)) || ~all(isfinite(X0(:))) || ~all(isfinite(X1(:)))
    error('orbit_network:partial_f_test:finite', 'y, X0 and X1 must be finite');
end
if p1 <= p0
    error('orbit_network:partial_f_test:nested', 'X1 must have more columns than X0 (got %d <= %d)', p1, p0);
end
df1 = p1 - p0;
df2 = n - p1;
if df2 <= 0
    error('orbit_network:partial_f_test:dof', 'no residual degrees of freedom: n = %d, p1 = %d', n, p1);
end
[Q1, R1] = qr(X1, 0);
if rank(R1) < p1 || rank(X0) < p0
    error('orbit_network:partial_f_test:rank', 'designs must have full column rank (rank X1 = %d of %d, rank X0 = %d of %d)', rank(R1), p1, rank(X0), p0);
end
proj = Q1*(Q1.'*X0);                                       % X0 projected onto span(X1)
if max(abs(proj(:) - X0(:))) > 1e-8*max(1, max(abs(X0(:))))
    error('orbit_network:partial_f_test:nested', 'X0 is not nested in X1 (its columns are not in the span of X1)');
end

b0 = X0 \ y;  r0 = y - X0*b0;  RSS0 = r0.'*r0;
b1 = R1 \ (Q1.'*y);  r1 = y - X1*b1;  RSS1 = r1.'*r1;
dRSS = RSS0 - RSS1;
if dRSS < 0
    if dRSS >= -1e-12*max(RSS0, eps)
        dRSS = 0;                                          % roundoff
    else
        error('orbit_network:partial_f_test:nested', 'RSS0 - RSS1 = %.3e < 0: the smaller model fits better, so the models are not nested', dRSS);
    end
end
if RSS1 <= 0 && dRSS <= 0
    F = 0;  p = 1;                                         % both exact: no improvement to claim
elseif RSS1 <= 0
    F = Inf;  p = 0;
else
    F = (dRSS/df1) / (RSS1/df2);
    p = betainc(df2/(df2 + df1*F), df2/2, df1/2);           % upper tail of F(df1, df2)
end
TSS = sum((y - mean(y)).^2);
s2  = RSS1/df2;                                            % residual variance of the larger model
Ri  = R1 \ eye(p1);                                        % (X1'X1)^-1 = Ri Ri'
se1 = sqrt(s2*sum(Ri.^2, 2));
t1  = b1./se1;
pt1 = betainc(df2./(df2 + t1.^2), df2/2, 0.5);             % two-sided p of each coefficient (t with df2)
xq  = betaincinv(0.05, df2/2, 0.5);                        % t_{0.975, df2} = sqrt(df2 (1/xq - 1))
tq  = sqrt(df2*(1/xq - 1));
ci1 = [b1 - tq*se1, b1 + tq*se1];
out = struct('b0', b0, 'b1', b1, 'RSS0', RSS0, 'RSS1', RSS1, 'df1', df1, 'df2', df2, ...
             'F', F, 'p', p, 'R2_0', 1 - RSS0/TSS, 'R2_1', 1 - RSS1/TSS, ...
             'se1', se1, 't1', t1, 'pt1', pt1, 'ci1', ci1);
end
