function out = partial_f_test(y, X0, X1)
%% Purpose:
%
%   Nested ordinary-least-squares models and the partial F test: does the
%   larger model X1 explain y better than the smaller model X0 by more than
%   chance? F = ((RSS0 - RSS1)/(p1 - p0)) / (RSS1/(n - p1)), and the
%   upper-tail p-value through the regularized incomplete beta function,
%   p = I_{d2/(d2 + d1 F)}(d2/2, d1/2), so no toolbox is needed.
%   Column spaces are compared, not columns: X0 nested in X1 means every
%   column of X0 lies in the span of X1's columns.
%
%% Inputs:
%
%  y                        [n x 1]                 Response
%
%  X0                       [n x p0]                Smaller model's design
%                                                   matrix (include the
%                                                   intercept column)
%
%  X1                       [n x p1]                Larger model's design
%                                                   matrix, p1 > p0
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
%                                                   .t1 and two-sided .pt1
%                                                   [p1 x 1] (t with df2)
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a 10-point line plus noise; the quadratic term should not help
     x  = (1:10).';
     y  = 2 + 0.5*x + [0.1; -0.2; 0.15; -0.05; 0.1; -0.1; 0.05; -0.15; 0.2; -0.1];
     out = partial_f_test(y, [ones(10, 1), x], [ones(10, 1), x, x.^2]);
     fprintf('F = %.4f  p = %.4f  (the quadratic term is not supported: p > 0.05)\n', out.F, out.p);
     return
end

y  = y(:);
n  = numel(y);
p0 = size(X0, 2);  p1 = size(X1, 2);
if size(X0, 1) ~= n || size(X1, 1) ~= n
    error('orbit_network:partial_f_test:shape', 'X0 and X1 must have numel(y) = %d rows (got %d, %d)', n, size(X0, 1), size(X1, 1));
end
if p1 <= p0
    error('orbit_network:partial_f_test:nested', 'X1 must have more columns than X0 (got %d <= %d)', p1, p0);
end
proj = X1*(X1 \ X0);                                       % X0 projected onto span(X1)
if max(abs(proj(:) - X0(:))) > 1e-8*max(1, max(abs(X0(:))))
    error('orbit_network:partial_f_test:nested', 'X0 is not nested in X1 (its columns are not in the span of X1)');
end

b0 = X0 \ y;  r0 = y - X0*b0;  RSS0 = r0.'*r0;
b1 = X1 \ y;  r1 = y - X1*b1;  RSS1 = r1.'*r1;
df1 = p1 - p0;
df2 = n - p1;
if df2 <= 0
    error('orbit_network:partial_f_test:dof', 'no residual degrees of freedom: n = %d, p1 = %d', n, p1);
end
F = ((RSS0 - RSS1)/df1) / (RSS1/df2);
if RSS1 <= 0
    p = 0;
else
    p = betainc(df2/(df2 + df1*F), df2/2, df1/2);           % upper tail of F(df1, df2)
end
TSS = sum((y - mean(y)).^2);
s2  = RSS1/df2;                                            % residual variance of the larger model
se1 = sqrt(s2*diag(inv(X1.'*X1)));
t1  = b1./se1;
pt1 = betainc(df2./(df2 + t1.^2), df2/2, 0.5);             % two-sided p of each coefficient (t with df2)
out = struct('b0', b0, 'b1', b1, 'RSS0', RSS0, 'RSS1', RSS1, 'df1', df1, 'df2', df2, ...
             'F', F, 'p', p, 'R2_0', 1 - RSS0/TSS, 'R2_1', 1 - RSS1/TSS, ...
             'se1', se1, 't1', t1, 'pt1', pt1);
end
