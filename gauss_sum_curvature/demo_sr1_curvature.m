% DEMO_SR1_CURVATURE  SR1 recovers indefinite measurement curvature; BFGS cannot.
%
% Core demonstration for the Gauss-sum / curvature idea:
%   (1) As the filter mean moves, the change in a measurement Jacobian row is a
%       SECANT observation of that component's Hessian:
%             grad h_k(x+s) - grad h_k(x) ~ H_k s.
%       SR1 recovers the (indefinite) H_3 = diag(1,-1) of the saddle exactly in
%       n=2 independent steps -- from Jacobians the EKF already computes.
%   (2) BFGS/DFP CANNOT: their curvature condition s'y>0 is violated in the
%       negative-curvature direction, so they reject the measurement. This is
%       why measurement-curvature estimation wants SR1, not BFGS.
%   (3) From the estimated curvature and the prior covariance, the R-weighted
%       max-(curvature x spread^2) direction is where a Gaussian-sum filter
%       should split; flat directions keep one wide component.
%
% REFERENCES:
%   Hennig & Kiefel, "Quasi-Newton Methods: A New Direction," JMLR 14 (2013).

clear; clc;

H3_true = [1 0; 0 -1];
fprintf('True saddle curvature H3 (indefinite): det = %+.1f\n', det(H3_true));

% ---- (1) SR1 from Jacobian-row differences as the mean moves ----
g3 = @(x) [x(1); -x(2)];          % grad of 3rd component = 3rd row of J(x) = H3*x
B  = zeros(2,2);                  % start with no curvature knowledge
x  = [0.3; -0.2];                 % arbitrary starting mean
steps = {[1;0], [0;1]};           % two independent moves of the operating point
fprintf('\n(1) SR1 curvature recovery (B0 = 0):\n');
for k = 1:numel(steps)
    s = steps{k}; xn = x + s;
    y = g3(xn) - g3(x);           % secant data = change in Jacobian row 3
    [B, d] = sr1_update(B, s, y);
    fprintf('   step %d: s=[%g;%g]  s''y=%+0.2f  denom=%+0.2f  ||B-H3||=%.2e\n', ...
            k, s(1), s(2), s'*y, d, norm(B - H3_true));
    x = xn;
end
disp('   recovered B ='); disp(B);

% ---- (2) Would BFGS accept these? curvature condition s'y > 0 ----
fprintf('(2) BFGS curvature condition s''y>0 on each move:\n');
x = [0.3; -0.2];
for k = 1:numel(steps)
    s = steps{k}; y = g3(x+s) - g3(x);
    if s'*y > 0, verdict = 'accepts';
    else,        verdict = 'REJECTS (cannot represent this curvature)'; end
    fprintf('   s=[%g;%g]: s''y=%+0.2f -> BFGS %s\n', s(1), s(2), s'*y, verdict);
end

% ---- (3) Curvature-driven split direction (R-weighted, prior spread P) ----
R = eye(3);  P = diag([1.0, 0.25]);     % elongated prior
L = chol(P,'lower');
errfun = @(v) (0.5*(v'*H3_true*v))^2;   % R^-1-weighted 2nd-order term (only k=3)
best_e = -inf; best_u = [1;0];
for th = linspace(0, pi, 3601)
    u = [cos(th); sin(th)]; v = L*u;    % direction scaled by prior spread
    e = errfun(v);
    if e > best_e, best_e = e; best_u = u; end
end
fprintf('\n(3) Split direction (max curvature x spread^2, R-metric):\n');
fprintf('   u = [%+0.3f; %+0.3f]  (err = %.4f) -> split the component along this\n', ...
        best_u(1), best_u(2), best_e);
fprintf('   flat directions keep one wide component.\n');
