function [Xc, Uc, tc, Jc, out] = collocation_transfer(x0, xf, tf, mu, N)
% COLLOCATION_TRANSFER  Min-energy transfer by direct trapezoidal collocation.
%
%   Transcribes the OCP into a nonlinear program: decision variables are the
%   states X(:,k) and controls U(:,k) at N+1 equally spaced nodes; the
%   dynamics become N trapezoidal defect equality constraints
%
%     d_k = x_{k+1} - x_k - (h/2) * ( f(x_k,u_k) + f(x_{k+1},u_{k+1}) ) = 0,
%
%   the cost is the trapezoidal quadrature of 0.5*||u||^2, and the boundary
%   states are pinned with equal lower/upper bounds. Solved with fmincon
%   (interior-point).
%
% INPUTS:
%   x0 - Initial state [4x1]
%   xf - Target state [4x1]
%   tf - Fixed final time [scalar]
%   mu - Gravitational parameter [scalar]
%   N  - Number of intervals (N+1 nodes) [scalar]
%
% OUTPUTS:
%   Xc  - State trajectory at the nodes [4 x N+1]
%   Uc  - Control trajectory at the nodes [2 x N+1]
%   tc  - Node times [1 x N+1]
%   Jc  - Achieved cost [scalar]
%   out - Struct: exitflag, max defect norm after solve, fmincon output,
%         and the equality-constraint Lagrange multipliers (lambda.eqnonlin)
%
% REFERENCES:
%   [1] Betts, "Practical Methods for Optimal Control", 3rd ed., Ch. 3.
%   [2] Kelly, "An Introduction to Trajectory Optimization", SIAM Review 2017.

nx = 4; nu = 2;
tc = linspace(0, tf, N+1);
h  = tf / N;

% --- pack/unpack helpers -----------------------------------------------
% Z = [X(:); U(:)], X is nx x (N+1), U is nu x (N+1)
nX = nx * (N+1);
unpackX = @(Z) reshape(Z(1:nX),      nx, N+1);
unpackU = @(Z) reshape(Z(nX+1:end),  nu, N+1);

% --- initial guess: polar interpolation, zero control --------------------
% A straight line in Cartesian coordinates from x0 to xf can pass through
% (or near) the origin -- where g(r) is singular -- making the constraint
% function undefined at the initial point. Interpolate radius and angle
% instead: a spiral sweep that stays well away from the primary.
r0g = norm(x0(1:2));  th0 = atan2(x0(2), x0(1));
rfg = norm(xf(1:2));  thf = atan2(xf(2), xf(1));
dth = mod(thf - th0, 2*pi);              % counter-clockwise sweep
rg  = r0g + (rfg - r0g) * (tc / tf);     % 1 x N+1
thg = th0 + dth * (tc / tf);
rdotg  = (rfg - r0g) / tf;
thdotg = dth / tf;
Xg = [ rg .* cos(thg);
       rg .* sin(thg);
       rdotg .* cos(thg) - rg .* thdotg .* sin(thg);
       rdotg .* sin(thg) + rg .* thdotg .* cos(thg) ];
Ug = zeros(nu, N+1);
Z0 = [Xg(:); Ug(:)];

% --- bounds: pin boundary states with lb = ub ---------------------------
lb = -inf(size(Z0));  ub = inf(size(Z0));
lb(1:nx)          = x0;  ub(1:nx)          = x0;   % X(:,1)   = x0
lb(nX-nx+1:nX)    = xf;  ub(nX-nx+1:nX)    = xf;   % X(:,N+1) = xf

% --- solve ---------------------------------------------------------------
opts = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
    'Display', 'off', 'MaxFunctionEvaluations', 4e5, 'MaxIterations', 3e3, ...
    'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10, ...
    'StepTolerance', 1e-12);
[Z, Jc, exitflag, output, mult] = fmincon(@cost, Z0, [], [], [], [], ...
                                          lb, ub, @defects, opts);

Xc = unpackX(Z);
Uc = unpackU(Z);
[~, ceq] = defects(Z);
out = struct('exitflag', exitflag, 'max_defect', max(abs(ceq)), ...
             'output', output, 'lam_eq', mult.eqnonlin);

% ======================= nested functions ================================
    function J = cost(Z)
        U = unpackU(Z);
        L = 0.5 * sum(U.^2, 1);              % 1 x N+1 running cost at nodes
        J = trapz(tc, L);
    end

    function [c, ceq] = defects(Z)
        X = unpackX(Z);
        U = unpackU(Z);
        F = zeros(nx, N+1);                  % dynamics at every node
        for k = 1:N+1
            r_vec  = X(1:2, k);
            F(:,k) = [ X(3:4, k);
                       two_body_accel(r_vec, mu) + U(:, k) ];
        end
        D   = X(:, 2:end) - X(:, 1:end-1) - (h/2) * (F(:, 2:end) + F(:, 1:end-1));
        ceq = D(:);                          % 4N equality constraints
        c   = [];                            % no inequalities
    end
end
