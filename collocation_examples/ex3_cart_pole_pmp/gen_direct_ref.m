function gen_direct_ref()
%% Purpose:
%
%   Produce the committed fixture the indirect demo seeds from and compares
%   against: one direct trapezoidal solution of the cart-pole swing-up at a
%   RELAXED force bound, so it approximates the unconstrained optimum the
%   PMP-BVP solves. Run once; the .mat is committed.
%
%   The formulation is the existing example's
%   (ex2_cart_pole_swing_up/try2), re-expressed as a function so the
%   multipliers can be harvested.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none (writes data/cartpole_direct_ref.mat)
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
p  = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
tf = 5;  N = 200;  nN = N + 1;  uMax = 2000;
tN = linspace(0, tf, nN);
h  = tN(2) - tN(1);

% seed: the example's "good initial guess"
q1  = 0.3*sin(2*pi*tN/tf);
q2  = pi*tN/tf;
q1d = (pi/tf)*cos(pi*tN/tf);
q2d = (pi/tf)*ones(1, nN);
u0  = 40*(sin(2*pi*tN/tf) + 0.5*sin(4*pi*tN/tf));
Z0  = [q1, q2, q1d, q2d, u0].';

obj = @(Z) trap_cost(Z, nN, h);
nlc = @(Z) deal([], defects(Z, nN, h, p));
lb  = -inf(5*nN, 1);  ub = inf(5*nN, 1);
lb(4*nN+1:end) = -uMax;  ub(4*nN+1:end) = uMax;

opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
    'MaxIterations', 5000, 'MaxFunctionEvaluations', 1e7, ...
    'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-8, ...
    'SpecifyObjectiveGradient', false);
[Z, J, exitflag, ~, lambda] = fmincon(obj, Z0, [], [], [], [], lb, ub, nlc, opts);
assert(exitflag > 0, 'gen_direct_ref:fmincon', 'fmincon exit flag %d', exitflag);

X = [Z(1:nN).'; Z(nN+1:2*nN).'; Z(2*nN+1:3*nN).'; Z(3*nN+1:4*nN).'];
U = Z(4*nN+1:5*nN).';
% the equality multipliers: 4 boundary rows at each end, then 4 per interval
mu = lambda.eqnonlin(:);
muDefect = reshape(mu(9:end), 4, N);

if ~isfolder(fullfile(here, 'data')), mkdir(fullfile(here, 'data')); end
save(fullfile(here, 'data', 'cartpole_direct_ref.mat'), ...
     'tN', 'X', 'U', 'muDefect', 'J', 'p', 'tf', 'uMax', 'N');
fprintf('gen_direct_ref: J = %.6f, max|u| = %.2f N, terminal miss %.2e\n', ...
        J, max(abs(U)), max(abs(X(:,end) - [0;pi;0;0])));
end

% ------------------------------------------------------------------------
function J = trap_cost(Z, nN, h)
%% Purpose:
%
%   Trapezoidal quadrature of the running cost u^2.
%
u = Z(4*nN+1:5*nN);
J = h*(sum(u.^2) - 0.5*(u(1)^2 + u(end)^2));
end

function ceq = defects(Z, nN, h, p)
%% Purpose:
%
%   Boundary conditions (8 rows) then the trapezoidal defects (4 per
%   interval), in that order -- the order gen_direct_ref harvests by.
%
X = [Z(1:nN).'; Z(nN+1:2*nN).'; Z(2*nN+1:3*nN).'; Z(3*nN+1:4*nN).'];
u = Z(4*nN+1:5*nN).';
F = zeros(4, nN);
for k = 1:nN
    F(:,k) = ref_field(X(:,k), u(k), p);
end
d = X(:,2:end) - X(:,1:end-1) - (h/2)*(F(:,2:end) + F(:,1:end-1));
ceq = [X(:,1) - [0;0;0;0]; X(:,end) - [0;pi;0;0]; d(:)];
end

function dx = ref_field(x, u, p)
%% Purpose:
%
%   The example's dynamics: cart and pendulum accelerations with the force.
%
q2 = x(2);  q1d = x(3);  q2d = x(4);
s = sin(q2);  c = cos(q2);
D1 = p.m1 + p.m2*(1 - c^2);
D2 = p.L*(p.m1 + p.m2)*(1 - (p.m2/(p.m1 + p.m2))*c^2);
dx = [q1d; q2d;
      (p.L*p.m2*s*q2d^2 + u + p.m2*p.g*c*s)/D1;
      (p.L*p.m2*c*s*q2d^2 + u*c + (p.m1 + p.m2)*p.g*s)/D2];
end
