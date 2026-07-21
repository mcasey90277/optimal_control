% block_main
% Block mover example from Matthew Kelly paper

%% intial guess
t0 = 0; tF = 1;
num_of_segments = 10; % number of segments
h = 1/num_of_segments; % time-step
num_of_nodes = num_of_segments + 1;% number of nodes
n_var = 3*num_of_nodes; % control;position;velocity
t_grid = (0:num_of_segments)'.*h; % discrete time
x0 = t_grid; % position
v0 = ones(length(x0),1);% velocity
u0 = zeros(length(x0),1);% control

X0 = [x0;v0;u0];

%% solve non-linear (quadratic) program

% Bounds (none needed, but can add if desired, e.g., lb = -inf, ub = inf)
lb = -inf(n_var, 1);
ub = inf(n_var, 1);

% Linear constraints: none (boundaries handled in nonlinear eq)
A = []; b = []; Aeq = []; beq = [];

% Options for fmincon
options = optimoptions('fmincon', ...
    'Display', 'iter', ...  % Show progress
    'MaxFunctionEvaluations', 1e5, ...
    'OptimalityTolerance', 1e-6, ...
    'ConstraintTolerance', 1e-6);


[X_opt, J_opt] = fmincon(@objective,X0,A,b,Aeq,beq,lb,ub,@constraints,options);



%% optimization functions

function J = objective(Z)
% Extract
n_nodes = (length(Z) / 3);
N = n_nodes - 1;
h = 1 / N;
u = Z(2*n_nodes+1:end);

% Trapezoidal quadrature for int u^2 dt
J = 0;
for k = 1:N
    J = J + (h/2) * (u(k)^2 + u(k+1)^2);
end
end

function [c, ceq] = constraints(Z)
    % Extract
    n_nodes = (length(Z) / 3);
    N = n_nodes - 1;
    h = 1 / N;
    x = Z(1:n_nodes);
    v = Z(n_nodes+1:2*n_nodes);
    u = Z(2*n_nodes+1:end);
    
    % Inequality constraints: none
    c = [];
    
    % Equality constraints: dynamics + boundaries
    ceq = zeros(2*N + 4, 1);  % 2 per segment + 4 boundaries
    
    % Boundaries
    ceq(1) = x(1) - 0;
    ceq(2) = v(1) - 0;
    ceq(3) = x(end) - 1;
    ceq(4) = v(end) - 0;
    
    % Dynamics defects
    idx = 5;
    for k = 1:N
        % x defect
        ceq(idx) = x(k+1) - x(k) - (h/2) * (v(k) + v(k+1));
        idx = idx + 1;
        % v defect
        ceq(idx) = v(k+1) - v(k) - (h/2) * (u(k) + u(k+1));
        idx = idx + 1;
    end
end