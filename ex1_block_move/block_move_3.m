% block_move_3                                                                                                     

%% set-up and time discretization
number_of_segments = 100;
number_of_nodes = number_of_segments + 1;
t0 = 0; tF = 1;
t_grid = linspace(t0, tF, number_of_nodes)';
h = 1/number_of_segments;% distance between evenly spaced nodes

%% intial guess
x0 = t_grid;
v0 = ones(length(x0),1);
u0 = zeros(length(x0),1);
X0 = [u0;x0;v0];

%% solve nlp

% Options for fmincon
options = optimoptions('fmincon', ...
    'Display', 'iter', ...  % Show progress
    'MaxFunctionEvaluations', 1e5, ...
    'OptimalityTolerance', 1e-6, ...
    'ConstraintTolerance', 1e-6);

% lower and upper bounds on X: none
lb = []; ub = [];
% Ineqaulity constraints: none
A = []; b=[];
% equality constraints: handled in non-linear constraint function
Aeq = []; beq = [];

[X_opt, J_opt] = fmincon(@objective,X0,A,b,Aeq,beq,lb,ub,@nonlcon,options);

%% nlp functions

function J = objective(X)
    num_of_nodes = length(X)/3;
    % extract control profile
    u = X(1:num_of_nodes);
    % length of each time-segment
    h = 1/(num_of_nodes-1);
    % compute the objective function using trapezoidal approximation
    u_sq = u.^2;
    u1 = u_sq(1:num_of_nodes-1);
    u2 = u_sq(2:num_of_nodes);
    u_3 = u1+u2;
    J = h/2 * sum(u_3,1);
end

function [c,ceq] = nonlcon(X)
    num_of_nodes = length(X)/3;
    h = 1/(num_of_nodes-1);
    % extract control and state (position and velocity)
    u = X(1:num_of_nodes);% control profile
    x = X(num_of_nodes+1:2*num_of_nodes);% position
    v = X(2*num_of_nodes+1:end);% velocity
    % inequality constraints
    c = [];
    % system dynamics (defects)
    xkp1 = x(2:end);
    xk = x(1:end-1);
    vkp1 = v(2:end);
    vk = v(1:end-1);
    ukp1 = u(2:end);
    uk = u(1:end-1);
    ceq_pos = xkp1 - xk - h/2.*(vkp1 + vk);
    ceq_vel = vkp1 - vk - h/2.*(ukp1 + uk);
    % boundary conditions
    ceq_boundary = [x(1);v(1);x(end)-1;v(end)];
    % set equality constraint vector
    ceq = [ceq_pos;ceq_vel;ceq_boundary];
end

%% plots

% extract control and state (position and velocity)
u = X_opt(1:number_of_nodes);% control profile
x = X_opt(number_of_nodes+1:2*number_of_nodes);% position
v = X_opt(2*number_of_nodes+1:end);% velocity


figure;
plot(t_grid, u,'.-');xlabel('Time [s]'); ylabel('Control (=force)');title('Control Profile');
figure;
plot(t_grid, x,'.-');xlabel('Time [s]'); ylabel('Position [m]');title('Position Profile');
figure;
plot(t_grid, v,'.-');xlabel('Time [s]'); ylabel('Velocity [m/s]');title('Velocity Profile');
















