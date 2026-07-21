% block_move_2                                                                                                            

%% set-up

number_of_segments = 100;
number_of_nodes = number_of_segments + 1;
t0 = 0; 
tF = 1;
t_grid = linspace(t0,tF, number_of_nodes)';
% initial guess
x0 = t_grid;
v0 = ones(length(t_grid),1);
u0 = zeros(length(t_grid),1);
X0 = [x0;v0;u0];
n_var = length(X0);

%% solve non-linear program

% Options for fmincon
options = optimoptions('fmincon', ...
    'Display', 'iter', ...  % Show progress
    'MaxFunctionEvaluations', 1e5, ...
    'OptimalityTolerance', 1e-6, ...
    'ConstraintTolerance', 1e-6);

% lower and upper bounds
lb = -inf(n_var,1);
ub = inf(n_var,1);

% constraints: linear inequality and equality
A = [];b = [];% no linear inequality constraints
Aeq = []; beq = [];% no linear equality constraints
% note: boundaries handled in non-linear contraints function


[X_opt,J_opt] = fmincon(@J_obj,X0,A,b,Aeq,beq,lb,ub,@nonlcon,options);
%extract state and control
x = X_opt(1:number_of_nodes,1);% position
v = X_opt(number_of_nodes+1:2*number_of_nodes);% velocity
u = X_opt(2*number_of_nodes+1:3*number_of_nodes);% control (= applied force)

%% plots

figure;
plot(t_grid,x,'.-');title('Position');xlabel('Time [sec]');ylabel('Postion');


figure;
plot(t_grid,v,'.-');title('Velocity');xlabel('Time [sec]');ylabel('Velocity');


figure;
plot(t_grid,u,'.-');title('Control');xlabel('Time [sec]');ylabel('Control');





%% functions: objective and non-linear constraints

function J = J_obj(X)
    % number of nodes
    num_of_nodes = length(X)/3;
    %extract state and control
    x = X(1:num_of_nodes,1);% position
    v = X(num_of_nodes+1:2*num_of_nodes);% velocity
    u = X(2*num_of_nodes+1:3*num_of_nodes);% control (= applied force)
    % 
    h = 1/(num_of_nodes-1);
    % compute the objective function
    J = 0;
    for knode = 1:(num_of_nodes-1)
        J = J + 0.5*h*(u(knode)^2 + u(knode+1)^2);
    end
end

function [c, ceq] = nonlcon(X)
    % number of nodes
    num_of_nodes = length(X)/3;
    %extract state and control
    x = X(1:num_of_nodes,1);% position
    v = X((num_of_nodes+1):(2*num_of_nodes));% velocity
    u = X((2*num_of_nodes+1):(3*num_of_nodes));% control (= applied force)
    % 
    h = 1/(num_of_nodes-1);
    % no inequality constriants
    c = [];
    % equality constraints
    % Boundary constraints
    ceq(1) = x(1) - 0;
    ceq(2) = v(1) - 0;
    ceq(3) = x(end)-1;
    ceq(4) = v(end)-0;
    % State dynamic errors
    idx = 5;
    for k = 1:(num_of_nodes-1)
        % derivative of position should be velocity
        ceq(idx) = x(k+1) - x(k) - 1/2 * h * (v(k+1)+v(k));
        % derivative of velocity should be force (=the control)
        idx = idx+1;
        ceq(idx) = v(k+1) - v(k) - 1/2 * h * (u(k+1)+u(k));
        idx = idx + 1;
    end
end
    








