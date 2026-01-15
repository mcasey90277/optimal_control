% cart_pole_swing_up_trap_1                                                                                                   

%% set-up
N = 200;% number of segments
num_nodes = N + 1;% number of discrete time points
d = 7.5;% meters
t0 = 0;% seconds, initial time
T = 5.0;% seconds, final time
L = 2;% length of pendulum, meters
m1 = 5;% mass of cart, kg
m2 = 1;% mass of pendulum bob, kg
g = 9.8;% gravitational acceleration, m/s^2
u_max = 2000;% Newtons, max control force

%% discretize time
t_grid = linspace(t0,T,num_nodes)';

%% initial guess

% horrid initial guess
% u0 = zeros(length(t_grid),1);% control profile initial guess
% q10 = zeros(length(t_grid),1);% cart poistion initial guess
% q20 = zeros(length(t_grid),1);% pendulum angle initial guess
% q1dot0 = zeros(length(t_grid),1);% cart velocity initial guess
% q2dot0 = zeros(length(t_grid),1);% cart angular velocity initial guess
% X0 = [q10;q20;q1dot0;q2dot0;u0];% initial state [cart position; pendulum angle; control profile]

%  possibly acceptable initial guess
% q10 = t_grid/T .* d/2;% poisiton goes from 0 to d
% q20 = t_grid/T .* pi;% angle goes from 0 to pi
% q1dot0 = d/T .* ones(num_nodes,1);% constant cart velocity
% q2dot0 = pi/T .* ones(num_nodes,1);% constant angular velocity
% u0 = zeros(length(t_grid),1);% control profile initial guess
% X0 = [q10;q20;q1dot0;q2dot0;u0];% initial state [cart position; pendulum angle; control profile]

% good initial guess
q10    = 0.3 * sin(2*pi*t_grid/T);% small oscillations
q20 = t_grid/T .* pi;% angle goes from 0 to pi
q1dot0 = (pi/T) * cos(pi*t_grid/T);% constant cart velocity
q2dot0 = pi/T .* ones(num_nodes,1);% constant angular velocity
u0 = 40 * (sin(2*pi*t_grid/T) + 0.5*sin(4*pi*t_grid/T));  % fundamental + harmonic
X0 = [q10;q20;q1dot0;q2dot0;u0];% initial state [cart position; pendulum angle; control profile]

%% NLP 

% Options for fmincon
% options = optimoptions('fmincon', ...
%     'Display', 'iter', ...  % Show progress
%     'MaxFunctionEvaluations', 1e5, ...
%     'OptimalityTolerance', 1e-6, ...
%     'ConstraintTolerance', 1e-6);

options = optimoptions('fmincon', ...
    'Algorithm',             'sqp', ...
    'EnableFeasibilityMode', true, ...
    'SubproblemAlgorithm',   'cg', ...
    'Display',               'iter', ...       % helpful to watch progress
    'ConstraintTolerance',   1e-6, ...
    'MaxIterations',         1e7);            % increase if needed

% lower and upper bounds on X
lbpos = -d.*ones(num_nodes,1);
lbang = -inf(size(lbpos));
lbcartvel = -inf(size(lbpos));
lbangvel = -inf(size(lbpos));
lbcontrol = -u_max.*ones(size(lbpos));
lb = [lbpos;lbang;lbcartvel;lbangvel;lbcontrol];

ubpos = d.*ones(num_nodes,1);
ubang = inf(size(lbpos));
ubcartvel = inf(size(lbpos));
ubangvel = inf(size(lbpos));
ubcontrol = u_max.*ones(size(lbpos));
ub = [ubpos;ubang;ubcartvel;ubangvel;ubcontrol];

% Inequality constraints: none
A = []; b=[];
% equality constraints: handled in non-linear constraint function
Aeq = []; beq = [];

[X_opt, J_opt, exitflag, output] = fmincon(@(X) objective(X,T),X0,A,b,Aeq,beq,lb,ub,@(X) nonlcon(X,T,L,m1,m2,g,d),options);

disp('=== Cart-Pole Results ===');
disp(['Final cost J = ' num2str(J_opt, '%.4f')]);
disp(['Exit flag = ' num2str(exitflag)]);
disp(['Iterations = ' num2str(output.iterations)]);
disp(['Func evals = ' num2str(output.funcCount)]);
disp(['Max constr violation = ' num2str(output.constrviolation, '%.2e')]);

% extract cart position, pendulum angle, and control
q1 = X_opt(1:num_nodes);% cart position
q2 = X_opt(num_nodes+1:2*num_nodes);% pendulum angle
q1dot = X_opt(2*num_nodes+1:3*num_nodes);% cart velocity
q2dot = X_opt(3*num_nodes+1:4*num_nodes);% pendulum angular velocity
u = X_opt(4*num_nodes+1:end);% control (= applied force) profile

%% warm start

Xwarm = X_opt;
u_max = 40;
% lower and upper bounds on X
lbpos = -d.*ones(num_nodes,1);
lbang = -inf(size(lbpos));
lbcartvel = -inf(size(lbpos));
lbangvel = -inf(size(lbpos));
lbcontrol = -u_max.*ones(size(lbpos));
lb = [lbpos;lbang;lbcartvel;lbangvel;lbcontrol];

ubpos = d.*ones(num_nodes,1);
ubang = inf(size(lbpos));
ubcartvel = inf(size(lbpos));
ubangvel = inf(size(lbpos));
ubcontrol = u_max.*ones(size(lbpos));
ub = [ubpos;ubang;ubcartvel;ubangvel;ubcontrol];

[X_opt2, J_opt2, exitflag2, output2] = fmincon(@(X) objective(X,T),Xwarm,A,b,Aeq,beq,lb,ub,@(X) nonlcon(X,T,L,m1,m2,g,d),options);

disp('=== Cart-Pole Results ===');
disp(['Final cost J2 = ' num2str(J_opt2, '%.4f')]);
disp(['Exit flag = ' num2str(exitflag2)]);
disp(['Iterations = ' num2str(output2.iterations)]);
disp(['Func evals = ' num2str(output2.funcCount)]);
disp(['Max constr violation = ' num2str(output2.constrviolation, '%.2e')]);

% extract cart position, pendulum angle, and control
q1warm = X_opt2(1:num_nodes);% cart position
q2warm = X_opt2(num_nodes+1:2*num_nodes);% pendulum angle
q1dotwarm = X_opt2(2*num_nodes+1:3*num_nodes);% cart velocity
q2dotwarm = X_opt2(3*num_nodes+1:4*num_nodes);% pendulum angular velocity
uwarm = X_opt2(4*num_nodes+1:end);% control (= applied force) profile


%% NLP functions

function J = objective(X,T)
    % number of nodes (discrete-time points)
    num_nodes = length(X)/5;
    num_of_segments = num_nodes - 1;
    % length of each segment
    h = T/num_of_segments;
    % extract cart position, pendulum angle, and control
    q1 = X(1:num_nodes);% cart position
    q2 = X(num_nodes+1:2*num_nodes);% pendulum angle
    q1dot = X(2*num_nodes+1:3*num_nodes);% cart velocity
    q2dot = X(3*num_nodes+1:4*num_nodes);% pendulum angular velocity
    u = X(4*num_nodes+1:end);% control (= applied force) profile
    % integrate to get total force effort
    J = h/2 .* sum (u(2:end).^2 + u(1:end-1).^2,1);% trapezoidal approximation
end

function [c,ceq] = nonlcon(X,T,L,m1,m2,g,d)
    % number of nodes (discrete-time points)
    num_nodes = length(X)/5;
    num_of_segments = num_nodes - 1;
    % length of each segment
    h = T/num_of_segments;
    % extract cart position, pendulum angle, and control
    q1 = X(1:num_nodes);% cart position
    q2 = X(num_nodes+1:2*num_nodes);% pendulum angle
    q1dot = X(2*num_nodes+1:3*num_nodes);% cart velocity
    q2dot = X(3*num_nodes+1:4*num_nodes);% pendulum angular velocity
    u = X(4*num_nodes+1:end);% control (= applied force) profile
    % inequality constraints: none
    c = [];
    % equality constraints
    % initial time boundary constraints
    ceq(1) = q1(1);% cart begins at position 0
    ceq(2) = q2(1);% pendulum begins straight-down
    ceq(3) = q1dot(1);% cart initial velocity is 0
    ceq(4) = q2dot(1);% pendulum initial angular velocity is 0
    % final time boundary conditions
    ceq(5) = q1(end) - 0;% cart is at position d
    ceq(6) = q2(end) - pi;% pendulum is straight-up
    ceq(7) = q1dot(end);% cart is stationary
    ceq(8) = q2dot(end);% pendulum is stationary
    % system dynamics defects
    ceq1 = q1(2:end) - q1(1:end-1) - h/2.*(q1dot(2:end) + q1dot(1:end-1));
    ceq2 = q2(2:end) - q2(1:end-1) - h/2.*(q2dot(2:end) + q2dot(1:end-1));

    term_kp1_cart_accel = L*m2.*sin(q2(2:end)).*q2dot(2:end).^2 + u(2:end) + m2*g.*cos(q2(2:end)).*sin(q2(2:end));
    term_kp1_cart_accel = term_kp1_cart_accel./(m1 + m2.*(1-cos(q2(2:end)).^2));
    term_k_cart_accel = L*m2.*sin(q2(1:end-1)).*q2dot(1:end-1).^2 + u(1:end-1) + m2*g.*cos(q2(1:end-1)).*sin(q2(1:end-1));
    term_k_cart_accel = term_k_cart_accel./(m1 + m2.*(1-cos(q2(1:end-1)).^2));
    ceq3 = q1dot(2:end) - q1dot(1:end-1) - h/2.*(term_kp1_cart_accel + term_k_cart_accel);

    term_kp1_pend_accel = (L*m2.*cos(q2(2:end)).*sin(q2(2:end)).*q2dot(2:end).^2 + u(2:end).*cos(q2(2:end)) + (m1+m2)*g.*sin(q2(2:end)));
    term_kp1_pend_accel = term_kp1_pend_accel./(L*(m1+m2).*(1-m2/(m1+m2).*cos(q2(2:end)).^2));
    term_k_pend_accel = (L*m2.*cos(q2(1:end-1)).*sin(q2(1:end-1)).*q2dot(1:end-1).^2 + u(1:end-1).*cos(q2(1:end-1)) + (m1+m2)*g.*sin(q2(1:end-1)));
    term_k_pend_accel = term_k_pend_accel./(L*(m1+m2).*(1-m2/(m1+m2).*cos(q2(1:end-1)).^2));
    ceq4 = q2dot(2:end) - q2dot(1:end-1) - h/2 .* (term_kp1_pend_accel + term_k_pend_accel);

    ceq = [ceq';ceq1;ceq2;ceq3;ceq4];
end


%% plots

% cart position
figure;
plot(t_grid, q1,'.-'); xlabel('Time [s]'); ylabel('Cart Position [m]');
title('Cart Position Profile');

figure;
plot(t_grid, q1dot,'.-'); xlabel('Time [s]'); ylabel('Cart Velocity [m]');
title('Cart Velocity Profile');


% pendulum angle
figure;
plot(t_grid, q2.*180/pi,'.-'); xlabel('Time [s]'); ylabel('Pendulum Angle [deg]');
title('Pendulum Angle Profile');

figure;
plot(t_grid, q2dot.*180/pi,'.-'); xlabel('Time [s]'); ylabel('Pendulum Angular Velocity [deg]');
title('Pendulum Angular Velocity Profile');


% pendulum angle
figure;
plot(t_grid, u,'.-'); xlabel('Time [s]'); ylabel('Control [N]');
title('Control Applied Force Profile');



% cart position
figure;
plot(t_grid, q1warm,'.-'); xlabel('Time [s]'); ylabel('Cart Position [m]');
title('Warm Cart Position Profile');

figure;
plot(t_grid, q1dotwarm,'.-'); xlabel('Time [s]'); ylabel('Cart Velocity [m]');
title('Warm Cart Velocity Profile');


% pendulum angle
figure;
plot(t_grid, q2warm.*180/pi,'.-'); xlabel('Time [s]'); ylabel('Pendulum Angle [deg]');
title('Warm Pendulum Angle Profile');

figure;
plot(t_grid, q2dotwarm.*180/pi,'.-'); xlabel('Time [s]'); ylabel('Pendulum Angular Velocity [deg]');
title('Warm Pendulum Angular Velocity Profile');


% pendulum angle
figure;
plot(t_grid, uwarm,'.-'); xlabel('Time [s]'); ylabel('Control [N]');
title('Warm Control Applied Force Profile');





















