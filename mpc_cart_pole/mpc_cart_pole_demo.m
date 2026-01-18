function mpc_cart_pole_demo()
%% MPC Cart-Pole Controller Demo
% Simple implementation to demonstrate MPC concepts

clear; clc; close all;

%% System Parameters
params = struct();
params.m1 = 5.0;    % Cart mass (kg)
params.m2 = 1.0;    % Pendulum mass (kg)
params.L = 2.0;     % Pendulum length (m)
params.g = 9.81;    % Gravity (m/s^2)

%% MPC Parameters
mpc = struct();
mpc.N = 50;         % Prediction horizon (steps)
mpc.dt = 0.05;      % Time step (20 Hz control)
mpc.u_max = 1000;   % Max force (N)

% Cost weights
mpc.Q = diag([1, 100, 1, 10]);  % [pos, angle, vel_pos, vel_ang]
mpc.R = 0.01;                    % Control weight
mpc.Q_terminal = 2 * mpc.Q;      % Terminal cost weight

%% Simulation Setup
T_sim = 60.0;       % Total simulation time
t_vec = 0:mpc.dt:T_sim;
N_sim = length(t_vec);

% Initial and target states
x0 = [0; 0; 0; 0];              % Start: cart at 0, pendulum down, at rest
x_target = [3; pi; 0; 0];       % Goal: cart at 3m, pendulum up, at rest

%% Initialize Storage
x_history = zeros(4, N_sim);
u_history = zeros(1, N_sim-1);
solve_time = zeros(1, N_sim-1);
cost_history = zeros(1, N_sim-1);

%% Initialize State
x_current = x0;
x_history(:,1) = x_current;

fprintf('Starting MPC Cart-Pole Simulation...\n');
fprintf('Horizon: %d steps, Control rate: %.1f Hz\n', mpc.N, 1/mpc.dt);

%% Main MPC Loop
for k = 1:N_sim-1
    tic;

    %% Solve MPC Problem
    [u_opt, x_pred, cost] = solve_mpc_step(x_current, x_target, params, mpc);

    solve_time(k) = toc;

    %% Apply First Control Action
    u_apply = u_opt(1);
    u_history(k) = u_apply;
    cost_history(k) = cost;

    %% Simulate Real System (with disturbance)
    % Add some random disturbance
    if k > 100  % Add disturbance after 5 seconds
        disturbance = [0; 0; 0; 0.1*randn()];  % Random torque on pendulum
    else
        disturbance = zeros(4,1);
    end

    x_next = simulate_one_step(x_current, u_apply, params, mpc.dt, disturbance);

    %% Update State
    x_current = x_next;
    x_history(:,k+1) = x_current;

    %% Progress Display
    if mod(k, 40) == 0  % Every 2 seconds
        fprintf('t = %.1fs, Cost = %.2f, Solve time = %.1fms\n', ...
                k*mpc.dt, cost, solve_time(k)*1000);
    end
end

%% Results and Visualization
fprintf('\nSimulation Complete!\n');
fprintf('Average solve time: %.2f ms\n', mean(solve_time)*1000);
fprintf('Max solve time: %.2f ms\n', max(solve_time)*1000);

plot_results(t_vec, x_history, u_history, x_target, solve_time, cost_history);
animate_cart_pole(t_vec, x_history, params);

end

%% ========================================================================
%  MPC Optimization Function
%  ========================================================================

function [u_opt, x_pred, cost] = solve_mpc_step(x0, x_ref, params, mpc)
%% Solve single MPC optimization problem

% Decision variables: control sequence u = [u(0), u(1), ..., u(N-1)]
n_vars = mpc.N;

% Initial guess (warm start with previous solution if available)
persistent u_prev
if isempty(u_prev)
    u_guess = zeros(n_vars, 1);
else
    % Shift previous solution and add zero at end
    u_guess = [u_prev(2:end); 0];
end

% Bounds on control
lb = -mpc.u_max * ones(n_vars, 1);
ub =  mpc.u_max * ones(n_vars, 1);

% Optimization options (fast for real-time)
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'off', ...
    'MaxIterations', 50, ...
    'MaxFunctionEvaluations', 500, ...
    'OptimalityTolerance', 1e-4, ...
    'ConstraintTolerance', 1e-4);

% Solve optimization
[u_opt, cost] = fmincon(@(u) mpc_cost(u, x0, x_ref, params, mpc), ...
                        u_guess, [], [], [], [], lb, ub, ...
                        @(u) mpc_constraints(u, x0, params, mpc), options);

% Store for warm start next time
u_prev = u_opt;

% Compute predicted trajectory for visualization
x_pred = predict_trajectory(x0, u_opt, params, mpc);

end

%% ========================================================================
%  Cost Function
%  ========================================================================

function J = mpc_cost(u, x0, x_ref, params, mpc)
%% MPC Cost Function

% Predict trajectory
x_traj = predict_trajectory(x0, u, params, mpc);

J = 0;

% Running cost
for i = 1:mpc.N-1
    % State cost
    e = x_traj(:,i) - x_ref;
    J = J + e' * mpc.Q * e;

    % Control cost
    J = J + mpc.R * u(i)^2;
end

% Terminal cost
e_terminal = x_traj(:,mpc.N) - x_ref;
J = J + e_terminal' * mpc.Q_terminal * e_terminal;

end

%% ========================================================================
%  Constraints Function
%  ========================================================================

function [c, ceq] = mpc_constraints(u, x0, params, mpc)
%% MPC Constraints

% For this simple example, we only have bounds (handled by fmincon)
% Could add state constraints here if needed

c = [];   % No inequality constraints
ceq = []; % No equality constraints (dynamics handled in prediction)

% Example of how to add state constraints:
% x_traj = predict_trajectory(x0, u, params, mpc);
% c = [x_traj(1,:)' - 5;    % position < 5
%     -x_traj(1,:)' - 5];   % position > -5

end

%% ========================================================================
%  Trajectory Prediction
%  ========================================================================

function x_traj = predict_trajectory(x0, u, params, mpc)
%% Predict trajectory using control sequence

x_traj = zeros(4, mpc.N);
x_traj(:,1) = x0;

for i = 1:mpc.N-1
    % Simulate one step forward
    x_traj(:,i+1) = rk4_step(@(t,x) cart_pole_dynamics(t, x, u(i), params), ...
                             0, x_traj(:,i), mpc.dt);
end

end

%% ========================================================================
%  Dynamics Functions
%  ========================================================================

function xdot = cart_pole_dynamics(t, x, u, params)
%% Cart-pole dynamics

q1 = x(1);      % Cart position
q2 = x(2);      % Pendulum angle
q1dot = x(3);   % Cart velocity
q2dot = x(4);   % Pendulum angular velocity

% Extract parameters
m1 = params.m1; m2 = params.m2; L = params.L; g = params.g;

% Compute accelerations
q1ddot = cart_accel(q1, q2, q1dot, q2dot, u, m1, m2, L, g);
q2ddot = pendulum_accel(q1, q2, q1dot, q2dot, u, m1, m2, L, g);

% State derivative
xdot = [q1dot; q2dot; q1ddot; q2ddot];

end

function q1ddot = cart_accel(q1, q2, q1dot, q2dot, u, m1, m2, L, g)
%% Cart acceleration
% Derived from Euler-Lagrange equations for cart-pole system

den = m1 + m2 - m2*cos(q2)^2;
num = u + m2*sin(q2)*(L*q2dot^2 + g*cos(q2));
q1ddot = num / den;

end

function q2ddot = pendulum_accel(q1, q2, q1dot, q2dot, u, m1, m2, L, g)
%% Pendulum acceleration
% Derived from Euler-Lagrange equations for cart-pole system

den = L*(m1 + m2 - m2*cos(q2)^2);
num = -u*cos(q2) - m2*L*q2dot^2*cos(q2)*sin(q2) - (m1+m2)*g*sin(q2);
q2ddot = num / den;

end

%% ========================================================================
%  Simulation and Utility Functions
%  ========================================================================

function x_next = simulate_one_step(x, u, params, dt, disturbance)
%% Simulate one time step with disturbance

if nargin < 5
    disturbance = zeros(4,1);
end

% RK4 integration
x_next = rk4_step(@(t,x) cart_pole_dynamics(t, x, u, params), 0, x, dt);

% Add disturbance
x_next = x_next + disturbance * dt;

end

function x_next = rk4_step(f, t, x, dt)
%% Single RK4 integration step

k1 = f(t, x);
k2 = f(t + dt/2, x + dt*k1/2);
k3 = f(t + dt/2, x + dt*k2/2);
k4 = f(t + dt, x + dt*k3);

x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

end

%% ========================================================================
%  Visualization Functions
%  ========================================================================

function plot_results(t_vec, x_history, u_history, x_target, solve_time, cost_history)
%% Plot simulation results

figure('Position', [100, 100, 1200, 800]);

% State trajectories
subplot(3,2,1);
plot(t_vec, x_history(1,:), 'b-', 'LineWidth', 2);
hold on;
plot(t_vec, x_target(1)*ones(size(t_vec)), 'r--', 'LineWidth', 1);
ylabel('Cart Position (m)');
legend('Actual', 'Target', 'Location', 'best');
grid on;

subplot(3,2,2);
plot(t_vec, x_history(2,:)*180/pi, 'b-', 'LineWidth', 2);
hold on;
plot(t_vec, x_target(2)*180/pi*ones(size(t_vec)), 'r--', 'LineWidth', 1);
ylabel('Pendulum Angle (deg)');
legend('Actual', 'Target', 'Location', 'best');
grid on;

subplot(3,2,3);
plot(t_vec, x_history(3,:), 'b-', 'LineWidth', 2);
ylabel('Cart Velocity (m/s)');
grid on;

subplot(3,2,4);
plot(t_vec, x_history(4,:)*180/pi, 'b-', 'LineWidth', 2);
ylabel('Pendulum Angular Vel (deg/s)');
grid on;

% Control input
subplot(3,2,5);
plot(t_vec(1:end-1), u_history, 'g-', 'LineWidth', 2);
ylabel('Control Force (N)');
xlabel('Time (s)');
grid on;

% Solve time and cost
subplot(3,2,6);
yyaxis left;
plot(t_vec(1:end-1), solve_time*1000, 'r-', 'LineWidth', 1);
ylabel('Solve Time (ms)', 'Color', 'r');

yyaxis right;
plot(t_vec(1:end-1), cost_history, 'b-', 'LineWidth', 1);
ylabel('MPC Cost', 'Color', 'b');
xlabel('Time (s)');
grid on;

sgtitle('MPC Cart-Pole Controller Results');

end

function animate_cart_pole(t_vec, x_history, params)
%% Animate the cart-pole system

figure('Position', [200, 200, 800, 600]);

% Animation parameters
L = params.L;
cart_width = 0.4;
cart_height = 0.2;

for i = 1:5:length(t_vec)  % Every 5th frame for speed
    clf;

    % Current state
    cart_pos = x_history(1,i);
    pend_angle = x_history(2,i);

    % Cart position
    cart_x = [cart_pos - cart_width/2, cart_pos + cart_width/2, ...
              cart_pos + cart_width/2, cart_pos - cart_width/2, cart_pos - cart_width/2];
    cart_y = [-cart_height/2, -cart_height/2, cart_height/2, cart_height/2, -cart_height/2];

    % Pendulum position (angle=0 is down, angle=pi is up)
    pend_x = cart_pos + L * sin(pend_angle);
    pend_y = -L * cos(pend_angle);  % Negative so angle=0 hangs DOWN

    % Draw cart
    fill(cart_x, cart_y, 'b', 'EdgeColor', 'k', 'LineWidth', 2);
    hold on;

    % Draw pendulum
    plot([cart_pos, pend_x], [0, pend_y], 'r-', 'LineWidth', 3);
    plot(pend_x, pend_y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

    % Ground
    plot([-15, 15], [-cart_height/2-0.1, -cart_height/2-0.1], 'k-', 'LineWidth', 2);

    % Target position
    plot([3, 3], [-1, 4], 'g--', 'LineWidth', 2);
    text(3.2, 3.5, 'Target', 'Color', 'g', 'FontSize', 12);

    axis equal;
    xlim([-15, 15]);
    ylim([-3, 5]);
    title(sprintf('MPC Cart-Pole Control - t = %.2f s', t_vec(i)));
    xlabel('Position (m)');
    ylabel('Height (m)');
    grid on;

    drawnow;
    pause(0.01);
end

end
