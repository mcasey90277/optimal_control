%  cart_pole_swing_up_trap_2 

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

[X_opt, J_opt, exitflag, output] = fmincon(@(X) objective(X,T,num_nodes),X0,A,b,Aeq,beq,lb,ub,@(X) nonlcon_cart_pole(X,T,L,m1,m2,g,num_nodes),options);

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


%% plots

figure;
hold on; grid on;
xlabel('Position [m]'); ylabel('Y [m]');
xlim([-6 6]);ylim([-3 3]);
for k=1:length(q1)
    X0 = [q1(k);0];% cart position
    X1 = X0 + L.*[sin(q2(k));-cos(q2(k))];
    plot(X0(1),X0(2),'rs');
    plot(X1(1),X1(2),'bo','MarkerFaceColor','b',MarkerSize=10);
    plot([X0(1) X1(1)],[X0(2) X1(2)], 'r-')
end


%% make a movie

% create video writer object
vidFile = VideoWriter('cart_pole_swingup.mp4', 'MPEG-4');
vidFile.FrameRate = 30;
open(vidFile);

figure('Position', [100 100 800 400]);
for k = 1:length(q1)
    clf;
    hold on; grid on;
    xlabel('Position [m]'); ylabel('Y [m]');
    title(sprintf('Cart-Pole Swing-Up | t = %.2f s', t_grid(k)));
    xlim([-6 6]); ylim([-3 3]);
    axis equal;

    % cart position
    cart_pos = [q1(k); 0];
    % pendulum bob position
    bob_pos = cart_pos + L.*[sin(q2(k)); -cos(q2(k))];

    % draw cart as a rectangle
    cart_width = 0.8;
    cart_height = 0.4;
    rectangle('Position', [cart_pos(1)-cart_width/2, cart_pos(2)-cart_height/2, cart_width, cart_height], ...
              'FaceColor', [0.3 0.3 0.8], 'EdgeColor', 'k', 'LineWidth', 2);

    % draw pendulum rod
    plot([cart_pos(1) bob_pos(1)], [cart_pos(2) bob_pos(2)], 'r-', 'LineWidth', 3);

    % draw pendulum bob
    plot(bob_pos(1), bob_pos(2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 20);

    % draw ground line
    plot([-6 6], [-cart_height/2 -cart_height/2], 'k-', 'LineWidth', 2);

    drawnow;
    frame = getframe(gcf);
    writeVideo(vidFile, frame);
end

close(vidFile);
disp('Movie saved to cart_pole_swingup.mp4');


