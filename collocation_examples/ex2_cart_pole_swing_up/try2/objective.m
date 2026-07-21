function J = objective(X,T,num_nodes)

% number of nodes (discrete-time points)
num_of_segments = num_nodes - 1;
% length of each segment
h = T/num_of_segments;% assumes discrete-time nodes are evenly spaced
% extract cart position, pendulum angle, and control
q1 = X(1:num_nodes);% cart position
q2 = X(num_nodes+1:2*num_nodes);% pendulum angle
q1dot = X(2*num_nodes+1:3*num_nodes);% cart velocity
q2dot = X(3*num_nodes+1:4*num_nodes);% pendulum angular velocity
u = X(4*num_nodes+1:end);% control (= applied force) profile
% integrate to get total force effort
J = h/2 .* sum (u(2:end).^2 + u(1:end-1).^2,1);% trapezoidal approximation