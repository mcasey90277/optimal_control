function [c,ceq]=nonlcon_cart_pole(X,T,L,m1,m2,g,num_nodes)

% number of nodes (discrete-time points)
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
ceq(1,1) = q1(1);% cart begins at position 0
ceq(2,1) = q2(1);% pendulum begins straight-down
ceq(3,1) = q1dot(1);% cart initial velocity is 0
ceq(4,1) = q2dot(1);% pendulum initial angular velocity is 0
% final time boundary conditions
ceq(5,1) = q1(end) - 0;% cart is at position d
ceq(6,1) = q2(end) - pi;% pendulum is straight-up
ceq(7,1) = q1dot(end);% cart is stationary
ceq(8,1) = q2dot(end);% pendulum is stationary
% system dynamics defects
ceq1 = q1(2:end) - q1(1:end-1) - h/2.*(q1dot(2:end) + q1dot(1:end-1));% delta cart position should equal integral of velocity between nodes
ceq2 = q2(2:end) - q2(1:end-1) - h/2.*(q2dot(2:end) + q2dot(1:end-1));% delta pendulum angle should equal integral of pendulum angular velocity between nodes

term_kp1_cart_accel = cart_accel(q2(2:end), q2dot(2:end), u(2:end), L, m1, m2, g);
term_k_cart_accel = cart_accel(q2(1:end-1), q2dot(1:end-1), u(1:end-1), L, m1, m2, g);
ceq3 = q1dot(2:end) - q1dot(1:end-1) - h/2.*(term_kp1_cart_accel + term_k_cart_accel);

term_kp1_pend_accel = pendulum_accel(q2(2:end), q2dot(2:end),u(2:end),L,m1,m2,g);
term_k_pend_accel = pendulum_accel(q2(1:end-1), q2dot(1:end-1),u(1:end-1),L,m1,m2,g);
ceq4 = q2dot(2:end) - q2dot(1:end-1) - h/2 .* (term_kp1_pend_accel + term_k_pend_accel);

ceq = [ceq;ceq1;ceq2;ceq3;ceq4];


