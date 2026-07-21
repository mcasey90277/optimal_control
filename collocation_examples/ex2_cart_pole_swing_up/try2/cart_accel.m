function cart_acceleration = cart_accel(q_angle, q_angular_velocity,u,L,m1,m2,g)


numerator_cart_accel = L*m2.*sin(q_angle).*(q_angular_velocity.^2) + u + m2*g.*cos(q_angle).*sin(q_angle);
demoninator_cart_accel = m1 + m2.*(1-cos(q_angle).^2);

% cart acceleration (= a function of pendulum angle, pendulum anglular velocity, and control)
cart_acceleration = numerator_cart_accel./demoninator_cart_accel;