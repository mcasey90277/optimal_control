function pend_acceleration = pendulum_accel(q_angle, q_angular_velocity,u,L,m1,m2,g)

numerator_pend_accel = (L*m2.*cos(q_angle).*sin(q_angle).*q_angular_velocity.^2 + u.*cos(q_angle) + (m1+m2)*g.*sin(q_angle));
denominator_pend_accel = L*(m1+m2).*(1-m2/(m1+m2).*cos(q_angle).^2);
% pendulum accceleration 
pend_acceleration = numerator_pend_accel./denominator_pend_accel;