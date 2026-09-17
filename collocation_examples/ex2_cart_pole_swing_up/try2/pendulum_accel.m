function pend_acceleration = pendulum_accel(q_angle, q_angular_velocity,u,L,m1,m2,g)

% SIGN FIXED 2026-09-17: the numerator carries a leading MINUS (Kelly, SIAM
% Review 2017, App. E.1). Without it the equations conserve no energy at
% u = 0 and the hanging equilibrium q2 = 0 is UNSTABLE -- found by GPT-6
% Astra's review of ex3_cart_pole_pmp, and pinned since by
% ex3_cart_pole_pmp/tests/test_cartpole_physics (power balance
% dE/dt = u*q1dot from the geometry alone).
numerator_pend_accel = -(L*m2.*cos(q_angle).*sin(q_angle).*q_angular_velocity.^2 + u.*cos(q_angle) + (m1+m2)*g.*sin(q_angle));
denominator_pend_accel = L*(m1+m2).*(1-m2/(m1+m2).*cos(q_angle).^2);
% pendulum accceleration 
pend_acceleration = numerator_pend_accel./denominator_pend_accel;