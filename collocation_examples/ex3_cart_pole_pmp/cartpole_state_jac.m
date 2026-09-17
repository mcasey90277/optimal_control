function A = cartpole_state_jac(x, u, p)
%% Purpose:
%
%   The exact 4x4 state Jacobian A = d(F + G u)/dx of the cart-pole field
%   (cartpole_field) at fixed control u, derived symbolically with the
%   Symbolic Math Toolbox (gen_state_jac.m) and written out here by hand as
%   the house-styled wrapper around the generated body.
%
%   Generated, not complex-stepped: the costate equation that consumes this
%   Jacobian sits INSIDE the PMP field, and Task 5's state-transition matrix
%   takes a complex-step derivative THROUGH that field. A complex step
%   inside a complex step shares the imaginary channel and would silently
%   corrupt the inner derivative. Deriving A symbolically breaks that
%   nesting.
%
%  ASSUMPTIONS / NOTES:
%
% • COMPLEX-STEP SAFE, and must stay so: no abs, no norm, no max/min, no
%   branch on a state value. This function itself is called from inside the
%   complex-stepped STM Jacobian (Task 5), so it must survive a
%   complex-perturbed x.
%
%% Inputs:
%
%  x                        [4 x 1]                 [q1; q2; q1dot; q2dot]:
%                                                   cart position (m),
%                                                   pendulum angle (rad, 0 =
%                                                   hanging), and their rates
%
%  u                        scalar                  Control force (N), held
%                                                   fixed for this Jacobian
%
%  p                        struct                  .m1 cart mass (kg), .m2
%                                                   bob mass (kg), .L
%                                                   pendulum length (m), .g
%                                                   gravity (m/s^2)
%
%% Outputs:
%
%  A                        [4 x 4]                 d(F + G u)/dx at (x, u)
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the state Jacobian at the hanging equilibrium, zero control:
     pd = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
     A = cartpole_state_jac([0; 0; 0; 0], 0, pd);
     fprintf('at rest, hanging:  A(3,2) = %.4f,  A(4,2) = %.4f\n', A(3,2), A(4,2));
     return
end

pv = [p.m1; p.m2; p.L; p.g];

L = pv(3,:);
g = pv(4,:);
m1 = pv(1,:);
m2 = pv(2,:);
q2 = x(2,:);
q2d = x(4,:);
t2 = cos(q2);
t3 = sin(q2);
t4 = m1.*2.0;
t5 = q2.*2.0;
t6 = q2d.^2;
t9 = 1.0./L;
t7 = t2.^2;
t8 = t3.^2;
t10 = t7.*2.0;
t11 = t7-1.0;
t12 = m2.*t11;
t13 = t10-1.0;
t14 = -t12;
t15 = m1+t14;
A = reshape([0.0,0.0,0.0,0.0,0.0,0.0,(g.*m2.*t7-g.*m2.*t8+L.*m2.*t2.*t6)./t15-m2.*t2.*t3.*1.0./t15.^2.*(u+L.*m2.*t3.*t6+g.*m2.*t2.*t3).*2.0,(t9.*(t3.*u.*-2.0+g.*m2.*t2.*2.0+g.*t2.*t4+L.*m2.*t6.*t13.*2.0))./(m2+t4-m2.*t13)-m2.*t2.*t3.*t9.*1.0./(m1+m2-m2.*t7).^2.*(t2.*u+g.*m1.*t3+g.*m2.*t3+L.*m2.*t2.*t3.*t6).*2.0,1.0,0.0,0.0,0.0,0.0,1.0,(L.*m2.*q2d.*t3.*2.0)./(m1+m2.*t8),(m2.*q2d.*sin(t5).*2.0)./(m2+t4-m2.*cos(t5))],[4,4]);
end
