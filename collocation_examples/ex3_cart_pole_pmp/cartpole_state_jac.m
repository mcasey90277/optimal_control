function A = cartpole_state_jac(x, u, p)
%% Purpose:
%
%   The exact 4x4 state Jacobian A = d(F + G u)/dx of the cart-pole dynamics
%   at FIXED control u -- what the costate equation lamdot = -A' lam needs.
%
%   The expression lives in cartpole_state_jac_gen.m, which gen_state_jac
%   derives symbolically FROM cartpole_field and writes with matlabFunction;
%   this wrapper only unpacks the parameter struct. Regenerate with
%   gen_state_jac whenever cartpole_field changes --
%   tests/test_cartpole_state_jac fails if the two drift apart.
%
%  ASSUMPTIONS / NOTES:
%
% • GENERATED, NOT COMPLEX-STEPPED, on purpose: this function is called from
%   inside the PMP field, and the propagator takes a complex step THROUGH
%   that field to build the STM. A complex step within a complex step shares
%   the imaginary channel and silently corrupts the inner derivative.
% • Complex-step clean itself: sines, cosines, powers and division only, so
%   it survives complex-perturbed x and u on that outer path.
%
%% Inputs:
%
%  x                        [4 x 1]                 [q1; q2; q1dot; q2dot]
%
%  u                        double                  Control force (N), held
%                                                   fixed in the derivative
%
%  p                        struct                  .m1 .m2 .L .g
%
%% Outputs:
%
%  A                        [4 x 4]                 d(F + G u)/dx
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  M. Casey  body moved to the generated file it now calls     (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the pendulum mode at both equilibria, hanging then upright:
     pd = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
     e0 = eig(cartpole_state_jac([0; 0; 0; 0], 0, pd));
     ep = eig(cartpole_state_jac([0; pi; 0; 0], 0, pd));
     fprintf('hanging: eig = [%s]\n', sprintf('%.4f%+.4fi ', [real(e0) imag(e0)].'));
     fprintf('upright: eig = [%s]\n', sprintf('%.4f%+.4fi ', [real(ep) imag(ep)].'));
     return
end

A = cartpole_state_jac_gen(x(:), u, [p.m1; p.m2; p.L; p.g]);
end
