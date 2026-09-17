function gen_state_jac()
%% Purpose:
%
%   Generate cartpole_state_jac.m: the exact 4x4 d(F + G u)/dx, derived with
%   the Symbolic Math Toolbox and written out with matlabFunction. Run once;
%   the generated file is committed.
%
%   It is generated rather than taken by complex step because the costate
%   equation that uses it sits INSIDE the PMP field, and the STM's Jacobian
%   complex-steps through that field. A complex step within a complex step
%   shares the imaginary channel and silently corrupts the inner derivative.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none (writes cartpole_state_jac.m beside this file)
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
syms q1 q2 q1d q2d u m1 m2 L g real
s = sin(q2);  c = cos(q2);
D1 = m1 + m2*(1 - c^2);
D2 = L*(m1 + m2)*(1 - (m2/(m1 + m2))*c^2);
f = [q1d;
     q2d;
     (L*m2*s*q2d^2 + u + m2*g*c*s)/D1;
     (L*m2*c*s*q2d^2 + u*c + (m1 + m2)*g*s)/D2];
A = simplify(jacobian(f, [q1; q2; q1d; q2d]));

out = fullfile(here, 'cartpole_state_jac_raw.m');
matlabFunction(A, 'File', out, 'Vars', {[q1; q2; q1d; q2d], u, [m1; m2; L; g]}, ...
               'Outputs', {'A'}, 'Optimize', true);
fprintf('gen_state_jac: wrote %s -- wrap it as cartpole_state_jac(x, u, p)\n', out);
end
