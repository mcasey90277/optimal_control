function gen_state_jac()
%% Purpose:
%
%   Generate cartpole_state_jac_gen.m: the exact 4x4 d(F + G u)/dx, derived
%   with the Symbolic Math Toolbox and written out with matlabFunction. The
%   generated file is committed, and the public cartpole_state_jac CALLS it,
%   so regenerating is this one command with no hand step in between (a hand
%   transfer of the generated body is where a stale derivative gets in).
%
%   The symbolic dynamics are obtained by calling cartpole_field ITSELF on
%   symbolic inputs, so this file holds NO copy of the physics: the Jacobian
%   is by construction the derivative of the field the solver integrates.
%   (Until 2026-09-17 it held its own copy, which carried the same pendulum
%   sign error as the field -- GPT-6 Astra's review.)
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
%  none (writes cartpole_state_jac_gen.m beside this file)
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  M. Casey  derive from cartpole_field; wrapper calls the gen (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(here);
syms q1 q2 q1d q2d u m1 m2 L g real
xs = [q1; q2; q1d; q2d];
ps = struct('m1', m1, 'm2', m2, 'L', L, 'g', g);
[F, G] = cartpole_field(xs, ps);
A = simplify(jacobian(F + G*u, xs));

out = fullfile(here, 'cartpole_state_jac_gen.m');
matlabFunction(A, 'File', out, 'Vars', {xs, u, [m1; m2; L; g]}, ...
               'Outputs', {'A'}, 'Optimize', true);
fprintf('gen_state_jac: wrote %s\n', out);
end
