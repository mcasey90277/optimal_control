function A = drone_state_jac(x, u, p)
%% Purpose:
%
%   The exact state Jacobian A = d(F(x) + B(x) u)/dx at fixed u, for each
%   model. kinematic and pointmass are constant and written here;
%   rigid6dof is the generated drone_state_jac_rigid6dof_gen (from
%   gen_drone_jacs, differentiated from drone_field itself).
%
%% Inputs:
%
%  x                        [nx x 1]                State
%
%  u                        [nu x 1]                Control
%
%  p                        struct                  drone_params(model)
%
%% Outputs:
%
%  A                        [nx x nx]               d(F + B u)/dx
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the rigid body's Jacobian at a hover state:
     pr = drone_params('rigid6dof');
     A = drone_state_jac(zeros(12,1), (pr.m*pr.g/4)*ones(4,1), pr);
     fprintf('A at hover: nnz = %d, eigenvalue magnitudes max %.3e\n', nnz(A), max(abs(eig(A))));
     return
end

switch p.model
    case 'kinematic'
        A = zeros(3);
    case 'pointmass'
        A = [zeros(3), eye(3); zeros(3, 6)];
    case 'rigid6dof'
        A = drone_state_jac_rigid6dof_gen(x(:), u(:), [p.m; p.I(:); p.l; p.kappa; p.g]);
    otherwise
        error('drone_state_jac:model', 'unknown model ''%s''', p.model);
end
end
