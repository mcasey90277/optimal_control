function [F, B] = drone_field(x, p)
%% Purpose:
%
%   The drone dynamics in control-affine form, xdot = F(x) + B(x) u, for
%   the three models drone_params knows. One home for the physics: the
%   generated Jacobians are differentiated FROM this function, and the
%   oracle (tests/test_drone_physics) pins it by mechanics, not by
%   agreement with a copy.
%
%   kinematic   x = p (3)                     u = v (3):  pdot = v
%   pointmass   x = [p; v] (6)                u = a (3):  vdot = a - g e3
%   rigid6dof   x = [p; v; Theta; omega] (12) u = f (4):
%                 pdot     = v
%                 vdot     = -g e3 + (T/m) R(Theta) e3
%                 Thetadot = W(Theta) omega
%                 omegadot = I^-1 (tau - omega x I omega)
%               T = sum f,  tau = [l (f2 - f4); l (f3 - f1); kappa (f1 - f2 + f3 - f4)]
%               rotor 1 on +x_b, 2 on +y_b, 3 on -x_b, 4 on -y_b ("+" layout)
%
%  ASSUMPTIONS / NOTES:
%
% • World frame: local-level, z up (ENU-compatible: x east, y north, z up);
%   gravity is -g e3. NED callers rotate at their boundary.
% • Conventions: Theta = [phi; theta; psi] is roll/pitch/yaw (rad) with
%   R = Rz(psi) Ry(theta) Rx(phi), body-to-world. omega is the BODY-frame
%   angular rate (rad/s), not Euler rates and not world-frame rates. Body
%   triad: x_b through rotor 1, y_b through rotor 2, z_b along thrust (up
%   at hover). Units: p in m, v in m/s, f in N; pointmass u = a in m/s^2;
%   kinematic u = v in m/s.
% • Complex-step safe: sin/cos/tan, products and divisions only; the cross
%   product is written out; no abs, norm, max, min, no branch on a value.
% • W is singular at theta = +/- pi/2; p.thetaMax is the box the solvers
%   impose. Nothing is clipped here.
%
%% Inputs:
%
%  x                        [nx x 1]                State, layout per model
%
%  p                        struct                  drone_params(model)
%
%% Outputs:
%
%  F                        [nx x 1]                Drift
%
%  B                        [nx x nu]               Control columns
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the rigid body at a 10-degree roll with hover thrust:
     pr = drone_params('rigid6dof');
     x0 = [zeros(6,1); 10*pi/180; 0; 0; zeros(3,1)];
     [F, B] = drone_field(x0, pr);
     xd = F + B*(pr.m*pr.g/4)*ones(4,1);
     fprintf('xdot at 10 deg roll, hover thrust: [%s]\n', sprintf('%.4f ', xd));
     return
end

switch p.model
    case 'kinematic'
        F = zeros(3, 1);
        B = eye(3);
    case 'pointmass'
        v = x(4:6);
        F = [v; 0; 0; -p.g];
        B = [zeros(3); eye(3)];
    case 'rigid6dof'
        v  = x(4:6);
        Th = x(7:9);
        om = x(10:12);
        [R, W] = euler_zyx(Th);
        Iw = p.I.*om;                                  % I omega, I diagonal
        omxIw = [om(2)*Iw(3) - om(3)*Iw(2);
                 om(3)*Iw(1) - om(1)*Iw(3);
                 om(1)*Iw(2) - om(2)*Iw(1)];
        F = [v;
             0; 0; -p.g;
             W*om;
             -omxIw./p.I];
        Mix = [0, p.l, 0, -p.l;
               -p.l, 0, p.l, 0;
               p.kappa, -p.kappa, p.kappa, -p.kappa];
        B = [zeros(3, 4);
             (R(:,3)/p.m)*ones(1, 4);
             zeros(3, 4);
             Mix./p.I];
    otherwise
        error('drone_field:model', 'unknown model ''%s''', p.model);
end
end
