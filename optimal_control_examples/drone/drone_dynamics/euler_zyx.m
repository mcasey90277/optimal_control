function [R, W, dRdTh] = euler_zyx(Theta)
%% Purpose:
%
%   ZYX (yaw-pitch-roll) Euler attitude: the body-to-inertial rotation
%   R = Rz(psi) Ry(theta) Rx(phi), the map W from body rates to Euler rates
%   (Thetadot = W omega), and dR/dTheta. Complex-step safe (sin, cos, tan
%   and division only; no abs, no branch on a value).
%
%  ASSUMPTIONS / NOTES:
%
% • W is singular at theta = +/- pi/2 (division by cos theta). The models
%   that use it carry a thetaMax box; this function does not guard, so a
%   caller sees Inf/NaN rather than a silent clip.
%
%% Inputs:
%
%  Theta                    [3 x 1]                 [phi; theta; psi]
%                                                   roll, pitch, yaw (rad)
%
%% Outputs:
%
%  R                        [3 x 3]                 Body-to-inertial
%
%  W                        [3 x 3]                 Thetadot = W omega
%
%  dRdTh                    {1 x 3}                 dR/dphi, dR/dtheta,
%                                                   dR/dpsi
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a 30-degree roll, 20-degree pitch, 45-degree yaw:
     [R, W] = euler_zyx([30; 20; 45]*pi/180);
     fprintf('R =\n'); disp(R); fprintf('W =\n'); disp(W);
     return
end

ph = Theta(1);  th = Theta(2);  ps = Theta(3);
cph = cos(ph);  sph = sin(ph);
cth = cos(th);  sth = sin(th);
cps = cos(ps);  sps = sin(ps);

Rx = [1 0 0; 0 cph -sph; 0 sph cph];
Ry = [cth 0 sth; 0 1 0; -sth 0 cth];
Rz = [cps -sps 0; sps cps 0; 0 0 1];
R  = Rz*Ry*Rx;

if nargout > 1
    tth = sth/cth;
    W = [1, sph*tth, cph*tth;
         0, cph,     -sph;
         0, sph/cth, cph/cth];
end
if nargout > 2
    dRx = [0 0 0; 0 -sph -cph; 0 cph -sph];
    dRy = [-sth 0 cth; 0 0 0; -cth 0 -sth];
    dRz = [-sps -cps 0; cps -sps 0; 0 0 0];
    dRdTh = {Rz*Ry*dRx, Rz*dRy*Rx, dRz*Ry*Rx};
end
end
