function ok = test_euler_zyx()
%% Purpose:
%
%   euler_zyx against things that do not come from euler_zyx:
%     1. R is a rotation: R'R = I, det R = +1, R(0) = I;
%     2. R equals the three elementary rotations written out longhand here;
%     3. dRdTh matches a complex step through euler_zyx's own R (this one
%        IS self-referential in the field, but the complex step is an
%        independent differentiation of it);
%     4. W is the inverse of the body-rate map E(Theta) built here from the
%        kinematics of successive rotations:
%            omega = [phidot; 0; 0] + Rx' [0; thetadot; 0] + Rx' Ry' [0; 0; psidot]
%        so E = [1, 0, -sin th; 0, cos ph, sin ph cos th; 0, -sin ph, cos ph cos th]
%        and W E = I at every attitude away from the singularity;
%     5. W(0) = I.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
rngWas = rng(23);  restore = onCleanup(@() rng(rngWas));

worstOrth = 0;  worstDet = 0;  worstLong = 0;  worstD = 0;  worstW = 0;
for trial = 1:50
    Th = (rand(3,1) - 0.5)*2;  Th(2) = Th(2)*1.2;        % |theta| < 1.2 rad, away from 90 deg
    [R, W, dR] = euler_zyx(Th);
    worstOrth = max(worstOrth, max(max(abs(R.'*R - eye(3)))));
    worstDet  = max(worstDet, abs(det(R) - 1));
    ph = Th(1); th = Th(2); ps = Th(3);
    Rx = [1 0 0; 0 cos(ph) -sin(ph); 0 sin(ph) cos(ph)];
    Ry = [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];
    Rz = [cos(ps) -sin(ps) 0; sin(ps) cos(ps) 0; 0 0 1];
    worstLong = max(worstLong, max(max(abs(R - Rz*Ry*Rx))));
    hs = 1e-20;
    for k = 1:3
        Tc = Th;  Tc(k) = Tc(k) + 1i*hs;
        Rc = euler_zyx(Tc);
        worstD = max(worstD, max(max(abs(imag(Rc)/hs - dR{k}))));
    end
    E = [1, 0, -sin(th); 0, cos(ph), sin(ph)*cos(th); 0, -sin(ph), cos(ph)*cos(th)];
    worstW = max(worstW, max(max(abs(W*E - eye(3)))));
end
ok = chk(ok, worstOrth < 1e-14, sprintf('R''R = I to %.1e', worstOrth));
ok = chk(ok, worstDet < 1e-14, sprintf('det R = 1 to %.1e', worstDet));
ok = chk(ok, worstLong < 1e-15, sprintf('R = Rz Ry Rx longhand to %.1e', worstLong));
ok = chk(ok, worstD < 1e-13, sprintf('dR/dTheta vs complex step to %.1e', worstD));
ok = chk(ok, worstW < 1e-13, sprintf('W E = I (E from successive-rotation kinematics) to %.1e', worstW));
[R0, W0] = euler_zyx([0; 0; 0]);
ok = chk(ok, isequal(R0, eye(3)) && isequal(W0, eye(3)), 'R(0) = I and W(0) = I exactly');

if ok, fprintf('TEST_EULER_ZYX: ALL PASS\n'); else, fprintf('TEST_EULER_ZYX: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
