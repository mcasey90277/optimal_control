function Y = quad_frames(mode, X, frame)
%% Purpose:
%
%   The frame layer between a user's ENU or NED description of a quadrotor
%   state and the library's world frame (drone_field: local-level, z up,
%   ENU-compatible, body x_b through rotor 1, z_b along thrust).
%
%   ENU   the library frame itself. Attitude = ZYX Euler angles of the
%         library body relative to ENU (hover = zero angles, yaw from east
%         counter-clockwise). Identity map.
%
%   NED   the aerospace convention: position/velocity in north-east-down,
%         attitude = ZYX Euler angles of an FRD body (x forward, y right,
%         z down) relative to NED (hover = zero angles, yaw from north
%         clockwise), body rates [p; q; r] of the FRD body. The FRD body is
%         the library body rotated 180 degrees about x_b (D = diag(1,-1,-1)),
%         and NED is ENU with the axes permuted and z flipped (P). So
%
%             p_lib = P p_ned,   v_lib = P v_ned,
%             R_lib = P R_ned D  (body-to-world, a proper rotation),
%             omega_lib = D omega_frd.
%
%         Consequences worth knowing: NED hover with zero angles is the
%         library's hover with x_b pointing north (library yaw 90 degrees);
%         NED yaw 90 degrees (heading east) is library yaw 0.
%
%   The control (four rotor thrusts) is frame-free and is not touched.
%
%% Inputs:
%
%  mode                     char                    'to_lib' | 'from_lib'
%
%  X                        [12 x N]                States, one per column,
%                                                   [p; v; Theta; omega] in
%                                                   the source frame
%
%  frame                    char                    'ENU' | 'NED'
%
%% Outputs:
%
%  Y                        [12 x N]                States in the target
%                                                   frame
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: NED hover at 10 m altitude heading north, in library terms:
     Y = quad_frames('to_lib', [0; 0; -10; zeros(9, 1)], 'NED');
     fprintf('NED hover at 10 m, heading north -> library p = [%s], Theta = [%s] deg\n', ...
             sprintf('%.0f ', Y(1:3)), sprintf('%.0f ', Y(7:9)*180/pi));
     return
end

if size(X, 1) ~= 12
    error('quad_frames:shape', 'X must be [12 x N], got [%d x %d]', size(X, 1), size(X, 2));
end
if ~any(strcmp(mode, {'to_lib', 'from_lib'}))
    error('quad_frames:mode', 'mode must be ''to_lib'' or ''from_lib'', got ''%s''', mode);
end
switch upper(frame)
    case 'ENU'
        Y = X;
        return
    case 'NED'
        P = [0 1 0; 1 0 0; 0 0 -1];                  % NED <-> ENU, self-inverse
        D = diag([1 -1 -1]);                         % FRD <-> library body, self-inverse
    otherwise
        error('quad_frames:frame', 'frame must be ''ENU'' or ''NED'', got ''%s''', frame);
end

N = size(X, 2);
Y = zeros(12, N);
Y(1:3, :) = P*X(1:3, :);
Y(4:6, :) = P*X(4:6, :);
% R_lib = P R_ned D and R_ned = P R_lib D: P and D are their own inverses,
% so the same product serves both directions.
for k = 1:N
    R = euler_zyx(X(7:9, k));
    Y(7:9, k) = euler_from_R(P*R*D);
end
Y(10:12, :) = D*X(10:12, :);
end

% ------------------------------------------------------------------------
function Theta = euler_from_R(R)
%% Purpose:
%
%   ZYX Euler angles [phi; theta; psi] of R = Rz(psi) Ry(theta) Rx(phi):
%   R(3,1) = -sin(theta), R(3,2) = cos(theta) sin(phi), R(3,3) =
%   cos(theta) cos(phi), R(2,1) = cos(theta) sin(psi), R(1,1) =
%   cos(theta) cos(psi). Pitch is taken in (-pi/2, pi/2), the branch the
%   library's W is regular on.
%
Theta = [atan2(R(3, 2), R(3, 3));
         -asin(max(-1, min(1, R(3, 1))));
         atan2(R(2, 1), R(1, 1))];
end
