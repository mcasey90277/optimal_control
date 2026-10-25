function body = sleigh_body(p)
%% Purpose:
%
%   The sleigh's drawing in BLADE-point body coordinates for
%   mech.animate_planar_body: a pointed hull whose length scales with
%   max(|a|, 0.1) and reaches from behind the blade to ahead of the center
%   of mass, the blade as a thick segment along e1 through the blade
%   point, and two markers (center of mass at (a, 0), blade at the
%   origin). The reference point the animator carries is the blade.
%
%% Inputs:
%
%  p                        struct                  From sleigh_params
%
%% Outputs:
%
%  body                     struct                  .outline [2 x 5];
%                                                   .markers (com, blade);
%                                                   .segments (the blade)
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: draw the defaults
     body = sleigh_body(sleigh_params());
     figure;  O = [body.outline, body.outline(:, 1)];
     plot(O(1, :), O(2, :), 'k-', body.segments(1).pts(1, :), body.segments(1).pts(2, :), 'c-', 'LineWidth', 2);
     hold on;  plot(body.markers(1).pos(1), body.markers(1).pos(2), 'ro', body.markers(2).pos(1), body.markers(2).pos(2), 'bo');
     axis equal;  grid on;  title('sleigh body, blade at the origin, +x forward');
     return
end
Lb = max(abs(p.a), 0.1);
xf = max(p.a, 0) + 0.5*Lb;             % nose
xr = min(p.a, 0) - 0.5*Lb;             % tail
hw = 0.3*Lb;                           % half width
body.outline  = [xf, xf - 0.3*Lb, xr, xr, xf - 0.3*Lb; 0, hw, 0.7*hw, -0.7*hw, -hw];
body.markers  = struct('name', {'com', 'blade'}, 'pos', {[p.a; 0], [0; 0]}, ...
                       'color', {[1 0.45 0.2], [0.5 0.9 1]});
body.segments = struct('pts', [-0.25*Lb, 0.25*Lb; 0, 0], 'color', [0.5 0.9 1], 'width', 4);
end
