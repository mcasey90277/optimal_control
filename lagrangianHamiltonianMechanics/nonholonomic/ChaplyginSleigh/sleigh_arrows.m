function arrows = sleigh_arrows(t, Y, out, p)
%% Purpose:
%
%   The two arrows the movie draws, as functions of time for
%   mech.animate_planar_body: the center-of-mass velocity at the center of
%   mass and the blade's constraint force lambda e2 at the blade point.
%   Each is scaled so the largest arrow over the run is one body length,
%   max(|a|, 0.1), and interpolated onto the movie's time with pchip.
%
%% Inputs:
%
%  t                        [1 x N]                 Times (s)
%
%  Y                        [5 x N]                 Reduced solution
%
%  out                      struct                  From sleigh_quantities
%
%  p                        struct                  From sleigh_params
%
%% Outputs:
%
%  arrows                   struct array [1 x 2]    .name, .fun = @(tk) ->
%                                                   [base; vec] world, .color
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the defaults, print both arrows at t = 1
     p = sleigh_params();  t = 0:0.01:3;
     [~, yy] = ode113(@(a, b) sleigh_rhs(a, b, p), t, [0; 0; 0; 1; 2], odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
     out = sleigh_quantities(t, yy.', p);
     arrows = sleigh_arrows(t, yy.', out, p);
     fprintf('%s at t = 1: %s\n', arrows(1).name, mat2str(arrows(1).fun(1).', 4));
     fprintf('%s at t = 1: %s\n', arrows(2).name, mat2str(arrows(2).fun(1).', 4));
     return
end
Lb = max(abs(p.a), 0.1);
th = Y(3, :);
e2 = [-sin(th); cos(th)];
vc = out.vc;
F  = out.lambda.*e2;
sv = Lb/max(max(sqrt(sum(vc.^2, 1))), eps);
sf = Lb/max(max(abs(out.lambda)), eps);
vcI = @(tk) interp1(t, vc.', tk, 'pchip').';
FI  = @(tk) interp1(t, F.', tk, 'pchip').';
rcI = @(tk) interp1(t, out.rc.', tk, 'pchip').';
rbI = @(tk) interp1(t, out.rb.', tk, 'pchip').';
arrows = struct('name', {'velocity', 'constraint force'}, ...
                'fun',  {@(tk) [rcI(tk); sv*vcI(tk)], @(tk) [rbI(tk); sf*FI(tk)]}, ...
                'color', {[1 0.75 0.3], [0.9 0.4 0.9]});
end
