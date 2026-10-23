function files = plot_quad_flight(t, X, U, info, outDir)
%% Purpose:
%
%   Static plots of a quadrotor flight, drawn in the user's frame (ENU or
%   NED), saved as two PNGs:
%
%     <label>_path3d.png   the 3D path with start/end markers and the body
%                          triad drawn every second (x_b red, y_b green,
%                          z_b blue), axes labelled for the frame
%     <label>_states.png   position, velocity, Euler angles, body rates,
%                          rotor thrusts, and altitude above the WGS-84
%                          ellipsoid against time
%
%% Inputs:
%
%  t                        [1 x N]                 Times (s)
%
%  X                        [12 x N]                States, LIBRARY frame
%                                                   (simulate_quad output)
%
%  U                        [4 x N]                 Rotor thrusts (N)
%
%  info                     struct                  .frame 'ENU' | 'NED';
%                                                   .lla0 [lat lon alt]
%                                                   deg/deg/m; .label (file
%                                                   stem)
%
%  outDir                   char                    Folder for the PNGs
%
%% Outputs:
%
%  files                    {1 x 2} cell            Paths written
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a short hover-and-climb flight:
     p = drone_params('rigid6dof');  h = p.m*p.g/4;
     [t, X, U] = simulate_quad([0; 0; 1; zeros(9, 1)], struct('t', [0 1], 'f', h*[1 1 1 1; 1.1 1.1 1.1 1.1].'), 3, p);
     files = plot_quad_flight(t, X, U, struct('frame', 'ENU', 'lla0', [33.77 -118.39 50], 'label', 'demo'), tempdir);
     fprintf('wrote %s and %s\n', files{:});
     return
end

if ~exist(outDir, 'dir'), mkdir(outDir); end
frame = upper(info.frame);
Xu    = quad_frames('from_lib', X, frame);
[lat, lon, alt] = quad_georef(X(1:3, :), info.lla0);
[axName, zDir] = frameLabels(frame);
p = drone_params('rigid6dof');

%% 3D path
f1 = figure('Color', 'w', 'Position', [100 100 900 700]);
ax = axes(f1);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
plot3(ax, Xu(1, :), Xu(2, :), Xu(3, :), '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 1.8);
plot3(ax, Xu(1, 1), Xu(2, 1), Xu(3, 1), 'o', 'MarkerFaceColor', [0.1 0.7 0.2], 'MarkerEdgeColor', 'k', 'MarkerSize', 9);
plot3(ax, Xu(1, end), Xu(2, end), Xu(3, end), 's', 'MarkerFaceColor', [0.85 0.2 0.2], 'MarkerEdgeColor', 'k', 'MarkerSize', 9);
armLen = max(0.12*max(max(Xu(1:3, :), [], 2) - min(Xu(1:3, :), [], 2)), 6*p.l);
kTri = 1:max(1, round(1/max(t(2) - t(1), 1e-9))):numel(t);
for k = kTri
    R = bodyAxesInFrame(X(7:9, k), frame);
    c = Xu(1:3, k);
    cols = {[0.9 0.1 0.1], [0.1 0.6 0.1], [0.1 0.1 0.9]};
    for a = 1:3
        e = c + armLen*R(:, a);
        plot3(ax, [c(1) e(1)], [c(2) e(2)], [c(3) e(3)], '-', 'Color', cols{a}, 'LineWidth', 1.3);
    end
end
xlabel(ax, [axName{1} '  [m]']);  ylabel(ax, [axName{2} '  [m]']);  zlabel(ax, [axName{3} '  [m]']);
set(ax, 'ZDir', zDir{1}, 'YDir', zDir{2});
axis(ax, 'equal');  view(ax, 35, 25);
title(ax, sprintf('%s flight in %s   (origin %.5f, %.5f, %.1f m)   start o, end []', ...
      info.label, frame, info.lla0(1), info.lla0(2), info.lla0(3)), 'Interpreter', 'none');
f1file = fullfile(outDir, [info.label '_path3d.png']);
exportgraphics(f1, f1file, 'Resolution', 110);
close(f1);

%% states
f2 = figure('Color', 'w', 'Position', [100 100 1200 800]);
tl = tiledlayout(f2, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
panel(tl, t, Xu(1:3, :), axName, 'position  [m]');
panel(tl, t, Xu(4:6, :), strcat('v ', axName), 'velocity  [m/s]');
panel(tl, t, Xu(7:9, :)*180/pi, {'roll', 'pitch', 'yaw'}, 'Euler ZYX  [deg]');
panel(tl, t, Xu(10:12, :)*180/pi, {'p', 'q', 'r'}, 'body rates  [deg/s]');
panel(tl, t, U, {'f_1', 'f_2', 'f_3', 'f_4'}, 'rotor thrust  [N]');
ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');
plot(ax, t, alt, 'LineWidth', 1.4);
ylabel(ax, 'altitude WGS-84  [m]');  xlabel(ax, 'time  [s]');
title(ax, sprintf('lat %.6f -> %.6f,  lon %.6f -> %.6f', lat(1), lat(end), lon(1), lon(end)));
title(tl, sprintf('%s: states in %s', info.label, frame), 'Interpreter', 'none');
f2file = fullfile(outDir, [info.label '_states.png']);
exportgraphics(f2, f2file, 'Resolution', 110);
close(f2);

files = {f1file, f2file};
end

% ------------------------------------------------------------------------
function panel(tl, t, Y, names, ylab)
%% Purpose:
%
%   One tile: rows of Y against t with a legend.
%
ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');
for r = 1:size(Y, 1)
    plot(ax, t, Y(r, :), 'LineWidth', 1.4);
end
ylabel(ax, ylab);  xlabel(ax, 'time  [s]');
legend(ax, names, 'Location', 'best');
end

% ------------------------------------------------------------------------
function [axName, dirs] = frameLabels(frame)
%% Purpose:
%
%   Axis names and the axis directions that keep the display right-handed
%   with up on the screen: NED reverses both Z (down) and Y (east), which
%   together are a proper rotation of the view.
%
switch frame
    case 'ENU', axName = {'east', 'north', 'up'};    dirs = {'normal', 'normal'};
    case 'NED', axName = {'north', 'east', 'down'};  dirs = {'reverse', 'reverse'};
end
end

% ------------------------------------------------------------------------
function R = bodyAxesInFrame(ThetaLib, frame)
%% Purpose:
%
%   Columns = the LIBRARY body axes (x_b rotor 1, y_b rotor 2, z_b thrust)
%   expressed in the user's frame, so the triad drawn is the same physical
%   triad in either frame.
%
R = euler_zyx(ThetaLib);
if strcmp(frame, 'NED')
    R = [0 1 0; 1 0 0; 0 0 -1]*R;
end
end
