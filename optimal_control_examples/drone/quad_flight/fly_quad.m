%% FLY_QUAD  Fly a quadrotor from a user-set initial state under a user-set
%  control sequence; make a movie and static plots.
%
%  Purpose:
%
%   The first thing to do with the drone library before any optimal
%   control: set where the quadcopter starts, in ENU or NED at a
%   user-defined latitude/longitude/altitude origin, give it a rotor-thrust
%   sequence (held zero-order between breakpoints), fly it through
%   drone_field with ode113, and look at the result: a 1280 x 720 MP4, a 3D
%   path plot with the body triad, the state histories, and a CSV of the
%   georeferenced path. Everything a user changes is in section 1.
%
%   Frames (quad_frames): ENU is the library's own frame (z up, attitude =
%   ZYX Euler of the library body: x_b through rotor 1, z_b along thrust;
%   hover = zero angles, yaw from east). NED is the aerospace convention
%   (north-east-down, attitude = ZYX Euler of an FRD body, hover = zero
%   angles, yaw from north clockwise, body rates [p q r] of the FRD body).
%
%   Control: rotor thrusts f = [f1 f2 f3 f4] in N, rotor 1 on +x_b, 2 on
%   +y_b, 3 on -x_b, 4 on -y_b. Hover is m g / 4 each (0.0662 N for the
%   Crazyflie-class constants). tau = [l (f2 - f4); l (f3 - f1); kappa (f1 -
%   f2 + f3 - f4)], so f2 > f4 rolls positive about x_b, f3 > f1 pitches
%   positive about y_b, f1 + f3 > f2 + f4 yaws positive about z_b.
%
%  Outputs (in results/ beside this script):
%
%   <label>.mp4            the movie
%   <label>_path3d.png     3D path with body triads
%   <label>_states.png     state and thrust histories, altitude
%   <label>_georef.csv     t, position in the chosen frame, lat, lon, alt
%
%  Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
p = drone_params('rigid6dof');
hover = p.m*p.g/4;

%% 1. USER SETTINGS ------------------------------------------------------
label = 'quad_flight';
frame = 'ENU';                        % 'ENU' (east, north, up) or 'NED' (north, east, down)
lla0  = [33.7716, -118.3944, 60];     % origin: latitude deg, longitude deg, altitude m (WGS-84)

pos0  = [0; 0; 2];                    % m, in `frame` (NED: a positive z is BELOW the origin)
vel0  = [0; 0; 0];                    % m/s, in `frame`
att0  = [0; 0; 0]*pi/180;             % roll, pitch, yaw (rad) in `frame`'s convention
rate0 = [0; 0; 0];                    % body rates (rad/s)
tEnd  = 8;                            % s

% The control sequence: breakpoint times (s) and the four rotor thrusts (N)
% held from each breakpoint to the next. Mind the scale: the Crazyflie's
% inertia is so small that a 2 % thrust differential gives ~8.5 rad/s^2, so
% attitude is set by SHORT paired pulses (a pulse then its mirror leaves a
% tilt of tau/I * duration^2 and no rate). This one climbs and brakes,
% pitches 10 degrees nose-down to run forward, levels, rolls 10 degrees to
% drift sideways, levels, and ends with a touch of extra thrust to arrest
% the sink the tilts caused.
ctrl.t = [0,    1.0,  1.5,  2.0,  2.1,  2.2,  3.7,  3.8,  3.9,  4.0,  4.1,  4.2,  5.5,  5.6,  5.7];
ctrl.f = hover*[ ...
    1.00 1.00 1.00 1.00;              % 0.0  hover
    1.15 1.15 1.15 1.15;              % 1.0  climb
    0.85 0.85 0.85 0.85;              % 1.5  brake: vertical speed back to zero
    0.98 1.00 1.02 1.00;              % 2.0  f3 > f1: +pitch about y_b (nose x_b tips down)
    1.02 1.00 0.98 1.00;              % 2.1  mirror pulse: rate back to zero, 10 deg nose-down stays
    1.00 1.00 1.00 1.00;              % 2.2  hover thrust, leaning: runs forward at g tan(10 deg)
    1.02 1.00 0.98 1.00;              % 3.7  -pitch pulse
    0.98 1.00 1.02 1.00;              % 3.8  mirror: level again, keeps its forward speed
    1.00 1.00 1.00 1.00;              % 3.9  hover
    1.00 1.02 1.00 0.98;              % 4.0  f2 > f4: +roll about x_b
    1.00 0.98 1.00 1.02;              % 4.1  mirror: 10 deg roll stays
    1.00 1.00 1.00 1.00;              % 4.2  leaning sideways: drifts along -y_b
    1.00 0.98 1.00 1.02;              % 5.5  -roll pulse
    1.00 1.02 1.00 0.98;              % 5.6  mirror: level
    1.02 1.02 1.02 1.02].';           % 5.7  a touch more thrust to arrest the sink

movieOpts = struct('fps', 30, 'seconds', tEnd);   % real time
dtOut     = 0.01;                                 % output grid (s)
%% -------------------------------------------------------------------------

%% 2. fly
xUser = [pos0; vel0; att0; rate0];
x0    = quad_frames('to_lib', xUser, frame);
[t, X, U] = simulate_quad(x0, ctrl, tEnd, p, struct('dt', dtOut));
Xu = quad_frames('from_lib', X, frame);
[lat, lon, alt] = quad_georef(X(1:3, :), lla0);

%% 3. report
outDir = fullfile(here, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
info = struct('frame', frame, 'lla0', lla0, 'label', label);
fprintf('\n%s: %s frame, origin (%.5f, %.5f, %.1f m), %.1f s flown, %d samples\n', ...
        label, frame, lla0(1), lla0(2), lla0(3), t(end), numel(t));
fprintf('  start  pos [%8.3f %8.3f %8.3f] m   lat %.6f lon %.6f alt %.2f m\n', Xu(1:3, 1), lat(1), lon(1), alt(1));
fprintf('  end    pos [%8.3f %8.3f %8.3f] m   lat %.6f lon %.6f alt %.2f m\n', Xu(1:3, end), lat(end), lon(end), alt(end));
fprintf('  end    vel [%8.3f %8.3f %8.3f] m/s  roll/pitch/yaw [%6.1f %6.1f %6.1f] deg\n', Xu(4:6, end), Xu(7:9, end)*180/pi);
fprintf('  speed max %.2f m/s, altitude %.2f .. %.2f m, |pitch| max %.1f deg, |roll| max %.1f deg\n', ...
        max(sqrt(sum(Xu(4:6, :).^2))), min(alt), max(alt), max(abs(Xu(8, :)))*180/pi, max(abs(Xu(7, :)))*180/pi);

csvFile = fullfile(outDir, [label '_georef.csv']);
fid = fopen(csvFile, 'w');
cols = frameNames(frame);
fprintf(fid, 't_s,%s_m,%s_m,%s_m,lat_deg,lon_deg,alt_m\n', cols{:});
fprintf(fid, '%.4f,%.6f,%.6f,%.6f,%.9f,%.9f,%.4f\n', [t; Xu(1:3, :); lat; lon; alt]);
fclose(fid);

files = plot_quad_flight(t, X, U, info, outDir);
mov   = movie_quad(t, X, info, fullfile(outDir, label), movieOpts);
fprintf('  wrote %s\n        %s\n        %s\n        %s\n', files{1}, files{2}, mov, csvFile);

%% -------------------------------------------------------------------------
function c = frameNames(frame)
%% Purpose:
%
%   CSV column names for the chosen frame.
%
switch upper(frame)
    case 'ENU', c = {'east', 'north', 'up'};
    case 'NED', c = {'north', 'east', 'down'};
end
end
