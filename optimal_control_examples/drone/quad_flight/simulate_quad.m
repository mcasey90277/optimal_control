function [t, X, U] = simulate_quad(x0, ctrl, tEnd, p, opts)
%% Purpose:
%
%   Fly the rigid6dof quadrotor (drone_field) from a library-frame state
%   under a user-given rotor-thrust sequence held zero-order between its
%   breakpoints, with ode113 restarted at every breakpoint so the thrust
%   steps are exact, and the state reported on a uniform grid that includes
%   the breakpoints. The integration stops early, with a warning, if the
%   pitch approaches the Euler singularity (89 degrees); what was flown is
%   returned. Thrusts outside [0, fMax] are flown as given, with one
%   warning: the script is for looking, not for certifying.
%
%% Inputs:
%
%  x0                       [12 x 1]                Initial state in the
%                                                   library frame
%
%  ctrl                     struct                  .t [1 x K] breakpoints
%                                                   (s), non-decreasing,
%                                                   t(1) = 0; .f [4 x K]
%                                                   rotor thrusts (N) held
%                                                   from t(k) to t(k+1)
%
%  tEnd                     scalar                  End time (s)
%
%  p                        struct                  drone_params('rigid6dof')
%
%  opts                     struct (optional)       .dt [0.01] output step
%                                                   (s); .RelTol [1e-10];
%                                                   .AbsTol [1e-12]
%
%% Outputs:
%
%  t                        [1 x N]                 Output times (s)
%
%  X                        [12 x N]                States, library frame
%
%  U                        [4 x N]                 Thrust in force at each
%                                                   time (right-continuous)
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: hover 1 s, then 10 % more thrust for 1 s:
     p  = drone_params('rigid6dof');
     h  = p.m*p.g/4;
     [t, X] = simulate_quad([0; 0; 1; zeros(9, 1)], struct('t', [0 1], 'f', h*[1 1 1 1; 1.1 1.1 1.1 1.1].'), 2, p);
     fprintf('after hover + 10%% thrust for 1 s: z = %.4f m, vz = %.4f m/s at t = %.1f s\n', X(3, end), X(6, end), t(end));
     return
end

if nargin < 5 || isempty(opts), opts = struct(); end
dt     = getf(opts, 'dt', 0.01);
relTol = getf(opts, 'RelTol', 1e-10);
absTol = getf(opts, 'AbsTol', 1e-12);

tb = ctrl.t(:).';
fb = ctrl.f;
if numel(tb) ~= size(fb, 2) || size(fb, 1) ~= 4
    error('simulate_quad:ctrl', 'ctrl.t must be [1 x K] and ctrl.f [4 x K]');
end
if tb(1) ~= 0 || any(diff(tb) < 0)
    error('simulate_quad:ctrl', 'ctrl.t must start at 0 and be non-decreasing');
end
if any(fb(:) < 0) || any(fb(:) > p.fMax)
    warning('simulate_quad:thrust', 'thrusts outside [0, fMax = %.3f N] are flown as given', p.fMax);
end

% Segment boundaries: the breakpoints inside (0, tEnd), plus 0 and tEnd.
edges = unique([0, tb(tb > 0 & tb < tEnd), tEnd]);
grid  = unique([0:dt:tEnd, tEnd]);
thetaStop = 89*pi/180;
odeOpts = odeset('RelTol', relTol, 'AbsTol', absTol, ...
                 'Events', @(tt, x) pitchEvent(x, thetaStop));

t = [];  X = [];  U = [];
x = x0(:);
for s = 1:numel(edges) - 1
    ta = edges(s);  tbnd = edges(s + 1);
    f  = fb(:, find(tb <= ta + 1e-12, 1, 'last'));
    tq = grid(grid > ta + 1e-12 & grid < tbnd - 1e-12);
    tspan = [ta, tq, tbnd];
    [ts, xs, te] = ode113(@(tt, xx) rhs(xx, f, p), tspan, x, odeOpts);
    ts = ts.';  xs = xs.';
    if numel(tspan) == 2 && isempty(te)
        ts = ts([1 end]);  xs = xs(:, [1 end]);     % a 2-point tspan returns every step
    end
    if ~isempty(te)
        warning('simulate_quad:pitch', 'pitch reached %.0f deg at t = %.3f s; stopping there', ...
                thetaStop*180/pi, te(1));
        t = [t, ts];  X = [X, xs];  U = [U, repmat(f, 1, numel(ts))];
        return
    end
    % the segment's first point duplicates the previous end: keep one, and
    % give it this segment's thrust (the hold is right-continuous)
    if ~isempty(t), ts = ts(2:end);  xs = xs(:, 2:end);  U(:, end) = f; end
    t = [t, ts];  X = [X, xs];  U = [U, repmat(f, 1, numel(ts))];
    x = X(:, end);
end
end

% ------------------------------------------------------------------------
function xd = rhs(x, f, p)
%% Purpose:
%
%   The field with the held thrust.
%
[F, B] = drone_field(x, p);
xd = F + B*f;
end

% ------------------------------------------------------------------------
function [value, isterminal, direction] = pitchEvent(x, thetaStop)
%% Purpose:
%
%   Terminal event: |pitch| reaches thetaStop.
%
value      = thetaStop - abs(x(8));
isterminal = 1;
direction  = -1;
end

% ------------------------------------------------------------------------
function v = getf(s, name, default)
%% Purpose:
%
%   Field with a default.
%
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
