%% Purpose:
%
%   THE front door for the Chaplygin sleigh. Set the geometry and the
%   initial velocity in section 1, run, and get: the integrated motion
%   (ode113 at 1e-12), five printed diagnostics with PASS/FAIL (hand-
%   derived equations against the symbolic Lagrange-d'Alembert builder,
%   the exact tanh solution, energy, the constraint, the stability
%   statement), two PNG figures and a 1280 x 720 MP4 in results/ beside
%   this script. The derivation and the conventions are in
%   doc/chaplygin_sleigh_note.pdf; the sign of a v is the whole story:
%   straight motion with the center of mass leading the blade in the
%   direction of travel (a v > 0) is stable to yaw -- one of a line of
%   equilibria, transversely attracting; a v < 0 (center of mass trailing)
%   is unstable and the sleigh reverses.
%
%% What the user sets (section 1):
%
%   label, m, J, a, v0, w0, tEnd, dtOut, movieOpts.
%
%% What it writes:
%
%   results/<label>_path.png, results/<label>_diag.png, results/<label>.mp4
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% 1. USER SETTINGS
label = 'sleigh_forward';
m     = 1.0;     % mass (kg)
J     = 0.1;     % inertia about the center of mass (kg m^2)
a     = 0.2;     % blade to center of mass (m); a > 0 = center of mass ahead of the blade
v0    = 1.0;     % initial forward speed along the blade (m/s)
w0    = 2.0;     % initial yaw rate (rad/s); positive = counterclockwise = a left turn; below ~1e-6 rad/s a trailing launch is ill-conditioned at the script's tolerances and D2 can fail
tEnd  = 10;      % final time (s)
dtOut = 0.01;    % output step (s)
movieOpts = struct('fps', 30, 'seconds', tEnd, 'follow', 1.0);   % seconds = the movie's length; frames are streamed, memory does not grow with fps*seconds; follow = camera half-height (m) tracking the sleigh; [] = fixed whole-path camera

%% 2. SETUP
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'mechlib'));
outDir = fullfile(here, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
tiny = 1e-300;        % floor for divisions (a script runs in the base workspace, where a user variable can shadow eps)
assert(tEnd > 0 && dtOut > 0 && dtOut <= tEnd/2, 'run_sleigh:time', 'need 0 < dtOut <= tEnd/2 (at least two output intervals)');
p = sleigh_params(m, J, a);
C = v0^2 + (p.Jb/p.m)*w0^2;
fprintf('sleigh: m %.3g kg, J %.3g kg m^2, a %+.3g m -> Jb %.4g, k %.4g 1/m; C %.6g, v -> %+.6g m/s\n', ...
        m, J, a, p.Jb, p.k, C, sign(a)*sqrt(C));

%% 3. INTEGRATE (blade-point state [x_b; y_b; th; v; w], blade at the origin heading east)
t  = unique([0:dtOut:tEnd, tEnd]);
y0 = [0; 0; 0; v0; w0];
[~, yy] = ode113(@(tt, y) sleigh_rhs(tt, y, p), t, y0, odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
Y   = yy.';
out = sleigh_quantities(t, Y, p);
fprintf('integrated %d steps to %.2f s: v(end) %+.6f m/s, w(end) %+.3e rad/s, path length %.3f m\n', ...
        numel(t) - 1, t(end), Y(4, end), Y(5, end), sum(sqrt(sum(diff(out.rb, 1, 2).^2, 1))));

%% 4. THE SYMBOLIC BUILDER (the independent derivation)
[q, qd, L, A] = sleigh_lagrangian(p);
sys = mech.lagrange_dalembert(q, qd, L, A);
fprintf('builder: reduced accelerations and lambda formed symbolically (n = %d, k = %d)\n', sys.n, sys.k);

%% 5. DIAGNOSTICS
d1 = 0;
for kk = 1:numel(t)
    d1 = max(d1, max(abs(sys.f(t(kk), out.Yfull(:, kk)) - out.Yfulld(:, kk))));
end
verdict('D1 hand-derived vs builder', d1, 1e-9);
d2 = max([abs(Y(4, :) - out.vExact), abs(Y(5, :) - out.wExact)]);
verdict('D2 exact solution (v, w)', d2, 1e-8);
chk = mech.conservation_check(t, out.Yfull, sys);
verdict('D3 energy drift (rel)', chk.driftRel, 1e-9);
verdict('D4 constraint residual |A qd|', chk.resid, 1e-9);
leads = a*Y(4, end) > 0;
settled = abs(Y(5, end)) <= 1e-3*max(max(abs(Y(5, :))), tiny);   % against the run's own peak |w|
S = sqrt(C);  b = sqrt(p.Jb/p.m);
if v0 == 0 && w0 == 0
    fprintf('  NOTE  D5 at rest: nothing moves\n');
elseif a == 0 && w0 ~= 0
    fprintf('  NOTE  D5 a = 0: the knife edge; v and w stay constant, the blade runs a circle of radius %.4g m\n', abs(v0/w0));
elseif w0 == 0 && a*v0 > 0
    fprintf('  NOTE  D5 w0 = 0: a straight ride at v = %.3g m/s; with the center of mass leading (a v > 0) it is stable to yaw perturbations\n', v0);
elseif w0 == 0 && a*v0 < 0
    fprintf('  NOTE  D5 w0 = 0: a straight ride at v = %.3g m/s with the center of mass TRAILING (a v < 0): an exact but unstable equilibrium, any yaw would grow and reverse it\n', v0);
elseif w0 == 0
    fprintf('  NOTE  D5 a = 0, w0 = 0: a straight ride at v = %.3g m/s; every straight ride is an equilibrium of the knife edge\n', v0);
elseif leads && settled && abs((Y(3, end) - Y(3, 1)) - (out.thetaInf - Y(3, 1))) < 1e-2   % 0.6 deg: the settle bar (1e-3 of peak |w|) leaves ~1e-3 rad of heading still to come
    fprintf('  PASS  D5 asymptotics: v -> sign(a) sqrt(C) = %+.6f (v(end) %+.6f), |w(end)| = %.1e, heading change %.4f rad vs exact %.4f; straight motion with the center of mass leading is stable to yaw (one of a line of equilibria: it attracts transversely, on its energy ellipse)\n', ...
            out.asymptote, Y(4, end), abs(Y(5, end)), Y(3, end) - Y(3, 1), out.thetaInf - Y(3, 1));
elseif leads && settled
    fprintf('  NOTE  D5 asymptotics: center of mass leads and w has settled, but the heading change %.4f rad has not converged to the exact total %.4f rad yet (lengthen tEnd)\n', Y(3, end) - Y(3, 1), out.thetaInf - Y(3, 1));
elseif leads
    fprintf('  NOTE  D5 asymptotics: center of mass leads at the end but |w(end)| = %.1e has not settled below 1e-3 of the run''s peak |w| by t = %.3g s (lengthen tEnd)\n', abs(Y(5, end)), tEnd);
else
    tRev = -asinh(v0/(b*abs(w0)))/(p.k*S);
    fprintf('  NOTE  D5 asymptotics: the center of mass still trails (a v < 0) at t = %.3g s; this is the unstable branch, |w| grows until the signed speed crosses zero (at t = %.3g s) and the sleigh then departs along sign(a) (lengthen tEnd to see it)\n', ...
            tEnd, tRev);
end

%% 6. PLOTS AND MOVIE
files = plot_sleigh(t, Y, out, p, struct('label', label, 'visible', false, 'resid', chk.residSeries), outDir);
fprintf('wrote %s\n       %s\n', files{1}, files{2});
body = sleigh_body(p);
body.arrows   = sleigh_arrows(t, Y, out, p);
lamI = @(tk) interp1(t, out.lambda, tk, 'pchip');
vI   = @(tk) interp1(t, Y(4, :), tk, 'pchip');
wI   = @(tk) interp1(t, Y(5, :), tk, 'pchip');
EI   = @(tk) interp1(t, (out.E - out.E(1))/max(abs(out.E(1)), tiny), tk, 'pchip');
body.subtitle = @(tk) sprintf('t =%6.2f s   v =%7.3f m/s   w =%7.3f rad/s   lambda =%7.3f N   dE/E0 = %+.1e', tk, vI(tk), wI(tk), lamI(tk), EI(tk));
mfile = mech.animate_planar_body(t, [out.rb; Y(3, :)], body, fullfile(outDir, label), ...
                                 struct('fps', movieOpts.fps, 'seconds', movieOpts.seconds, 'follow', movieOpts.follow, ...
                                        'traces', cat(3, out.rc, out.rb), 'traceColors', [0.85 0.33 0.1; 0.45 0.70 1.00], ...
                                        'title', sprintf('Chaplygin sleigh: a = %+.3g m, v0 = %.3g m/s, w0 = %.3g rad/s', a, v0, w0)));
fprintf('wrote %s\n', mfile);

% ------------------------------------------------------------------------
function verdict(name, value, bar)
%% Purpose:
%
%   One PASS/FAIL line with the measured value and its bar.
%
if value <= bar, s = 'PASS'; else, s = 'FAIL'; end
fprintf('  %s  %s: %.2e (bar %.0e)\n', s, name, value, bar);
end
