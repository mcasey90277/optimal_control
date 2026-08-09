function movie_landing(out, sol, P, outfile, opts)
% MOVIE_LANDING  Landing movie: booster + thrust vector | throttle trace.
%
% Left panel: altitude-view booster marker with a thrust-vector arrow
% scaled by throttle, ground line, pad. Right panel: closed-loop throttle
% trace ||Tcmd||/Tmax with a moving time cursor and the three bound lines
% that matter on this campaign -- the engine floor (Tmin/Tmax), the
% GUIDANCE de-rated ceiling (etaT, what sol was solved against), and the
% engine's true ceiling (1.0, what the tracker in out.Tcmd is allowed to
% use for feedback authority -- see booster_params.m's P.etaT note).
%
% Frame size LOCKED to 1280x720 (divisible by 16) via imresize on every
% frame -- an unlocked figure size produces a diagonal-color-streak H.264
% shear artifact in the exported MP4/GIF, not a MATLAB render bug (house
% lesson, see matlab-polished-graphics skill / MEMORY.md).
%
% COLOR SCHEME (fix report, 2026-08-09): the binding colorblind-safe
% blue/orange pair applies here too, not only to the two static plots --
% booster body/marker/throttle-trace use the blue triplet
% [0 0.447 0.741], thrust-plume arrow and the moving time cursor use the
% orange triplet [0.85 0.325 0.098]. Bound lines (Tmin/Tmax, etaT, Tmax)
% and the ground-track path sketch stay neutral black/gray -- they are
% reference geometry, not data series, and adding them to the two-color
% scheme would dilute it. The pad marker is the one deliberate exception:
% it stays a distinct green rather than folding into blue/orange, because
% a landing pad/helipad has its own real-world conventional color
% (aviation safety green, e.g. helipad "H" markings) that a viewer
% already associates with "this is the target," independent of this
% campaign's solver-comparison scheme.
%
% INPUTS:
%   out     - sim_closed_loop trace (.t .X .Tcmd)
%   sol     - guidance solution, either solver (.tf used for the throttle
%             panel's x-axis extent only)
%   P       - booster_params (.Tmax .Tmin .etaT .pad_radius)
%   outfile - output .mp4 path
%   opts    - (optional) .duration [s of playback, def 12], .fps [def 30]
% OUTPUTS: none (writes outfile)
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-08-08-booster-landing-design.md
%   [2] matlab-polished-graphics skill -- 1280x720 divide-by-16 lock
if nargin < 5, opts = struct(); end
if ~isfield(opts,'duration'), opts.duration = 12; end
if ~isfield(opts,'fps'),      opts.fps = 30;      end

colBooster = [0 0.447 0.741];      % colorblind-safe blue   -- booster/throttle
colAccent  = [0.85 0.325 0.098];   % colorblind-safe orange -- plume/cursor

fig = figure('Visible','off','Position',[100 100 1280 720],'Color','w');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
axT = nexttile(tl);  axR = nexttile(tl);

%% Fixed limits computed ONCE (no autoscale jitter):
xmax = max([abs(out.X(:,1)); abs(out.X(:,2)); 50]);
zmax = max(out.X(:,3)) * 1.05;
Tmag = sqrt(sum(out.Tcmd.^2, 2));

nf = round(opts.duration * opts.fps);
ts = interp1(linspace(0,1,numel(out.t)), out.t, linspace(0,1,nf));
vw = VideoWriter(outfile, 'MPEG-4');  vw.FrameRate = opts.fps;  open(vw);
for k = 1:nf
    xk = interp1(out.t, out.X, ts(k)).';
    Tk = interp1(out.t, out.Tcmd, ts(k)).';
    % left: altitude view (downrange = sqrt(x^2+y^2) signed by x), booster
    % marker, thrust arrow DOWN-scaled by throttle, pad + ground line
    cla(axT);  hold(axT, 'on');
    plot(axT, sqrt(sum(out.X(:,1:2).^2,2)).*sign(out.X(:,1)+eps), ...
         out.X(:,3), '-', 'Color', [0.7 0.7 0.7]);
    dk = sqrt(sum(xk(1:2).^2)) * sign(xk(1)+eps);
    arr = -0.15 * zmax * Tk / P.Tmax;    % plume points opposite thrust
    quiver(axT, dk, xk(3), arr(1), arr(3), 0, 'Color', colAccent, ...
           'LineWidth', 2, 'MaxHeadSize', 0.5);
    plot(axT, dk, xk(3), 's', 'Color', colBooster, ...
         'MarkerFaceColor', colBooster, 'MarkerSize', 9);
    plot(axT, [-P.pad_radius P.pad_radius], [0 0], 'g-', 'LineWidth', 4);  % pad: conventional green, see header note
    xlim(axT, [-xmax xmax]*1.1);  ylim(axT, [-0.02*zmax zmax]);
    title(axT, sprintf('t = %6.2f s   alt = %7.1f m', ts(k), xk(3)));
    xlabel(axT, 'downrange [m]');  ylabel(axT, 'altitude [m]');
    % right: throttle trace + cursor, with the three bound lines
    cla(axR);  hold(axR, 'on');
    plot(axR, out.t, Tmag/P.Tmax, '-', 'Color', colBooster, 'LineWidth', 1.5);
    yline(axR, P.Tmin/P.Tmax, 'k:',  'LineWidth', 1.25, ...
          'Label', sprintf('Tmin/Tmax=%.2f', P.Tmin/P.Tmax));
    yline(axR, P.etaT, 'k--', 'LineWidth', 1.25, ...
          'Label', sprintf('etaT=%.2f (guidance)', P.etaT));
    yline(axR, 1, 'k-',  'LineWidth', 1.0, 'Label', 'Tmax (tracker)');
    xline(axR, ts(k), '-', 'Color', colAccent, 'LineWidth', 1.5);
    ylim(axR, [0 1.1]);  xlim(axR, [0 out.t(end)]);
    xlabel(axR, 'time [s]');  ylabel(axR, 'throttle T/Tmax');
    title(axR, 'throttle');
    drawnow;
    frame = getframe(fig);
    frame.cdata = imresize(frame.cdata, [720 1280]);   % divide-by-16 lock
    writeVideo(vw, frame);
end
close(vw);  close(fig);
end
