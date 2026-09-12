function out = movie_phase_sweep(opts)
%% Purpose:
%
%   MOVIE OF THE PHASE SWEEP: the study script's 3D plot, one frame per
%   (departure phase, arrival phase) pair, so the whole trajectory is seen
%   to deform as the endpoints walk around their orbits.
%
%   The frames ARE transfer_study's picture -- this calls plot_transfer_3d,
%   it does not reimplement it -- with the four things that must not move
%   between frames pinned afterwards: the axis box, the view angles, the
%   colour axis, and the title's field widths. A movie whose axes or colour
%   scale re-fit every frame shows the fitting, not the physics.
%
%   THE COLOUR AXIS is fixed at [0, max t_f] over the SELECTED pairs, not
%   per frame. The arc is coloured by elapsed time, so with one scale for
%   the whole movie a colour means the same number of days in every frame,
%   and a short transfer visibly stops partway up the bar: the arc's end
%   colour reads as its duration.
%
%   A PAIR WITH NO CERTIFIED ENTRY still gets a frame -- both orbits, the
%   same axes, the same title layout with the numbers dashed out, and a
%   note. The sweep's rhythm stays uniform and the library's gaps become
%   visible, which is itself worth seeing (two whole arrival columns are
%   empty).
%
%  ASSUMPTIONS / NOTES:
%
% • Frames come from the catalog's CERTIFIED costates, flown once each
%   through the shared flight function -- about a second a frame. Nothing is
%   re-solved: these entries are already certified, and re-solving would
%   return the same numbers an order of magnitude slower.
% • Frame size is forced to exactly 1280 x 720. A size not divisible by 16
%   makes H.264 shear the picture into diagonal colour streaks, which reads
%   as a render bug and is not one (memory: matlab-movie-diagonal-streaks).
% • The R2026a Aerospace Toolbox puts planetEphemeris on the path without
%   its data package, which blocks pumpkyn's analytic fallback; it is hidden
%   for the render and the path restored on exit (as movie_70mN does).
%
%% Inputs:
%
%  opts                     struct (optional)
%   PHASES -- give VALUES or INDICES, never both for one axis:
%   .sD    'all' | [values]   departure phases (grid: k/12)
%   .sA    'all' | [values]   arrival phases (grid: 0.0754 + j/12)
%   .sDIdx 'all' | [1..12]    departure phases by grid index
%   .sAIdx 'all' | [1..12]    arrival phases by grid index
%   (default: sD = the first grid value, sA = 'all')
%   .order [1]                1 = outer sD, inner sA; 2 = outer sA, inner sD
%   .view  [-35 22]           [az el], held for every frame
%   .slow [1]                 hold every frame this many times longer: 2 is
%                             half speed, 0.5 is double speed. It scales BOTH
%                             outputs together -- the video frame rate is
%                             divided by it and the gif delay multiplied --
%                             so the mp4 and the gif always run at the same
%                             pace. Prefer this to .fps when all you want is
%                             a slower movie: .fps is the base rate, .slow is
%                             how long each frame is held at that rate.
%   .fps [6] .gifDelay [1/fps] .outStem '' .catMat '' (default the shipped
%   70 mN catalog) .dark [true] .nThrust [40] .pad [0.06] axis padding
%   .snapTol [1e-6] how close a requested value must be to a grid value
%
%% Outputs:
%
%  out                      struct                  .mp4 .gif .nFrames
%                                                   .nCertified .nGaps
%                                                   .clim .axLim .pairs
%                                                   [n x 2] (sD, sA)
%                                                   .tfDays [n x 1] (NaN at
%                                                   a gap) .titles {n x 1}
%                                                   (the two title lines of
%                                                   each frame -- returned so
%                                                   the CONSTANT-WIDTH claim
%                                                   is testable) .slow .fps
%                                                   (effective) .gifDelay
%                                                   .runSec (length of the
%                                                   movie) .wallSec
%
%% Revision History:
%  M. Casey                                                   (c) 09/12/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));

order = d('order', 1);
assert(isscalar(order) && ismember(order, [1 2]), 'movie_phase_sweep:order', ...
       'opts.order must be 1 (outer sD, inner sA) or 2 (outer sA, inner sD)');
viewAng = d('view', [-35 22]);
fps = d('fps', 6);
slow = d('slow', 1);
assert(isscalar(slow) && isreal(slow) && isfinite(slow) && slow > 0, ...
       'movie_phase_sweep:slow', 'opts.slow must be a positive finite scalar (2 = half speed)');
gifDelay = d('gifDelay', 1/fps) * slow;     % both outputs scale together
fpsEff = fps / slow;
pad = d('pad', 0.06);  snapTol = d('snapTol', 1e-6);
outStem = d('outStem', fullfile(here, 'results', 'phase_sweep_movie'));

%% The catalog and its grid:
catMat = d('catMat', fullfile(here, 'results', 'costate_catalog_dro_tulip_70mN.mat'));
assert(isfile(catMat), 'movie_phase_sweep:catalog', 'catalog not found: %s', catMat);
L = load(catMat);  fn = fieldnames(L);  cat_ = L.(fn{1});  sh = cat_.sheets(1);
gridD = sh.sD_frac(:).';  gridA = sh.sA_frac(:).';
[B, ~] = arclength_arrival('setup');
pr = B.problem;
phys = struct('Tnd', B.Tnd, 'cnd', B.cnd, 'muStar', B.mu, ...
              'lStar', pr.lStar, 'tStar', pr.tStar, 'm0kg', pr.m0kg);

%% Which pairs, in which order:
iD = selectAxis(opts, 'sD', 'sDIdx', gridD, snapTol, 'departure', 1);
iA = selectAxis(opts, 'sA', 'sAIdx', gridA, snapTol, 'arrival', []);
% written as the two loops it is, because a meshgrid here is one transpose
% away from silently sweeping the other axis
pairIdx = zeros(numel(iD)*numel(iA), 2);  k = 0;
if order == 1                               % outer sD, inner sA
    for a = iD, for b = iA, k = k + 1;  pairIdx(k,:) = [a, b];  end, end
else                                        % outer sA, inner sD
    for b = iA, for a = iD, k = k + 1;  pairIdx(k,:) = [a, b];  end, end
end
nF = size(pairIdx, 1);
fprintf('PHASE SWEEP MOVIE: %d frame(s), outer %s, inner %s\n', nF, ...
        pick(order == 1, 'sD', 'sA'), pick(order == 1, 'sA', 'sD'));
if slow ~= 1
    fprintf('  slow x%.3g: %.3g fps, gif delay %.3f s per frame\n', slow, fpsEff, gifDelay);
end

%% PRE-PASS -- fly every certified pair once, so the axis box and the colour
%  axis are known BEFORE the first frame is drawn. This is the whole reason
%  the movie can hold them fixed.
F = cell(1, nF);  tfDays = nan(nF, 1);
lo = inf(1,3);  hi = -inf(1,3);             % extents only: a growing point list
for k = 1:nF                                % would carry 4000 rows per frame
    i = pairIdx(k,1);  j = pairIdx(k,2);
    if sh.has_solution(i,j,1) ~= 1, continue, end
    z8 = sh.z8(:, sh.entry_index(i,j,1));
    rv0 = B.stateD(gridD(i));  rvf = B.stateA(gridA(j));
    fl = fly_transfer(z8, rv0(1:6), rvf(1:6), phys);
    F{k} = fl;  tfDays(k) = fl.tfDays;
    lo = min(lo, min(fl.Y(:,1:3), [], 1));  hi = max(hi, max(fl.Y(:,1:3), [], 1));
end
nCert = nnz(~isnan(tfDays));
assert(nCert > 0, 'movie_phase_sweep:empty', ...
       'none of the %d requested pairs is certified in %s', nF, catMat);

% the orbits belong in the box too: they are drawn in every frame
ss = linspace(0, 1, 600);
Dall = B.stateD(ss);  Aall = B.stateA(ss);
lo = min([lo; Dall(1:3,:).'; Aall(1:3,:).'], [], 1);
hi = max([hi; Dall(1:3,:).'; Aall(1:3,:).'], [], 1);
axLim = boxOf(lo, hi, pad);
clim_ = [0, max(tfDays, [], 'omitnan')];
fprintf('  %d certified, %d gap(s); t_f %.3f .. %.3f d; colour axis [0 %.3f] d\n', ...
        nCert, nF - nCert, min(tfDays, [], 'omitnan'), max(tfDays, [], 'omitnan'), clim_(2));

%% Ephemeris guard (see the header):
restoreAero = hideEphemeris();
cleaner = onCleanup(restoreAero);

%% RENDER:
t0 = tic;
titles = cell(nF, 1);
vw = VideoWriter([outStem '.mp4'], 'MPEG-4');
vw.FrameRate = fpsEff;  vw.Quality = 95;  open(vw);
gifFile = [outStem '.gif'];
for k = 1:nF
    i = pairIdx(k,1);  j = pairIdx(k,2);
    if isempty(F{k})
        titles{k} = titleGap(pr, gridD(i), gridA(j));
        figh = gapFrame(B, Dall, Aall, gridD(i), gridA(j), viewAng, axLim, clim_, d('dark', true));
    else
        T = struct('z', F{k}.z8, 'sD', gridD(i), 'sA', gridA(j));
        P = plot_transfer_3d(T, B, struct('visible', false, 'view', viewAng, ...
                                          'dark', d('dark', true), 'flight', F{k}, ...
                                          'nThrust', d('nThrust', 40)));
        figh = P.fig;
        % PIN what must not move between frames
        set(P.ax, 'XLim', axLim(1,:), 'YLim', axLim(2,:), 'ZLim', axLim(3,:), ...
                  'CLim', clim_);
        titles{k} = titleText(pr, gridD(i), gridA(j), P.tfDays, P.dvKms, P.propellantKg);
        title(P.ax, titles{k}, ...
              'Color', pick(d('dark', true), 'w', 'k'), 'FontWeight', 'normal', ...
              'FontName', get(0, 'FixedWidthFontName'), 'Interpreter', 'tex');
    end
    figh.Position = [60 60 1280 720];
    drawnow;
    img = frame1280(figh);
    writeVideo(vw, img);
    [Aq, map] = rgb2ind(img, 256);
    if k == 1
        imwrite(Aq, map, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', gifDelay);
    else
        imwrite(Aq, map, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', gifDelay);
    end
    close(figh);
    if mod(k, 5) == 0 || k == nF, fprintf('    frame %3d / %3d\n', k, nF); end
end
close(vw);

out = struct('mp4', [outStem '.mp4'], 'gif', gifFile, 'nFrames', nF, ...
             'nCertified', nCert, 'nGaps', nF - nCert, 'clim', clim_, ...
             'axLim', axLim, 'pairs', [gridD(pairIdx(:,1)).', gridA(pairIdx(:,2)).'], ...
             'tfDays', tfDays, 'titles', {titles}, 'slow', slow, 'fps', fpsEff, ...
             'gifDelay', gifDelay, 'runSec', nF/fpsEff, 'wallSec', toc(t0));
fprintf('  -> %s.mp4 / .gif  (%d frames, 1280x720, %.3g fps = %.1f s of video, built in %.0f s)\n', ...
        outStem, nF, fpsEff, out.runSec, out.wallSec);
end

% ------------------------------------------------------------------------
function idx = selectAxis(opts, fVal, fIdx, gridv, tol, name, defaultIdx)
% SELECTAXIS  Grid indices from a VALUE list (fVal) or an INDEX list (fIdx),
% 'all' in either, refusing anything off the grid by name.
% INPUTS: opts; fVal; fIdx; gridv; tol; name; defaultIdx.  OUTPUTS: idx.
hasV = isfield(opts, fVal) && ~isempty(opts.(fVal));
hasI = isfield(opts, fIdx) && ~isempty(opts.(fIdx));
assert(~(hasV && hasI), 'movie_phase_sweep:bothForms', ...
       'give opts.%s (values) or opts.%s (indices) for the %s phase, not both', fVal, fIdx, name);
if ~hasV && ~hasI
    if isempty(defaultIdx), idx = 1:numel(gridv); else, idx = defaultIdx; end
    return
end
if hasI
    v = opts.(fIdx);
    if ischar(v) || isstring(v)
        assert(strcmpi(v, 'all'), 'movie_phase_sweep:badIndex', 'opts.%s must be ''all'' or indices', fIdx);
        idx = 1:numel(gridv);  return
    end
    idx = v(:).';
    assert(all(idx == round(idx)) && all(idx >= 1 & idx <= numel(gridv)), ...
           'movie_phase_sweep:badIndex', 'opts.%s must be integers in 1..%d', fIdx, numel(gridv));
    return
end
v = opts.(fVal);
if ischar(v) || isstring(v)
    assert(strcmpi(v, 'all'), 'movie_phase_sweep:badValue', 'opts.%s must be ''all'' or phase values', fVal);
    idx = 1:numel(gridv);  return
end
v = mod(v(:).', 1);                         % a phase is a fraction of a period
idx = zeros(1, numel(v));
for k = 1:numel(v)
    [dmin, m] = min(abs(gridv - v(k)));
    assert(dmin <= tol, 'movie_phase_sweep:offGrid', ...
           ['%s phase %.6f is not a certified grid value (nearest %.6f, %.2e away).\n' ...
            'The grid is: %s'], name, v(k), gridv(m), dmin, ...
           strjoin(compose('%.4f', gridv), ' '));
    idx(k) = m;
end
end

% ------------------------------------------------------------------------
function s = titleText(pr, sD, sA, tfDays, dvKms, propKg)
% TITLETEXT  Two lines, every numeric field at CONSTANT WIDTH so the title
% geometry is identical frame to frame and only the digits move.
% INPUTS: pr; sD; sA; tfDays; dvKms; propKg.  OUTPUTS: s {1 x 2}.
s = {sprintf('Minimum-time DRO \\rightarrow %d-petal tulip   |   %.0f mN, I_{sp} %g s, %g kg', ...
             pr.NpTulip, pr.thrustN*1000, pr.ispS, pr.m0kg), ...
     sprintf(['s_D = %6.4f   s_A = %6.4f   |   t_f = %6.3f d   ' ...
              '\\DeltaV = %6.4f km/s   propellant %5.2f kg'], sD, sA, tfDays, dvKms, propKg)};
end

% ------------------------------------------------------------------------
function s = titleGap(pr, sD, sA)
% TITLEGAP  The same two lines with the numbers dashed out at the SAME
% widths, so a gap frame does not move the title.
% INPUTS: pr; sD; sA.  OUTPUTS: s {1 x 2}.
s = {sprintf('Minimum-time DRO \\rightarrow %d-petal tulip   |   %.0f mN, I_{sp} %g s, %g kg', ...
             pr.NpTulip, pr.thrustN*1000, pr.ispS, pr.m0kg), ...
     sprintf(['s_D = %6.4f   s_A = %6.4f   |   t_f = %6s d   ' ...
              '\\DeltaV = %6s km/s   propellant %5s kg'], sD, sA, '-----', '------', '----')};
end

% ------------------------------------------------------------------------
function fig = gapFrame(B, Dall, Aall, sD, sA, viewAng, axLim, clim_, dark)
% GAPFRAME  A pair with no certified entry: the orbits, the same box, view,
% colour axis and title layout, and a note where the arc would be.
% INPUTS: B; Dall; Aall; sD; sA; viewAng; axLim; clim_; dark.  OUTPUTS: fig.
bg = pick(dark, 'k', 'w');  fg = pick(dark, 'w', 'k');
fig = figure('Color', bg, 'Position', [60 60 1280 720], 'Visible', 'off');
ax = axes(fig);  hold(ax, 'on');
set(ax, 'Color', bg, 'XColor', fg, 'YColor', fg, 'ZColor', fg, 'GridAlpha', 0.25);
plot3(ax, Dall(1,:), Dall(2,:), Dall(3,:), '-', 'Color', [0.25 0.60 0.35], ...
      'LineWidth', 1.4, 'DisplayName', sprintf('DRO (\\tau = %.2f)', B.problem.tauDRO));
plot3(ax, Aall(1,:), Aall(2,:), Aall(3,:), '-', 'Color', [0.75 0.28 0.28], ...
      'LineWidth', 1.4, 'DisplayName', sprintf('%d-petal tulip', B.problem.NpTulip));
plot3(ax, 1 - B.mu, 0, 0, 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.7 0.7 0.72], ...
      'MarkerEdgeColor', 'none', 'DisplayName', 'Moon');
colormap(ax, parula);
cb = colorbar(ax);  cb.Label.String = 'elapsed time [days]';
cb.Color = fg;  cb.Label.Color = fg;
grid(ax, 'on');  box(ax, 'off');
axis(ax, 'equal');  ax.DataAspectRatio = [1 1 1];
set(ax, 'XLim', axLim(1,:), 'YLim', axLim(2,:), 'ZLim', axLim(3,:), 'CLim', clim_);
view(ax, viewAng);
xlabel(ax, 'x [ND, rotating frame]');  ylabel(ax, 'y [ND]');  zlabel(ax, 'z [ND]');
lg = legend(ax, 'Location', 'northeast');
lg.TextColor = fg;  lg.Color = pick(dark, [0.1 0.1 0.1], 'w');
lg.EdgeColor = pick(dark, [0.3 0.3 0.3], [0.7 0.7 0.7]);
text(ax, mean(axLim(1,:)), mean(axLim(2,:)), mean(axLim(3,:)), ...
     'no certified solution at this phase pair', 'Color', [1 0.75 0.35], ...
     'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
title(ax, titleGap(B.problem, sD, sA), 'Color', fg, 'FontWeight', 'normal', ...
      'FontName', get(0, 'FixedWidthFontName'), 'Interpreter', 'tex');
end

% ------------------------------------------------------------------------
function L = boxOf(lo, hi, pad)
% BOXOF  A common EQUAL-ASPECT box around the extents, padded.  A shared box
% that is not equal-aspect would distort differently at different view
% angles.  INPUTS: lo [1x3]; hi [1x3]; pad.  OUTPUTS: L [3 x 2].
c = (lo + hi)/2;  half = max((hi - lo)/2);
half = half*(1 + pad);
L = [c(:) - half, c(:) + half];
end

% ------------------------------------------------------------------------
function img = frame1280(fig)
% FRAME1280  getframe forced to exactly 720 x 1280 by index resampling: a
% frame size not divisible by 16 shears under H.264.
% INPUTS: fig.  OUTPUTS: img [720 x 1280 x 3 uint8].
F = getframe(fig);
ri = max(1, min(size(F.cdata,1), round(linspace(1, size(F.cdata,1), 720))));
ci = max(1, min(size(F.cdata,2), round(linspace(1, size(F.cdata,2), 1280))));
img = F.cdata(ri, ci, :);
end

% ------------------------------------------------------------------------
function restore = hideEphemeris()
% HIDEEPHEMERIS  Hide a data-less planetEphemeris for the render, returning
% the restorer.  INPUTS: none.  OUTPUTS: restore (function handle).
restore = @() [];
if ~exist('planetEphemeris', 'file'), return, end
try
    planetEphemeris(juliandate(datetime(2030,1,1)), 'Earth', 'Moon');
catch
    aeroDir = fileparts(which('planetEphemeris'));
    if ~isempty(aeroDir)
        rmpath(aeroDir);
        restore = @() addpath(aeroDir);
    end
end
end

% ------------------------------------------------------------------------
function v = pick(c, a, b)
% PICK  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
