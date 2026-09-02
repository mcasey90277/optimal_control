function outStem = minfuel_race_movie(outStem)
%% Purpose:
%
%   Renders the race's deepest min-fuel solution (run_minfuel_race, eps
%   arm) as a two-panel movie: LEFT the lit pumpkynPie CR3BP scene (DRO
%   departure, tulip arrival, flown transfer with trail), RIGHT the
%   min-fuel STORY -- the near-bang throttle profile s(t) with its coast
%   arcs (full curve as a dim static sketch, drawn bright up to the moving
%   cursor, a BURN/COAST indicator colored by the current throttle) above
%   the mass trace with the final-mass readout. Both panels run on one
%   clock. House mechanics: 1280 x 720 exact via index resampling, .mp4 +
%   .gif, ephemeris guard.
%
%  ASSUMPTIONS / NOTES:
%
% • Reads direct/results/minfuel_race.mat (the eps arm's deepest step) and
%   the matching minenergy_pilot record (endpoints, tf, thruster).
% • The trajectory shown is FLOWN from the solution's own initial PMP
%   state with the smoothed field at the deepest eps -- what you watch is
%   the actual extremal, coast arcs included.
%
%% Inputs:
%
%  outStem                  char                    Output stem [default
%                                                   direct/results/
%                                                   minfuel_25_g12_movie]
%
%% Outputs:
%
%  outStem                  char                    The stem written ('' if
%                                                   the animator is
%                                                   unavailable)
%
%% Revision History:
%  M. Casey                                                   (c) 09/02/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
if nargin < 1
    outStem = fullfile(here, '..', 'direct', 'results', 'minfuel_25_g12_movie');
end
if isempty(which('pumpkynPie.plot.SatelliteAnimator'))
    warning('minfuel_race_movie:noAnimator', 'SatelliteAnimator missing -- no movie.');
    outStem = '';  return
end

%% Ephemeris guard (aero-toolbox trap; see costate_catalog_extremes_movies):
if exist('planetEphemeris', 'file')
    try
        planetEphemeris(juliandate(datetime(2030,1,1)), 'Earth', 'Moon');
    catch
        aeroDir = fileparts(which('planetEphemeris'));
        rmpath(aeroDir);
        cleanupPath = onCleanup(@() addpath(aeroDir));
    end
end

%% The solution, its endpoints, and the flight:
L  = load(fullfile(here, '..', 'direct', 'results', 'minfuel_race.mat'));
out = L.out;
A  = out.arms.eps;
Lp = load(fullfile(here, '..', 'direct', 'results', 'minenergy_pilot.mat'));
ki = find(arrayfun(@(r) isequal([r.iD r.iA], out.cellIdx) && ...
                        abs(r.gam - out.gamma) < 1e-9, Lp.R), 1);
rec = Lp.R(ki);
tf = rec.tf;  Tmax = rec.Tmax;  c = rec.c;  muStar = rec.muStar;
tStar = 382981.289129055;  d2day = tStar/86400;
sm = struct('family', 'eps', 'p', A.p(end));
y1 = A.Y{end}(:,1);
[~, ~, Tt, Yt] = cr3bp_minfuel_prop(tf, y1, false, Tmax, c, muStar, sm);

% throttle along the flight (vectorized from the eps law):
nlv = sqrt(sum(Yt(:,11:13).^2, 2) + 1e-300);
Q   = Tmax*(nlv./Yt(:,7) + Yt(:,14)/c);
sT  = min(1, max(0, (Q - (1 - sm.p))/(2*sm.p)));
mT  = Yt(:,7);
mfE = rec.ms.z;  %#ok<NASGU>  % (min-energy lambda0, unused; kept for provenance)
mfEnergy = rec.direct.mf;
coastPct = 100*nnz(sT < 1e-3)/numel(sT);

% endpoint orbits from the fine sheet's meta recipe:
Qs = load(fullfile(here, '..', 'direct', 'results', 'thrust_ladder_12x12.mat'), 'meta');
[tD, rvD, tT2, rvT2] = ladder_endpoints(Qs.meta);

%% LEFT panel: the lit scene:
 epoch = juliandate(datetime(2030,1,1,0,0,0));
params = struct('muStar', muStar, 'lStar', 389703.264829278, 'tStar', tStar);
figL = figure('Color','k', 'Position',[40 40 800 720], ...
              'Visible','off', 'InvertHardcopy','off');
anim = pumpkynPie.plot.SatelliteAnimator(figL, ...
    'Epoch', epoch, 'CR3BPParameters', params, 'ShowEarth', true, ...
    'FigureTheme', 'dark', 'BodyLighting', 'lambertian', ...
    'FitMode', 'trajectory', 'TrailMode', 'distance', ...
    'EnableTrailColor', true, 'EnableSparkles', false, ...
    'EnableSatelliteMarker', true, 'CameraUpVector', [0 0 1], ...
    'View', [-35 22]);
anim.addOrbit(epoch + tD*tStar/86400, rvD, ...
    'Name', 'DRO (depart)', 'Color', [0.25 0.60 0.35], 'LineWidth', 0.9);
anim.addOrbit(epoch + tT2*tStar/86400, rvT2, ...
    'Name', 'tulip (arrive)', 'Color', [0.75 0.28 0.28], 'LineWidth', 0.9);
[tu, iu] = unique(Tt);
anim.addSatellite(epoch + tu*d2day, Yt(iu,1:6), ...
    'Name', 'min-fuel transfer', ...
    'MarkerColor', [1.00 0.92 0.55], 'TrailColor', [0.35 0.72 1.00]);
anim.initialize();
camzoom(1.35);
annotation(figL, 'textbox', [0.02 0.90 0.96 0.09], 'String', ...
    {'Coasting is free: the catalog line''s first MIN-FUEL transfer', ...
     sprintf('DRO(\\tau=1) \\rightarrow tulip(N_p=7)  \\bullet  t_f fixed at 1.2\\timest_f^{min}  \\bullet  \\epsilon \\rightarrow %.4g', sm.p)}, ...
    'Color', 'w', 'EdgeColor', 'none', 'FontSize', 11, ...
    'HorizontalAlignment', 'center');
clockL = annotation(figL, 'textbox', [0.02 0.02 0.96 0.05], 'String', '', ...
    'Color', [1 0.92 0.55], 'EdgeColor', 'none', 'FontSize', 13, ...
    'FontWeight', 'bold', 'FontName', 'FixedWidth', 'HorizontalAlignment', 'center');

%% RIGHT panel: throttle + mass story:
figR = figure('Color','k', 'Position',[860 40 474 720], 'Visible','off', ...
              'InvertHardcopy','off');
tDay = Tt*d2day;
axS = subplot('Position', [0.16 0.56 0.80 0.30]);
plot(axS, tDay, sT, '-', 'Color', [1 1 1 0.25], 'LineWidth', 0.8); hold(axS, 'on');
hS = plot(axS, NaN, NaN, '-', 'Color', [0.35 0.72 1.00], 'LineWidth', 1.8);
hSc = plot(axS, NaN, NaN, 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', [1 0.92 0.55], 'MarkerEdgeColor', 'none');
set(axS, 'Color', 'k', 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], ...
    'XLim', [0 tDay(end)], 'YLim', [-0.05 1.1], 'GridAlpha', 0.3, ...
    'FontSize', 10);
grid(axS, 'on');
ylabel(axS, 'throttle s', 'Color', 'w');
title(axS, 'near-bang throttle: burn \bullet coast \bullet burn', ...
      'Color', 'w', 'FontSize', 11);
hMode = annotation(figR, 'textbox', [0.16 0.88 0.80 0.06], 'String', '', ...
    'Color', 'w', 'EdgeColor', 'none', 'FontSize', 22, 'FontWeight', 'bold', ...
    'FontName', 'FixedWidth', 'HorizontalAlignment', 'center');

axM = subplot('Position', [0.16 0.12 0.80 0.30]);
plot(axM, tDay, mT, '-', 'Color', [1 1 1 0.25], 'LineWidth', 0.8); hold(axM, 'on');
hM = plot(axM, NaN, NaN, '-', 'Color', [0.55 1.0 0.55], 'LineWidth', 1.8);
hMc = plot(axM, NaN, NaN, 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', [1 0.92 0.55], 'MarkerEdgeColor', 'none');
yl = [min(mT) - 0.002, 1.001];
set(axM, 'Color', 'k', 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], ...
    'XLim', [0 tDay(end)], 'YLim', yl, 'GridAlpha', 0.3, 'FontSize', 10);
grid(axM, 'on');
xlabel(axM, 'time [days]', 'Color', 'w');
ylabel(axM, 'mass fraction m/m_0', 'Color', 'w');
annotation(figR, 'textbox', [0.16 0.43 0.80 0.06], 'String', ...
    sprintf('m_f = %.5f   (min-energy: %.5f)  \\bullet  %.0f%% coast', ...
            mT(end), mfEnergy, coastPct), ...
    'Color', [0.55 1.0 0.55], 'EdgeColor', 'none', 'FontSize', 13, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center');
annotation(figR, 'textbox', [0.10 0.005 0.86 0.03], ...
    'String', 'Coorbital | Casey & Koblick', 'Color', [0.5 0.5 0.5], ...
    'FontSize', 10, 'HorizontalAlignment', 'right', 'EdgeColor', 'none');

%% Frames:
 nFrames = 240;
     fps = 30;
     pwL = 800;  pwR = 474;                 % 800 + 6 + 474 = 1280
 divider = uint8(90*ones(720, 1280 - pwL - pwR, 3));
vw = VideoWriter([outStem '.mp4'], 'MPEG-4');
vw.FrameRate = fps;  vw.Quality = 95;  open(vw);
cleanupVW = onCleanup(@() close(vw));
gifFile = [outStem '.gif'];
for kf = 1:nFrames
    tNow = tf*(kf-1)/(nFrames-1);
    anim.update(epoch + tNow*d2day);
    set(clockL, 'String', sprintf('t = %5.2f of %.2f days', tNow*d2day, tf*d2day));
    msk = Tt <= tNow;
    kNow = max(1, nnz(msk));
    set(hS, 'XData', tDay(msk), 'YData', sT(msk));
    set(hSc, 'XData', tDay(kNow), 'YData', sT(kNow));
    set(hM, 'XData', tDay(msk), 'YData', mT(msk));
    set(hMc, 'XData', tDay(kNow), 'YData', mT(kNow));
    if sT(kNow) > 0.5
        set(hMode, 'String', 'BURN',  'Color', [1.0 0.55 0.35]);
    else
        set(hMode, 'String', 'COAST', 'Color', [0.35 0.72 1.00]);
    end
    drawnow limitrate
    img = [panel_frame(figL, pwL), divider, panel_frame(figR, pwR)];
    writeVideo(vw, img);
    if mod(kf, 3) == 1
        [Ai, cmap] = rgb2ind(img, 256);
        if kf == 1
            imwrite(Ai, cmap, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
        else
            imwrite(Ai, cmap, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
end
close(figL);  close(figR);
fprintf('  movie -> %s.mp4 / .gif  (%d frames, 1280x720)\n', outStem, nFrames);
end

% ------------------------------------------------------------------------
function img = panel_frame(fig, w)
% PANEL_FRAME  getframe forced to 720 x w by index resampling (the H.264
% shear cure).  INPUTS: fig; w.  OUTPUTS: img [720 x w x 3 uint8].
F = getframe(fig);
ri = max(1, min(size(F.cdata,1), round(linspace(1, size(F.cdata,1), 720))));
ci = max(1, min(size(F.cdata,2), round(linspace(1, size(F.cdata,2), w))));
img = F.cdata(ri, ci, :);
end
