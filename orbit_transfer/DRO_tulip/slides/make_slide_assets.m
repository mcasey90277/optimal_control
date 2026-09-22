function make_slide_assets()
%% Purpose:
%
%   Render the assets for the DRO -> tulip overview deck from the 70 mN
%   library of record (24 x 24): the FASTEST and SLOWEST certified transfers
%   side by side on one clock (.mp4 + .gif, 1280 x 720), and a CR3BP
%   geometry figure (rotating frame, Earth, Moon, L1/L2, the DRO and the
%   7-petal tulip). Endpoints use the catalog's own rule (ladder_endpoints +
%   phase_state), exactly as audit_phase_catalog does.
%
%% Inputs:
%
%  (none)
%
%% Outputs:
%
%  (files)                  slides/assets/extremes_70mN.{mp4,gif},
%                           slides/assets/cr3bp_geometry.png,
%                           slides/assets/extremes_70mN.txt (numbers)
%
%% Revision History:
%  M. Casey                                                   (c) 09/22/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here   = fileparts(mfilename('fullpath'));
outDir = fullfile(here, 'assets');
addpath(fullfile(here, '..', '..', 'costate_common'));
catMat = fullfile(here, '..', 'indirect', 'results', 'library_70mN_24x24_final', ...
                  'costate_catalog_dro_tulip_70mN.mat');

L = load(catMat);  fn = fieldnames(L);  c = L.(fn{1});
s = c.sheets(1);
lStar = c.constants.lStar_km;  tStar = c.constants.tStar_s;  mu = c.constants.muStar;
m0 = c.thruster.m0_kg;  isp = c.thruster.isp_s;
ndp = nd_propulsion(c.rungs_N(1), isp, m0, lStar, tStar);
cnd = ndp.cnd;  Tnd = ndp.Tnd;
ob = struct('muStar', mu, 'lStar', lStar, 'tStar', tStar, 'tauDRO', s.tauDRO, ...
            'NpTulip', s.Np, 'tauTulip', s.period_tulip_nd, 'pmTulip', s.pm, ...
            'ispS', isp, 'm0kg', m0);
[tD, rvD, tT, rvT] = ladder_endpoints(ob);
stD = phase_state(tD, rvD);  stA = phase_state(tT, rvT);

%% The two extremes of the torus:
tf = s.tf_nd(:, :, 1);  tf(~s.has_solution(:, :, 1)) = NaN;
[~, kMin] = min(tf(:));  [~, kMax] = max(tf(:));
K = [kMin kMax];  lab = {'FASTEST', 'SLOWEST'};
fid = fopen(fullfile(outDir, 'extremes_70mN.txt'), 'w');
for k = 1:2
    [iD, iA] = ind2sub(size(tf), K(k));
    z8  = s.z8(:, s.entry_index(iD, iA, 1));
    rv0 = stD(s.sD_frac(iD));  rvf = stA(s.sA_frac(iA));
    [tj, yj] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(1:6); 1; z8(1:7)], Tnd, cnd, mu);
    [~, iu] = unique(tj);  yj = yj(iu, :);
    % dense re-flight of the same all-burn field (the movie trail is polygonal
    % on the integrator's own steps); endpoint checked against tfMinProp
    y0 = [rv0(1:6); 1; z8(1:7)];
    [tu, yd] = ode45(@pumpkyn.cr3bp.tfMinEoM, linspace(0, z8(8), 3000), y0, ...
                     odeset('RelTol', 1e-10, 'AbsTol', 1e-12), Tnd, cnd, mu);
    assert(max(abs(yd(end,1:7) - yj(end,1:7))) < 1e-8, 'dense re-flight disagrees with tfMinProp');
    yj  = yd;
    mf  = yj(end, 7);
    E(k) = struct('sD', s.sD_frac(iD), 'sA', s.sA_frac(iA), 'tf', z8(8), ...
                  'tfDays', z8(8)*tStar/86400, 'dV', cnd*log(1/mf)*lStar/tStar, ...
                  'prop', m0*(1 - mf), 'missKm', sqrt(sum((yj(end,1:3) - rvf(1:3)').^2))*lStar, ...
                  't', tu, 'y', yj); %#ok<AGROW>
    fprintf(fid, '%s sD=%.4f sA=%.4f tf=%.3f d dV=%.4f km/s prop=%.2f kg miss=%.3g km\n', ...
            lab{k}, E(k).sD, E(k).sA, E(k).tfDays, E(k).dV, E(k).prop, E(k).missKm);
    fprintf('%s sD=%.4f sA=%.4f tf=%.3f d dV=%.4f km/s prop=%.2f kg miss=%.3g km\n', ...
            lab{k}, E(k).sD, E(k).sA, E(k).tfDays, E(k).dV, E(k).prop, E(k).missKm);
end
fclose(fid);

cr3bp_figure(outDir, mu, lStar, tD, rvD, tT, rvT);

%% Ephemeris guard (memory: aero-toolbox-ephemeris-trap), as in movie_70mN:
if exist('planetEphemeris', 'file')
    try
        planetEphemeris(juliandate(datetime(2030,1,1)), 'Earth', 'Moon');
    catch
        aeroDir = fileparts(which('planetEphemeris'));
        rmpath(aeroDir);
        cleanupPath = onCleanup(@() addpath(aeroDir)); %#ok<NASGU>
    end
end

%% One lit pumpkyn scene per panel, shared clock (costate_catalog_extremes_movies):
params = struct('muStar', mu, 'lStar', lStar, 'tStar', tStar);
epoch  = juliandate(datetime(2030,1,1,0,0,0));
for k = 1:2
    S(k).fig = figure('Color','k', 'Position',[60 60 640 720], ...
                      'Visible','off', 'InvertHardcopy','off'); %#ok<AGROW>
    S(k).anim = pumpkynPie.plot.SatelliteAnimator(S(k).fig, ...
        'Epoch', epoch, 'CR3BPParameters', params, 'ShowEarth', true, ...
        'FigureTheme', 'dark', 'BodyLighting', 'lambertian', ...
        'FitMode', 'trajectory', 'TrailMode', 'distance', 'TrailPointCount', 1500, ...
        'EnableTrailColor', true, 'EnableSparkles', false, ...
        'EnableSatelliteMarker', true, 'CameraUpVector', [0 0 1], 'View', [-35 22]);
    S(k).anim.addOrbit(epoch + tD*tStar/86400, rvD, ...
        'Name', 'DRO (\tau = 1)', 'Color', [0.25 0.60 0.35], 'LineWidth', 0.9);
    S(k).anim.addOrbit(epoch + tT*tStar/86400, rvT, ...
        'Name', '7-petal tulip', 'Color', [0.75 0.28 0.28], 'LineWidth', 0.9);
    S(k).anim.addSatellite(epoch + E(k).t*tStar/86400, E(k).y(:,1:6), ...
        'Name', sprintf('%s transfer', lower(lab{k})), ...
        'MarkerColor', [1.00 0.92 0.55], 'TrailColor', [0.35 0.72 1.00]);
    S(k).anim.initialize();
    camzoom(1.4);
    annotation(S(k).fig, 'textbox', [0.02 0.89 0.96 0.10], 'String', ...
        {sprintf('%s:  t_f = %.2f d', lab{k}, E(k).tfDays), ...
         sprintf('\\DeltaV %.3f km/s,  %.2f kg,  s_D %.3f / s_A %.3f', ...
                 E(k).dV, E(k).prop, E(k).sD, E(k).sA)}, ...
        'Color', 'w', 'EdgeColor', 'none', 'FontSize', 13, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
    S(k).clock = annotation(S(k).fig, 'textbox', [0.02 0.02 0.96 0.05], ...
        'String', '', 'Color', [1 0.92 0.55], 'EdgeColor', 'none', ...
        'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    S(k).tf = E(k).tf;
end
tEnd = max(S(1).tf, S(2).tf);

nFrames = 240;  fps = 30;  pw = 637;
divider = uint8(90*ones(720, 1280 - 2*pw, 3));
stem = fullfile(outDir, 'extremes_70mN');
vw = VideoWriter([stem '.mp4'], 'MPEG-4');  vw.FrameRate = fps;  vw.Quality = 95;  open(vw);
gifFile = [stem '.gif'];
for kf = 1:nFrames
    tNow = tEnd*(kf - 1)/(nFrames - 1);
    for k = 1:2
        S(k).anim.update(epoch + min(tNow, S(k).tf)*tStar/86400);
        if tNow >= S(k).tf
            set(S(k).clock, 'String', sprintf('ARRIVED  t = %.2f days', S(k).tf*tStar/86400));
        else
            set(S(k).clock, 'String', sprintf('t = %.2f days', tNow*tStar/86400));
        end
    end
    drawnow limitrate
    img = [panel_frame(S(1).fig, pw), divider, panel_frame(S(2).fig, pw)];
    writeVideo(vw, img);
    if mod(kf, 3) == 1
        [A, cmap] = rgb2ind(img, 256);
        if kf == 1
            imwrite(A, cmap, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
        else
            imwrite(A, cmap, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
end
close(vw);
imwrite(img, fullfile(outDir, 'extremes_70mN_poster.png'));
close(S(1).fig);  close(S(2).fig);
fprintf('movies -> %s.{mp4,gif}\n', stem);
end

% ---------------------------------------------------------------------------
function cr3bp_figure(outDir, mu, lStar, tD, rvD, tT, rvT)
% CR3BP_FIGURE  Rotating-frame geometry: Earth, Moon, L1/L2, DRO, tulip.
% INPUTS: outDir; mu; lStar [km]; tD,rvD / tT,rvT orbit samples [Nx1, Nx6].
% OUTPUTS: (file) cr3bp_geometry.png
f  = @(x) x - (1-mu)*(x+mu)./abs(x+mu).^3 - mu*(x-1+mu)./abs(x-1+mu).^3;
xL1 = fzero(f, [0.5 1-mu-1e-3]);  xL2 = fzero(f, [1-mu+1e-3 1.5]);
fig = figure('Color', 'k', 'Position', [60 60 1200 800], 'Visible', 'off', ...
             'InvertHardcopy', 'off');
ax = axes(fig, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', ...
          'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.2);
hold(ax, 'on');  grid(ax, 'on');  ax.GridColor = [0.4 0.4 0.4];
kk = lStar/1e3;                                     % axes in 10^3 km
xm = 1 - mu;
plot3(ax, (rvD(:,1)-xm)*kk, rvD(:,2)*kk, rvD(:,3)*kk, 'Color', [0.35 0.80 0.45], 'LineWidth', 2.5);
plot3(ax, (rvT(:,1)-xm)*kk, rvT(:,2)*kk, rvT(:,3)*kk, 'Color', [0.95 0.40 0.40], 'LineWidth', 2.0);
plot3(ax, 0, 0, 0, 'o', 'MarkerSize', 14, 'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', 'w');
plot3(ax, (xL1-xm)*kk, 0, 0, 'd', 'MarkerSize', 10, 'MarkerFaceColor', [1 0.85 0.3], 'MarkerEdgeColor', 'none');
plot3(ax, (xL2-xm)*kk, 0, 0, 'd', 'MarkerSize', 10, 'MarkerFaceColor', [1 0.85 0.3], 'MarkerEdgeColor', 'none');
text(ax, 0, 0, -6, 'Moon', 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(ax, (xL1-xm)*kk, 0, -8, '\leftarrow Earth (384,400 km)', 'Color', [0.6 0.8 1], 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(ax, (xL1-xm)*kk, 0, 8, 'L_1', 'Color', [1 0.85 0.3], 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(ax, (xL2-xm)*kk, 0, 8, 'L_2', 'Color', [1 0.85 0.3], 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
xlabel(ax, 'x  [10^3 km, Moon-centred]');  ylabel(ax, 'y  [10^3 km]');  zlabel(ax, 'z  [10^3 km]');
legend(ax, {sprintf('DRO (period %.1f d)', tD(end)*382981.289129055/86400), ...
        sprintf('7-petal tulip (period %.1f d)', tT(end)*382981.289129055/86400), 'Moon', 'L_1', 'L_2'}, ...
       'TextColor', 'w', 'Color', 'k', 'EdgeColor', [0.5 0.5 0.5], 'Location', 'northeast', 'FontSize', 12);
axis(ax, 'equal');  view(ax, [-30 25]);
exportgraphics(fig, fullfile(outDir, 'cr3bp_geometry.png'), 'Resolution', 150, 'BackgroundColor', 'k');
close(fig);
end

% ---------------------------------------------------------------------------
function img = panel_frame(fig, w)
% PANEL_FRAME  getframe forced to 720 x w by index resampling (H.264 shear cure).
% INPUTS: fig; w target width.  OUTPUTS: img [720 x w x 3 uint8].
F = getframe(fig);
ri = max(1, min(size(F.cdata,1), round(linspace(1, size(F.cdata,1), 720))));
ci = max(1, min(size(F.cdata,2), round(linspace(1, size(F.cdata,2), w))));
img = F.cdata(ri, ci, :);
end
