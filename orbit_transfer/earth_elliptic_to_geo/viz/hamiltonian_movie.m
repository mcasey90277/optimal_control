function hamiltonian_movie(H, outStem, cfg)
% HAMILTONIAN_MOVIE  Two-panel synced animation of the Hamiltonians along a
% min-fuel MEE transfer: top = L-domain Hamiltonian H_L(t) (breathes once per
% orbit), bottom = time-domain Hamiltonian H_t(t) (conserved -> flat). Both are
% drawn on the SAME vertical scale, so the pedagogical contrast is immediate:
% the same physical extremal has a Hamiltonian that oscillates when parametrized
% by true longitude but is a first integral in physical time.
%
% Playback is UNIFORM IN PHYSICAL TIME (days). A bright "revealed" trace grows
% left-to-right with a lead marker and a shared time cursor.
%
% INPUTS:
%   H       - struct from hamiltonian_along_traj (uses .tdays .HL .Ht .HtMean
%             .HtCoV .revs) [struct]
%   outStem - output basename WITHOUT extension; writes <outStem>.mp4 + .gif [char]
%   cfg     - OPTIONAL struct: .thrustN [N, for the title], .nFrames [240],
%             .fps [20], .titleStr [override headline]
%
% OUTPUTS: none (files written; paths printed)
%
% REFERENCES:
%   [1] verify/hamiltonian_along_traj.m (producer of H).
%   [2] skill matlab-polished-graphics (dark theme, fixed limits, div-by-16 frame).
%   [3] memory matlab-movie-diagonal-streaks (resize frames to a /16 size).
if nargin < 3, cfg = struct(); end
d = @(f,v) optdef(cfg, f, v);
thrustN  = d('thrustN', NaN);
nFrames  = d('nFrames', 240);
fps      = d('fps', 20);

td  = H.tdays(:).';
tdense = linspace(td(1), td(end), 3000);
HLd = interp1(td, H.HL, tdense, 'pchip');
Htd = interp1(td, H.Ht, tdense, 'pchip');
tf  = linspace(td(1), td(end), nFrames);

ylo = min(H.HL);  yhi = max(H.HL);  pad = 0.10*(yhi - ylo);
ylims = [ylo - pad, yhi + pad];
xlims = [td(1), td(end)];

dim  = [0.30 0.34 0.42];      % unrevealed trace
briL = [0.25 0.72 1.00];      % revealed H_L (cyan)
briT = [0.35 1.00 0.55];      % revealed H_t (green)
lead = [1.00 0.90 0.30];      % lead marker (amber)

if isnan(thrustN), thrTxt = ''; else, thrTxt = sprintf('%g N  \\bullet  ', thrustN); end
headline = d('titleStr', sprintf(['A min-fuel spiral''s Hamiltonian:  ' ...
    'H_L breathes each orbit  \\bullet  H_t is conserved']));

fig = figure('Color','k','Position',[100 100 1280 720],'Visible','off');
tl  = tiledlayout(fig, 2, 1, 'TileSpacing','compact','Padding','compact');
title(tl, headline, 'Color','w','FontSize',16,'FontWeight','bold');
subtitle(tl, sprintf('%s%.1f revolutions  \\bullet  %.0f-day transfer  \\bullet  GTO \\rightarrow GEO', ...
    thrTxt, H.revs(end), td(end)), 'Color',[0.7 0.75 0.8],'FontSize',12);

axT = nexttile;   axB = nexttile;
for ax = [axT axB]
    hold(ax,'on');  grid(ax,'on');  box(ax,'on');
    set(ax,'Color','k','XColor',[0.7 0.7 0.7],'YColor',[0.7 0.7 0.7], ...
        'GridColor',[0.4 0.4 0.4],'GridAlpha',0.35,'XLim',xlims,'YLim',ylims, ...
        'FontSize',11);
end
ylabel(axT, 'H_L  (L-domain Hamiltonian)', 'Color','w','FontSize',13);
ylabel(axB, 'H_t  (time-domain Hamiltonian)', 'Color','w','FontSize',13);
xlabel(axB, 'physical time  [days]', 'Color','w','FontSize',13);

annotation(fig,'textbox',[0.68 0.005 0.31 0.03],'String', ...
    'Coorbital | Casey & Koblick','Color',[0.5 0.5 0.5],'FontSize',10, ...
    'HorizontalAlignment','right','EdgeColor','none');

mp4 = [outStem '.mp4'];
vw  = VideoWriter(mp4,'MPEG-4');  vw.FrameRate = fps;  vw.Quality = 98;  open(vw);
gifFrames = cell(1, nFrames);

for fi = 1:nFrames
    tc   = tf(fi);
    mask = tdense <= tc;
    HLnow = interp1(td, H.HL, tc);
    Htnow = interp1(td, H.Ht, tc);

    cla(axT);
    plot(axT, tdense, HLd, '-', 'Color',dim, 'LineWidth',1.0);
    plot(axT, tdense(mask), HLd(mask), '-', 'Color',briL, 'LineWidth',2.2);
    xline(axT, tc, ':', 'Color',[0.6 0.6 0.6], 'LineWidth',1.0);
    plot(axT, tc, HLnow, 'o', 'MarkerFaceColor',lead, 'MarkerEdgeColor','k', 'MarkerSize',9);
    title(axT, sprintf('H_L = %+7.3f   (varies)', HLnow), ...
        'Color',briL, 'FontName','FixedWidth', 'FontSize',15, 'FontWeight','bold');

    cla(axB);
    plot(axB, tdense, Htd, '-', 'Color',dim, 'LineWidth',1.0);
    plot(axB, tdense(mask), Htd(mask), '-', 'Color',briT, 'LineWidth',2.2);
    xline(axB, tc, ':', 'Color',[0.6 0.6 0.6], 'LineWidth',1.0);
    plot(axB, tc, Htnow, 'o', 'MarkerFaceColor',lead, 'MarkerEdgeColor','k', 'MarkerSize',9);
    title(axB, sprintf('H_t = %+7.3f   (conserved: CoV %.0e)', Htnow, H.HtCoV), ...
        'Color',briT, 'FontName','FixedWidth', 'FontSize',15, 'FontWeight','bold');

    drawnow;
    fr = getframe(fig);
    im = imresize(fr.cdata, [720 1280]);        % force /16 size (no H.264 shear)
    writeVideo(vw, im);
    gifFrames{fi} = im;
end
close(vw);

gif = [outStem '.gif'];
for fi = 1:nFrames
    [A,map] = rgb2ind(gifFrames{fi}, 256);
    if fi == 1, imwrite(A,map,gif,'gif','LoopCount',Inf,'DelayTime',1/fps);
    else,       imwrite(A,map,gif,'gif','WriteMode','append','DelayTime',1/fps); end
end
close(fig);
fprintf('hamiltonian_movie -> %s\n                     %s\n', mp4, gif);
end
