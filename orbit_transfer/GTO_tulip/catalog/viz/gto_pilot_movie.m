function gto_pilot_movie(entryMat, outStem)
% GTO_PILOT_MOVIE  Rotating-frame movie of one GTO->tulip catalog entry.
%
%   Flies the entry's z8 with pumpkyn tfMinProp and renders the transfer in
%   the Earth-Moon rotating frame: GTO departure locus (dashed), tulip
%   arrival orbit (dashed), transfer trail colored by time, moving
%   spacecraft marker. House movie mechanics per the catalog extremes
%   movies: black background, frames EXACTLY 1280x720 (H.264 shear lesson),
%   MPEG-4 + gif (every 3rd frame, 10 fps).
%
% INPUTS:
%   entryMat - .mat with fields z8 [8x1], tf, sD, sA, rungN, meta
%              (meta carries muStar/lStar/tStar/ispS/m0kg + families) [char]
%   outStem  - output path stem; writes <stem>.mp4 and <stem>.gif [char]
%
% OUTPUTS:
%   none (files written)
%
% REFERENCES:
%   [1] costate_catalog_extremes_movies.m (the frame/codec conventions)
%   [2] costate_common/get_family_orbit.m ('gto' + 'tulip' construction)

E = load(entryMat);
mu = E.meta.muStar;  lStar = E.meta.lStar;  tStar = E.meta.tStar;
g0 = 9.80665*tStar^2/(1000*lStar);
cnd = (E.meta.ispS/tStar)*g0;
Tnd = (E.rungN/E.meta.m0kg)*tStar^2/(lStar*1000);

% Endpoint orbits from the family recipes:
[tD, rvD] = get_family_orbit('gto',   struct('orientDeg', E.meta.tauDRO));
[tA, rvA] = get_family_orbit('tulip', struct('Np', E.meta.NpTulip, 'pm', E.meta.pmTulip));
rv0 = interp1(tD, rvD, mod(E.sD,1)*tD(end), 'spline');
rvf = interp1(tA, rvA, mod(E.sA,1)*tA(end), 'spline');

% Fly the entry:
[tj, yj] = pumpkyn.cr3bp.tfMinProp(E.z8(8), [rv0(:); 1; E.z8(1:7)], Tnd, cnd, mu);
[tu, iu] = unique(tj);
nF = 96;
tq = linspace(0, E.z8(8), nF);
Y = interp1(tu, yj(iu,1:7), tq, 'pchip');

fig = figure('Color','k', 'Position',[60 60 1280 720], 'MenuBar','none');
ax = axes(fig, 'Color','k', 'Position',[0 0 1 1]); hold(ax,'on');
axis(ax,'equal','off');
plot3(ax, rvD(:,1), rvD(:,2), rvD(:,3), '--', 'Color',[0.35 0.55 0.95], 'LineWidth',0.8);
plot3(ax, rvA(:,1), rvA(:,2), rvA(:,3), '--', 'Color',[0.55 0.95 0.55], 'LineWidth',0.8);
plot3(ax, -mu, 0, 0, 'o', 'MarkerSize',10, 'MarkerFaceColor',[0.2 0.45 1], 'MarkerEdgeColor','none');
plot3(ax, 1-mu, 0, 0, 'o', 'MarkerSize',5, 'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerEdgeColor','none');
text(-mu, -0.09, 0, 'Earth', 'Color',[0.6 0.7 1], 'FontSize',10, 'Horiz','center', 'Parent',ax);
text(1-mu, -0.09, 0, 'Moon', 'Color',[0.75 0.75 0.75], 'FontSize',10, 'Horiz','center', 'Parent',ax);
cmap = turbo(nF);
ttl = title(ax, '', 'Color','w', 'FontSize',13, 'FontWeight','normal');
view(ax, 3); axis(ax, 'vis3d');
xlim(ax, [-0.35 1.15]); ylim(ax, [-0.6 0.6]); zlim(ax, [-0.35 0.35]);

vw = VideoWriter([outStem '.mp4'], 'MPEG-4');  vw.FrameRate = 30;  open(vw);
gifFile = [outStem '.gif'];
sat = []; trail = [];
for kf = 1:nF
    if ~isempty(sat), delete(sat); end
    if ~isempty(trail) && all(isvalid(trail)), delete(trail); end
    if kf >= 2
        trail = surface(ax, [Y(1:kf,1) Y(1:kf,1)]', [Y(1:kf,2) Y(1:kf,2)]', ...
                        [Y(1:kf,3) Y(1:kf,3)]', [tq(1:kf); tq(1:kf)], ...
                        'FaceColor','none', 'EdgeColor','interp', 'LineWidth',2.2);
        colormap(ax, cmap);
    else
        trail = [];
    end
    sat = plot3(ax, Y(kf,1), Y(kf,2), Y(kf,3), 'o', 'MarkerSize',7, ...
                'MarkerFaceColor','w', 'MarkerEdgeColor','w');
    set(ttl, 'String', sprintf(['GTO \\rightarrow tulip (Np=%d)   %g N   ', ...
        't = %.2f / %.2f d   m = %.1f kg'], E.meta.NpTulip, E.rungN, ...
        tq(kf)*tStar/86400, E.z8(8)*tStar/86400, Y(kf,7)*E.meta.m0kg));
    drawnow;
    F = getframe(fig);
    img = frame_720p(F.cdata);
    writeVideo(vw, img);
    if mod(kf,3) == 1
        [A, map] = rgb2ind(img, 128);
        if kf == 1
            imwrite(A, map, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
        else
            imwrite(A, map, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
end
% hold the final frame:
for kf = 1:20, writeVideo(vw, img); end
close(vw); close(fig);
fprintf('movie -> %s.mp4 / .gif (%d frames, 1280x720, tf %.2f d)\n', ...
        outStem, nF, E.z8(8)*tStar/86400);
end

% ------------------------------------------------------------------------
function img = frame_720p(cdata)
% FRAME_720P  Force exactly 1280x720 by index resampling (H.264 shear fix).
% INPUTS: cdata [HxWx3 uint8]  OUTPUTS: img [720x1280x3 uint8]
ri = max(1, min(size(cdata,1), round(linspace(1, size(cdata,1), 720))));
ci = max(1, min(size(cdata,2), round(linspace(1, size(cdata,2), 1280))));
img = cdata(ri, ci, :);
end
