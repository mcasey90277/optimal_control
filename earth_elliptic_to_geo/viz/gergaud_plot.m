function gergaud_plot(res, outPng, titleStr)
% GERGAUD_PLOT  Static throttle-colored 3D trajectory PNG for a LEO-ellipse
% ->GEO min-fuel transfer (named for the Haberkorn-Martinon-Gergaud 2004
% elements convention that mee_res_to_cart_res reconstructs from).
%
% A single-frame, publication-style companion to transfer_movie.m: same
% burn/coast color law and GEO-ring/Earth backdrop, but one static image
% (no animation, no video I/O) via exportgraphics. Takes a Cartesian `res`
% as produced by mee_res_to_cart_res -- it does not itself run
% elements_to_cart, so callers already holding an MEE/L-domain solution
% should reconstruct via mee_res_to_cart_res first and pass the result in.
%
% INPUTS:
%   res     - Cartesian results struct (mee_res_to_cart_res.m output, or a
%             run_transfer_mee.m Cartesian res): res.cfg (.thrustN .ctf),
%             res.fuel.X (9x(N+1) = [r(3);v(3);m;t;cScale], row 9 unused),
%             res.fuel.U (4x(N+1) = [alpha_ECI(3);throttle s]) [struct]
%   outPng  - output PNG path (with extension) [char/string]
%   titleStr- optional plot title; defaults to a thrust/ctf-derived title
%             built from res.cfg if omitted or empty [char/string]
%
% OUTPUTS: none (PNG file written; path printed)
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/transfer_movie.m (animated sibling this plot
%       reuses the burn/coast color law and GEO-ring/Earth backdrop from).
%   [2] earth_elliptic_to_geo/mee_res_to_cart_res.m (the reconstruction this
%       function consumes -- not duplicated here).
if nargin < 3 || isempty(titleStr)
    titleStr = sprintf('LEO ellipse \\rightarrow GEO min-fuel  (T=%g N, c_{tf}=%.2f)', ...
        res.cfg.thrustN, res.cfg.ctf);
end

X = res.fuel.X;  U = res.fuel.U;
r = X(1:3,:);
s = U(4,:);
burnTol = 0.05;                          % throttle above this = BURNING (red).
burn = s > burnTol;

% --- GEO ring backdrop (target orbit: radius-1 circle, equatorial plane) ---
thG = linspace(0, 2*pi, 361);
geoRing = [cos(thG); sin(thG); 0*thG].';

% --- fixed axis limits ------------------------------------------------------
allP = [r.'; geoRing];
pad = 0.05*(max(allP)-min(allP)+eps);
xl = [min(allP(:,1))-pad(1), max(allP(:,1))+pad(1)];
yl = [min(allP(:,2))-pad(2), max(allP(:,2))+pad(2)];
zl = [min(allP(:,3))-pad(3), max(allP(:,3))+pad(3)];
if diff(zl) < eps, zl = [-0.1 0.1]; end  % near-coplanar case: avoid degenerate zlim

% --- figure -----------------------------------------------------------------
fig = figure('Color','w','Position',[100 100 900 760],'Visible','off');
try, theme(fig,'light'); catch, end
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

xb = r(1,:); yb = r(2,:); zb = r(3,:);
plot3(ax, maskv(xb,~burn), maskv(yb,~burn), maskv(zb,~burn), '-', ...
    'Color',[0.15 0.35 0.85], 'LineWidth',2.0);               % coast (blue)
plot3(ax, maskv(xb,burn),  maskv(yb,burn),  maskv(zb,burn),  '-', ...
    'Color',[0.85 0.15 0.15], 'LineWidth',2.2);               % burn (red)
plot3(ax, geoRing(:,1), geoRing(:,2), geoRing(:,3), '-', ...
    'Color',[0.45 0.72 0.45], 'LineWidth',1.0);               % GEO ring
plot3(ax, 0,0,0,'o','MarkerFaceColor',[0.10 0.35 0.80], ...
    'MarkerEdgeColor','k','MarkerSize',11);                   % Earth
text(ax, 0,0,0.05,'Earth','FontSize',9);

xlim(ax,xl); ylim(ax,yl); zlim(ax,zl);
xlabel(ax,'x (ND)'); ylabel(ax,'y'); zlabel(ax,'z');
view(ax,-37,24); daspect(ax,[1 1 1]);
title(ax, titleStr);

exportgraphics(fig, outPng, 'ContentType','image', 'Resolution',200);
close(fig);
fprintf('WROTE %s\n', outPng);
end

function v = maskv(x,m)
% MASKV  NaN-mask a vector so plotted segments break where the mask is false.
v = x; v(~m) = nan;
end
