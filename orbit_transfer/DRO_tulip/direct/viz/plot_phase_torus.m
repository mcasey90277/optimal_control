function plot_phase_torus(ckptFile, outStem)
% PLOT_PHASE_TORUS  The phasing torus, colored by sweep status.
%
% The (s_D, s_A) domain is a genuine torus -- departure phase around the DRO is
% the big circle, arrival phase around the tulip the small one, both wrapping.
% Two figures:
%   1. the 3-D torus (saved as .fig so it can be rotated in MATLAB), and
%   2. the flat unwrap on [0,1) x [0,1).
% Colour scheme (Mike's): GREY cells untried, GREEN tried-and-passed the map
% gate, RED tried-and-failed. Near-misses (converged but unresolved at map
% resolution, global error 1-100x the gate) draw ORANGE, between the two.
%
% INPUTS:
%   ckptFile - path to a sweep checkpoint/result .mat holding sD, sA, PASS,
%              GLOBKM (and optionally TF) as written by sweep_phasing_direct
%   outStem  - path stem for '<stem>_torus.fig/.png' and '<stem>_flat.png'
%
% OUTPUTS: none (files written)
%
% REFERENCES:
%   [1] sweep_phasing_direct.m -- the data source and gate definitions.

D = load(ckptFile);
sDdat = D.sD(:);  sAdat = D.sA(:);
tried0 = ~isnan(D.TF);
pass0  = D.PASS;
nearm0 = tried0 & ~pass0 & D.GLOBKM < 100*100;  % within 100x of the 100-km gate

% DISPLAY GRID, decoupled from the data grid: a 1xN transect must not paint
% the whole departure circle -- untried rows stay grey. Data points map to
% their nearest display cell (circular distance).
% Use the data's own resolution when it is a genuine 2-D grid; only inflate
% the display grid for transect-style data (a single row or column), where
% painting full rings would misrepresent coverage.
if numel(sDdat) >= 4 && numel(sAdat) >= 4
    nD = numel(sDdat);  nA = numel(sAdat);
else
    nD = max(numel(sDdat), 18);  nA = max(numel(sAdat), 18);
end
sD = (0:nD-1).'/nD;  sA = (0:nA-1).'/nA;
tried = false(nD,nA);  pass = tried;  nearm = tried;  TFdisp = nan(nD,nA);
for a = 1:numel(sDdat)
    for b = 1:numel(sAdat)
        if ~tried0(a,b), continue, end
        [~,ia] = min(min(abs(sD - mod(sDdat(a),1)), 1-abs(sD - mod(sDdat(a),1))));
        [~,ib] = min(min(abs(sA - mod(sAdat(b),1)), 1-abs(sA - mod(sAdat(b),1))));
        tried(ia,ib) = true;
        pass(ia,ib)  = pass0(a,b);
        nearm(ia,ib) = nearm0(a,b);
        TFdisp(ia,ib) = D.TF(a,b);
    end
end

cGrey  = [0.80 0.80 0.80];
cGreen = [0.15 0.65 0.25];
cRed   = [0.85 0.20 0.15];
cOrng  = [0.95 0.60 0.10];

% per-cell colors [nD x nA x 3]
C = repmat(reshape(cGrey,1,1,3), nD, nA);
for kD = 1:nD
    for kA = 1:nA
        if ~tried(kD,kA), continue, end
        if pass(kD,kA),      col = cGreen;
        elseif nearm(kD,kA), col = cOrng;
        else,                col = cRed;
        end
        C(kD,kA,:) = reshape(col,1,1,3);
    end
end

%% ---- figure 1: the 3-D torus --------------------------------------------
% Cell-centred double sampling so each grid cell is one flat-colored tile.
R = 2.0;  r = 0.8;
uEdge = 2*pi*[sD - 1/(2*nD); sD(end) + 1/(2*nD)];   % departure = big circle
vEdge = 2*pi*[sA - 1/(2*nA); sA(end) + 1/(2*nA)];   % arrival   = small circle
f1 = figure('Color','w','Position',[60 60 900 750]);
hold on
for kD = 1:nD
    for kA = 1:nA
        uu = linspace(uEdge(kD), uEdge(kD+1), 8);
        vv = linspace(vEdge(kA), vEdge(kA+1), 8);
        [U,V] = meshgrid(uu, vv);
        X = (R + r*cos(V)).*cos(U);
        Y = (R + r*cos(V)).*sin(U);
        Z = r*sin(V);
        surf(X, Y, Z, 'FaceColor', squeeze(C(kD,kA,:)).', ...
             'EdgeColor','none','FaceAlpha',1);
    end
end
% cell borders for readability
for kD = 1:nD+1
    vv = linspace(0,2*pi,60);  u0 = uEdge(min(kD,nD+1));
    plot3((R+r*cos(vv))*cos(u0),(R+r*cos(vv))*sin(u0),r*sin(vv),'-','Color',[1 1 1]*0.55,'LineWidth',0.3);
end
for kA = 1:nA+1
    uu = linspace(0,2*pi,120); v0 = vEdge(min(kA,nA+1));
    plot3((R+r*cos(v0))*cos(uu),(R+r*cos(v0))*sin(uu),r*sin(v0)*ones(size(uu)),'-','Color',[1 1 1]*0.55,'LineWidth',0.3);
end
axis equal off;  view(-35, 35);   % flat colors -- no lighting, it washes cells out
title(sprintf(['DRO \\rightarrow tulip phasing torus  |  big circle = departure s_D, ' ...
    'small = arrival s_A\ngreen = CONTINUOUSLY VERIFIED (flown control lands <100 km, map res) / ' ...
    'orange near-miss / red = no valid trajectory / grey untried  (%d tried, %d verified)'], ...
    nnz(tried), nnz(pass)));
savefig(f1, [outStem '_torus.fig']);
exportgraphics(f1, [outStem '_torus.png'], 'Resolution', 140);
fprintf('  torus  -> %s_torus.fig (rotatable) + .png\n', outStem);

%% ---- figure 2: the flat unwrap ------------------------------------------
f2 = figure('Color','w','Position',[60 60 820 700]);
hold on
for kD = 1:nD
    for kA = 1:nA
        x = mod(sD(kD),1);  y = mod(sA(kA),1);
        rectangle('Position',[x-1/(2*nD), y-1/(2*nA), 1/nD, 1/nA], ...
            'FaceColor', squeeze(C(kD,kA,:)).', 'EdgeColor',[1 1 1]*0.5, 'LineWidth',0.3);
        if tried(kD,kA) && pass(kD,kA) && isfinite(TFdisp(kD,kA))
            text(x, y, sprintf('%.2f', TFdisp(kD,kA)), 'HorizontalAlignment','center', ...
                'FontSize',8, 'Color','w', 'FontWeight','bold');
        end
    end
end
axis([-0.5/nD 1-0.5/nD+1/nD -0.5/nA 1-0.5/nA+1/nA]);  axis square
xlabel('departure phase s_D (wraps)');  ylabel('arrival phase s_A (wraps)');
title({'flat torus: green = continuously verified, flown control lands <100 km (t_f shown)', ...
    'orange near-miss / red = solver certificate but NO valid trajectory / grey untried; axes periodic'});
exportgraphics(f2, [outStem '_flat.png'], 'Resolution', 140);
fprintf('  flat   -> %s_flat.png\n', outStem);

%% ---- figure 3: the t_f-VALUE map ----------------------------------------
% The actual product: t_f over the torus, colormapped on the PASSED cells.
% Untried cells light grey; tried-but-failed cells white with a red cross, so
% holes in the map read as "unresolved", never as data.
f3 = figure('Color','w','Position',[60 60 880 720]);
hold on
tfOK = TFdisp(pass & tried);
if isempty(tfOK), tfLo = 0; tfHi = 1; else
    tfLo = min(tfOK);  tfHi = max(tfOK);
    if tfHi <= tfLo, tfHi = tfLo + eps(tfLo); end
end
cmap = turbo(256);
for kD = 1:nD
    for kA = 1:nA
        x = mod(sD(kD),1);  y = mod(sA(kA),1);
        px = [x-1/(2*nD), y-1/(2*nA), 1/nD, 1/nA];
        if tried(kD,kA) && pass(kD,kA)
            w = (TFdisp(kD,kA) - tfLo)/(tfHi - tfLo);
            col = cmap(1 + round(w*255), :);
            rectangle('Position',px,'FaceColor',col,'EdgeColor',[1 1 1]*0.5,'LineWidth',0.3);
            text(x, y, sprintf('%.2f', TFdisp(kD,kA)), 'HorizontalAlignment','center', ...
                'FontSize',8, 'Color','w', 'FontWeight','bold');
        elseif tried(kD,kA)
            rectangle('Position',px,'FaceColor',[1 1 1],'EdgeColor',[0.85 0.2 0.15],'LineWidth',0.8);
            plot(x, y, 'x', 'Color',[0.85 0.2 0.15], 'MarkerSize',7, 'LineWidth',1.2);
        else
            rectangle('Position',px,'FaceColor',[0.88 0.88 0.88],'EdgeColor',[1 1 1]*0.6,'LineWidth',0.3);
        end
    end
end
axis([-0.5/nD 1-0.5/nD+1/nD -0.5/nA 1-0.5/nA+1/nA]);  axis square
xlabel('departure phase s_D (wraps)');  ylabel('arrival phase s_A (wraps)');
colormap(f3, cmap);
cb = colorbar;  caxis([tfLo tfHi]);
cb.Label.String = sprintf('t_f [ND]   (%.1f - %.1f days)', ...
    tfLo*4.4327, tfHi*4.4327);
title({'minimum transfer time over the phasing torus (map resolution, pre-certification)', ...
       'colored where the CONTINUOUS trajectory is verified (<100 km flown);  red x = discrete-only, rejected;  grey = untried'});
exportgraphics(f3, [outStem '_tfmap.png'], 'Resolution', 140);
fprintf('  tf map -> %s_tfmap.png\n', outStem);
end
