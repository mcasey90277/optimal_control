function plot_costate_torus(lib, thrustN, selDepDays, selArrDays, outStem)
%% Purpose:
%
%   Draws the phasing torus of a costate library at ONE thrust level, colored
%   by what the library holds there, with the selected phase pair marked by a
%   yellow star. Two views of the same information:
%
%     Figure 1  the flat unwrap -- departure phase across, arrival phase up,
%               each cell a tile, flight time printed where a solution exists
%     Figure 2  the genuine torus in 3-D, rotatable (departure phase runs
%               around the big circle, arrival phase around the tube), saved
%               as a .fig so it can be spun in MATLAB
%
%   Colour scheme:
%     GREEN   the library has a verified entry at THIS thrust
%     RED     this phase pair solves at some other thrust, but not here --
%             the ladder reached it and then stopped
%     GREY    this phase pair has no entry at any thrust
%     YELLOW STAR  the pair selected by the caller
%
%  ASSUMPTIONS / NOTES:
%
% • Both phase axes are periodic, which is why the domain is a torus and not
%   a rectangle: the left and right edges of the flat plot are the same
%   points, as are the top and bottom.
% • thrustN must be one of lib.thruster.thrust_rungs_N (the rung actually
%   stored); interpolated thrusts have no map of their own.
%
%% Inputs:
%
%  lib                      struct                  costate_lib_dro_tulip_v2
%                                                   structure
%
%  thrustN                  double                  Thrust rung to display
%                                                   (N)
%
%  selDepDays               double                  Selected departure phase
%                                                   (days) -- marked with a
%                                                   star; [] for none
%
%  selArrDays               double                  Selected arrival phase
%                                                   (days); [] for none
%
%  outStem                  char                    Optional path stem for
%                                                   '<stem>_flat.png' and
%                                                   '<stem>_torus.fig/.png';
%                                                   '' or omitted = no files
%
%% Outputs:
%
%  none (two figures drawn; files written when outStem is given)
%
%% Revision History:
%  M. Casey                                                   (c) 08/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the torus at 5 N with an arbitrary phase pair selected:
       L = load('costate_lib_dro_tulip_v2.mat');
     lib = L.costate_lib_dro_tulip_v2;
     plot_costate_torus(lib, 5, 1.0, 10.0);
     return;
end
if ~exist('outStem','var'), outStem = ''; end
if ~exist('selDepDays','var'), selDepDays = []; end
if ~exist('selArrDays','var'), selArrDays = []; end

%% Locate the requested rung:
      kr = find(abs(lib.grid.thrust_N - thrustN) < 1e-9, 1);
if isempty(kr)
    error('plot_costate_torus:noRung', ...
          'thrust %.3f N is not a stored rung; available: %s', ...
          thrustN, num2str(lib.grid.thrust_N));
end

%% Status of every cell at this rung:
    here = lib.grid.has_solution(:,:,kr);      % green
 anyRung = any(lib.grid.has_solution, 3);      % solves somewhere
    else_ = anyRung & ~here;                   % red
   never  = ~anyRung;                          % grey
      TFd = lib.grid.tf_days(:,:,kr);
[nD, nA] = size(here);

  cGreen = [0.15 0.65 0.25];
    cRed = [0.85 0.20 0.15];
   cGrey = [0.80 0.80 0.80];
       C = zeros(nD, nA, 3);
for kD = 1:nD
    for kA = 1:nA
        if here(kD,kA),      col = cGreen;
        elseif else_(kD,kA), col = cRed;
        else,                col = cGrey;
        end
        C(kD,kA,:) = reshape(col,1,1,3);
    end
end

%% Where the selection falls (nearest cell, and its exact phase fractions):
      Pd = lib.departure_params.period_days;
      Pa = lib.arrival_params.period_days;
   selOn = ~isempty(selDepDays) && ~isempty(selArrDays);
if selOn
      fdS = mod(selDepDays/Pd, 1);
      faS = mod(selArrDays/Pa, 1);
end

      sD = lib.grid.departure_phase_frac(:);
      sA = lib.grid.arrival_phase_frac(:);
   sDday = lib.grid.departure_phase_days(:);
   sAday = lib.grid.arrival_phase_days(:);

%% ---- Figure 1: the flat unwrap ----------------------------------------
figure('Color','w','Position',[60 60 860 720]); hold on
for kD = 1:nD
    for kA = 1:nA
        rectangle('Position', [sDday(kD)-Pd/(2*nD), sAday(kA)-Pa/(2*nA), ...
                               Pd/nD, Pa/nA], ...
                  'FaceColor', squeeze(C(kD,kA,:))', ...
                  'EdgeColor', [1 1 1]*0.5, 'LineWidth', 0.3);
        if here(kD,kA)
            text(sDday(kD), sAday(kA), sprintf('%.2f', TFd(kD,kA)), ...
                 'HorizontalAlignment','center', 'FontSize',7, ...
                 'Color','w', 'FontWeight','bold');
        end
    end
end
if selOn
    plot(fdS*Pd, faS*Pa, 'p', 'MarkerSize',26, ...
         'MarkerFaceColor',[1 0.85 0], 'MarkerEdgeColor','k', 'LineWidth',1.2);
end
axis([-Pd/(2*nD), Pd+Pd/(2*nD), -Pa/(2*nA), Pa+Pa/(2*nA)]);
xlabel('departure phase along the DRO [days]   (wraps)');
ylabel('arrival phase along the tulip [days]   (wraps)');
title({sprintf('costate library at %.2f N: %d of %d phase pairs solved', ...
                thrustN, nnz(here), numel(here)), ...
       'green = entry here / red = solves at another thrust / grey = nowhere;  numbers are t_f [days]'});
grid off
if ~isempty(outStem)
    exportgraphics(gcf, [outStem '_flat.png'], 'Resolution', 140);
end

%% ---- Figure 2: the torus in 3-D, rotatable ----------------------------
       R = 2.0;                                  % big radius (departure)
       r = 0.8;                                  % tube radius (arrival)
   uEdge = 2*pi*[sD - 1/(2*nD); sD(end) + 1/(2*nD)];
   vEdge = 2*pi*[sA - 1/(2*nA); sA(end) + 1/(2*nA)];
f2 = figure('Color','w','Position',[60 60 920 780]); hold on
for kD = 1:nD
    for kA = 1:nA
        uu = linspace(uEdge(kD), uEdge(kD+1), 8);
        vv = linspace(vEdge(kA), vEdge(kA+1), 8);
        [U,V] = meshgrid(uu, vv);
        X = (R + r*cos(V)).*cos(U);
        Y = (R + r*cos(V)).*sin(U);
        Z = r*sin(V);
        surf(X, Y, Z, 'FaceColor', squeeze(C(kD,kA,:))', ...
             'EdgeColor','none', 'FaceAlpha',1);
    end
end
% faint cell borders so the grid reads
for kD = 1:nD+1
    vv = linspace(0, 2*pi, 60);  u0 = uEdge(min(kD,nD+1));
    plot3((R+r*cos(vv))*cos(u0), (R+r*cos(vv))*sin(u0), r*sin(vv), ...
          '-', 'Color', [1 1 1]*0.55, 'LineWidth', 0.3);
end
for kA = 1:nA+1
    uu = linspace(0, 2*pi, 120);  v0 = vEdge(min(kA,nA+1));
    plot3((R+r*cos(v0))*cos(uu), (R+r*cos(v0))*sin(uu), ...
          r*sin(v0)*ones(size(uu)), '-', 'Color', [1 1 1]*0.55, 'LineWidth', 0.3);
end
if selOn
    uS = 2*pi*fdS;  vS = 2*pi*faS;
    xS = (R + 1.06*r*cos(vS))*cos(uS);           % lifted clear of the surface
    yS = (R + 1.06*r*cos(vS))*sin(uS);
    zS = 1.06*r*sin(vS);
    plot3(xS, yS, zS, 'p', 'MarkerSize',26, ...
          'MarkerFaceColor',[1 0.85 0], 'MarkerEdgeColor','k', 'LineWidth',1.2);
end
axis equal off;  view(-35, 35);                  % flat colors: no lighting
title({sprintf('phasing torus at %.2f N  (big circle = departure phase, tube = arrival phase)', thrustN), ...
       sprintf('green = %d solved / red = %d solve at another thrust / grey = %d nowhere', ...
               nnz(here), nnz(else_), nnz(never))});
rotate3d on
if ~isempty(outStem)
    savefig(f2, [outStem '_torus.fig']);
    exportgraphics(f2, [outStem '_torus.png'], 'Resolution', 140);
end

end
