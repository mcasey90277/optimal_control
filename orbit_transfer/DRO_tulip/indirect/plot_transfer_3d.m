function P = plot_transfer_3d(T, B, opts)
%% Purpose:
%
%   The interactive 3D view of one transfer: the movie's last frame, as a
%   figure you can rotate. Departure orbit, target orbit, the transfer arc
%   coloured by elapsed time, the thrust direction along the way, and
%   the Moon, in the rotating Earth-Moon frame.
%
%   Everything drawn comes from ONE flight of the certificate's own costates.
%   A caller that already owns that flight (the study script) passes it in
%   as opts.flight, so the figure and the printed numbers describe the same
%   trajectory; otherwise it is flown here. The annotations are recomputed
%   from whichever flight is drawn, never copied from the certificate.
%
%% Inputs:
%
%  T                        struct                  a certify_root output
%                                                   (needs .z, and .sD/.sA)
%  B                        struct                  arclength_arrival setup
%                                                   (.Tnd .cnd .mu .stateD
%                                                   .stateA .problem)
%  opts                     struct (optional)
%   .visible [true] .nThrust [40] arrows .outPng '' .view [-35 22]
%   .dark [true] .flight (struct .t [N x 1], .Y [N x 14]) a flight of T.z
%   to draw instead of re-flying
%
%% Outputs:
%
%  P                        struct                  .fig .ax .hDep .hArr
%                                                   .hTx .startMissKm
%                                                   .endMissKm .tfDays
%                                                   .dvKms .nThrust
%                                                   .flightSupplied
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
dark = d('dark', true);
pr = B.problem;  lStar = pr.lStar;  tStar = pr.tStar;  mu = pr.muStar;

% ---- the flight: supplied, or flown here -------------------------------
rv0 = B.stateD(T.sD);  rvf = B.stateA(T.sA);
z = T.z(:);
fl = d('flight', []);
P = struct('flightSupplied', ~isempty(fl));
if P.flightSupplied
    tu = fl.t(:);  Y = fl.Y;
    assert(size(Y, 2) >= 14 && numel(tu) == size(Y, 1), 'opts.flight needs .t [N x 1] and .Y [N x 14]');
else
    [tu, Y] = pumpkyn.cr3bp.tfMinProp(z(8), [rv0(1:6); 1; z(1:7)], B.Tnd, B.cnd, mu);
end
r = Y(:, 1:3);
lamv = Y(:, 11:13);
alpha = -lamv ./ max(vecnorm(lamv, 2, 2), realmin);   % PMP thrust direction

% ---- the two periodic orbits, as the residual sees them ------------------
% vectorised: the closures are ppval-based, which takes a whole query vector
ss = linspace(0, 1, 600);
Dall = B.stateD(ss);   Aall = B.stateA(ss);
D = Dall(1:3, :).';    A = Aall(1:3, :).';

fig = figure('Color', pick(dark, 'k', 'w'), 'Position', [80 80 1100 800], ...
             'Visible', pick(d('visible', true), 'on', 'off'), 'InvertHardcopy', 'off');
ax = axes(fig);  hold(ax, 'on');
set(ax, 'Color', pick(dark, 'k', 'w'), 'XColor', pick(dark, 'w', 'k'), ...
        'YColor', pick(dark, 'w', 'k'), 'ZColor', pick(dark, 'w', 'k'), 'GridAlpha', 0.25);

P.hDep = plot3(ax, D(:,1), D(:,2), D(:,3), '-', 'Color', [0.25 0.60 0.35], ...
               'LineWidth', 1.4, 'DisplayName', sprintf('DRO (\\tau = %.2f)', pr.tauDRO));
P.hArr = plot3(ax, A(:,1), A(:,2), A(:,3), '-', 'Color', [0.75 0.28 0.28], ...
               'LineWidth', 1.4, 'DisplayName', sprintf('%d-petal tulip', pr.NpTulip));
% the arc, coloured by elapsed time
P.hTx = patch(ax, 'XData', [r(:,1); nan], 'YData', [r(:,2); nan], 'ZData', [r(:,3); nan], ...
    'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2.0, ...
    'FaceVertexCData', [tu(:); tu(end)]*tStar/86400, 'DisplayName', 'minimum-time transfer');
colormap(ax, parula);
cb = colorbar(ax);  cb.Label.String = 'elapsed time [days]';
cb.Color = pick(dark, 'w', 'k');  cb.Label.Color = pick(dark, 'w', 'k');

% thrust direction, on a subsample
nT = min(d('nThrust', 40), size(r, 1));
ix = round(linspace(1, size(r, 1), nT));
sc = 0.06*max(range(r(:,1)), range(r(:,2)));
quiver3(ax, r(ix,1), r(ix,2), r(ix,3), sc*alpha(ix,1), sc*alpha(ix,2), sc*alpha(ix,3), 0, ...
    'Color', [1.00 0.92 0.55], 'LineWidth', 0.9, 'MaxHeadSize', 0.5, ...
    'DisplayName', 'thrust direction');
P.nThrust = nT;

plot3(ax, r(1,1), r(1,2), r(1,3), 'o', 'MarkerSize', 9, 'LineWidth', 1.6, ...
      'MarkerEdgeColor', [0.35 1 0.55], 'DisplayName', 'departure');
plot3(ax, r(end,1), r(end,2), r(end,3), 'p', 'MarkerSize', 14, 'LineWidth', 1.4, ...
      'MarkerFaceColor', [1 0.85 0.3], 'MarkerEdgeColor', 'k', 'DisplayName', 'arrival');
% The Moon only. Both orbits and the whole transfer sit around the Moon, so
% the Earth marker at x = -mu just stretches the equal-aspect box and shrinks
% everything worth looking at.
plot3(ax, 1-mu, 0, 0, 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.7 0.7 0.72], ...
      'MarkerEdgeColor', 'none', 'DisplayName', 'Moon');

grid(ax, 'on');  box(ax, 'off');
axis(ax, 'equal');  ax.DataAspectRatio = [1 1 1];
view(ax, d('view', [-35 22]));
rotate3d(fig, 'on');                      % the point of this figure
xlabel(ax, 'x [ND, rotating frame]');  ylabel(ax, 'y [ND]');  zlabel(ax, 'z [ND]');
lg = legend(ax, 'Location', 'northeast');
lg.TextColor = pick(dark, 'w', 'k');  lg.Color = pick(dark, [0.1 0.1 0.1], 'w');
lg.EdgeColor = pick(dark, [0.3 0.3 0.3], [0.7 0.7 0.7]);

% RECOMPUTED from the flight that is drawn, not copied from the certificate,
% so the title cannot disagree with the picture.
P.tfDays = tu(end)*tStar/86400;
P.dvKms  = B.cnd*log(1/Y(end,7))*lStar/tStar;
P.propellantKg = pr.m0kg*(1 - Y(end,7));
P.finalMassKg  = pr.m0kg*Y(end,7);
title(ax, {sprintf(['Minimum-time DRO \\rightarrow %d-petal tulip   |   ' ...
                    '%.0f mN, I_{sp} %g s, %g kg'], pr.NpTulip, pr.thrustN*1000, pr.ispS, pr.m0kg), ...
           sprintf(['s_D = %.4f, s_A = %.4f   |   t_f = %.4f d, \\DeltaV = %.4f km/s, ' ...
                    'propellant %.2f kg'], T.sD, T.sA, P.tfDays, P.dvKms, P.propellantKg)}, ...
      'Color', pick(dark, 'w', 'k'), 'FontWeight', 'normal');

P.fig = fig;  P.ax = ax;
P.startMissKm = norm(r(1,:) - rv0(1:3)')*lStar;
P.endMissKm   = norm(r(end,:) - rvf(1:3)')*lStar;
if ~isempty(d('outPng', '')), exportgraphics(fig, d('outPng', ''), 'Resolution', 150); end
end

function v = pick(c, a, b)
% PICK  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
