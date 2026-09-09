function fig = plot_arrival_arcs(files, S, outPng)
%% Purpose:
%
%   Draw the arrival-phase branch map: t_f against arrival phase for every
%   arclength arc (one colour each, phase wrapped to [0,1)), folds as
%   filled markers, grid crossings as open circles, and -- when a sheet is
%   given -- the certified minimum per grid point as a black diamond.
%
%% Inputs:
%
%  files                    cellstr | struct array  arc .mat files (each
%                                                   holding A) or the A
%                                                   structs themselves
%  S                        struct (optional)       sheet (build_arrival_sheet)
%  outPng                   char (optional)         save path
%
%% Outputs:
%
%  fig                      figure handle
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, S = []; end
if nargin < 3, outPng = ''; end
tStar = 382981.289129055;
fig = figure('Color', 'w', 'Position', [100 100 1100 600]);  hold on
cols = lines(numel(files));
for k = 1:numel(files)
    if iscell(files), L = load(files{k});  A = L.A;  nm = files{k};
    else,             A = files(k);  nm = sprintf('arc %d', k); end
    ctf = A.anc.ctf;
    tf = cellfun(@(p) p(ctf), A.p)*tStar/86400;
    q = mod(A.q, 1);
    % break the line where the phase wraps
    br = [false, abs(diff(q)) > 0.5];  qq = q;  tt = tf;  qq(br) = NaN;  tt(br) = NaN;
    [~, nm2] = fileparts(nm);
    plot(qq, tt, '-', 'Color', cols(k,:), 'LineWidth', 1.2, 'DisplayName', strrep(nm2, '_', '\_'));
    for f = A.folds
        plot(mod(f.q, 1), interp1(1:numel(tf), tf, f.index), 'v', 'MarkerFaceColor', cols(k,:), ...
             'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'HandleVisibility', 'off');
    end
    for c = A.crossings
        if c.converged
            plot(mod(c.level, 1), c.p(ctf)*tStar/86400, 'o', 'Color', cols(k,:), 'MarkerSize', 7, 'HandleVisibility', 'off');
        end
    end
end
if ~isempty(S)
    m = isfinite(S.TF);
    plot(S.sA(m), S.TF(m), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'certified minimum');
end
grid on;  xlabel('arrival phase s_A (fraction of tulip period)');  ylabel('t_f [days]');
title('DRO \rightarrow tulip min-time, 70 mN / Isp 900 s / s_D = 0: arrival-phase branches');
legend('Location', 'best');
if ~isempty(outPng), exportgraphics(fig, outPng, 'Resolution', 150); end
end
