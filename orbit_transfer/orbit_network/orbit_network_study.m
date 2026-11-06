%% ORBIT_NETWORK_STUDY  The low-thrust orbital network over the costate catalogs
%
%   Braik & Ross (arXiv 2605.31543, 2026) organize periodic orbits into a
%   weighted graph whose edges come from reachable-set overlap under
%   heading-change impulses, then read hubs, relays and budget regimes off
%   it. Here the nodes are the 18 pumpkyn orbits of our four min-time
%   catalogs and the edges are the 56 sheets, every number over ACCEPTED
%   cells (certified witnesses). No proxies, no new solves. The edge
%   attributes are catalog coverage across thrust: the lowest accepted
%   rung and a phase-resolved time of flight; their (dV cap, time cap)
%   budget plane becomes our (thrust capability, time cap) plane, shown
%   as raw coverage and as the known-feasible envelope.
%
%   What this data cannot say, said once: an unaccepted cell or rung is
%   unresolved, not infeasible; no pair has both a direct edge and a
%   two-leg route, so relay-versus-direct cost is not tested (the lemma in
%   the note: the unrestricted infimum is never above the relay infimum);
%   relay CLOSURE is tested, with the phase-compatible assembled time, on
%   the only multileg structure the catalogs hold, L1 halo -> L2 halo ->
%   tulip.
%
%     0  settings                          (one place for every choice)
%     1  load                              (G1 periods; G2 edges, keys, 1 N / 0.5 N nesting)
%     2  graph and edge tables             (digraph; per-rung accepted fraction, t_f, lowest accepted rung)
%     3  planes                            (raw coverage and envelope; G3 monotone in the cap)
%     4  strengths                         (in/out strength at 1 N; a topology caveat)
%     5  relay closure                     (phase-compatible assembled times; G4)
%     6  hypothetical reverse graph        (G5 symmetry at an off-plane point; HYPOTHETICAL edges)
%     7  figures                           (results/*.png)
%
%   G gates are code gates: a FAIL of G1-G4 means the data or the
%   computation is not what this script assumes, and the summary throws;
%   G5 throws only when the hypothetical graph is on. Every printed number
%   from section 6 carries the word HYPOTHETICAL and never enters sections
%   3-5 or the figures.
%
%   A pipeline test runs this headless: if a struct `study_opts` exists in
%   the base workspace, its .outDir, .visible and .hypothetical are
%   honoured.
%
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.

%% 0  settings
clearvars -except study_opts; clc
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'lib'), fullfile(fileparts(here), 'costate_common'));

F_ref        = 1;                 % N, the rung every sheet has an accepted cell at: strengths and the graph figure use it
F_deep       = 0.5;               % N, the catalogs' deepest rung
nCaps        = 25;                % time caps, from the smallest accepted cell time to the largest (or the largest relay time)
hypothetical = false;             % the reverse graph by the x-axis symmetry: geometry only, off by default
relay        = struct('src', 'halo_L1', 'via', 'halo_L2', 'dst', 'tulip');   % key prefixes of the only multileg structure in the data
outDir       = fullfile(here, 'results');
visible      = 'on';
if exist('study_opts', 'var')
    if isfield(study_opts, 'outDir'),       outDir       = study_opts.outDir;       end
    if isfield(study_opts, 'visible'),      visible      = study_opts.visible;      end
    if isfield(study_opts, 'hypothetical'), hypothetical = study_opts.hypothetical; end
end
if ~exist(outDir, 'dir'), mkdir(outDir); end
gates = struct('id', {}, 'ok', {}, 'msg', {}, 'throws', {});
lines = {};                                                 % the summary file, one printed fact per line

%% 1  load
T     = load_catalog_sheets([]);
N     = orbit_node_table(T, struct('verbose', false));
nS    = numel(T);  nN = numel(N);
rungs = T(1).rungs_N;  nR = numel(rungs);
rC    = find(rungs == F_ref, 1);  rD = find(rungs == F_deep, 1);
tStar = T(1).constants.tStar_s;
keyN  = {N.key};
fprintf('\n== 1  load: %d sheets, %d orbits, rungs %s N\n', nS, nN, mat2str(rungs));
ePer = max([abs(cellfun(@(k) N(strcmp(keyN, k)).period_nd, {T.dep_key}) - [T.T_dep]), ...
            abs(cellfun(@(k) N(strcmp(keyN, k)).period_nd, {T.arr_key}) - [T.T_arr])]);
nNest = 0;  nViolAll = 0;  nAccX = 0;
for kk = 1:nS
    a = T(kk).accepted;  h = T(kk).has_solution;
    nNest    = nNest    + nnz(h(:, :, rD) & ~h(:, :, rC));   % SOLVED sets nest (the ladder descends from 1 N)
    nAccX    = nAccX    + nnz(a(:, :, rD) & ~a(:, :, rC));   % accepted at 0.5 N but refuted/gated at 1 N: a fact, not a gate
    nViolAll = nViolAll + nnz(h(:, :, 2:end) & ~h(:, :, 1:end-1));
end
lines{end+1} = sprintf('cells accepted at %g N but not at %g N (refuted or gated at the shallower rung): %d', F_deep, F_ref, nAccX);
fprintf('%s\n', lines{end});
perSolved   = sum(cell2mat(arrayfun(@(r) reshape(sum(sum(r.has_solution, 1), 2), 1, []), T, 'UniformOutput', false)), 1);
perAccepted = sum(cell2mat(arrayfun(@(r) reshape(sum(sum(r.accepted, 1), 2), 1, []), T, 'UniformOutput', false)), 1);
lines{end+1} = sprintf('solved cells per rung %s: %s', mat2str(rungs), mat2str(perSolved));
lines{end+1} = sprintf('accepted cells per rung: %s (certified witnesses; every count below uses these)', mat2str(perAccepted));
lines{end+1} = sprintf('whole-ladder non-monotone solved cells: %d (the 15 N starting rung closes fewer than 12 N: solver coverage, not gated)', nViolAll);
fprintf('%s\n%s\n%s\n', lines{end-2}, lines{end-1}, lines{end});
gates(end+1) = struct('id', 'G1', 'ok', ePer < 1e-9, 'msg', sprintf('node periods match the catalogs (max |dT| = %.1e)', ePer), 'throws', true);

%% 2  graph and edge tables
E = edge_tables(T);
G = network_graph(E, N);
src = findnode(G, {T.dep_key});  dst = findnode(G, {T.arr_key});          % node indices in SHEET order (used only for key checks)
fprintf('\n== 2  graph: %d nodes, %d edges (every attribute read from G.Edges; MATLAB sorts edges)\n', numnodes(G), numedges(G));
for rr = [rC rD]
    fprintf('\n  edge table at %g N (accepted cells)\n', rungs(rr));
    fprintf('  %-18s %-16s %7s %6s %8s %8s %20s\n', 'departure', 'arrival', 'frac', 'solved', 'tf_min', 'tf_med', 'lowest accepted');
    for kk = 1:nS
        fprintf('  %-18s %-16s %7.3f %6d %8.3f %8.3f %20s\n', E(kk).dep_key, E(kk).arr_key, E(kk).closed_frac(rr), E(kk).n_solved(rr), E(kk).tf_min_days(rr), E(kk).tf_med_days(rr), E(kk).F_low_str);
    end
end
nLowest = nnz(contains({E.F_low_str}, 'lowest tested'));
lines{end+1} = sprintf('edges accepted at the lowest tested rung (%g N): %d of %d; edges with at least one accepted cell at %g N: %d of %d', F_deep, nLowest, nS, F_ref, nnz(arrayfun(@(e) e.n_accepted(rC) > 0, E)), nS);
fprintf('%s\n', lines{end});
gates(end+1) = struct('id', 'G2', 'ok', nS == 56 && numedges(G) == 56 && all(src > 0) && all(dst > 0) && nNest == 0, ...
                      'msg', sprintf('56 edges on named nodes (findnode round-trips every key); every cell solved at %g N is solved at %g N (%d violations)', F_deep, F_ref, nNest), 'throws', true);

%% 3  planes: raw coverage and the known-feasible envelope
tfCells = [];
for kk = 1:nS
    tf = T(kk).tf_nd(T(kk).accepted)*tStar/86400;
    tfCells = [tfCells; tf(:)];
end
capLo = min(tfCells);  capHi = max(tfCells);                % section 5 may extend capHi for the relay times
caps  = linspace(capLo, capHi, nCaps);
P     = accessibility_plane(T, caps);
fprintf('\n== 3  planes: %d rungs x %d caps (%.3f .. %.3f days)\n', nR, nCaps, caps(1), caps(end));
fprintf('  raw coverage, pairs (rows = rungs, cols = caps):\n');  disp(P.pairs);
fprintf('  known-feasible envelope, pairs:\n');  disp(P.pairs_env);
monoCap  = all(all(diff(P.pairs, 1, 2) >= 0)) && all(all(diff(P.cells, 1, 2) >= 0)) && all(all(diff(P.pairs_env, 1, 2) >= 0)) && all(all(diff(P.cells_env, 1, 2) >= 0));
monoEnv  = all(all(diff(P.pairs_env, 1, 1) <= 0)) && all(all(diff(P.cells_env, 1, 1) <= 0));
cFull = find(P.pairs(rC, :) == nS, 1, 'first');
if isempty(cFull)
    lines{end+1} = sprintf('at %g N no sampled cap gives all %d direct pairs (max %d of %d under %.3f days)', F_ref, nS, max(P.pairs(rC, :)), nS, caps(end));
else
    lines{end+1} = sprintf('at %g N all %d direct pairs are covered from the first sampled cap of %.3f days', F_ref, nS, caps(cFull));
end
lines{end+1} = sprintf('at %g N, raw coverage %d of %d pairs under the largest cap; envelope at the %g N bound %d of %d', F_deep, P.pairs(rD, end), nS, rungs(1), P.pairs_env(1, end), nS);
fprintf('%s\n%s\n', lines{end-1}, lines{end});
gates(end+1) = struct('id', 'G3', 'ok', monoCap && monoEnv, 'msg', 'raw and envelope counts non-decreasing along the cap; envelope non-increasing toward deeper rungs', 'throws', true);

%% 4  strengths at the reference rung (descriptive; the between-family ranking is the catalog's topology)
tfRef = G.Edges.tf_min_days(:, rC);
eSrc  = findnode(G, G.Edges.EndNodes(:, 1));  eDst = findnode(G, G.Edges.EndNodes(:, 2));
outS  = accumarray(eSrc, 1./tfRef, [nN 1]);
inS   = accumarray(eDst, 1./tfRef, [nN 1]);
fprintf('\n== 4  strengths at %g N (weights 1 / t_f[days]); each tulip has 12 in-edges, each L2 halo 2 in and 4 out: rankings reflect the catalog, not transport\n', F_ref);
hdr = sprintf('  %-20s %4s %4s %9s %9s', 'orbit', 'in', 'out', 'in-str', 'out-str');
fprintf('%s\n', hdr);  lines{end+1} = ['strengths at the reference rung:' hdr];
for kk = 1:nN
    row = sprintf('  %-20s %4d %4d %9.3f %9.3f', G.Nodes.Name{kk}, indegree(G, kk), outdegree(G, kk), inS(kk), outS(kk));
    fprintf('%s\n', row);  lines{end+1} = row;
end

%% 5  relay closure (L1 halo -> L2 halo -> tulip: the only multileg structure in the data)
R = relay_closure(T, N, relay.src, relay.via, relay.dst);
fprintf('\n== 5  relay closure: %d routes %s -> %s -> %s, phase-compatible assembled times\n', numel(R), relay.src, relay.via, relay.dst);
pairKeys = unique(strcat({R.src_key}, ' -> ', {R.dst_key}), 'stable');
best = struct('pair', pairKeys, 't_route_days', [], 't_sum_days', [], 'relay_key', [], 'deepest_N', []);
for pp = 1:numel(pairKeys)
    sel = find(strcmp(strcat({R.src_key}, ' -> ', {R.dst_key}), pairKeys{pp}));
    tr  = vertcat(R(sel).t_route_days);                     % [nRelays x nR]
    ts  = vertcat(R(sel).t_sum_days);
    [tb, ib] = min(tr, [], 1, 'omitnan');                   % best relay per rung by the ASSEMBLED time
    tb(all(isnan(tr), 1)) = NaN;
    best(pp).t_route_days = tb;
    best(pp).t_sum_days   = min(ts, [], 1, 'omitnan');
    best(pp).relay_key    = arrayfun(@(rr) R(sel(ib(rr))).relay_key, 1:nR, 'UniformOutput', false);
    best(pp).relay_key(isnan(tb)) = {''};                  % an open rung has no relay (min over all-NaN returns index 1)
    deep = find(~isnan(tb), 1, 'last');
    if isempty(deep), best(pp).deepest_N = NaN; else, best(pp).deepest_N = rungs(deep); end
    if isnan(tb(rC)), rk = 'none'; else, rk = best(pp).relay_key{rC}; end
    row = sprintf('  %-32s deepest rung closed %5g N | at %g N: best relay %-16s t_route %.3f d (t_sum %.3f d)', ...
                  pairKeys{pp}, best(pp).deepest_N, F_ref, rk, tb(rC), best(pp).t_sum_days(rC));
    fprintf('%s\n', row);  lines{end+1} = row;
end
% extend the caps if the relay times exceed the direct ones, then the directed-reachability overlay
tRel = [best.t_route_days];
if max(tRel, [], 'omitnan') > caps(end)
    caps = linspace(capLo, max(tRel, [], 'omitnan'), nCaps);
    P = accessibility_plane(T, caps);
    lines{end+1} = sprintf('caps extended to %.3f days to cover the relay routes', caps(end));
    fprintf('%s\n', lines{end});
end
P.multileg = zeros(nR, nCaps);
for rr = 1:nR
    for cc = 1:nCaps
        P.multileg(rr, cc) = nnz(arrayfun(@(b) b.t_route_days(rr) <= caps(cc), best));
    end
end
P.reach = P.pairs + P.multileg;                            % directed reachable ordered pairs under the cap (direct + two-leg)
usage = zeros(nN, 1);
for pp = 1:numel(best)
    if ~isnan(best(pp).t_route_days(rC))
        kU = find(strcmp(keyN, best(pp).relay_key{rC}), 1);
        usage(kU) = usage(kU) + 1;
    end
end
for kk = find(usage > 0).'
    row = sprintf('  relay usage at %g N: %-16s serves %d of %d (L1 halo, tulip) pairs as the best assembled relay', F_ref, keyN{kk}, usage(kk), numel(best));
    fprintf('%s\n', row);  lines{end+1} = row;
end
lines{end+1} = sprintf('%d relay routes over %d (src, dst) pairs; pairs reaching %g N by relay: %d of %d; directed reachable ordered pairs at %g N under the largest cap: %d (direct %d + two-leg %d) of %d candidate pairs (direct sheets + relay pairs)', ...
                       numel(R), numel(best), F_deep, nnz([best.deepest_N] == F_deep), numel(best), F_ref, P.reach(rC, end), P.pairs(rC, end), P.multileg(rC, end), nS + numel(best));
fprintf('%s\n', lines{end});
% G4: enumeration and consistency with the sheets (not a physics check)
g4 = numel(R) == 32 && numel(best) == 8;
for kk = 1:numel(R)
    L1 = T(R(kk).s1);  L2 = T(R(kk).s2);
    legs = squeeze(any(any(L1.accepted, 1), 2)).' & squeeze(any(any(L2.accepted, 1), 2)).';
    g4 = g4 && strcmp(L1.arr_key, L2.dep_key) && strncmp(L1.dep_key, relay.src, numel(relay.src)) && strncmp(L2.arr_key, relay.dst, numel(relay.dst)) ...
            && isequal(R(kk).closed, legs) && all(isnan(R(kk).t_route_days(~legs))) && all(~isnan(R(kk).t_route_days(legs))) ...
            && all(R(kk).t_route_days(legs) >= R(kk).t_sum_days(legs) - 1e-9) && isequal(sort(L1.sA_frac), sort(L2.sD_frac));
end
gates(end+1) = struct('id', 'G4', 'ok', g4, 'msg', '32 routes over 8 pairs, legs share the relay key and its phase grid, closed iff both legs accepted, t_route >= t_sum and NaN exactly when open', 'throws', true);

%% 6  hypothetical reverse graph (x-axis symmetry: geometry only; the mirrored arc belongs to a mass-GAIN model)
mu    = T(1).constants.muStar;
resP  = zeros(nN, 1);  resS = zeros(nN, 1);
for kk = 1:nN
    x0 = N(kk).x0;  Tp = N(kk).period_nd;
    [~, X] = pumpkyn.cr3bp.prop([0 0.37*Tp Tp], x0.', mu);
    resP(kk) = max(abs(X(end, 1:6).' - x0));                % periodic closure, independent of the symmetry
    resS(kk) = symmetry_residual(X(2, 1:6).', Tp, mu);      % mirror image of an OFF-PLANE point lies on the orbit
end
fprintf('\n== 6  symmetry: periodic closure max %.1e; off-plane mirror residual max %.1e (every stored x0 is on the plane, so the check is made at phi(x0, 0.37 T))\n', max(resP), max(resS));
barG5 = 1e-6;                                               % pumpkyn.cr3bp.prop runs at 1e-10 tolerances: the longest tulip closes to 8e-8; a non-symmetric state fails by 1e-3
gates(end+1) = struct('id', 'G5', 'ok', all(resP < barG5) && all(resS < barG5), ...
                      'msg', sprintf('every node closes (%.1e) and is x-axis symmetric at an off-plane point (%.1e), bar %.0e (propagator floor 8e-8 on tulip_Np12); throws only when the hypothetical graph is on', max(resP), max(resS), barG5), 'throws', hypothetical);
if hypothetical
    Eh = G.Edges;
    Eh.EndNodes = fliplr(Eh.EndNodes);
    Eh.hypothetical = true(height(Eh), 1);
    Eh.n_accepted   = zeros(size(Eh.n_accepted));          % no certificate travels with a mirrored arc
    Gh = digraph([G.Edges; Eh], G.Nodes);
    [~, sccSize] = conncomp(Gh, 'Type', 'strong');
    D  = distances(Gh);  D(logical(eye(nN))) = Inf;
    nReach = nnz(isfinite(D));
    row = sprintf('  HYPOTHETICAL: %d edges (%d + %d reversed); strongly connected components %d (largest %d nodes); ordered node pairs joined by a directed path %d of %d', ...
                  numedges(Gh), nS, nS, numel(sccSize), max(sccSize), nReach, nN*(nN - 1));
    fprintf('%s\n', row);  lines{end+1} = row;
    lines{end+1} = 'HYPOTHETICAL numbers never enter sections 3-5 or the figures';
else
    fprintf('  hypothetical reverse graph off (hypothetical = false); G5 is still reported\n');
end

%% 7  figures (never the hypothetical edges)
axStyle = @(ax) set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
rungLab = arrayfun(@(r) sprintf('%g', r), rungs, 'UniformOutput', false);
famOf   = @(k) regexprep(k, '_tau.*|_Np.*', '');          % dro, dpo, halo_L1, halo_L2, tulip
famList = {'dro', 'dpo', 'halo_L1', 'halo_L2'};

% F1: the graph at the reference rung and at the deepest rung; width = accepted fraction, color = t_f; built from the edge table
fig = figure('Color', 'w', 'Position', [60 60 1400 620], 'Visible', visible);
for pp = 1:2
    rr = [rC rD];  rr = rr(pp);
    Er = G.Edges(G.Edges.closed_frac(:, rr) > 0, :);        % accepted edges at this rung, attributes travel with the rows
    Gr = digraph(Er, G.Nodes);
    ax = subplot(1, 2, pp);  axStyle(ax);
    h  = plot(ax, Gr, 'Layout', 'layered', 'Direction', 'right', 'NodeLabel', Gr.Nodes.Name, 'Interpreter', 'none', ...
              'MarkerSize', 6, 'NodeColor', [0.1 0.3 0.7], 'ArrowSize', 8);
    h.NodeLabelColor = 'k';
    h.LineWidth = 0.5 + 4*Gr.Edges.closed_frac(:, rr);
    h.EdgeCData = Gr.Edges.tf_min_days(:, rr);
    axStyle(ax);  xl = xlim(ax);  xlim(ax, [xl(1) xl(2) + 0.2*(xl(2) - xl(1))]);       % plot() resets the axes colors; room for the right-hand labels
    colormap(ax, 'parula');  cb = colorbar(ax);  cb.Color = 'k';  cb.Label.String = 'min t_f over accepted phases [days]';
    title(ax, sprintf('%g N: %d edges, width = accepted fraction of the 6 x 6 grid', rungs(rr), numedges(Gr)), 'Color', 'k', 'Interpreter', 'none');
end
f1 = fullfile(outDir, 'network_graph.png');  exportgraphics(fig, f1, 'Resolution', 150);  if strcmp(visible, 'off'), close(fig); end

% F2: raw coverage and the envelope, with boolean-mask contours (all direct pairs; first relay route)
fig = figure('Color', 'w', 'Position', [60 60 1400 560], 'Visible', visible);
planes = {P.pairs, 'raw coverage: pairs with an accepted cell AT this thrust under the cap'; P.pairs_env, 'known-feasible envelope: pairs with an accepted cell at this thrust OR LOWER under the cap'};
for pp = 1:2
    ax = subplot(1, 2, pp);  axStyle(ax);  hold(ax, 'on');
    imagesc(ax, caps, 1:nR, planes{pp, 1});  axis(ax, 'tight');  set(ax, 'YDir', 'normal', 'YTick', 1:nR, 'YTickLabel', rungLab);
    colormap(ax, 'hot');  cb = colorbar(ax);  cb.Color = 'k';  cb.Label.String = sprintf('pairs (of %d)', nS);
    mAll = double(planes{pp, 1} == nS);
    if any(mAll(:)) && ~all(mAll(:)), contour(ax, caps, 1:nR, mAll, [0.5 0.5], 'c-', 'LineWidth', 2); end
    mRel = double(P.multileg >= 1);
    if any(mRel(:)) && ~all(mRel(:)), contour(ax, caps, 1:nR, mRel, [0.5 0.5], 'g--', 'LineWidth', 2); end
    xlabel(ax, 'time cap [days]', 'Color', 'k');  ylabel(ax, 'thrust bound [N]', 'Color', 'k');
    title(ax, {planes{pp, 2}, 'cyan: all 56 direct pairs; green dashed: first relay route (assembled time)'}, 'Color', 'k', 'FontSize', 9);
end
f2 = fullfile(outDir, 'accessibility_plane.png');  exportgraphics(fig, f2, 'Resolution', 150);  if strcmp(visible, 'off'), close(fig); end

% F3: catalog closure down the ladder, by departure family (accepted fraction; the 15 N dip is ladder coverage)
fig = figure('Color', 'w', 'Position', [60 60 900 480], 'Visible', visible);
ax  = axes(fig);  axStyle(ax);  hold(ax, 'on');  grid(ax, 'on');
cf  = vertcat(E.closed_frac);  depFam = cellfun(famOf, {E.dep_key}, 'UniformOutput', false);
for ff = 1:numel(famList)
    sel = strcmp(depFam, famList{ff});
    plot(ax, 1:nR, mean(cf(sel, :), 1), '-o', 'LineWidth', 1.8, 'DisplayName', sprintf('%s departures (%d edges)', famList{ff}, nnz(sel)));
end
set(ax, 'XTick', 1:nR, 'XTickLabel', rungLab);
xlabel(ax, 'thrust rung [N] (shallow to deep)', 'Color', 'k');  ylabel(ax, 'mean accepted fraction of the phasing grid', 'Color', 'k');
legend(ax, 'Location', 'southwest', 'Interpreter', 'none', 'TextColor', 'k', 'Color', 'w');
title(ax, 'catalog closure down the ladder (an unaccepted cell is unresolved, not infeasible; the 15 N dip is ladder coverage)', 'Color', 'k', 'FontSize', 9);
f3 = fullfile(outDir, 'ladder_closure.png');  exportgraphics(fig, f3, 'Resolution', 150);  if strcmp(visible, 'off'), close(fig); end

% F4: relay closure, best route per (L1 halo, tulip) pair: assembled time (solid) and the optimistic sum (dotted)
fig = figure('Color', 'w', 'Position', [60 60 900 520], 'Visible', visible);
ax  = axes(fig);  axStyle(ax);  hold(ax, 'on');  grid(ax, 'on');
cols = parula(numel(best));                                % NOT lines(): that name is the script's summary cellstr
for pp = 1:numel(best)
    plot(ax, 1:nR, best(pp).t_route_days, '-o', 'Color', cols(pp, :), 'LineWidth', 1.6, 'DisplayName', best(pp).pair);
    plot(ax, 1:nR, best(pp).t_sum_days,   ':',  'Color', cols(pp, :), 'LineWidth', 1.0, 'HandleVisibility', 'off');
end
set(ax, 'XTick', 1:nR, 'XTickLabel', rungLab);
xlabel(ax, 'thrust rung [N]', 'Color', 'k');  ylabel(ax, 'two-leg time [days] (solid assembled, dotted sum)', 'Color', 'k');
legend(ax, 'Location', 'northwest', 'Interpreter', 'none', 'TextColor', 'k', 'Color', 'w', 'FontSize', 8);
title(ax, 'relay closure L1 halo -> L2 halo -> tulip: best assembled relay per pair and rung', 'Color', 'k', 'Interpreter', 'none');
f4 = fullfile(outDir, 'relay_closure.png');  exportgraphics(fig, f4, 'Resolution', 150);  if strcmp(visible, 'off'), close(fig); end
lines{end+1} = sprintf('figures: %s, %s, %s, %s', f1, f2, f3, f4);
fprintf('%s\n', lines{end});

%% summary
fprintf('\n== SUMMARY\n');
for kk = 1:numel(gates)
    if gates(kk).ok, s = 'PASS'; else, s = 'FAIL'; end
    fprintf('  %s  %s  %s\n', gates(kk).id, s, gates(kk).msg);
end
fid = fopen(fullfile(outDir, 'network_summary.txt'), 'w');
fprintf(fid, 'orbit_network_study summary, %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
for kk = 1:numel(lines), fprintf(fid, '%s\n', lines{kk}); end
for kk = 1:numel(gates)
    if gates(kk).ok, s = 'PASS'; else, s = 'FAIL'; end
    fprintf(fid, '%s %s %s\n', gates(kk).id, s, gates(kk).msg);
end
fclose(fid);
fprintf('summary: %s\n', fullfile(outDir, 'network_summary.txt'));
if any(~[gates.ok] & [gates.throws])
    error('orbit_network:orbit_network_study:gate', 'a G gate failed; see the summary');
end
