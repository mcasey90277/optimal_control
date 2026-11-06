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

%% 5  relay closure

%% 6  hypothetical reverse graph

%% 7  figures

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
