%% STABILITY_COST_STUDY  Does orbit stability predict low-thrust transfer cost?
%
%   Braik & Ross (arXiv 2605.31543, 2026) read off their impulsive orbital
%   network that linearly STABLE periodic orbits are the hard nodes (they
%   sit in tori) while unstable ones share a chaotic sea and connect
%   cheaply. This script asks whether that survives in our setting: a
%   mass-varying MIN-TIME optimal control problem with certified extremals,
%   56 (departure orbit, arrival orbit) sheets across the four costate
%   catalogs, each sheet a 6 x 6 phasing grid walked down a thrust ladder
%   from 15 N to 0.5 N.
%
%   Each sheet is one point. Covariates: the Jacobi gap dCJ = C_J(arr) -
%   C_J(dep) (the energy the transfer must buy), and the instability rates
%   sigma_dep, sigma_arr = ln|lambda|max / T of the two orbits (exactly 0
%   when |lambda|max < 1 + 1e-3). Responses: tf1, the min-over-phases time
%   of flight at the shared 1 N rung, and surv, the fraction of cells
%   accepted at 1 N that are still accepted at the 0.5 N rung (catalog
%   retention down the ladder). The question is whether sigma_dep carries
%   information BEYOND energy: a nested regression per response, dCJ alone
%   against dCJ + sigma_dep + sigma_arr, with the verdict on the sigma_dep
%   coefficient's own t-test (the partial F and sigma_arr beside it).
%
%     0  settings                          (one place for every choice)
%     1  load and census                   (P1-P4: the data is what we think)
%     2  orbit properties                  (P5-P7: the spectra are trustworthy)
%     3  covariates and responses          (the 56-row table, printed)
%     4  regression                        (H1-H2 verdicts)
%     5  figure                            (results/stability_cost.png)
%
%   Every count and minimum is over ACCEPTED cells (solved, conjugate-pass,
%   hypothesis gates); an unaccepted cell is unresolved, not infeasible.
%   P gates are code gates: a FAIL means the data or the computation is
%   not what this script assumes, and the summary throws. H verdicts are
%   HYPOTHESES: DETECTED or NOT DETECTED (nominal OLS) are both results,
%   and NOT DETECTED is never "no effect".
%
%   A pipeline test runs this headless: if a struct `study_opts` exists in
%   the base workspace, its .outDir and .visible are honoured.
%
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.

%% 0  settings
clearvars -except study_opts; clc
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'lib'), fullfile(fileparts(here), 'costate_common'));

F_cost   = 1;                     % N, the shared cost rung (every sheet closes somewhere here)
F_wall   = 0.5;                   % N, the deepest rung the catalogs hold
alpha    = 0.05;                  % significance level for the partial F test
outDir   = fullfile(here, 'results');
visible  = 'on';
if exist('study_opts', 'var')
    if isfield(study_opts, 'outDir'),  outDir  = study_opts.outDir;  end
    if isfield(study_opts, 'visible'), visible = study_opts.visible; end
end
if ~exist(outDir, 'dir'), mkdir(outDir); end
gates = struct('id', {}, 'ok', {}, 'msg', {});

%% 1  load and census
T  = load_catalog_sheets([]);
nS = numel(T);
rungs = T(1).rungs_N;
rC = find(rungs == F_cost, 1);
rW = find(rungs == F_wall, 1);
tStar = T(1).constants.tStar_s;
days  = @(tnd) tnd*tStar/86400;

fprintf('\n== 1  census: %d sheets, rungs %s N (counts are ACCEPTED cells: solved, conjugate-pass, hypothesis gates; solved beside)\n', nS, mat2str(rungs));
fprintf('%-3s %-11s %-18s %-14s %6s %6s %6s %6s %9s\n', '#', 'catalog', 'departure', 'arrival', 'a(1N)', 's(1N)', 'a(.5N)', 's(.5N)', 'tf1 [d]');
nC = zeros(nS, 1);  nW = zeros(nS, 1);  tf1 = nan(nS, 1);  sC = zeros(nS, 1);  sW = zeros(nS, 1);  nX = 0;
for kk = 1:nS
    aC = T(kk).accepted(:, :, rC);
    aW = T(kk).accepted(:, :, rW);
    tC = T(kk).tf_nd(:, :, rC);
    nC(kk) = nnz(aC);  nW(kk) = nnz(aW & aC);                % retained = accepted at BOTH rungs (nested by construction)
    nX = nX + nnz(aW & ~aC);                                  % accepted at the wall rung but refuted/gated at the cost rung: printed, not counted
    sC(kk) = nnz(T(kk).has_solution(:, :, rC));  sW(kk) = nnz(T(kk).has_solution(:, :, rW));
    if nC(kk) > 0, tf1(kk) = days(min(tC(aC))); end
    fprintf('%-3d %-11s %-18s %-14s %6d %6d %6d %6d %9.3f\n', kk, T(kk).catalog, T(kk).dep_key, T(kk).arr_key, nC(kk), sC(kk), nW(kk), sW(kk), tf1(kk));
end

% P1 the sheet count; P2 one rung vector (the loader refuses otherwise; re-stated here);
% P3 the SOLVED sets are nested between the two rungs used (the ladder walks
%    down from the cost rung, so a wall-rung solution descends from a cost-
%    rung one); the ACCEPTED sets need not be, since certification can fail
%    at the shallower rung and pass at the deeper one, so surv counts cells
%    accepted at BOTH rungs over cells accepted at the cost rung and the
%    exceptions are printed; across the WHOLE ladder closure is not
%    monotone -- the 15 N starting rung closes fewer cells than 12 N --
%    printed as a fact;
% P4 every sheet has an accepted cell at the cost rung
nViolPair = 0;  nViolAll = 0;
for kk = 1:nS
    h = T(kk).has_solution;
    nViolPair = nViolPair + nnz(h(:, :, rW) & ~h(:, :, rC));
    nViolAll  = nViolAll  + nnz(h(:, :, 2:end) & ~h(:, :, 1:end-1));   % solved deeper but open one rung shallower
end
fprintf('solved cells per rung, summed over sheets:   %s\n', mat2str(sum(cell2mat(arrayfun(@(r) reshape(sum(sum(r.has_solution, 1), 2), 1, []), T, 'UniformOutput', false)), 1)));
fprintf('accepted cells per rung, summed over sheets: %s\n', mat2str(sum(cell2mat(arrayfun(@(r) reshape(sum(sum(r.accepted, 1), 2), 1, []), T, 'UniformOutput', false)), 1)));
fprintf('whole-ladder non-monotone solved cells (closed at a rung, open one rung shallower): %d -- solver coverage, not a data error; an unclosed cell is unresolved, not infeasible\n', nViolAll);
fprintf('cells accepted at %g N but not at %g N (refuted or gated at the shallower rung): %d -- excluded from surv, which counts cells accepted at both rungs\n', F_wall, F_cost, nX);
gates(end+1) = struct('id', 'P1', 'ok', nS == 56, 'msg', sprintf('56 sheets loaded (%d)', nS));
gates(end+1) = struct('id', 'P2', 'ok', all(arrayfun(@(r) isequal(r.rungs_N, rungs), T)) && ~isempty(rC) && ~isempty(rW), 'msg', sprintf('one rung vector, holding %g N and %g N', F_cost, F_wall));
gates(end+1) = struct('id', 'P3', 'ok', nViolPair == 0, 'msg', sprintf('every cell solved at %g N is solved at %g N (%d violations); surv is the nested ratio accepted-at-both / accepted-at-%g N', F_wall, F_cost, nViolPair, F_cost));
gates(end+1) = struct('id', 'P4', 'ok', all(nC > 0), 'msg', sprintf('every sheet has an accepted cell at %g N (min %d cells)', F_cost, min(nC)));

%% 2  orbit properties
fprintf('\n== 2  orbit properties (sigma = ln|lambda|max / T, 0 when stable)\n');
N = orbit_node_table(T, struct('verbose', false));
fprintf('%-20s %8s %10s %12s %8s %7s\n', 'orbit', 'T [ND]', 'C_J', '|lam|max', 'sigma', 'stable');
for kk = 1:numel(N)
    fprintf('%-20s %8.4f %10.5f %12.4e %8.4f %7d\n', N(kk).key, N(kk).period_nd, N(kk).CJ, N(kk).lam_max, N(kk).sigma, N(kk).stable);
end
keyN = {N.key};
iDep = cellfun(@(k) find(strcmp(keyN, k), 1), {T.dep_key});
iArr = cellfun(@(k) find(strcmp(keyN, k), 1), {T.arr_key});

ePer  = max([abs([N(iDep).period_nd] - [T.T_dep]), abs([N(iArr).period_nd] - [T.T_arr])]);
eSym  = max(abs([N.lam_max].*[N.lam_min] - 1));
eDet  = max(abs([N.detM] - 1));
isTul = strcmp({T.arr_family}, 'tulip');
gates(end+1) = struct('id', 'P5', 'ok', eSym < 1e-6 && eDet < 1e-6, 'msg', sprintf('symplectic pairing |lam|max |lam|min = 1 (%.1e) and det M = 1 (%.1e) on all %d orbits', eSym, eDet, numel(N)));
gates(end+1) = struct('id', 'P6', 'ok', ePer < 1e-9, 'msg', sprintf('node periods match the catalogs (max |dT| = %.1e)', ePer));
gates(end+1) = struct('id', 'P7', 'ok', all([N(iArr(isTul)).stable]), 'msg', sprintf('every tulip arrival is linearly stable (%d sheets)', nnz(isTul)));

%% 3  covariates and responses
dCJ   = [N(iArr).CJ].' - [N(iDep).CJ].';
sDep  = [N(iDep).sigma].';
sArr  = [N(iArr).sigma].';
stDep = [N(iDep).stable].';
surv  = nW./nC;
fprintf('\n== 3  the sheet table (%d rows)\n', nS);
fprintf('%-3s %-18s %-14s %8s %8s %8s %9s %7s\n', '#', 'departure', 'arrival', 'dCJ', 'sig_dep', 'sig_arr', 'tf1 [d]', 'surv');
for kk = 1:nS
    fprintf('%-3d %-18s %-14s %8.4f %8.4f %8.4f %9.3f %7.3f\n', kk, T(kk).dep_key, T(kk).arr_key, dCJ(kk), sDep(kk), sArr(kk), tf1(kk), surv(kk));
end
fprintf('stable-departure sheets: %d, unstable-departure sheets: %d; distinct departure orbits %d, distinct arrivals %d (sheets share nodes and are not independent observations)\n', ...
        nnz(stDep), nnz(~stDep), numel(unique({T.dep_key})), numel(unique({T.arr_key})));

%% 4  regression
fprintf('\n== 4  nested regression, M0: y ~ 1 + dCJ   M1: y ~ 1 + dCJ + sigma_dep + sigma_arr\n');
X0 = [ones(nS, 1), dCJ];
X1 = [ones(nS, 1), dCJ, sDep, sArr];
resp = struct('name', {'tf1 [days] at 1 N', 'retention to 0.5 N'}, 'y', {tf1, surv}, 'sign', {-1, +1}, 'id', {'H1', 'H2'});
verdict = struct('id', cell(1, numel(resp)), 'supported', [], 'msg', '');
fits    = cell(1, numel(resp));
for rr = 1:numel(resp)
    o = partial_f_test(resp(rr).y, X0, X1);
    gS = mean(resp(rr).y(stDep));  gU = mean(resp(rr).y(~stDep));
    fprintf('\n  %s\n', resp(rr).name);
    fprintf('    M0: b = [%.4f %.4f]              R2 = %.3f\n', o.b0, o.R2_0);
    fprintf('    M1: b = [%.4f %.4f %.4f %.4f]  R2 = %.3f\n', o.b1, o.R2_1);
    fprintf('    partial F(%d, %d) = %.3f, p = %.4f   (both sigma terms jointly)\n', o.df1, o.df2, o.F, o.p);
    fprintf('    M1 coefficient t-tests (nominal OLS, %d dof): dCJ t = %6.2f p = %.4f | sigma_dep t = %6.2f p = %.4f | sigma_arr t = %6.2f p = %.4f\n', ...
            o.df2, o.t1(2), o.pt1(2), o.t1(3), o.pt1(3), o.t1(4), o.pt1(4));
    fprintf('    95%% CI: dCJ [%.4f, %.4f] | sigma_dep [%.4f, %.4f] | sigma_arr [%.4f, %.4f]\n', o.ci1(2, :), o.ci1(3, :), o.ci1(4, :));
    fprintf('    group means (descriptive): stable departure %.4f, unstable departure %.4f\n', gS, gU);
    % the paper's claim is about DEPARTURE stability beyond energy, so the
    % verdict rests on the sigma_dep coefficient's own test and sign; the
    % joint F and the sigma_arr term are reported beside it (sigma_arr only
    % varies on the 8 halo-to-halo sheets, so it is confounded with catalog).
    % NOT DETECTED means not detected by this analysis at the nominal
    % threshold, never "no effect": the p-values are nominal (sheets share
    % nodes) and 1 N / 0.5 N sit far above the 70 mN regime of the DRO wall.
    signOK = sign(o.b1(3)) == resp(rr).sign;
    sup = o.pt1(3) < alpha && signOK;
    if resp(rr).sign < 0, want = 'negative (unstable departures faster)'; else, want = 'positive (unstable departures retained better down the ladder)'; end
    verdict(rr) = struct('id', resp(rr).id, 'supported', sup, 'msg', sprintf('%s: sigma_dep coefficient %.4f, 95%% CI [%.4f, %.4f], t-test p = %.4f (nominal; alpha %.2f; joint F p = %.4f), paper predicts %s', ...
                         resp(rr).name, o.b1(3), o.ci1(3, 1), o.ci1(3, 2), o.pt1(3), alpha, o.p, want));
    fits{rr}    = o;
end

%% 5  figure
fig = figure('Color', 'w', 'Position', [100 100 1100 440], 'Visible', visible);
cols = [0.15 0.45 0.85; 0.85 0.30 0.10];
for rr = 1:numel(resp)
    ax = subplot(1, 2, rr);  hold(ax, 'on');  grid(ax, 'on');
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);   % the shared session may default axes to dark
    y  = resp(rr).y;
    scatter(ax, dCJ(stDep),  y(stDep),  55, cols(1, :), 'filled', 'DisplayName', 'stable departure');
    scatter(ax, dCJ(~stDep), y(~stDep), 55, sDep(~stDep), 'filled', 'MarkerEdgeColor', cols(2, :), 'LineWidth', 1.2, 'DisplayName', 'unstable departure (color = \sigma_{dep})');
    colormap(ax, 'autumn');
    xx = linspace(min(dCJ), max(dCJ), 50).';
    plot(ax, xx, [ones(50, 1), xx]*fits{rr}.b0, 'k-', 'LineWidth', 1.2, 'DisplayName', 'M0 fit (dCJ only)');
    xlabel(ax, '\Delta C_J = C_J(arr) - C_J(dep)', 'Color', 'k');
    ylabel(ax, resp(rr).name, 'Color', 'k');
    title(ax, sprintf('%s: sigma_{dep} t-test p = %.3f, nominal (joint F p = %.3f)', resp(rr).id, fits{rr}.pt1(3), fits{rr}.p), 'Color', 'k');
    legend(ax, 'Location', 'best', 'TextColor', 'k', 'Color', 'w');
end
figFile = fullfile(outDir, 'stability_cost.png');
exportgraphics(fig, figFile, 'Resolution', 150);
if strcmp(visible, 'off'), close(fig); end
fprintf('\nfigure: %s\n', figFile);

%% summary
fprintf('\n== SUMMARY\n');
for kk = 1:numel(gates)
    if gates(kk).ok, s = 'PASS'; else, s = 'FAIL'; end
    fprintf('  %s  %s  %s\n', gates(kk).id, s, gates(kk).msg);
end
for kk = 1:numel(verdict)
    if verdict(kk).supported, s = 'DETECTED    '; else, s = 'NOT DETECTED'; end
    fprintf('  %s  %s  %s\n', verdict(kk).id, s, verdict(kk).msg);
end
if ~all([gates.ok])
    error('orbit_network:stability_cost_study:gate', 'a P gate failed; see the summary');
end
