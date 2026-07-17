function run_ctf_sweep()
% RUN_CTF_SWEEP  M3 front: m_f vs c_tf at 10 N + thrust law across {10,5,2.5} N.
%
% Each point is one run_transfer call, saved individually and skipped when its
% results file exists (resume-after-crash pattern; the sporadic CasADi MEX
% fatal kills the process ~1 in 10 solves, and the first opti.solve() of a
% fresh MATLAB process after an idle gap has repeatedly died at IPOPT plugin
% init -- just rerun this script, appending the log; completed points skip).
%
% Leg 1: c_tf front at 10 N, neighbor-chain seeded upward from M2_manifold.mat
%        (c_tf = 1.2, 1.5, 2.0, 2.5, 3.0). run_transfer now only caches
%        CERTIFIED results (campaign rule, Task 14 controller triage): an
%        uncertified chain point warns and leaves no file, and the neighbor
%        pointer only advances past a point that actually saved -- so a bad
%        point can never poison the rest of the chain as a silent seed.
%        CONTROLLER-AUTHORIZED deviation from the brief's "reseed from the
%        other neighbor" rule: basin scatter was found across the WHOLE
%        upper tail (c_tf=1.2 uncertified; c_tf=2.5/3.0 certified but
%        dominated -- both neighbors of a scattered point can themselves be
%        scattered, so a same-strategy reseed cannot fix it). For c_tf in the
%        REDO set {1.2, 2.5, 3.0}, a second FRESH full-schedule solve (no
%        seedMat, cold tangential seed, tag '..._fresh') is also attempted;
%        the collection step keeps the better CERTIFIED m_f of {chain,
%        fresh} per point (best-of envelope, campaign pattern) with
%        provenance printed. This matches the paper's own Fig 18, which
%        shows exactly this kind of basin scatter in the same c_tf range.
%
% Dominated-basin recovery (Task 14, 2026-07-17): the densify pass exposed
% two monotonicity violations invisible on the sparser 5-point front:
% c_tf=1.35 (1360.21 kg) below c_tf=1.20 (1360.37 kg), and c_tf=2.00
% (1377.10 kg, an ORIGINAL chain point that predates the densify targets)
% below c_tf=1.75 (1379.86 kg). fix_dominated_basins.m (separate driver)
% attempted alternate-neighbor-seeded retries, tags '..._alt1'/'..._alt2':
% c_tf=1.35 alt1 (seeded UP from c120_fresh) recovered to 1368.20 kg,
% restoring monotonicity there. c_tf=2.00 alt1 (seeded from c175) reached
% 1378.71 kg and alt2 (seeded from c225) reached 1375.60 kg -- best-of is
% 1378.71 kg, which STILL sits 1.15 kg below c175's 1379.86 kg: a residual
% dominated point that survived both retry directions (reported, not
% hidden -- see task-14-report.md). The collection loop below folds
% alt1/alt2 into the same max-certified-m_f best-of rule as chain/fresh.
%
% Leg 2: thrust law at c_tf = 1.5, fresh fuel pipelines at T_max = 10, 5, 2.5 N
%        (N = 600, 1200, 2400). Each thrust needs its own min-time anchor
%        (run_mintime). Two earlier fix attempts (thrust-continuation
%        stretch-seeds, warm-starting a lower thrust's anchor from a
%        neighboring-thrust anchor's time-stretched trajectory, with and
%        without a stage-1 warm-up) both re-stalled or were found to have a
%        TOPOLOGY FLAW: stretching only rescales physical time, it does not
%        add revolutions, so a 4.5-rev shape was being forced toward an
%        ~8.4-rev min-time. CURRENT RECIPE (Task 14 controller triage round
%        3): every thrust solves its OWN cold tangential-seed min-time
%        anchor (run_mintime's default recipe, no seedAnchor) at
%        nodes-per-rev parity -- mesh scaled with that thrust's own rev
%        count (600 @ 10 N/~4.2 revs, 1200 @ 5 N/~8.4 revs, 2400 @ 2.5
%        N/~17 revs) -- with a recalibrated continuation guard (up to 24
%        rounds, 0.15-decade stall floor; ~0.24 decades/round observed at
%        5 N, ~5 min/round, so a full anchor solve can take ~2 h). The 10,
%        5, 2.5 N loop order plus immediate error propagation from
%        run_mintime already enforces "2.5 N only after 5 N succeeds"
%        without extra bookkeeping.
%
% After both legs: prints the collection table (T, c_tf, m_f, switches,
% t_f,min, provenance), checks law R0 (T_max*t_f,min spread across the three
% thrusts), and saves a front figure of m_f vs c_tf at 10 N.
%
% CLOSURE (Task 14 controller triage, 2026-07-17): the thrust-law leg is
% BLOCKED at the 5 N min-time anchor after six documented solver attempts
% (see run_mintime.m header + task-14-report.md attempt table), the last of
% which was Route-B (energy warm-start, this cell's documented cure
% elsewhere in the campaign) -- the energy stage itself hit
% Maximum_Iterations_Exceeded (defect 3.0e-2) at N=1200/tf=89 ND. Pattern
% finding: every 5 N attempt sits at <=70 nodes/rev vs >=130 for every 10 N
% success; nodes-per-rev parity (N~2400 for the 5 N ENERGY stage, with its
% own continuation rounds) is the recommended FUTURE attempt, out of scope
% for this task. blockedThrusts (leg 2) stops the sweep from re-attempting
% this known-doomed solve on every future rerun -- the c_tf FRONT (leg 1) is
% unaffected and fully delivered; the law leg reports blocked with the 10 N
% point as the sole available law-R0 value.
%
% INPUTS:  none   OUTPUTS: none (results/sweep_*.mat + front figure + printed table)
% REFERENCES: [1] DESIGN.md sec 5 milestone M3. [2] paper Figs 18/23, law R0.
%   [3] Task 14 controller triage (this task): certified-only caching,
%       best-of chain/fresh envelope for basin scatter, thrust-continuation
%       anchor seeding (run_mintime seedAnchor, run_transfer certified-gate),
%       Route-B, closure with the law leg blocked at 5 N.
here = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');

% --- leg 1: c_tf front at 10 N (neighbor-chain upward from M2) ---------------
% ctfs includes the 4 densify_front.m points (1.35, 1.75, 2.25, 2.75; Task 14
% cleanup pass) alongside the original 5 -- all 9 are already cached, so
% adding them here only widens the collection table/figure below, it does
% not change the main sweep loop's solve behavior (every point skips).
ctfs = [1.2 1.35 1.5 1.75 2.0 2.25 2.5 2.75 3.0];
redoSet = [1.2 2.5 3.0];              % controller-authorized best-of-envelope set
noChainRetry = 1.2;                   % Task 14 M3 closure: the chain-seeded solve at
                                       % c_tf=1.2 deterministically fails to certify (same
                                       % seed, same solve -- reproducibly uncertified across
                                       % every rerun this task, confirmed at least 4x). The
                                       % redo-set FRESH solve below is authoritative for this
                                       % point (mf=1360.37, certified), so stop reattempting
                                       % the doomed chain solve on every future rerun.
prev = fullfile(resDir, 'M2_manifold.mat');
for cf = ctfs
    tag = sprintf('sweep_T100_c%03d', round(100*cf));
    fn  = fullfile(resDir, [tag '.mat']);
    if isfile(fn)
        fprintf('skip %s\n', tag);  prev = fn;  continue;
    end
    if ismember(cf, noChainRetry)
        fprintf('skip %s (chain seed known to deterministically fail to certify -- see header; using fresh-only)\n', tag);
        continue;    % do NOT advance prev: the next point still seeds from the last GOOD one
    end
    run_transfer(struct('thrustN',10, 'ctf',cf, 'hx0',0.0612, 'term','manifold', ...
                        'tag',tag, 'seedMat',prev));
    if isfile(fn), prev = fn; end   % only advance the chain past a CERTIFIED point
end
% redo set: second fresh full-schedule solve for points needing another basin
for cf = redoSet
    tag = sprintf('sweep_T100_c%03d_fresh', round(100*cf));
    fn  = fullfile(resDir, [tag '.mat']);
    if isfile(fn), fprintf('skip %s\n', tag); continue; end
    run_transfer(struct('thrustN',10, 'ctf',cf, 'hx0',0.0612, 'term','manifold', ...
                        'tag',tag, 'N',600));
end

% --- leg 2: thrust law at c_tf = 1.5 (fresh pipelines; N scaled by revs) -----
% Anchors: 10 N is already cached (M0-M2). 5 N and 2.5 N use ROUTE-B (Task 14
% controller triage round 5): a plain cold tangential-seed min-time
% continuation genuinely stalls at 5 N (defect floor ~5e-3, IPOPT declares
% local infeasibility -- a documented failure signature for this cell, cured
% in the parent tulip campaign by warm-starting min-time from a converged
% smooth ENERGY solution instead of a raw thrust-propagation seed; see
% run_mintime.m header and NLP_lowThrust_GTO_tulip/elfo/gen_elfo_mintime.m).
% routeB.tDur is picked at ~2x each thrust's own min-time estimate from the
% C-law (T*tf ~ const, C ~ 846.6 N.h from the 10 N anchor): 5 N -> tDur=89 ND
% (169 h est. x2); 2.5 N -> tDur=176 ND (338 h est. x2), chained from a
% converged 5 N anchor (only attempted once 5 N succeeds -- an error halts
% the whole script before reaching it). N is nodes-per-rev parity (~4.2/8.4/
% 17 revs => 600/1200/2400).
% BLOCKED (Task 14 M3 closure): 5 N failed all six documented solver
% attempts, the last being Route-B (see header). 2.5 N was never reached
% (chained after 5 N). blockedThrusts stops the sweep from re-attempting a
% known-doomed multi-hour solve on every future rerun; remove an entry once
% a working recipe (nodes-per-rev-parity energy stage, see header) is found.
blockedThrusts = [5 2.5];
thr = [10 5 2.5];  Ns = [600 1200 2400];  mintimeNs = [600 1200 2400];
routeBOpts = {[], struct('sbar',0.6,'tDur',89), struct('sbar',0.6,'tDur',176)};
for kt = 1:numel(thr)
    tag = sprintf('sweep_T%03d_c150', round(10*thr(kt)));
    fn  = fullfile(resDir, [tag '.mat']);
    anchorFn = fullfile(resDir, sprintf('mintime_T%d_i7.mat', round(10*thr(kt))));
    if ~isfile(anchorFn)
        if ismember(thr(kt), blockedThrusts)
            fprintf('BLOCKED %g N anchor (known-doomed after 6 solver attempts, see header) -- skip, not reattempting\n', thr(kt));
            continue;
        end
        run_mintime(thr(kt), 0.0612, mintimeNs(kt), '', routeBOpts{kt});
    end
    if isfile(fn), fprintf('skip %s\n', tag); continue; end
    run_transfer(struct('thrustN',thr(kt), 'ctf',1.5, 'hx0',0.0612, ...
                        'term','manifold', 'tag',tag, 'N',Ns(kt)));
end

% --- collect + gates ----------------------------------------------------------
% Candidate tags per c_tf: base chain, cold-fresh redo, and (Task 14 dominated-
% basin recovery pass) alt1/alt2 alternate-neighbor-seeded retries. Same
% max-certified-m_f best-of rule and provenance recording as chain/fresh.
fprintf('\n%-8s %-6s %-9s %-6s %-9s %-6s\n', 'T [N]', 'c_tf', 'mf [kg]', 'sw', 'tfmin [h]', 'src');
mfC = zeros(1, numel(ctfs));
srcTags = {'', '_fresh', '_alt1', '_alt2'};
srcNames = {'chain', 'fresh', 'alt1', 'alt2'};
for kc = 1:numel(ctfs)
    cf = ctfs(kc);
    cand = struct('mf', {}, 'S', {}, 'src', {});
    for kt = 1:numel(srcTags)
        fn = fullfile(resDir, sprintf('sweep_T100_c%03d%s.mat', round(100*cf), srcTags{kt}));
        if isfile(fn)
            S = load(fn);
            cand(end+1) = struct('mf', S.res.report.m_f_kg, 'S', S, 'src', srcNames{kt}); %#ok<AGROW>
        end
    end
    if isempty(cand)
        error('run_ctf_sweep:missingPoint', ...
              'c_tf=%.2f has no certified result (chain/fresh/alt1/alt2) -- rerun the sweep', cf);
    end
    [~, ibest] = max([cand.mf]);
    mfC(kc) = cand(ibest).mf;
    S = cand(ibest).S;
    fprintf('%-8g %-6.2f %-9.2f %-6d %-9.1f %-6s\n', 10, cf, mfC(kc), ...
            S.res.report.switches, S.res.mintime.tfmin_h, cand(ibest).src);
end
fprintf('\n--- thrust law (c_tf=1.5) ---\n');
C = [];
for kt = 1:numel(thr)
    chainFn = fullfile(resDir, sprintf('sweep_T%03d_c150.mat', round(10*thr(kt))));
    freshFn = fullfile(resDir, sprintf('sweep_T%03d_c150_fresh.mat', round(10*thr(kt))));
    cand = struct('mf', {}, 'S', {}, 'src', {});
    if isfile(chainFn)
        Sc = load(chainFn);
        cand(end+1) = struct('mf', Sc.res.report.m_f_kg, 'S', Sc, 'src', 'chain'); %#ok<AGROW>
    end
    if isfile(freshFn)
        Sf = load(freshFn);
        cand(end+1) = struct('mf', Sf.res.report.m_f_kg, 'S', Sf, 'src', 'fresh'); %#ok<AGROW>
    end
    if isempty(cand)
        fprintf('%-8g %-6.2f BLOCKED (min-time anchor did not converge -- see task-14-report.md)\n', ...
                thr(kt), 1.5);
        continue;
    end
    [~, ibest] = max([cand.mf]);
    S = cand(ibest).S;
    C(end+1) = thr(kt) * S.res.mintime.tfmin_h; %#ok<AGROW>
    fprintf('%-8g %-6.2f %-9.2f %-6d %-9.1f %-6s\n', thr(kt), 1.5, ...
            S.res.report.m_f_kg, S.res.report.switches, S.res.mintime.tfmin_h, cand(ibest).src);
end
if numel(C) == numel(thr)
    fprintf('law R0: T*tfmin = %s N.h  (spread %.1f%%)\n', mat2str(round(C)), ...
            100*(max(C)-min(C))/mean(C));
else
    fprintf(['law R0: only %d/%d thrust point(s) available -- T*tfmin = %s N.h ' ...
             '(spread NOT evaluable; thrust-law leg BLOCKED, see task-14-report.md)\n'], ...
            numel(C), numel(thr), mat2str(round(C)));
end
fig = figure('Visible','off');
try, theme(fig,'light'); catch, end
set(fig, 'Color', 'w');
plot(ctfs, mfC, 'o-', 'Color',[0.10 0.35 0.85], 'MarkerFaceColor',[0.10 0.35 0.85], ...
     'MarkerEdgeColor','k', 'LineWidth', 1.8, 'MarkerSize', 7);
grid on
xlabel('c_{tf}');  ylabel('m_f [kg]');
title('GEO transfer: final mass vs transfer-time multiplier (T_{max}=10 N)');
exportgraphics(fig, fullfile(resDir, 'front_mf_ctf.png'), 'Resolution', 150);
close(fig);
end
