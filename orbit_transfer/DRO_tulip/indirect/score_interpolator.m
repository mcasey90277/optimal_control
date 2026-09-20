function R = score_interpolator(cat_, opts)
%% Purpose:
%
%   THE HIT RATE OF THE COSTATE LIBRARY AS A SOURCE OF GUESSES. Draw phase
%   pairs at random, off the grid; for each, blend the entries around it into
%   a guess (catalog_blend_guess), polish it under a hard cap
%   (polish_costate_guess), and ask: did it converge, and onto the branch of
%   the corners it was blended from? The same solve from the NEAREST entry
%   alone is the baseline -- what a user without an interpolator would do.
%
%   One number answers "is this library fine enough": the fraction of random
%   queries that are USABLE (converged AND on its corners' branch). It is the
%   yardstick for every resolution decision (24 x 24, 24 x 48, ...).
%
%  ASSUMPTIONS / NOTES:
%
% • Campaign discipline: rows are saved after EVERY query, a rerun resumes
%   (same seed -> same queries; asking for more keeps the ones done), every
%   solve is fenced, and the verdict is a FILE. Launch long runs detached.
% • A solve the fence stopped, or one that did not converge, COUNTS AGAINST
%   the hit rate. A query whose cell holds nothing (tier 'none') does too:
%   the library had no guess to offer there.
% • "Same branch": the root's costates within opts.tolBranch (relative) of
%   the guess, and its t_f within opts.tolBranchDays of the used corners'
%   range. A converged solve that fails this found A transfer, not the
%   library's.
% • No certification here (a minute per query); interp_study certifies one.
%   opts.solver is the TEST SEAM: a handle (zGuess, rv0, rvf, phys) ->
%   [finished, z8, info] replacing the fenced polish.
%
%% Inputs:
%
%  cat_                     char | struct           a catalog file, a catalog
%                                                   struct, or a
%                                                   catalog_blend_setup
%                                                   output (has .GD)
%
%  opts                     struct (optional)       .nQuery [60] .seed [1]
%                                                   .outDir [<catalog's
%                                                   folder>/interp_score]
%                                                   .doBaseline [true]
%                                                   .capSec [420] per solve
%                                                   .pool [capped_pool(2)]
%                                                   .allowUnfenced [false]
%                                                   .tolBranch [0.5]
%                                                   .tolBranchDays [1]
%                                                   .sDfixed [] (set for a
%                                                   one-row library; taken
%                                                   from it automatically)
%                                                   .print [true] .solver []
%
%% Outputs:
%
%  R                        struct                  .rows (per query: .sD .sA
%                                                   .tier .reason .outcome
%                                                   .converged .sameBranch
%                                                   .usable .guessMissKm
%                                                   .iters .sec .q2 (guess vs
%                                                   root, relative) .q2Min
%                                                   .tfDays .base (the
%                                                   baseline's .converged
%                                                   .sameRoot .guessMissKm
%                                                   .iters .sec)),
%                                                   .summary (.nQuery
%                                                   .nUsable .hitRate
%                                                   .byTier(.tier .n .usable)
%                                                   .baseHitRate
%                                                   .medianIters
%                                                   .baseMedianIters ...),
%                                                   .outDir
%
%% Revision History:
%  M. Casey                                                   (c) 09/20/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
d = @(f, v) pick(opts, f, v);
if isstruct(cat_) && isfield(cat_, 'GD'), K = cat_; else, K = catalog_blend_setup(cat_); end
nQuery = d('nQuery', 60);  seed = d('seed', 1);  doBaseline = d('doBaseline', true);  capSec = d('capSec', 420);
tolBranch = d('tolBranch', 0.5);  tolBranchDays = d('tolBranchDays', 1);  doPrint = d('print', true);
defaultOut = fullfile(tempdir, 'interp_score');  if ~isempty(K.catMat), defaultOut = fullfile(fileparts(K.catMat), 'interp_score'); end
outDir = d('outDir', defaultOut);  if ~isfolder(outDir), mkdir(outDir); end
rowsF = fullfile(outDir, 'score_rows.mat');  logF = fullfile(outDir, 'score.log');  verdictF = fullfile(outDir, 'SCORE_VERDICT.txt');
day = @(t) t*K.phys.tStar/86400;

% ---- the solver: the fenced polish, or the test's stand-in -----------------
solver = d('solver', []);
if isempty(solver)
    pool = d('pool', []);  if ~isfield(opts, 'pool'), pool = capped_pool(2); end
    assert(~isempty(pool) || d('allowUnfenced', false), 'score_interpolator:noFence', ...
           'no worker pool: the solves would run unfenced. Pass allowUnfenced = true to accept that');
    solver = @(zG, rv0, rvf, phys) fencedPolish(pool, capSec, zG, rv0, rvf, phys);
end

% ---- the queries: seeded, so a rerun asks the same ones ----------------------
stream = RandStream('mt19937ar', 'Seed', seed);
Q = rand(stream, nQuery, 2);
oneRowD = isscalar(K.sheet.sD_frac);  oneRowA = isscalar(K.sheet.sA_frac);
if oneRowD, Q(:, 1) = K.sheet.sD_frac; end
if oneRowA, Q(:, 2) = K.sheet.sA_frac; end

rows = emptyRows();
if isfile(rowsF)
    L = load(rowsF);
    if isfield(L, 'meta') && L.meta.seed == seed && strcmp(L.meta.catMat, K.catMat), rows = L.rows; end
end
meta = struct('seed', seed, 'catMat', K.catMat, 'when', datestr(now)); %#ok<NASGU>

for kq = numel(rows) + 1:nQuery
    sD = Q(kq, 1);  sA = Q(kq, 2);
    r = newRow(sD, sA);
    [zG, blend] = catalog_blend_guess(K, sD, sA);
    r.tier = blend.tier;  r.reason = blend.reason;
    if strcmp(blend.tier, 'none')
        r.outcome = 'no guess';
    else
        x0 = K.stateD(sD);  xf = K.stateA(sA);  rv0 = x0(1:6);  rvf = xf(1:6);
        t0 = tic;  [finished, z, info] = solver(zG, rv0, rvf, K.phys);  r.sec = toc(t0);
        [r, zRoot] = judge(r, finished, z, info, zG, blend, K, tolBranch, tolBranchDays, day);
        if doBaseline
            [zN, ~] = catalog_blend_guess(K, sD, sA, struct('nearestOnly', true));
            t0 = tic;  [finB, zB, infoB] = solver(zN, rv0, rvf, K.phys);
            r.base.sec = toc(t0);  r.base.ran = true;
            if finB && isstruct(infoB)
                r.base.converged = logical(infoB.converged) && all(isfinite(zB));
                r.base.iters = infoB.iters;  r.base.guessMissKm = infoB.guessMissKm;
            end
            if r.base.converged && r.converged
                r.base.sameRoot = sqrt(sum((zB(1:7) - zRoot(1:7)).^2))/sqrt(sum(zRoot(1:7).^2)) < 1e-6;
            end
        end
    end
    rows(end+1) = r; %#ok<AGROW>
    save(rowsF, 'rows', 'meta');
    logLine(logF, sprintf('%3d (%.4f, %.4f) %-8s -> %-22s miss %8.0f km, %4g it, %4.0f s | baseline conv %d, %4g it', kq, sD, sA, r.tier, ...
                          r.outcome, r.guessMissKm, r.iters, r.sec, r.base.converged, r.base.iters));
end
rows = rows(1:min(numel(rows), nQuery));

% ---- the summary -------------------------------------------------------------
S = struct('nQuery', numel(rows), 'nUsable', nnz([rows.usable]), 'hitRate', NaN, 'nConverged', nnz([rows.converged]), ...
           'nOtherBranch', nnz([rows.converged] & ~[rows.sameBranch]), 'byTier', [], 'baseHitRate', NaN, 'nBaseConverged', NaN, ...
           'medianIters', median([rows([rows.usable]).iters]), 'baseMedianIters', NaN, 'medianMissKm', median([rows.guessMissKm], 'omitnan'), ...
           'baseMedianMissKm', NaN, 'medianQ2', median([rows.q2], 'omitnan'));
S.hitRate = S.nUsable/max(S.nQuery, 1);
tiers = {'entry', 'bilinear', 'linear', 'nearest', 'none'};
S.byTier = struct('tier', tiers, 'n', 0, 'usable', 0);
for kt = 1:numel(tiers)
    in = strcmp({rows.tier}, tiers{kt});
    S.byTier(kt).n = nnz(in);  S.byTier(kt).usable = nnz(in & [rows.usable]);
end
if doBaseline && ~isempty(rows)
    base = [rows.base];
    S.nBaseConverged = nnz([base.converged]);  S.baseHitRate = S.nBaseConverged/S.nQuery;
    S.baseMedianIters = median([base([base.converged]).iters]);  S.baseMedianMissKm = median([base.guessMissKm], 'omitnan');
end
R = struct('rows', rows, 'summary', S, 'outDir', outDir);

txt = sprintf(['INTERPOLATOR SCORE  %s\n  %d random off-grid queries (seed %d); grid steps 1/%.0f (departure), 1/%.0f (arrival)\n' ...
               '  USABLE (converged on its corners'' branch): %d of %d = %.0f%%\n  converged %d, of which on ANOTHER branch %d\n'], ...
              K.catMat, S.nQuery, seed, 1/max(K.gapD, eps), 1/max(K.gapA, eps), S.nUsable, S.nQuery, 100*S.hitRate, S.nConverged, S.nOtherBranch);
for kt = 1:numel(tiers)
    if S.byTier(kt).n > 0, txt = [txt sprintf('    tier %-8s: %3d queries, %3d usable\n', tiers{kt}, S.byTier(kt).n, S.byTier(kt).usable)]; end %#ok<AGROW>
end
txt = [txt sprintf('  blend   : median unsolved miss %.0f km, median %g iterations when usable, guess off by %.2f%% (median)\n', ...
                   S.medianMissKm, S.medianIters, 100*S.medianQ2)];
if doBaseline
    txt = [txt sprintf('  baseline (nearest entry alone): converged %d of %d = %.0f%%, median unsolved miss %.0f km, median %g iterations\n', ...
                       S.nBaseConverged, S.nQuery, 100*S.baseHitRate, S.baseMedianMissKm, S.baseMedianIters)];
end
fid = fopen(verdictF, 'w');  fprintf(fid, '%sDONE\n', txt);  fclose(fid);
if doPrint, fprintf('%s', txt); end
end

% ==========================================================================
function [finished, z, info] = fencedPolish(pool, capSec, zG, rv0, rvf, phys)
% FENCEDPOLISH  polish_costate_guess under run_capped (direct if no pool).
% INPUTS: pool; capSec; zG [8x1]; rv0, rvf [6x1]; phys.  OUTPUTS: finished; z; info.
if isempty(pool)
    [z, info] = polish_costate_guess(zG, rv0, rvf, phys);  finished = true;
else
    [finished, z, info] = run_capped(pool, @polish_costate_guess, 2, capSec, zG, rv0, rvf, phys);
end
end

function [r, zRoot] = judge(r, finished, z, info, zG, blend, K, tolBranch, tolBranchDays, day)
% JUDGE  One solve's outcome, in the row's words.
% INPUTS: r row; finished; z; info (solver outputs); zG guess; blend; K setup;
% tolBranch; tolBranchDays; day (ND -> days).  OUTPUTS: r; zRoot [8x1] | NaN.
zRoot = NaN(8, 1);
if ~finished || ~isstruct(info), r.outcome = 'stopped by the cap';  return, end
r.guessMissKm = info.guessMissKm;  r.iters = info.iters;
if ~(logical(info.converged) && isnumeric(z) && numel(z) == 8 && all(isfinite(z))), r.outcome = 'did not converge';  return, end
zRoot = z(:);  r.converged = true;  r.tfDays = day(zRoot(8));
r.q2 = sqrt(sum((zRoot(1:7) - zG(1:7)).^2))/sqrt(sum(zRoot(1:7).^2));
r.q2Min = (zG(8) - zRoot(8))*K.phys.tStar/60;
usedTf = arrayfun(@(a, b) K.sheet.tf_nd(a, b), blend.iD, blend.iA);
outside = max([0, day(min(usedTf) - zRoot(8)), day(zRoot(8) - max(usedTf))]);
r.sameBranch = r.q2 <= tolBranch && outside <= tolBranchDays;
r.usable = r.sameBranch;
if r.usable, r.outcome = 'usable'; else, r.outcome = 'converged, other branch'; end
end

function r = newRow(sD, sA)
% NEWROW  A row with every field born (rows are appended into one array).
% INPUTS: sD, sA.  OUTPUTS: r.
r = struct('sD', sD, 'sA', sA, 'tier', '', 'reason', '', 'outcome', '', 'converged', false, 'sameBranch', false, 'usable', false, ...
           'guessMissKm', NaN, 'iters', NaN, 'sec', NaN, 'q2', NaN, 'q2Min', NaN, 'tfDays', NaN, ...
           'base', struct('ran', false, 'converged', false, 'sameRoot', false, 'guessMissKm', NaN, 'iters', NaN, 'sec', NaN));
end

function rows = emptyRows()
% EMPTYROWS  A 0 x 0 array with newRow's fields.  INPUTS: none.  OUTPUTS: rows.
rows = newRow(0, 0);  rows = rows([]);
end

function logLine(f, s)
% LOGLINE  Append one timestamped line to the log FILE.  INPUTS: f; s.  OUTPUTS: none.
fid = fopen(f, 'a');  fprintf(fid, '%s %s\n', datestr(now, 'HH:MM:SS'), s);  fclose(fid);
end

function v = pick(s, f, d_)
% PICK  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
