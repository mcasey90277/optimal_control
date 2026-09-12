function ok = test_crossings_from_arc()
%% Purpose:
%
%   Refining the arrival grid rests on one claim: a saved arc already
%   contains every crossing there is, so new levels can be extracted from
%   its stored path instead of walking it again (hours each). This tests
%   that claim rather than trusting it.
%
%     1. re-scanning at the levels the walk was GIVEN reproduces every
%        crossing it recorded, to machine zero in arrival phase;
%     2. the only extra is the arc's own START point -- a walk does not
%        cross the level it begins on -- which is why the sheet seeds
%        itself from the library at those phases;
%     3. an interpolated seed POLISHES TO THE SAME ROOT as the recorded one.
%        This is the claim that matters: adjacent continuation steps are
%        ~2e-4 apart in phase, and the test measures that this is inside the
%        Newton basin rather than assuming it. Note it compares the ROOT,
%        not the verdict: most crossings legitimately fail a gate (the sheet
%        keeps every candidate with its reason), and two seeds that fail the
%        SAME gate at the SAME value have demonstrably converged to the same
%        place. Crossings are matched by .afterIndex, not by phase: a folded
%        arc crosses one level several times and matching on phase alone
%        picks an arbitrary one of them;
%     4. folds are kept: a level the arc crosses more than once yields more
%        than one candidate;
%     5. doubling the levels roughly doubles the crossings and keeps every
%        coarse one.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/12/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
ind = fullfile(fileparts(here), 'DRO_tulip', 'indirect');
addpath(here, ind);
arcFile = fullfile(ind, 'results', 'arrival_arc_cell11_dn_long.mat');
if ~isfile(arcFile), fprintf('  SKIP  no saved arc (%s)\n', arcFile);  return, end
L = load(arcFile);  A = L.A;
sA0 = 0.0754;  coarse = sA0 + (0:11)/12;  fine = sA0 + (0:23)/24;

% ---- 1, 2. the re-scan reproduces the walk, plus the start point --------
X = crossings_from_arc(A, coarse);
rec = sort([A.crossings.q]);  got = sort([X.q]);
dmin = arrayfun(@(r) min(abs(got - r)), rec);
ok = chk(ok, max(dmin) < 1e-12, ...
         sprintf('every one of the %d recorded crossings is recovered (worst %.1e in phase)', ...
                 numel(rec), max(dmin)));
extra = setdiff(round(got, 10), round(rec, 10));
ok = chk(ok, numel(extra) <= 1, sprintf('at most one extra crossing (%d)', numel(extra)));
if ~isempty(extra)
    ok = chk(ok, abs(extra(1) - A.q(1)) < 1e-9, ...
             sprintf('and it is the arc''s own start point (%.4f vs q(1) = %.4f)', extra(1), A.q(1)));
end

% ---- 4. folds keep repeat crossings -------------------------------------
lv = round([X.level], 6);
[u, ~, g] = unique(lv);
cnt = accumarray(g, 1);
ok = chk(ok, any(cnt > 1), ...
         sprintf('a folded arc gives a level more than one candidate (max %d at phase %.4f)', ...
                 max(cnt), u(find(cnt == max(cnt), 1))));

% ---- 5. doubling the grid ------------------------------------------------
XF = crossings_from_arc(A, fine);
qc = round(sort([X.q]), 9);  qf = round(sort([XF.q]), 9);
ok = chk(ok, all(ismember(qc, qf)), 'every coarse crossing survives at double resolution');
ok = chk(ok, numel(XF) >= 1.5*numel(X), ...
         sprintf('and the count roughly doubles: %d -> %d', numel(X), numel(XF)));

% ---- 3. THE claim: an interpolated seed certifies to the same root -------
[B, anc] = arclength_arrival('setup');
pool = gcp('nocreate');
copts = struct('pool', pool, 'wallSec', 600);
nTry = min(3, numel(A.crossings));
worst = 0;  nPaired = 0;  nVerdict = 0;
for k = 1:nTry
    cr = A.crossings(k);
    ix = find([X.afterIndex] == cr.afterIndex, 1);   % THE same crossing, fold or not
    if isempty(ix)
        fprintf('      crossing %d: no re-scanned crossing at step %d\n', k, cr.afterIndex);
        continue
    end
    Crec = certify_crossing(cr.p,    cr.level, B, anc, copts);
    Cint = certify_crossing(X(ix).p, cr.level, B, anc, copts);
    if all(isfinite(Crec.z)) && all(isfinite(Cint.z))
        dz = norm(Crec.z(:) - Cint.z(:))/max(norm(Crec.z), realmin);
        worst = max(worst, dz);  nPaired = nPaired + 1;
        nVerdict = nVerdict + (Crec.ok == Cint.ok);
        fprintf('      crossing %d at phase %.4f (step %d): t_f %.6f vs %.6f d, rel |dz| %.2e, ok %d/%d\n', ...
                k, cr.level, cr.afterIndex, Crec.tfDays, Cint.tfDays, dz, Crec.ok, Cint.ok);
    end
end
ok = chk(ok, nPaired == nTry, sprintf('%d of %d sampled crossings polished from BOTH seeds', nPaired, nTry));
ok = chk(ok, nPaired == nTry && worst < 1e-6, ...
         sprintf('and landed on the SAME root (worst relative |dz| %.2e)', worst));
ok = chk(ok, nVerdict == nTry, ...
         sprintf('with the same verdict from either seed (%d of %d)', nVerdict, nTry));

if ok, fprintf('TEST_CROSSINGS_FROM_ARC: ALL PASS\n'); else, fprintf('TEST_CROSSINGS_FROM_ARC: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
