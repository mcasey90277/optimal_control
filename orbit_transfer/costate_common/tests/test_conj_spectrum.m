function ok = test_conj_spectrum()
%% Purpose:
%
%   Tests conj_spectrum -- the DENSE singular-spectrum scan that closes the
%   two blind spots of the sampled sign test: two conjugate times inside one
%   segment, and an even-multiplicity crossing that produces no sign change.
%
%   POSITIVE AND NEGATIVE CONTROL in one test. It runs on a certified entry
%   (the conjugate test says PASS) and on an entry the conjugate test
%   REFUTED, and requires the scan to agree with both verdicts. A detector
%   that fires on everything, or on nothing, fails here.
%
%   It also pins the lesson from the 2026-09-10 adjudication: the smallest
%   singular value AT t_f does not discriminate -- a refuted entry and a
%   certified one share it to two percent -- so the scan must report interior
%   structure, not endpoint smallness.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));

[B, anc] = arclength_arrival('setup');
C = certify_crossing(anc.p, anc.sA, B, anc);
assert(C.ok, 'fixture: the anchor must certify (%s)', C.reason);

% ---- the CERTIFIED entry: no interior collapse -------------------------
S = conj_spectrum(C.z, B.rv0(1:6), B.Tnd, B.cnd, B.mu, struct('K', 24, 'nSub', 8));
ok = chk(ok, numel(S.t) > 150, sprintf('dense scan: %d samples over %d segments', numel(S.t), 24));
ok = chk(ok, isequal(size(S.sv), [6 numel(S.t)]), 'the full 6-value spectrum is returned, not a determinant');
ok = chk(ok, S.nInterior == 0, sprintf('certified entry: %d interior crossing(s)', S.nInterior));
ok = chk(ok, S.multiplicity == 0, sprintf('and no multiplicity event (%d)', S.multiplicity));
% Candidates are LOCATED and CLASSIFIED, not counted. A dip in the first
% samples is the start-up transient (Phi_rv -> 0 at t = 0); a dip in the
% last samples is the graded endpoint collapse; only an INTERIOR dip is a
% candidate conjugate point, and for those the scan says whether the sampled
% minimum FALLS under 4x refinement (a zero) or PLATEAUS (a near-miss).
ok = chk(ok, isfield(S, 'candidates') && isstruct(S.candidates), 'candidates are returned as a struct array');
cl = {S.candidates.class};
ok = chk(ok, all(ismember(cl, {'start', 'endpoint', 'interior'})), ...
         sprintf('every candidate is classified: %s', strjoin(cl, ' ')));
ok = chk(ok, ~any(strcmp(cl, 'interior')), 'the certified anchor has NO interior candidate');
ok = chk(ok, all(cellfun(@isscalar, {S.candidates.tOverTf})), 'each candidate carries its t/t_f');
ok = chk(ok, isfield(S, 'nInteriorCand') && S.nInteriorCand == 0 && isfield(S, 'nNearMiss') && isfield(S, 'nZero'), ...
         'interior candidates are summarised as near-miss / zero counts');

% ---- a REFUTED entry: the scan must find what the test found ------------
Sh = load(fullfile(fileparts(here), 'DRO_tulip', 'indirect', 'results', ...
                   'arrival_sheet_70mN_pass1.mat'));
bad = [];
for j = 1:numel(Sh.S.sA)
    c = Sh.S.cand{j};
    for k = 1:numel(c)
        if ~c(k).ok && contains(c(k).reason, 'conjugate') && isfinite(c(k).tfDays)
            bad = c(k);  badSA = Sh.S.sA(j);  break
        end
    end
    if ~isempty(bad), break, end
end
if isempty(bad)
    fprintf('  SKIP  no refuted candidate on disk to use as a positive control\n');
else
    rv0b = B.rv0(1:6);
    Sb = conj_spectrum(bad.z, rv0b, B.Tnd, B.cnd, B.mu, struct('K', 24, 'nSub', 8));
    ok = chk(ok, Sb.nInterior >= 1, ...
             sprintf('refuted entry (%.2f d, sA %.4f): %d interior crossing(s) found', ...
                     bad.tfDays, badSA, Sb.nInterior));
    ok = chk(ok, Sb.tFirst > 0 && Sb.tFirst < Sb.tf, ...
             sprintf('and it is INTERIOR, at t/t_f = %.4f', Sb.tFirst/Sb.tf));
    % the crossing must also appear as an INTERIOR candidate whose sampled
    % minimum FALLS under refinement -- a zero, not a near-miss
    ic = Sb.candidates(strcmp({Sb.candidates.class}, 'interior'));
    ok = chk(ok, ~isempty(ic), sprintf('%d interior candidate(s) located', numel(ic)));
    if ~isempty(ic)
        ok = chk(ok, any(strcmp({ic.kind}, 'zero')), ...
                 sprintf('refined: kinds = %s (min fell by %s, then %s)', strjoin({ic.kind}, ' '), ...
                         strjoin(compose('%.2g', [ic.refineRatio]), ' '), ...
                         strjoin(compose('%.2g', [ic.refineRatio2]), ' ')));
        ok = chk(ok, all(isfinite([ic.refineRatio2])), 'a SECOND refinement level is reported');
    end
end

% ---- the endpoint value must NOT be used as the verdict -----------------
ok = chk(ok, isfield(S, 'svEnd') && ~isfield(S, 'passFromEndpoint'), ...
         'the endpoint spectrum is reported but is not the verdict');

if ok, fprintf('TEST_CONJ_SPECTRUM: ALL PASS\n'); else, fprintf('TEST_CONJ_SPECTRUM: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
