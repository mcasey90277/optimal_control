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
