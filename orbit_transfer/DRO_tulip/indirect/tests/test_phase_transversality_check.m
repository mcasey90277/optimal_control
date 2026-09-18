function ok = test_phase_transversality_check()
% TEST_PHASE_TRANSVERSALITY_CHECK  The catalog-level cross-check X3 on the
% library of record: the phase sensitivities exist at every entry, ONE entry
% is put through the exact test (re-solve at s +/- delta, which never looks
% at a costate), and the departure-phase edge map is self-consistent.
% About a minute: four polishes.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));            % DRO_tulip/indirect
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
catMat = fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat');
X = phase_transversality_check(catMat, struct('print', false, 'exactCells', [6 17]));

ok = chk(ok, isequal(size(X.dTf_dsD), [24 24]) && isequal(size(X.dTf_dsA), [24 24]) && ...
             nnz(isfinite(X.dTf_dsD)) == 576 && nnz(isfinite(X.dTf_dsA)) == 576, 'a sensitivity in each phase at every one of the 576 entries');
e = X.exact(1);
ok = chk(ok, e.iD == 6 && e.iA == 17 && e.relErrA < 1e-3 && e.relErrD < 1e-3, ...
         sprintf('exact test at (6,17): formula vs re-solve, relative error %.1e (arrival) %.1e (departure)', e.relErrA, e.relErrD));
ok = chk(ok, X.okExact, 'the exact gate passes');
ok = chk(ok, isequal(size(X.edgeD), [24 24]) && X.nEdgeD == nnz(isfinite(X.edgeD)) && X.nEdgeDJump <= X.nEdgeD, 'the departure-phase edge map is counted consistently');
ok = chk(ok, isfield(X, 'stationary') && height(X.stationary) >= 1 && all(diff(X.stationary.gradNorm) >= 0), ...
         'the cells nearest first-order stationarity in BOTH phases are listed, smallest gradient first');
% ---- the exact test's verdict is THREE-valued --------------------------------
% A re-solve that does not converge says nothing about a costate: it is
% UNRESOLVED, not FAIL. A resolved disagreement is FAIL whatever else happened.
H = phase_transversality_check('localfunctions');
mk = @(eD, eA) struct('relErrD', eD, 'relErrA', eA);
ok = chk(ok, strcmp(H.exactVerdict([mk(1e-6, 1e-5), mk(1e-7, 2e-4)], 1e-3), 'PASS'), 'every component resolved and within tolerance -> PASS');
ok = chk(ok, strcmp(H.exactVerdict([mk(1e-6, NaN), mk(1e-7, 2e-4)], 1e-3), 'UNRESOLVED'), 'one re-solve did not converge, none disagrees -> UNRESOLVED');
ok = chk(ok, strcmp(H.exactVerdict([mk(1e-6, NaN), mk(5e-3, 2e-4)], 1e-3), 'FAIL'), 'a resolved disagreement -> FAIL, even beside an unresolved one');
ok = chk(ok, strcmp(H.exactVerdict(mk([], []), 1e-3), 'UNRESOLVED') || strcmp(H.exactVerdict(struct('relErrD', {}, 'relErrA', {}), 1e-3), 'UNRESOLVED'), 'nothing sampled -> UNRESOLVED, never a pass');
ok = chk(ok, isfield(X, 'verdictExact') && strcmp(X.verdictExact, 'PASS') && X.okExact, 'the run above reports its verdict by name');

if ok, fprintf('test_phase_transversality_check: ALL PASS\n'); else, fprintf('test_phase_transversality_check: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
