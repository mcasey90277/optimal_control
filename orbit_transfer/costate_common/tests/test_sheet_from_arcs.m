function ok = test_sheet_from_arcs()
%% Purpose:
%
%   Tests sheet_from_arcs -- the assembler that turns arclength_ms
%   crossings (from any number of arcs) into a certified arrival-phase
%   sheet: unwrapped level -> grid index (mod 1), duplicate candidates
%   from different arcs merged, EVERY candidate certified and kept with its
%   verdict, S.TF = min t_f over the CERTIFIED candidates only, failure
%   reasons stored. Uses a stub certifier so the logic is tested alone.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
sA0 = 0.0754;  nA = 12;

% p = [tfDays; tag]; the stub certifies when tag > 0, t_f = p(1)
mk = @(level, tfD, tag) struct('level', level, 'q', level, 'p', [tfD; tag], ...
                               'converged', true, 'normR', 1e-12, 'afterIndex', 1);
arc1.crossings = [mk(sA0 + 1/12, 17.0, 1), mk(sA0 + 2/12, 16.5, 1), mk(sA0 + 2/12, 15.9, -1)];
arc2.crossings = [mk(sA0 + 1/12 - 1, 17.0, 1), ...    % same root as arc1's first, unwrapped by -1
                  mk(sA0 + 3/12, 20.0, 1), mk(sA0 + 3/12, 19.0, 1)];
arc3.crossings = mk(sA0 + 4/12, 0, 0);  arc3.crossings(1).converged = false;

certFn = @(p, sA) stubCert(p, sA);
S = sheet_from_arcs({arc1, arc2, arc3}, struct('sA0', sA0, 'nA', nA, 'certFn', certFn));

ok = chk(ok, numel(S.sA) == nA && abs(S.sA(2) - (sA0 + 1/12)) < 1e-12, 'grid built from sA0 and nA');
ok = chk(ok, numel(S.cand{2}) == 1, sprintf('duplicate root across arcs merged: %d candidate at j = 2', numel(S.cand{2})));
ok = chk(ok, numel(S.cand{3}) == 2 && abs(S.TF(3) - 16.5) < 1e-12, ...
         sprintf('refused candidate kept but not chosen: TF(3) = %.2f from %d candidates', S.TF(3), numel(S.cand{3})));
ok = chk(ok, abs(S.TF(4) - 19.0) < 1e-12, sprintf('minimum over certified candidates: TF(4) = %.1f', S.TF(4)));
ok = chk(ok, isnan(S.TF(5)) && numel(S.cand{5}) == 1 && ~S.cand{5}(1).ok && contains(S.cand{5}(1).reason, 'not converged'), ...
         'unconverged crossing recorded with a reason, TF = NaN');
ok = chk(ok, all(isnan(S.TF([1 6:12]))), 'grid points with no crossing stay NaN');
ok = chk(ok, S.nCand == 6 && S.nCert == 4, sprintf('counts: %d candidates, %d certified', S.nCand, S.nCert));

if ok, fprintf('TEST_SHEET_FROM_ARCS: ALL PASS\n'); else, fprintf('TEST_SHEET_FROM_ARCS: FAIL\n'); end
end

function C = stubCert(p, sA)
% STUBCERT  Stand-in certifier: tag > 0 certifies.  INPUTS: p; sA.
% OUTPUTS: C.
C = struct('ok', p(2) > 0, 'reason', 'certified', 'z', [zeros(7,1); p(1)], ...
           'tfDays', p(1), 'sA', sA);
if ~C.ok, C.reason = 'stub refused'; end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
