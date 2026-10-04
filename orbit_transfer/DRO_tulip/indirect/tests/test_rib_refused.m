function ok = test_rib_refused()
% TEST_RIB_REFUSED  A rib keeps every refused candidate above the floor,
% including the point it stopped at; with walkPastConjugate a status-2
% point advances the walk and is recorded, never as a certified point; off
% (the default) behaves as before.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
% nD = 48: the 1/12 and 1/24 steps out of the anchor fail the polish (status 0, nothing
% to keep -- see test_rib_from_crossing); the 1/48 step solves, so a forced
% coarse FAIL is what refuses it (status 2).
ot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(fullfile(ot, 'costate_common'), fullfile(ot, 'DRO_tulip', 'indirect'));
[B, anc] = arclength_arrival('setup');
C0 = certify_crossing(anc.p, anc.sA, B, anc);
assert(C0.ok, 'fixture: the anchor must certify (%s)', C0.reason);
forced = struct('pool', capped_pool(), 'override', struct('conj', struct('pass', false, 'verdict', 'FAIL')));

R = rib_from_crossing(C0, B, anc, struct('nD', 48, 'direction', -1, 'nPts', 1, 'maxBisect', 1, 'copts', forced));
ok = chk(ok, isempty(R.pts) && startsWith(R.stop, 'stalled') && numel(R.refused) == 2, ...
         sprintf('default: stops, keeps both refused trials (%d)', numel(R.refused)));
ok = chk(ok, all([R.refused.status] == 2) && ~any([R.refused.walked]) && all([R.refused.overridden]), ...
         'refusals carry status 2, not walked, flagged overridden');

R = rib_from_crossing(C0, B, anc, struct('nD', 48, 'direction', -1, 'nPts', 1, 'maxBisect', 1, 'copts', forced, ...
                                         'walkPastConjugate', true));
ok = chk(ok, isempty(R.pts) && strcmp(R.stop, 'complete') && numel(R.refused) == 1 && R.refused(1).walked ...
         && abs(R.refused(1).sD - 47/48) < 1e-9, 'walkPastConjugate: advances through the status-2 point and records it');
if ok, fprintf('test_rib_refused: ALL PASS\n'); else, fprintf('test_rib_refused: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
