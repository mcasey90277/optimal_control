function ok = test_rib_from_crossing()
%% Purpose:
%
%   Tests rib_from_crossing -- the DEPARTURE-phase walk that hangs a rib off
%   one certified arrival-phase crossing. Departure steps are cheap (the
%   departure state moves the first junction only), so the rib is a fixed
%   walker with bisection, not an arclength arc; every point must pass the
%   full gate stack (certify_root) before it is kept.
%
%   Reference: the certified departure ring measured 2026-09-09 at the
%   anchor's arrival phase, walked in the NEGATIVE departure direction --
%   t_f = 17.877 d at sD = 11/12 and 18.069 d at 10/12. (The positive
%   direction failed out of the anchor in that sweep, and does here too:
%   the rib is not symmetric in sD.)
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

[B, anc] = arclength_arrival('setup');
C0 = certify_crossing(anc.p, anc.sA, B, anc);
assert(C0.ok, 'fixture: the anchor must certify first (%s)', C0.reason);

R = rib_from_crossing(C0, B, anc, struct('nD', 12, 'direction', -1, 'nPts', 2, 'maxBisect', 8));

ok = chk(ok, numel(R.pts) == 2, sprintf('two departure points returned (%d)', numel(R.pts)));
ok = chk(ok, all([R.pts.ok]), sprintf('both certified (reasons: %s)', strjoin({R.pts.reason}, '; ')));
ok = chk(ok, abs(R.pts(1).tfDays - 17.877) < 5e-3, sprintf('sD = 11/12: t_f = %.4f d (ring 17.877)', R.pts(1).tfDays));
ok = chk(ok, abs(R.pts(2).tfDays - 18.069) < 5e-3, sprintf('sD = 10/12: t_f = %.4f d (ring 18.069)', R.pts(2).tfDays));
ok = chk(ok, abs(R.pts(1).sD - 11/12) < 1e-9 && abs(R.pts(2).sD - 10/12) < 1e-9, 'departure phases land on the grid');
ok = chk(ok, all([R.pts.sA] == anc.sA), 'arrival phase held fixed along the rib');

% EXPLICIT TARGETS. Deriving a step count by rounding 1/|delta| does not
% reproduce an arbitrary requested phase -- asking for 5/12 gives round(2.4)
% = 2 and walks a HALF period instead. A caller that needs an exact phase
% passes it. (Astra chain review 2026-09-10.)
Rt = rib_from_crossing(C0, B, anc, struct('targets', -1/24, 'maxBisect', 8));
ok = chk(ok, numel(Rt.pts) == 1 && abs(Rt.pts(1).sD - (1 - 1/24)) < 1e-12, ...
         sprintf('walks to an explicit off-grid target: sD = %.6f (want %.6f)', ...
                 Rt.pts(1).sD, 1 - 1/24));
ok = chk(ok, Rt.pts(1).ok, sprintf('and it is certified there (%s)', Rt.pts(1).reason));

% a rib that cannot advance must report the failure, not return a short
% success: one point, zero bisections allowed, a full HALF period per step
R2 = rib_from_crossing(C0, B, anc, struct('nD', 2, 'direction', -1, 'nPts', 1, 'maxBisect', 0, 'wallSec', 60));
ok = chk(ok, isempty(R2.pts) && ~isempty(R2.stop), sprintf('impossible step reported: %s', R2.stop));

if ok, fprintf('TEST_RIB_FROM_CROSSING: ALL PASS\n'); else, fprintf('TEST_RIB_FROM_CROSSING: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
