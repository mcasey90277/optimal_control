function ok = test_verify_with_pumpkyn()
%% Purpose:
%
%   Tests verify_with_pumpkyn -- the INDEPENDENT verification step, made
%   visible: hand our converged costates to pumpkyn's own minimum-time
%   solver and show, component by component, that they do not move.
%
%   Two things this must get right. It must reproduce the agreement the
%   certifier measured (same solver, same question). And it must actually
%   MOVE when given something wrong -- a check that returns "no change"
%   whatever you feed it is not a check, and that was the live concern in
%   the 2026-09-10 review.
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
T = certify_crossing(anc.p, anc.sA, B, anc);
assert(T.ok, 'fixture: the anchor must certify (%s)', T.reason);

V = verify_with_pumpkyn(T, B, struct('quiet', true));
ok = chk(ok, V.converged, sprintf('the independent solve ran (%s)', V.note));
ok = chk(ok, abs(V.dz - T.dz) < 1e-12, ...
         sprintf('reproduces the certifier''s agreement: %.2e vs %.2e', V.dz, T.dz));
ok = chk(ok, numel(V.z8) == 8 && numel(V.z8Pumpkyn) == 8 && numel(V.dComponent) == 8, ...
         'both costate vectors and their component differences are returned');
ok = chk(ok, V.flyKm < 1, sprintf('and pumpkyn''s own solution flies to the target (%.4f km)', V.flyKm));
ok = chk(ok, V.moved == false, 'verdict: the costates did NOT move');

% THE CONTROL EXPERIMENT: perturb the costates and require the solver to move
Tb = T;  Tb.z(1:7) = 1.5*Tb.z(1:7);
Vb = verify_with_pumpkyn(Tb, B, struct('quiet', true));
ok = chk(ok, Vb.dz > 1, sprintf('a perturbed seed MOVES: |dz| = %.3g', Vb.dz));
ok = chk(ok, Vb.moved == true, 'verdict: the costates moved, so the check can fail');

if ok, fprintf('TEST_VERIFY_WITH_PUMPKYN: ALL PASS\n'); else, fprintf('TEST_VERIFY_WITH_PUMPKYN: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
