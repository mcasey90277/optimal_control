function ok = test_report_optimality()
%% Purpose:
%
%   Tests report_optimality -- the two-section optimality report for one
%   transfer. Section 1 is the NECESSARY conditions (Pontryagin: the
%   boundary-value residual, the Hamiltonian along the arc, the flown
%   arrival, the independent witness). Section 2 is the SUFFICIENCY
%   hypotheses (strengthened Legendre, the all-burn switching function, no
%   abnormal lift, no conjugate time).
%
%   The report must never overstate: a missing diagnostic reads NOT CHECKED,
%   never PASS, and the overall claim is made only when every line passes.
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

good = struct('ok', true, 'reason', 'certified', 'normR', 2.1e-11, 'flyKm', 0.0004, ...
    'flyVms', 0.0002, 'dz', 0, 'flyKmWitness', 0.0004, 'conj', 1, 'tfDays', 17.7976, ...
    'dvKms', 0.7485, 'propellantKg', 12.20, 'sD', 0, 'sA', 0.0754, ...
    'g', struct('minLamV', 3.24, 'minQmt', 3.82, 'dimS', 1, 'Hresid', 4e-12, ...
                'nullResid', 1e-9, 'svRatio', 2e-7));

good.h6Ok = true;  good.h6Margin = 9.0;  good.adjErr = 6e-9;  good.dirGap = 0;
good.lamMf = 3e-12;  good.flyVmsWitness = 0.0002;
R = report_optimality(good, struct('quiet', true));
% cross-checks live OUTSIDE the PMP conjunction: their failure must not read
% as "not an extremal"
bad_x = good;  bad_x.dz = 1;  bad_x.flyKmWitness = 5e4;
Rx = report_optimality(bad_x, struct('quiet', true));
ok = chk(ok, Rx.necessary && ~Rx.crossChecks && ~Rx.claim, ...
         sprintf('a failed cross-check leaves necessary intact but blocks the claim (%s)', Rx.verdict));
% three-way status survives aggregation: an UNRESOLVED conjugate verdict is
% not a failed hypothesis
unr = good;  unr.conj = -1;  unr.conjVerdict = 'ENDPOINT';
Ru = report_optimality(unr, struct('quiet', true));
ok = chk(ok, ~Ru.claim && contains(lower(Ru.verdict), 'unresolved'), ...
         sprintf('ENDPOINT propagates as UNRESOLVED, not as a failure (%s)', Ru.verdict));
% H6 is an instrument-validity prerequisite; missing means NOT CHECKED
noh = good;  noh = rmfield(noh, {'h6Ok', 'h6Margin'});
Rh = report_optimality(noh, struct('quiet', true));
ok = chk(ok, ~Rh.claim && any(contains(Rh.lines, 'H6')) && any(contains(Rh.lines, 'NOT CHECKED')), ...
         'missing H6 reads NOT CHECKED and blocks the claim');
ok = chk(ok, R.necessary && R.sufficient && R.claim, ...
         sprintf('a complete certificate passes both sections (%s)', R.verdict));
ok = chk(ok, numel(R.lines) >= 12 && all(cellfun(@(s) ischar(s), R.lines)), ...
         sprintf('%d report lines produced', numel(R.lines)));

% A MISSING diagnostic reads NOT CHECKED and blocks the section it is in --
% never PASS. The Hamiltonian residual is produced by the gate computation,
% so with no gates it is genuinely unverified, and a report claiming the
% necessary conditions "still pass" would itself be the overstatement this
% function exists to prevent.
noG = good;  noG.g = [];
R2 = report_optimality(noG, struct('quiet', true));
ok = chk(ok, ~R2.necessary && ~R2.sufficient && ~R2.claim, ...
         sprintf('missing gates block BOTH sections (%s)', R2.verdict));
ok = chk(ok, any(contains(R2.lines, 'NOT CHECKED')), 'and the missing lines say NOT CHECKED');
ok = chk(ok, contains(R2.verdict, 'NOT CHECKED') || contains(R2.verdict, 'not checked'), ...
         'the verdict says what was not checked rather than what failed');

% a conjugate failure must sink sufficiency but leave the first-order truth
badC = good;  badC.conj = 0;  badC.ok = false;  badC.reason = 'conjugate test verdict 0';
R3 = report_optimality(badC, struct('quiet', true));
ok = chk(ok, R3.necessary && ~R3.sufficient && ~R3.claim, ...
         'a conjugate failure sinks sufficiency and leaves the first-order truth');

% a bad residual must sink BOTH -- the second-order test is meaningless off
% an extremal, which is the whole reason the first section exists
badR = good;  badR.normR = 1e-3;
R4 = report_optimality(badR, struct('quiet', true));
ok = chk(ok, ~R4.necessary && ~R4.claim, 'a bad residual sinks the whole report');
ok = chk(ok, any(contains(R4.lines, 'meaningless')) || ~R4.sufficient, ...
         'and sufficiency is not asserted off an extremal');

if ok, fprintf('TEST_REPORT_OPTIMALITY: ALL PASS\n'); else, fprintf('TEST_REPORT_OPTIMALITY: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
