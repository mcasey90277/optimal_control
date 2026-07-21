% test_sosc_gate_wiring  Unit test of apply_sosc_gate's tiered-gate mapping
% (process/DESIGN_sosc.md sec 11.6): only FAIL (a proven saddle) demotes a
% feasibility-certified report; PASS/WEAK_MIN/INCONCLUSIVE/ERROR all keep
% report.certified true (WEAK_MIN is a POSITIVE certificate, not merely
% non-demoting). NO SOLVE -- pure struct-mapping test.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;

rep = struct('certified', true);

gFail = apply_sosc_gate(rep, struct('verdict','FAIL','status','feasible-only','reason','x'));
assert(~gFail.certified, 'FAIL must demote');
assert(strcmp(gFail.sosc.verdict,'FAIL'), 'sosc struct must be attached');

gPass = apply_sosc_gate(rep, struct('verdict','PASS','status','certified-sosc','reason','y'));
assert(gPass.certified, 'PASS keeps certified');

gWeak = apply_sosc_gate(rep, struct('verdict','WEAK_MIN','status','certified-weak-min','reason','w'));
assert(gWeak.certified, 'WEAK_MIN keeps certified (positive certificate)');

gInc = apply_sosc_gate(rep, struct('verdict','INCONCLUSIVE', ...
    'status','certified-feasibility+sosc-inconclusive','reason','z'));
assert(gInc.certified, 'INCONCLUSIVE keeps certified (annotated)');

gErr = apply_sosc_gate(rep, struct('verdict','ERROR', ...
    'status','certified-feasibility+sosc-inconclusive','reason','e'));
assert(gErr.certified, 'ERROR keeps certified (annotated, non-demoting)');

fprintf('test_sosc_gate_wiring PASSED\n');
