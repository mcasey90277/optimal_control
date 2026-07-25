function test_foc_report
% TEST_FOC_REPORT  Test suite for foc_report — first-order optimality print+save.
%
% Tests the full house header, status reporting, advisory line, and sidecar save.

import matlab.unittest.TestCase
import matlab.unittest.TestSuite

% Setup path (tests dir is in verify_common, so parent is verify_common)
addpath(fileparts(fileparts(mfilename('fullpath'))));

testCase = TestCase.forInteractiveUse;

% --- Test 1: Full header and content with all fields populated --------
fprintf('\n[Test 1] Full report with all fields, sidecar save...\n');
rep = build_fake_rep();
out_text = evalc('foc_report(rep, ''unittest'', tempdir)');

assert(contains(out_text, 'FIRST-ORDER OPTIMALITY REPORT'), ...
    'Test 1 FAIL: output missing FIRST-ORDER OPTIMALITY REPORT header');
assert(contains(out_text, 'unittest'), ...
    'Test 1 FAIL: output missing tag name');
assert(contains(out_text, 'ADVISORY'), ...
    'Test 1 FAIL: output missing ADVISORY line');
assert(contains(out_text, '(report-only burn-in: does not alter certified status)'), ...
    'Test 1 FAIL: output missing advisory disclaimer');

% Check sidecar exists
sidecar_path = fullfile(tempdir, 'foc_unittest.mat');
assert(isfile(sidecar_path), ...
    sprintf('Test 1 FAIL: sidecar not created at %s', sidecar_path));

% Verify rep is in the sidecar
saved_data = load(sidecar_path);
assert(isfield(saved_data, 'rep'), ...
    'Test 1 FAIL: sidecar missing rep field');
assert(isequal(saved_data.rep, rep), ...
    'Test 1 FAIL: saved rep does not match input');

fprintf('  PASS: report generated, sidecar saved\n');

% --- Test 2: Empty resDir — print only, no save --------
fprintf('\n[Test 2] Empty resDir — print only, no save...\n');
rep2 = build_fake_rep();
out_text2 = evalc('foc_report(rep2, ''printonly'', '''')');

assert(contains(out_text2, 'FIRST-ORDER OPTIMALITY REPORT'), ...
    'Test 2 FAIL: print-only missing header');

sidecar_path2 = fullfile('', 'foc_printonly.mat');
% should not exist in current dir
[~,current_dir] = fileparts(pwd);
assert(~isfile(fullfile(pwd, 'foc_printonly.mat')), ...
    'Test 2 FAIL: sidecar should not be created when resDir empty');

fprintf('  PASS: print-only mode works\n');

% --- Test 3: Omitted resDir (default) — print only --------
fprintf('\n[Test 3] Omitted resDir — print only...\n');
rep3 = build_fake_rep();
out_text3 = evalc('foc_report(rep3, ''omitted'')');

assert(contains(out_text3, 'FIRST-ORDER OPTIMALITY REPORT'), ...
    'Test 3 FAIL: omitted resDir missing header');

fprintf('  PASS: omitted resDir defaults to print-only\n');

% --- Test 4: NaN values print as '--' --------
fprintf('\n[Test 4] NaN values print as ''--''...\n');
rep4 = build_fake_rep();
rep4.lamMassEndRel = NaN;
rep4.lamTimeCoV = NaN;
out_text4 = evalc('foc_report(rep4, ''nantest'', '''')');

assert(contains(out_text4, '--'), ...
    'Test 4 FAIL: NaN not rendered as --');

fprintf('  PASS: NaN values render as --\n');

% --- Test 5: PASS/FAIL status lines --------
fprintf('\n[Test 5] Status lines (PASS/FAIL)...\n');
rep5a = build_fake_rep();
rep5a.kktStatInf = 1e-7;  % << 1e-6 threshold
rep5a.dirTanMax = 1e-7;   % << 1e-6 threshold
rep5a.signPct = 99.5;     % >= 99 threshold
rep5a.lamMassEndRel = 5e-4;  % < 1e-3 threshold
rep5a.singularArcNodes = 0;  % == 0
rep5a.sdotMinRel = 2e-3;  % > 1e-3 threshold
out_text5a = evalc('foc_report(rep5a, ''pass_test'', '''')');
assert(contains(out_text5a, 'PASS'), ...
    'Test 5 FAIL: good values should produce PASS statuses');

rep5b = build_fake_rep();
rep5b.kktStatInf = 2e-6;  % >> 1e-6 threshold
out_text5b = evalc('foc_report(rep5b, ''fail_test'', '''')');
assert(contains(out_text5b, 'FAIL'), ...
    'Test 5 FAIL: bad kktStatInf should produce FAIL');

fprintf('  PASS: status lines work correctly\n');

% --- Test 6: ipopt verdict (if present) --------
fprintf('\n[Test 6] ipopt verdict line (optional)...\n');
rep6 = build_fake_rep();
rep6.ipopt = struct('verdict', 'locally optimal');
out_text6 = evalc('foc_report(rep6, ''ipopt_test'', '''')');
assert(contains(out_text6, 'IPOPT inertia'), ...
    'Test 6 FAIL: ipopt.verdict not printed');
assert(contains(out_text6, 'locally optimal'), ...
    'Test 6 FAIL: ipopt verdict value not printed');

fprintf('  PASS: ipopt verdict line printed\n');

% --- Test 7: ipopt empty/missing does not crash --------
fprintf('\n[Test 7] ipopt missing or empty...\n');
rep7a = build_fake_rep();
% no ipopt field
out_text7a = evalc('foc_report(rep7a, ''no_ipopt'', '''')');
assert(contains(out_text7a, 'FIRST-ORDER OPTIMALITY REPORT'), ...
    'Test 7a FAIL: missing ipopt crashed');

rep7b = build_fake_rep();
rep7b.ipopt = struct('verdict', '');  % empty
out_text7b = evalc('foc_report(rep7b, ''empty_ipopt'', '''')');
assert(contains(out_text7b, 'FIRST-ORDER OPTIMALITY REPORT'), ...
    'Test 7b FAIL: empty ipopt crashed');
assert(~contains(out_text7b, 'IPOPT inertia'), ...
    'Test 7b FAIL: empty ipopt should not print verdict line');

fprintf('  PASS: ipopt absent/empty handled gracefully\n');

% --- Test 8: horizonNote included --------
fprintf('\n[Test 8] horizonNote printed...\n');
rep8 = build_fake_rep();
out_text8 = evalc('foc_report(rep8, ''horizon_test'', '''')');
assert(contains(out_text8, 'Horizon:'), ...
    'Test 8 FAIL: horizonNote not printed');

fprintf('  PASS: horizonNote printed\n');

fprintf('\n========== ALL TESTS PASSED ==========\n');

end

%--------------------------------------------------------------------------
function rep = build_fake_rep()
% Build a fully populated fake rep struct for testing.

rep = struct();

% Numeric fields
rep.kktStatInf = 1e-7;
rep.sLag = +1;
rep.dirTanMax = 1e-8;
rep.dirTanMed = 5e-9;
rep.signPct = 99.8;
rep.Sd = randn(1, 10);  % [1 x (N+1)]
rep.lam = randn(4, 10);  % fake [nx x (N+1)]
rep.lamTimeCoV = 0.002;
rep.lamTimeEnd = +0.567;
rep.lamMassEndRel = 2e-4;
rep.singularArcNodes = 0;
rep.sdotMinRel = 5e-3;
rep.nSwitches = 8;

% Text fields
rep.horizonNote = 'fixed t_f: H=const generally nonzero; constancy via lamTimeCoV';
rep.checksRun = {'kktStat'; 'signLaw'; 'dirTangential'; 'costates'; 'lamTime'; 'transversality'; 'singularArc'; 'sdotRegular'};

% Logical field
rep.pass = true;

end
