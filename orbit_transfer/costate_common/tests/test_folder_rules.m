function ok = test_folder_rules()
%% Purpose:
%
%   THE FOLDER'S OWN RULES, enforced (README -> "Rules: what belongs here").
%   The cleanup of 2026-09-16 removed job control, campaign-only tests, dead
%   code and a campaign dependency that had accumulated over months without
%   anyone having to justify them. This test is what makes the next drift
%   fail in seconds instead:
%
%     1. NO CAMPAIGN ON THE PATH: no file here may addpath a campaign
%        folder. A library that reaches into a campaign is not a library.
%     2. NO CAMPAIGN DATA: no executable line may name a campaign folder
%        (comments may; a self-demo that loads a campaign artifact must be
%        an EXEMPTION with a reason).
%     3. HEADERS: every file opens with the pumpkyn %% Purpose block.
%     4. NO SUPPRESSIONS: no %#ok pragmas, here or in tests -- a warning
%        worth silencing is worth fixing, and one of them hid dead code.
%     5. NO i/j LOOP VARIABLES (they are sqrt(-1) in MATLAB).
%     6. EVERY FILE IS EXERCISED: some test in tests/ names it, or it is
%        an EXEMPTION carrying a written reason.
%     7. EXEMPTIONS DO NOT ROT: every exemption must name a file that still
%        exists and still needs it (one that has since acquired a test is a
%        failure, not a pass) -- the list may only shrink by accident.
%
%   Rule 6 is the one with teeth: an exception has to be written down, in
%   this file, with its reason, where the next reader sees it.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every rule holds
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
cc = fileparts(fileparts(mfilename('fullpath')));       % costate_common
files = dir(fullfile(cc, '*.m'));
tests = dir(fullfile(cc, 'tests', '*.m'));
campaigns = {'DRO_tulip', 'HALO_tulip', 'DPO_tulip', 'HALO_HALO', 'GTO_tulip'};

%% The exemptions, each with the reason it is one:
%  (kind: 'entry'    = no calling code, run by hand
%         'elsewhere'= a test outside this folder exercises it, and must
%         'indirect' = another test here drives it end to end
%         'gap'      = genuinely untested, with the TODO that says so
%         'data'     = a self-demo may read a named campaign artifact)
EX = { ...
 'conj_catalog_pass',    'entry',     'entry point: the catalog-scale conjugate sweep, run by hand over a shipped catalog';
 'gates_catalog_pass',   'entry',     'entry point: the catalog-scale hypothesis gates, run by hand';
 'golden_cells',         'entry',     'entry point: it IS the regression others are run against';
 'casadi_mintime_dro',   'elsewhere', 'DRO_tulip/direct/certify/tests/test_minenergy_objective -- a bitwise regression against a stored campaign fixture';
 'assert_periodic_orbit','indirect',  'driven by tests/test_survey_family_bounds (the survey refuses a non-closing orbit)';
 'current_pool',         'indirect',  'driven by tests/test_run_capped and the second-order parallel test';
 'certify_dro_mintime',  'gap',       'NO test of its own; arrived 2026-09-16 with the ladder engine and is covered only through campaign re-solves (TODO.md)';
 'ladder_endpoints',     'data',      'its nargin==0 demo rebuilds the DRO fine sheet endpoints from that campaign results file' };

%% 1-5: the mechanical rules, file by file:
badPath = {};  badData = {};  badHead = {};  badPragma = {};  badLoop = {};
for k = 1:numel(files)
    nm = files(k).name;
    txt = fileread(fullfile(cc, nm));
    lines = strsplit(txt, newline);
    for L = 1:numel(lines)
        raw = lines{L};
        code = strip_comment(raw);
        if contains(code, 'addpath') && any(cellfun(@(c) contains(code, c), campaigns))
            badPath{end+1} = sprintf('%s:%d', nm, L);
        end
        if any(cellfun(@(c) contains(code, c), campaigns)) && ~isExempt(EX, nm, 'data')
            badData{end+1} = sprintf('%s:%d  %s', nm, L, strtrim(code));
        end
        if ~isempty(regexp(code, 'for\s+[ij]\s*=', 'once'))
            badLoop{end+1} = sprintf('%s:%d', nm, L);
        end
    end
    if ~contains(txt, '%% Purpose:'), badHead{end+1} = nm; end
    if contains(txt, '%#ok'), badPragma{end+1} = nm; end
end
for k = 1:numel(tests)
    txt = fileread(fullfile(cc, 'tests', tests(k).name));
    % this file is excepted because it carries the pattern as DATA
    if contains(txt, '%#ok') && ~strcmp(tests(k).name, 'test_folder_rules.m')
        badPragma{end+1} = ['tests/' tests(k).name];
    end
    lines = strsplit(txt, newline);
    for L = 1:numel(lines)
        if ~isempty(regexp(strip_comment(lines{L}), 'for\s+[ij]\s*=', 'once'))
            badLoop{end+1} = sprintf('tests/%s:%d', tests(k).name, L);
        end
    end
end
ok = report(ok, 'no file adds a campaign folder to the path', badPath);
ok = report(ok, 'no executable line names a campaign folder', badData);
ok = report(ok, 'every file opens with %% Purpose', badHead);
ok = report(ok, 'no %#ok pragmas (this file excepted: it holds the pattern as data)', badPragma);
ok = report(ok, 'no i/j loop variables', badLoop);

%% 6: every file is exercised by a test here, or exempt with a reason:
% THIS file is left out of the corpus deliberately: it holds every exempted
% name as a STRING, and counting those would make each exemption look tested.
allTests = '';
for k = 1:numel(tests)
    if strcmp(tests(k).name, 'test_folder_rules.m'), continue, end
    allTests = [allTests newline fileread(fullfile(cc, 'tests', tests(k).name))];
end
untested = {};
for k = 1:numel(files)
    stem = files(k).name(1:end-2);
    if named(allTests, stem), continue, end
    if isExemptAny(EX, stem), continue, end
    untested{end+1} = stem;
end
ok = report(ok, 'every file is exercised by a test, or exempt with a reason', untested);

%% 7: the exemptions still apply:
stale = {};
for k = 1:size(EX, 1)
    stem = EX{k,1};
    if ~isfile(fullfile(cc, [stem '.m']))
        stale{end+1} = sprintf('%s (no such file)', stem);
    elseif ~strcmp(EX{k,2}, 'data') && named(allTests, stem)
        stale{end+1} = sprintf('%s (a test names it now -- drop the exemption)', stem);
    end
end
ok = report(ok, 'no stale exemptions', stale);

fprintf('  %d files, %d tests, %d exemption(s): %s\n', numel(files), numel(tests), size(EX,1), ...
        strjoin(EX(:,1).', ', '));
if ok, fprintf('TEST_FOLDER_RULES: ALL PASS\n');
else,  fprintf('TEST_FOLDER_RULES: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = report(ok, label, offenders)
%% Purpose:
%
%   One rule's verdict, naming every offender.
%
if isempty(offenders)
    fprintf('  PASS  %s\n', label);
else
    fprintf('  FAIL  %s\n', label);
    fprintf('          %s\n', offenders{:});
    ok = false;
end
end

function tf = named(txt, stem)
%% Purpose:
%
%   True when the tests CALL this function -- 'stem(' or '@stem' -- rather
%   than merely mentioning it in prose. A test that names a function in a
%   comment does not exercise it.
%
e = regexptranslate('escape', stem);
tf = ~isempty(regexp(txt, ['(\<' e '\s*\()|(@' e '\>)'], 'once'));
end

function tf = isExempt(EX, nameOrFile, kind)
%% Purpose:
%
%   True when this file carries an exemption of the given kind.
%
stem = nameOrFile;
if endsWith(stem, '.m'), stem = stem(1:end-2); end
tf = any(strcmp(EX(:,1), stem) & strcmp(EX(:,2), kind));
end

function tf = isExemptAny(EX, stem)
%% Purpose:
%
%   True when this file carries any exemption.
%
tf = any(strcmp(EX(:,1), stem));
end

function code = strip_comment(line)
%% Purpose:
%
%   The executable part of a MATLAB line: everything before the first % that
%   is not inside a single-quoted string.
%
code = line;
inStr = false;
for k = 1:numel(line)
    ch = line(k);
    if ch == '''', inStr = ~inStr; end
    if ch == '%' && ~inStr, code = line(1:k-1); return, end
end
end
