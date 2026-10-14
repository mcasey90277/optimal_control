function ok = test_examples_run()
% TEST_EXAMPLES_RUN  Rehearsal of the recipient self-test: a scratch folder
% holds ONLY the exported catalog and the deliverable file list; a CHILD
% MATLAB (matlab -batch, pumpkynPie startup, then cd to the scratch folder)
% runs both example scripts. Asserts exit 0, '|dz|' printed, two
% 'ALL SECTIONS DONE', and that every listed function resolves INSIDE the
% scratch folder (nothing from the repo was used). Takes several minutes.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
[L, roots] = deliverable_file_list();
record = fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat');
scratch = tempname;  mkdir(scratch);
export_deliverable_catalog(record, fullfile(scratch, 'costate_catalog_dro_tulip_70mN.mat'));
for k = 1:numel(L)
    src = fullfile(roots.(L(k).from), L(k).file);
    ok = chk(ok, isfile(src), sprintf('listed file exists: %s/%s', L(k).from, L(k).file));
    if isfile(src), copyfile(src, scratch); end
end
ok = chk(ok, numel(dir(fullfile(scratch, '*.m'))) == numel(L), 'scratch holds exactly the listed files');

names = erase({L.file}, '.m');
whichCmd = sprintf('which(''%s'')\n', names{:});
childScript = fullfile(scratch, 'run_child_.txt');
fid = fopen(childScript, 'w');
fprintf(fid, ['cd(''/Users/msc/Desktop/proj7/external/pumpkynPie''); startup(); cd(''%s'');\n' ...
              'run(''costate_catalog_example.m''); run(''status_layer_example.m'');\n' ...
              'fprintf(''WHICHCHECK\\n'');\n'], scratch);
for k = 1:numel(names)
    fprintf(fid, 'fprintf(''WHICH %%s\\n'', which(''%s''));\n', names{k});
end
fclose(fid);
logFile = fullfile(scratch, 'child.log');
tic;
cmd = sprintf('/Applications/MATLAB_R2026a.app/bin/matlab -batch "run(''%s'')" > %s 2>&1', ...
              strrep(childScript, '.txt', '.m'), logFile);
copyfile(childScript, strrep(childScript, '.txt', '.m'));
[status, ~] = system(cmd);
dt = toc;
out = fileread(logFile);
fprintf('child run: %.0f s, exit %d, log %s\n', dt, status, logFile);
ok = chk(ok, status == 0, 'child MATLAB exit code 0');
ok = chk(ok, contains(out, '|dz|'), 'child output contains |dz|');
ok = chk(ok, numel(strfind(out, 'ALL SECTIONS DONE')) == 2, 'two ALL SECTIONS DONE lines');
tok = regexp(out, 'WHICH ([^\n]*)', 'tokens');
sc = fileparts(fullfile(scratch, 'x'));
allIn = numel(tok) == numel(names);
for k = 1:numel(tok)
    p = strtrim(tok{k}{1});
    inside = startsWith(p, [sc filesep]) || startsWith(p, ['/private' sc filesep]);
    if ~inside, fprintf('  not in scratch: %s\n', p); end
    allIn = allIn && inside;
end
ok = chk(ok, allIn, 'every listed function resolves inside the scratch folder');
lines = splitlines(out);  lines = lines(~cellfun(@isempty, lines));
fprintf('--- child log tail ---\n%s\n', strjoin(lines(max(1, end-60):end), newline));
if ok, fprintf('test_examples_run: ALL PASS\n'); else, fprintf('test_examples_run: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
