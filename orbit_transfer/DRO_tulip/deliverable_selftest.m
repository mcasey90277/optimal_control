function R = deliverable_selftest(outDir, opts)
%% Purpose:
%
%   The RECIPIENT self-test of a costate-library deliverable folder: a
%   fresh child MATLAB (matlab -batch, launched with system) does what the
%   recipient would -- pumpkynPie startup, cd into the folder, run both
%   example scripts, call every shipped helper's bare self-demo -- and then
%   proves the folder is complete and clean:
%
%     1. every .m in the folder which-resolves INSIDE the folder (the child
%        inherits this machine's startup path, where repo copies of the same
%        functions may sit, so resolution has to be checked, not assumed);
%     2. matlab.codetools.requiredFilesAndProducts over every .m returns
%        only files under the folder, matlabroot, or the pumpkyn /
%        pumpkynPie / movie2gifLite roots -- anything else is named;
%     3. no shipped .m and not the README contains 'addpath(' or '/Users/'
%        (ruling R7). Ruling R11: an addpath whose argument is a folder
%        under matlabroot (a MathWorks toolbox restore, the ShowEarth
%        guard) is allowed. Statically: the argument must be a single
%        variable v, and the same file must bind v to matlabroot -- an
%        assignment of v that calls toolboxdir( or matlabroot, or a guard
%        startsWith(v, matlabroot). Every other addpath fails.
%
%   The child prints 'SELFTEST OK' or one 'SELFTEST FAIL: <reason>' line
%   per finding; the parent parses them. Repo-only: NOT shipped.
%
%  ASSUMPTIONS / NOTES:
%
% • The log and the child script live in a work folder under tempdir,
%   never inside outDir. Files the child creates inside outDir (a demo's
%   movie, say) are deleted afterwards and listed in R.removed, so the
%   folder is zipped exactly as built.
% • The self-demos are found, not typed: every FUNCTION file in the folder
%   whose code has an 'if nargin == 0' branch (R.demos).
% • A missing root (e.g. movie2gifLite not installed) is simply not
%   allowed; it never widens the allow-list.
%
%% Inputs:
%
%  outDir                   char                    the deliverable folder
%
%  opts                     struct   (optional)     .matlab [R2026a binary]
%                                                   .pumpkynPie [proj7/external/
%                                                   pumpkynPie]
%                                                   .demos {cellstr} override
%                                                   the detected self-demos
%                                                   .examples [{'costate_catalog_example.m',
%                                                   'status_layer_example.m'}]
%
%% Outputs:
%
%  R                        struct                  .ok (logical) .log (char,
%                                                   log path) .reasons
%                                                   {cellstr} .demos {cellstr}
%                                                   .removed {cellstr} .status
%                                                   (child exit) .seconds
%
%% Revision History:
%  M. Casey                                                   (c) 10/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
mlab   = getf(opts, 'matlab', '/Applications/MATLAB_R2026a.app/bin/matlab');
pieDir = getf(opts, 'pumpkynPie', '/Users/msc/Desktop/proj7/external/pumpkynPie');
exs    = {'costate_catalog_example.m', 'status_layer_example.m'};
if isfield(opts, 'examples'), exs = opts.examples; end      % {} = run none
assert(isfolder(outDir), 'deliverable_selftest:outDir', 'no folder %s', outDir);
D0 = dir(outDir);  outDir = D0(1).folder;                  % absolute
demos = findDemos(outDir);
if isfield(opts, 'demos'), demos = opts.demos; end           % {} = run none

before = listFiles(outDir);
work = tempname;  mkdir(work);
child = fullfile(work, 'selftest_child_.m');
logFile = fullfile(work, 'selftest.log');
writeChild(child, outDir, pieDir, exs, demos);

fprintf('[selftest] child MATLAB on %s (%d self-demos: %s); log %s\n', ...
        outDir, numel(demos), strjoin(demos, ', '), logFile);
t0 = tic;
% the user's startup changes folder, so the child is found by path, not by cwd
[status, ~] = system(sprintf('''%s'' -batch "addpath(''%s''); selftest_child_" > ''%s'' 2>&1', ...
                             mlab, work, logFile));
secs = toc(t0);
txt = '';
if isfile(logFile), txt = fileread(logFile); end

% byproducts the child left inside the deliverable
after = listFiles(outDir);
removed = setdiff(after, before);
for k = 1:numel(removed), delete(fullfile(outDir, removed{k})); end

tok = regexp(txt, 'SELFTEST FAIL: ([^\n]*)', 'tokens');
reasons = cellfun(@(x) strtrim(x{1}), tok, 'UniformOutput', false);
if status ~= 0 && isempty(reasons)
    reasons = {sprintf('child MATLAB exited %d without a verdict (see the log)', status)};
end
okLine = ~isempty(regexp(txt, '(^|\n)SELFTEST OK', 'once'));
R = struct('ok', status == 0 && okLine && isempty(reasons), 'log', logFile, ...
           'reasons', {reasons}, 'demos', {demos}, 'removed', {removed}, ...
           'status', status, 'seconds', secs);
if R.ok
    fprintf('[selftest] SELFTEST OK (%.0f s)\n', secs);
else
    fprintf('[selftest] SELFTEST FAIL (%.0f s, exit %d):\n  %s\n', secs, status, strjoin(reasons, '\n  '));
end
if ~isempty(removed)
    fprintf('[selftest] removed child byproducts from the folder: %s\n', strjoin(removed, ', '));
end
end

% ------------------------------------------------------------------------
function demos = findDemos(outDir)
% FINDDEMOS  Function files in outDir with an 'if nargin == 0' branch.
% INPUTS: outDir (char).  OUTPUTS: demos {cellstr} (names, no .m).
D = dir(fullfile(outDir, '*.m'));
has = false(1, numel(D));
for k = 1:numel(D)
    t = fileread(fullfile(outDir, D(k).name));
    isFun = ~isempty(regexp(t, '^\s*function\s', 'once', 'lineanchors'));
    has(k) = isFun && ~isempty(regexp(t, '^\s*if\s+nargin\s*==\s*0', 'once', 'lineanchors'));
end
demos = erase({D(has).name}, '.m');
end

% ------------------------------------------------------------------------
function names = listFiles(d)
% LISTFILES  Bare names of the files (not folders) in d.  INPUTS: d (char).
% OUTPUTS: names {cellstr}.
D = dir(d);  D = D(~[D.isdir]);  names = {D.name};
end

% ------------------------------------------------------------------------
function writeChild(child, outDir, pieDir, exs, demos)
% WRITECHILD  Write the child function file the -batch MATLAB runs.
% INPUTS: child (path); outDir; pieDir (char); exs, demos {cellstr}.
% OUTPUTS: none.
q = @(s) ['''' strrep(s, '''', '''''') ''''];
cellLit = @(c) ['{' strjoin(cellfun(q, c, 'UniformOutput', false), ', ') '}'];
L = {
'function selftest_child_()'
'% SELFTEST_CHILD_  Generated by deliverable_selftest; runs in a child MATLAB.'
'% INPUTS: none.  OUTPUTS: none (prints SELFTEST OK / SELFTEST FAIL lines).'
'bad = {};'
['outDir = ' q(outDir) ';']
['cd(' q(pieDir) ');  startup();  cd(outDir);']
'od = pwd;'
['exs = ' cellLit(exs) ';']
'for k = 1:numel(exs)'
'    try, runIsolated(exs{k}); fprintf(''SELFTEST step: example %s ran\n'', exs{k});'
'    catch e, bad{end+1} = sprintf(''example %s threw: %s'', exs{k}, e.message); end'
'    close all force;  cd(od);'
'end'
['demos = ' cellLit(demos) ';']
'for k = 1:numel(demos)'
'    try, feval(demos{k}); fprintf(''SELFTEST step: self-demo %s ran\n'', demos{k});'
'    catch e, bad{end+1} = sprintf(''self-demo %s threw: %s'', demos{k}, e.message); end'
'    close all force;  cd(od);'
'end'
'% 1. which-resolution inside the folder'
'M = dir(fullfile(od, ''*.m''));'
'for k = 1:numel(M)'
'    w = which(erase(M(k).name, ''.m''));'
'    if ~inside(w, {od, outDir}), bad{end+1} = sprintf(''%s resolves outside the folder: %s'', M(k).name, w); end'
'end'
'% 2. dependency scan'
'roots = {od, outDir, matlabroot};'
'r = pkgRoot(which(''pumpkyn.cr3bp.tfMin''), ''+pumpkyn'');  if ~isempty(r), roots{end+1} = r; end'
'r = pkgRoot(which(''pumpkynPie.plot.SatelliteAnimator''), ''+pumpkynPie'');  if ~isempty(r), roots{end+1} = r; end'
'w = which(''movie2gifLite'');  if ~isempty(w), roots{end+1} = fileparts(w); end'
'fprintf(''SELFTEST step: allowed roots: %s\n'', strjoin(roots, '' ; ''));'
'files = fullfile(od, {M.name});'
'try'
'    fl = matlab.codetools.requiredFilesAndProducts(files);'
'    for k = 1:numel(fl)'
'        if ~inside(fl{k}, roots), bad{end+1} = sprintf(''dependency outside the folder and pumpkyn: %s'', fl{k}); end'
'    end'
'    fprintf(''SELFTEST step: dependency scan of %d files -> %d required files\n'', numel(files), numel(fl));'
'catch e, bad{end+1} = sprintf(''dependency scan threw: %s'', e.message); end'
'% 3. no addpath( and no /Users/ in shipped code or the README'
'G = [files, {fullfile(od, ''README.md'')}];'
'for k = 1:numel(G)'
'    if ~isfile(G{k}), continue, end'
'    t = splitlines(fileread(G{k}));'
'    for n = 1:numel(t)'
'        ln = t{n};'
'        if contains(ln, ''addpath('') && ~mathworksRestore(ln, t), bad{end+1} = sprintf(''%s line %d has addpath( of a folder not bound to matlabroot'', G{k}, n); end'
'        if contains(ln, ''/Users/''), bad{end+1} = sprintf(''%s line %d has /Users/'', G{k}, n); end'
'    end'
'end'
'if isempty(bad)'
'    fprintf(''SELFTEST OK\n'');'
'else'
'    for k = 1:numel(bad), fprintf(''SELFTEST FAIL: %s\n'', strrep(bad{k}, newline, '' '')); end'
'    exit(1);'
'end'
'end'
''
'function tf = mathworksRestore(ln, lines)'
'% MATHWORKSRESTORE  Ruling R11: addpath(v) with v a single variable that the'
'% same file binds to matlabroot (v = ...toolboxdir(/matlabroot..., or a guard'
'% startsWith(v, matlabroot)).  INPUTS: ln (char); lines {cellstr}.  OUTPUTS: tf.'
'tf = false;'
'tk = regexp(ln, ''addpath\(\s*(\w+)\s*\)'', ''tokens'', ''once'');'
'if isempty(tk) || numel(regexp(ln, ''addpath\('')) > 1, return, end'
'v = tk{1};  body = strjoin(lines, newline);'
'assigned = ~isempty(regexp(body, [''(^|\n)\s*'' v ''\s*=[^\n]*(toolboxdir\(|matlabroot)''], ''once''));'
'guarded = ~isempty(regexp(body, [''startsWith\(\s*'' v ''\s*,\s*matlabroot''], ''once''));'
'tf = assigned || guarded;'
'end'
''
'function runIsolated(f)'
'% RUNISOLATED  run a script in its own workspace.  INPUTS: f (char).  OUTPUTS: none.'
'run(f);'
'end'
''
'function tf = inside(p, roots)'
'% INSIDE  True when path p sits under one of roots (with or without /private).'
'% INPUTS: p (char); roots {cellstr}.  OUTPUTS: tf (logical).'
'tf = false;'
'if isempty(p), return, end'
'p = regexprep(p, ''^/private'', '''');'
'for k = 1:numel(roots)'
'    r = regexprep(roots{k}, ''^/private'', '''');'
'    if startsWith(p, [r filesep]), tf = true; return, end'
'end'
'end'
''
'function r = pkgRoot(w, pkg)'
'% PKGROOT  The folder above the package folder of w (its repo root).'
'% INPUTS: w (which result); pkg (''+name'').  OUTPUTS: r (char, '''' when absent).'
'r = '''';'
'if isempty(w), return, end'
'k = strfind(w, [filesep pkg filesep]);'
'if isempty(k), return, end'
'r = fileparts(w(1:k(1)-1));'
'end'
};
fid = fopen(child, 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
end

% ------------------------------------------------------------------------
function v = getf(s, fld, d)
% GETF  Field with default.  INPUTS: s; fld; d.  OUTPUTS: v.
if isfield(s, fld) && ~isempty(s.(fld)), v = s.(fld); else, v = d; end
end
