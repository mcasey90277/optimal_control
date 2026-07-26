% TEST_ARTIFACT_PATHS  Every hard-coded .mat path in the campaign must resolve.
%
% WHY THIS EXISTS. The 2026-07-26 flatten moved two data artifacts
% (sundman_minfuel_certified.mat, minfuel_from_energy_seed.mat) into lib/ and
% relocated the results tree. Most references were updated; ONE was missed --
% run_certified_minfuel.m still loaded the seed from the campaign root, so
% reproducing the flagship 1.15x result died instantly on "Unable to find file
% or directory".
%
% Nothing caught it. Every fast test passed, because no test loads that seed,
% and the only thing that exercises run_certified_minfuel is a ~15 minute solve
% nobody runs in a test sweep. The bug was found by a human asking "did you
% actually try to reproduce the 1.15x result?"
%
% This test closes that gap WITHOUT running a solve. It reads every .m file in
% the campaign, finds each literal fullfile(...) expression that names a .mat,
% substitutes the owning file's own directory for `here`, evaluates it, and
% asserts the file exists. A moved artifact then fails in seconds, in a test,
% instead of fifteen minutes into a reproduction run.
%
% SCOPE AND LIMITS, stated plainly:
%   * only LITERAL fullfile(here, ...) expressions -- paths built from
%     variables, cfg.dirs, or concatenation are not checked here (cfg.dirs is
%     covered by test_run_gto_tulip check 4);
%   * a path inside a branch that never executes is still checked, which is the
%     conservative direction;
%   * files under results/ and data/ are SKIPPED: those are generated products,
%     legitimately absent on a clean checkout.
%
% INPUTS:  none
% OUTPUTS: none (throws naming every unresolved path)
%
% REFERENCES:
%   [1] ../README.md (the flatten this guards).
here = fileparts(mfilename('fullpath'));
root = fileparts(here);                       % GTO_tulip/direct

files = dir(fullfile(root, '**', '*.m'));
files = files(~contains({files.folder}, [filesep 'attic']));

% fullfile( ... ) expressions that mention a .mat literal
pat = 'fullfile\(\s*here\s*,[^;]*?\.mat''\s*\)';

nChecked = 0;  bad = {};
for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    txt   = fileread(fpath);

    % Strip comment-only lines so documentation examples are not "references".
    lines = regexp(txt, '\r?\n', 'split');
    lines = lines(~startsWith(strtrim(lines), '%'));
    code  = strjoin(lines, newline);

    exprs = regexp(code, pat, 'match');
    for m = 1:numel(exprs)
        e = exprs{m};
        % `here` in the source means that file's own folder.
        try
            resolved = eval(strrep(e, 'here', sprintf('''%s''', files(k).folder)));
        catch
            continue    % not a pure literal (built from variables) -- out of scope
        end
        % Generated products are legitimately absent.
        if contains(resolved, [filesep 'results' filesep]) || ...
           contains(resolved, [filesep 'data' filesep]) || ...
           contains(resolved, tempdir)
            continue
        end
        nChecked = nChecked + 1;
        if ~isfile(resolved)
            bad{end+1} = sprintf('%s\n      %s\n      -> %s', ...
                files(k).name, strtrim(e), resolved); %#ok<AGROW>
        end
    end
end

if ~isempty(bad)
    error('test_artifact_paths:unresolved', ...
        ['%d hard-coded artifact path(s) do not resolve. An artifact moved and a ' ...
         'reference was left behind:\n\n   %s\n'], numel(bad), strjoin(bad, sprintf('\n   ')));
end

fprintf('\ntest_artifact_paths: ALL PASS (%d literal .mat paths across %d files resolve)\n', ...
        nChecked, numel(files));
