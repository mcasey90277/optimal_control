function ok = test_pipeline()
%% Purpose:
%
%   The study scripts, headless: stability_cost_study.m runs with figures
%   invisible into a temp folder (through the base-workspace `study_opts`
%   hook the scripts honour), its captured output carries the summary
%   block with every P gate PASS and both H verdicts present, the PNG is
%   written, and the script parses clean (checkcode, errors only). A
%   script runs in the base workspace, so it is driven with evalin.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every check passed
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok   = true;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root, fullfile(root, 'lib'), fullfile(fileparts(root), 'costate_common'));

outDir = fullfile(tempdir, 'test_orbit_network_pipeline');
if ~exist(outDir, 'dir'), mkdir(outDir); end
delete(fullfile(outDir, '*.png'));

scripts = {'stability_cost_study', {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7'}, {'H1', 'H2'}, {'stability_cost.png'}
           'orbit_network_study',  {'G1', 'G2', 'G3', 'G4', 'G5'},                      {},           {'network_graph.png', 'accessibility_plane.png', 'ladder_closure.png', 'relay_closure.png'}};
for ss = 1:size(scripts, 1)
    name = scripts{ss, 1};
    file = fullfile(root, [name '.m']);
    msgs = checkcode(file, '-id');
    bad  = msgs(arrayfun(@(m) ~isempty(regexp(m.id, '^(SYNER|PARSE)', 'once')) || contains(m.message, 'Parse error'), msgs));
    ok = check(ok, isempty(bad), sprintf('%s parses clean (%d analyzer messages, %d parse errors)', name, numel(msgs), numel(bad)));

    assignin('base', 'study_opts', struct('outDir', outDir, 'visible', 'off'));
    t0 = tic;
    try
        txt = evalc(sprintf('evalin(''base'', ''run(''''%s'''')'')', file));
        ran = true;
    catch err
        txt = err.message;
        ran = false;
    end
    evalin('base', 'clear study_opts');
    ok = check(ok, ran, sprintf('%s ran without error (%.0f s)', name, toc(t0)));
    if ~ran, fprintf('      %s\n', strrep(txt, newline, [newline '      '])); continue; end
    ok = check(ok, contains(txt, '== SUMMARY'), sprintf('%s printed its summary block', name));
    gates = scripts{ss, 2};
    for gg = 1:numel(gates)
        ok = check(ok, ~isempty(regexp(txt, ['\n\s*' gates{gg} '\s+PASS'], 'once')), sprintf('%s: gate %s PASS', name, gates{gg}));
    end
    verd = scripts{ss, 3};
    for vv = 1:numel(verd)
        ok = check(ok, ~isempty(regexp(txt, ['\n\s*' verd{vv} '\s+(DETECTED|NOT DETECTED)'], 'once')), sprintf('%s: verdict %s present', name, verd{vv}));
    end
    figs = scripts{ss, 4};
    for ff = 1:numel(figs)
        ok = check(ok, exist(fullfile(outDir, figs{ff}), 'file') == 2, sprintf('%s wrote %s', name, figs{ff}));
    end
end

%% the hypothetical reverse graph, which the default run never exercises
assignin('base', 'study_opts', struct('outDir', outDir, 'visible', 'off', 'hypothetical', true));
try
    txt = evalc(sprintf('evalin(''base'', ''run(''''%s'''')'')', fullfile(root, 'orbit_network_study.m')));
    ran = true;
catch err
    txt = err.message;  ran = false;
end
evalin('base', 'clear study_opts');
ok = check(ok, ran, 'orbit_network_study ran with hypothetical = true');
ok = check(ok, ran && ~isempty(regexp(txt, '\n\s*G5\s+PASS', 'once')), 'hypothetical run: G5 PASS');
ok = check(ok, ran && contains(txt, 'HYPOTHETICAL'), 'hypothetical run labels its numbers HYPOTHETICAL');
ok = check(ok, ran && contains(txt, '112 edges'), 'hypothetical run reports 112 edges (56 + 56 reversed)');

fprintf('test_pipeline: %s\n', verdict(ok));
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
function s = verdict(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
