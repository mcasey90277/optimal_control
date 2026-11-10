function ok = test_gto_dro_study()
%% Purpose:
%
%   The study script headless, three ways. This tests the INSTRUMENT, not
%   the physics verdict of one cell: (1) a full run at thrustN = 12 through
%   study_opts (the printed ladder must end at 12 N and the summary line
%   must say 12 N), strict self-check OFF; every V/N/X gate row must be
%   PASS or NOT APPLICABLE, every S gate row must be present with a status
%   in {PASS, FAIL, UNRESOLVED, NOT APPLICABLE}, the output must carry a
%   VERDICT line, and three PNGs plus the .mat must be written; (2) a REUSE
%   run on that .mat with tol.earthKm = 1e6 must print N7a FAIL and throw
%   gto_dro_study:brokenRoot (the gate can fail, and the script fails when
%   it should) and must save to gto_dro_study_reuse.mat, leaving the
%   original .mat in place with a finite run_meta.solveWallS; (3) a REUSE
%   run with a different DRO period must be refused (gto_dro_study:reuse).
%   Figures created by each run are closed. The script parses clean (checkcode, parse errors
%   only). Scripts run in the base workspace, so they are driven with evalin.
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
%  M. Casey                                                   (c) 11/09/2026
%  M. Casey  reuse provenance checks; close figures           (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok   = true;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
run(fullfile(root, 'setup_gto_dro_paths.m'));
file   = fullfile(root, 'gto_dro_study.m');
outDir = fullfile(tempdir, 'test_gto_dro_study');
if ~isfolder(outDir), mkdir(outDir); end
for f = {'transfer_3d.png', 'timeseries.png', 'ladder.png', 'gto_dro_study.mat', 'gto_dro_study_reuse.mat'}
    if isfile(fullfile(outDir, f{1})), delete(fullfile(outDir, f{1})); end
end

msgs = checkcode(file, '-id');
bad  = msgs(arrayfun(@(m) ~isempty(regexp(m.id, '^(SYNER|PARSE)', 'once')) || contains(m.message, 'Parse error'), msgs));
ok = check(ok, isempty(bad), sprintf('gto_dro_study parses clean (%d analyzer messages, %d parse errors)', numel(msgs), numel(bad)));

% (1) full run at 12 N, strict off
assignin('base', 'study_opts', struct('thrustN', 12, 'out', struct('dir', outDir, 'visible', 'off'), 'selfCheck', struct('strict', false)));
t0 = tic;  ran = true;
try
    txt = evalc(sprintf('evalin(''base'', ''run(''''%s'''')'')', file));
catch err
    txt = err.message;  ran = false;
end
evalin('base', 'clear study_opts');
close all;                                                   % figures are invisible; drop them
ok = check(ok, ran, sprintf('full run at 12 N completed (%.0f s)', toc(t0)));
if ~ran, fprintf('      %s\n', strrep(txt, newline, [newline '      '])); ok = false; return; end
fprintf('%s\n', txt);
ok = check(ok, ~isempty(regexp(txt, '== SUMMARY.*12 N', 'once')), 'summary line reports 12 N (study_opts reached the solve)');
ok = check(ok, ~isempty(regexp(txt, '\n\s+12\s+[\d.]+\s+[\d.]+.*\n\s+top-rung', 'once')), 'ladder table''s last row is the 12 N rung');
for id = {'V1 GTO geometry', 'V2 DRO closure', 'V3 ratio screen', 'N0 ladder', 'N1 BVP residual', 'N2 Hamiltonian', ...
          'N3 arrival', 'N4 transversality', 'N5 adjoint eqs', 'N6 min principle', 'N7a Earth floor', 'N7b Moon floor', ...
          'X1 tfMin acceptance', 'X2 field agreement'}
    pat = ['\n\s+' regexptranslate('escape', id{1}) '\s+(PASS|NOT APPLICABLE)'];
    ok = check(ok, ~isempty(regexp(txt, pat, 'once')), sprintf('gate %s is PASS or NOT APPLICABLE', id{1}));
end
for id = {'S1 strong Legendre', 'S2 Q_mt > 0', 'S3 dim S = 1', 'S4 conjugate'}
    pat = ['\n\s+' regexptranslate('escape', id{1}) '\s+(PASS|FAIL|UNRESOLVED|NOT APPLICABLE)'];
    ok = check(ok, ~isempty(regexp(txt, pat, 'once')), sprintf('gate %s is present with a status', id{1}));
end
ok = check(ok, contains(txt, 'VERDICT:'), 'output carries a VERDICT line');
for f = {'transfer_3d.png', 'timeseries.png', 'ladder.png', 'gto_dro_study.mat'}
    ok = check(ok, isfile(fullfile(outDir, f{1})), sprintf('%s written', f{1}));
end

% (2) reuse with an impossible Earth floor: N7a must FAIL and the script must throw brokenRoot
matFile = fullfile(outDir, 'gto_dro_study.mat');
assignin('base', 'study_opts', struct('thrustN', 12, 'out', struct('dir', outDir, 'visible', 'off', 'reuse', matFile), ...
                                       'tol', struct('earthKm', 1e6), 'selfCheck', struct('strict', false)));
threw = '';  txt2 = '';
try
    txt2 = evalc(sprintf('evalin(''base'', ''run(''''%s'''')'')', file));
catch err
    threw = err.identifier;  txt2 = err.message;
end
evalin('base', 'clear study_opts');
close all;
ok = check(ok, strcmp(threw, 'gto_dro_study:brokenRoot'), sprintf('impossible Earth floor throws gto_dro_study:brokenRoot (got ''%s'')', threw));
ok = check(ok, contains(txt2, 'N7a'), 'the failure names N7a');
ok = check(ok, isfile(fullfile(outDir, 'gto_dro_study_reuse.mat')), 'the reuse run saved to gto_dro_study_reuse.mat');

% (3) reuse with a different DRO period is refused
assignin('base', 'study_opts', struct('thrustN', 12, 'dro', struct('tau', 2), 'out', struct('dir', outDir, 'visible', 'off', 'reuse', matFile)));
threw = '';
try
    evalc(sprintf('evalin(''base'', ''run(''''%s'''')'')', file));
catch err
    threw = err.identifier;
end
evalin('base', 'clear study_opts');
close all;
ok = check(ok, strcmp(threw, 'gto_dro_study:reuse'), sprintf('reuse with tau = 2 refused with gto_dro_study:reuse (got ''%s'')', threw));

% the original (solved) .mat survives both reuse runs with its solve provenance
R = load(matFile, 'run_meta');
ok = check(ok, isfinite(R.run_meta.solveWallS), sprintf('original .mat unchanged after reuse (run_meta.solveWallS = %.0f s)', R.run_meta.solveWallS));

fprintf('test_gto_dro_study: %s\n', tern(ok, 'PASS', 'FAIL'));
end

function ok = check(ok, c, msg)
% CHECK  Print one line and fold the verdict.  INPUTS: ok; c; msg.  OUTPUTS: ok.
if c, s = 'ok  '; else, s = 'FAIL'; end
fprintf('  %s %s\n', s, msg);
ok = ok && c;
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end
