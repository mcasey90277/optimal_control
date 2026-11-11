function ok = test_lunar_pnt_mission_study()
%% Purpose:
%
%   The study script headless, five runs, each in its own outDir and each
%   run INSIDE a function workspace (run_study below), never the base
%   workspace (the script clears its workspace). This tests the
%   INSTRUMENT on reduced problems, not the mission verdict:
%   (A) default path (south-pole GDOP series), serviceDays = 3, the
%       schedule -> polish -> exclude loop run for real: gate rows I1 I2 C1
%       P1..P6 PF PR PP PH PC PA H1 H2 G0 present, I1 I2 C1 PF H1 H2 G0
%       PASS, no R rows; the saved pole struct has consistent lengths,
%       nOn rises from 0 to 6 exactly at the arrival days, GDOP is NaN
%       wherever fewer than four are on station, nVis <= nOn; the four
%       PNGs exist (no service_map.png); the VERDICT names the reduced knob;
%       no figure left open.
%   (B) validate = false (polishes reused from A): STOPPED AT VALIDATION
%       (G0 NOT RUN), no pole series saved.
%   (C) maxRounds = 1 (polishes reused from A): STOPPED AT POLISH, the
%       table labelled as the last attempted round, the six attempted
%       polishes saved.
%   (D) tol = 1e-7: C1 FAIL and a VERDICT naming it, quickly, no polishing.
%   (E) scoreRequirements = true, chunkSel 5, serviceDays 3 (polishes
%       reused from A): G0 and G0R PASS, R1-R3 rows present, R3 finite,
%       service_map.png written.
%   (F) reuse from a file with a different catalog fingerprint, and from a
%       file without one, is refused (lunar_pnt_mission_study:reuse).
%   (G) a tampered cached entry (sc 3's tArr shifted 1e-6 ND) is IGNORED and
%       re-polished fresh; the other five are reused.
%   (I) the G0-FAIL stop: G0's own series given a slot moved 0.3 ND
%       (g0SlotOverride seam) -> STOPPED AT VALIDATION (G0 FAIL), saved.
%   (E) also checks the REUSED disclosure: VERDICT text, row tags,
%       run_meta.polishSource.
%   The script parses clean (checkcode).
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
%  M. Casey                                                   (c) 11/10/2026
%  M. Casey  rescope: pole series default, stops, reuse         (c) 11/10/2026
%  M. Casey  guarded reuse (F, G), G0-FAIL stop (I)            (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok   = true;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
run(fullfile(root, 'setup_lunar_pnt_paths.m'));
file = fullfile(root, 'lunar_pnt_mission_study.m');
base = fullfile(tempdir, 'lunar_pnt_test');
if isfolder(base), rmdir(base, 's'); end
dirs = struct('A', fullfile(base, 'A'), 'B', fullfile(base, 'B'), 'C', fullfile(base, 'C'), ...
              'D', fullfile(base, 'D'), 'E', fullfile(base, 'E'), 'F', fullfile(base, 'F'), ...
              'G', fullfile(base, 'G'), 'I', fullfile(base, 'I'));
for f = fieldnames(dirs).', mkdir(dirs.(f{1})); end
cleanup = onCleanup(@() rmdir_if(base));

msgs = checkcode(file, '-id');
bad  = msgs(arrayfun(@(m) ~isempty(regexp(m.id, '^(SYNER|PARSE)', 'once')) || ...
                     contains(m.message, 'Parse error'), msgs));
ok = check(ok, isempty(bad), sprintf('study parses clean (%d analyzer messages, %d parse errors)', ...
                                     numel(msgs), numel(bad)));

%% (A) default path, 3 service days
figsBefore = numel(findall(groot, 'Type', 'figure'));
[txt, ran, wA] = run_study(file, struct('headless', true, 'outDir', dirs.A, 'serviceDays', 3));
ok = check(ok, ran, sprintf('(A) default-path run completed (%.0f s)', wA));
if ~ran, ok = false; return; end
ok = check(ok, numel(findall(groot, 'Type', 'figure')) == figsBefore, '(A) no figure left open');
for id = {'I1', 'I2', 'C1', 'PF', 'H1', 'H2', 'G0'}
    ok = check(ok, has_gate(txt, id{1}, 'PASS'), sprintf('(A) gate %s PASS', id{1}));
end
for id = {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'PR', 'PP', 'PH', 'PC', 'PA'}
    ok = check(ok, has_gate(txt, id{1}, '(PASS|FAIL|UNRESOLVED)'), sprintf('(A) gate %s present', id{1}));
end
ok = check(ok, ~any(cellfun(@(id) has_gate(txt, id, '\S+'), {'R1', 'R2', 'R3', 'G0R'})), ...
           '(A) no R1-R3 / G0R rows in the default path');
ok = check(ok, ~isempty(regexp(txt, 'VERDICT: [^\n]*\(reduced: serviceDays 3\)', 'once')), ...
           '(A) the VERDICT names the reduced knob');
ok = check(ok, ~isempty(regexp(txt, 'round\s+1:', 'once')), '(A) the loop printed its round lines');
for f = {'schedule.png', 'assembly_3d.png', 'pole_gdop.png', 'propellant.png'}
    ok = check(ok, isfile(fullfile(dirs.A, f{1})), sprintf('(A) %s written', f{1}));
end
ok = check(ok, ~isfile(fullfile(dirs.A, 'service_map.png')), '(A) no service_map.png in the default path');
matA = fullfile(dirs.A, 'lunar_pnt_mission_study.mat');
R = load(matA);
ok = check(ok, isfield(R.run_meta, 'date') && numel(R.rounds) >= 1, ...
           sprintf('(A) run_meta and %d round(s) saved', numel(R.rounds)));
ok = check(ok, numel(R.P) == 6 && all([R.P.ok]) && isa(R.H.r, 'single'), '(A) six accepted polishes, single positions');
pl = R.pole;
n  = numel(pl.tDays);
ok = check(ok, numel(pl.gdop) == n && numel(pl.nOn) == n && numel(pl.nVis) == n && numel(pl.arrDays) == 6, ...
           sprintf('(A) pole series lengths consistent (%d samples)', n));
onExp = sum(pl.tDays(:) >= pl.arrDays(:).', 2);
ok = check(ok, pl.nOn(1) == 0 && pl.nOn(end) == 6 && all(diff(pl.nOn) >= 0) && isequal(pl.nOn(:), onExp), ...
           '(A) nOn rises 0 -> 6, stepping exactly at the arrival days');
ok = check(ok, all(isnan(pl.gdop(pl.nOn < 4))) && all(pl.nVis(:) <= pl.nOn(:)), ...
           sprintf('(A) GDOP NaN wherever < 4 on station (%d samples), nVis <= nOn', nnz(pl.nOn < 4)));
ok = check(ok, isfield(pl, 'stats') && isfinite(pl.stats.fracLt6) && isfinite(pl.stats.longestOutageH), ...
           sprintf('(A) summary statistics saved (GDOP < 6 %.3f of the service window)', pl.stats.fracLt6));
ok = check(ok, isfield(R, 'fingerprint') && numel(R.fingerprint.z8md5) == 32 && ...
           all(strcmp({R.run_meta.polishSource.origin}, 'fresh')) && ~contains(txt, 'REUSED'), ...
           '(A) fingerprint saved, all six polishes fresh, no REUSED text');

%% (B) validate = false: stop before the GDOP section
[txt, ran, wB] = run_study(file, struct('headless', true, 'outDir', dirs.B, 'serviceDays', 3, ...
                                        'validate', false, 'reuse', matA));
ok = check(ok, ran, sprintf('(B) validate = false run completed (%.0f s)', wB));
ok = check(ok, has_gate(txt, 'G0', 'NOT RUN') && contains(txt, 'VERDICT: STOPPED AT VALIDATION (G0 NOT RUN)'), ...
           '(B) STOPPED AT VALIDATION (G0 NOT RUN)');
RB = load(fullfile(dirs.B, 'lunar_pnt_mission_study.mat'), 'pole');
ok = check(ok, isempty(RB.pole) && ~contains(txt, 'SOUTH-POLE GDOP SERIES'), '(B) no pole series computed');

%% (C) maxRounds = 1: stop with the attempted polishes saved
[txt, ran, wC] = run_study(file, struct('headless', true, 'outDir', dirs.C, 'serviceDays', 3, ...
                                        'maxRounds', 1, 'reuse', matA));
ok = check(ok, ran, sprintf('(C) maxRounds = 1 run completed (%.0f s)', wC));
ok = check(ok, contains(txt, 'VERDICT: STOPPED AT POLISH') && contains(txt, 'LAST ATTEMPTED ROUND'), ...
           '(C) STOPPED AT POLISH, table labelled as the last attempted round');
RC = load(fullfile(dirs.C, 'lunar_pnt_mission_study.mat'), 'Pattempted', 'P');
ok = check(ok, numel(RC.Pattempted) == 6 && ~all(cellfun(@(p) p.ok, RC.Pattempted)) && isempty(RC.P), ...
           sprintf('(C) six attempted polishes saved (%d accepted), no realised P', ...
                   nnz(cellfun(@(p) p.ok, RC.Pattempted))));

%% (D) impossible tolerance: C1 FAIL, quick, nothing polished
[txt, ran, wD] = run_study(file, struct('headless', true, 'outDir', dirs.D, 'tol', 1e-7));
ok = check(ok, ran, sprintf('(D) tol = 1e-7 run completed without a throw (%.0f s)', wD));
ok = check(ok, has_gate(txt, 'C1', 'FAIL'), '(D) gate C1 FAIL');
ok = check(ok, ~isempty(regexp(txt, 'VERDICT: [^\n]*C1', 'once')), '(D) the VERDICT names C1');
ok = check(ok, ~contains(txt, 'polish sc') && wD < 60, sprintf('(D) nothing polished, quick (%.0f s)', wD));

%% (E) requirements path, reduced: chunk 5, 3 service days
[txt, ran, wE] = run_study(file, struct('headless', true, 'outDir', dirs.E, 'serviceDays', 3, ...
                                        'scoreRequirements', true, 'chunkSel', 5, 'reuse', matA));
ok = check(ok, ran, sprintf('(E) scoreRequirements run completed (%.0f s)', wE));
for id = {'G0', 'G0R'}
    ok = check(ok, has_gate(txt, id{1}, 'PASS'), sprintf('(E) gate %s PASS', id{1}));
end
for id = {'R1', 'R2', 'R3'}
    ok = check(ok, has_gate(txt, id{1}, '(PASS|FAIL)'), sprintf('(E) gate %s present', id{1}));
end
RE = load(fullfile(dirs.E, 'lunar_pnt_mission_study.mat'), 'Gservice');
ok = check(ok, isfinite(RE.Gservice.r3Floor), sprintf('(E) R3 floor finite (%.3f h)', RE.Gservice.r3Floor));
ok = check(ok, isfile(fullfile(dirs.E, 'service_map.png')), '(E) service_map.png written');
ok = check(ok, ~isempty(regexp(txt, 'VERDICT: PASS \(REUSED polishes for sc 1 2 3 4 5 6 from [^\n]*certified [^\n]* at [^\n]*; not re-certified in this run\)', 'once')), ...
           '(E) a clean reuse run''s VERDICT carries the REUSED disclosure');
tags = cellfun(@(id) ~isempty(regexp(txt, ['\n\s+' id '\s+\S+\s+\S+\s+REUSED'], 'once')), ...
               {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'PF', 'PR', 'PP', 'PH', 'PC', 'PA'});
ok = check(ok, all(tags), sprintf('(E) P1..P6, PF, PR, PP, PH, PC, PA rows tagged REUSED (%d of 12)', nnz(tags)));
RE = load(fullfile(dirs.E, 'lunar_pnt_mission_study.mat'), 'run_meta');
ok = check(ok, all(strcmp({RE.run_meta.polishSource.origin}, 'reused')) && ...
           all(strcmp({RE.run_meta.polishSource.file}, matA)), '(E) run_meta.polishSource names the source file');

%% (F) reuse refused: a different catalog fingerprint, and a file without one
matF = fullfile(dirs.F, 'other.mat');
copyfile(matA, matF);
fingerprint = R.fingerprint;  fingerprint.z8md5 = repmat('0', 1, 32);
save(matF, 'fingerprint', '-append');
[txt, ran, ~, idF] = run_study(file, struct('headless', true, 'outDir', dirs.F, 'serviceDays', 3, 'reuse', matF));
ok = check(ok, ~ran && strcmp(idF, 'lunar_pnt_mission_study:reuse') && contains(txt, 'z8md5'), ...
           sprintf('(F) different fingerprint refused (%s), naming z8md5', idF));
matF2 = fullfile(dirs.F, 'nofp.mat');
L = rmfield(load(matA), 'fingerprint');
save(matF2, '-struct', 'L', '-v7.3');
[~, ran, ~, idF2] = run_study(file, struct('headless', true, 'outDir', dirs.F, 'serviceDays', 3, 'reuse', matF2));
ok = check(ok, ~ran && strcmp(idF2, 'lunar_pnt_mission_study:reuse'), ...
           sprintf('(F) a file without a fingerprint refused (%s)', idF2));

%% (G) a tampered cached entry (sc 3's tArr shifted) is ignored and re-polished
matG = fullfile(dirs.G, 'tampered.mat');
copyfile(matA, matG);
polishCache = R.polishCache;
key3 = sprintf('%d_%d_%d_%d_%d', 3, R.S.assign(3), R.S.iD(3), R.S.iA(3), R.S.nRev(3));
ik = find(strcmp(polishCache.key, key3));
polishCache.P{ik}.tArr = polishCache.P{ik}.tArr + 1e-6;
save(matG, 'polishCache', '-append');
[txt, ran, wG] = run_study(file, struct('headless', true, 'outDir', dirs.G, 'serviceDays', 3, ...
                                        'makeFigures', false, 'reuse', matG));
ok = check(ok, ran && contains(txt, 'IGNORED cached sc 3') && contains(txt, 'polish sc 3 ->') && ...
           has_gate(txt, 'P3', 'PASS') && isempty(regexp(txt, '\n\s+P3\s+PASS\s+\S+\s+REUSED', 'once')), ...
           sprintf('(G) tampered entry ignored, sc 3 re-polished fresh, P3 not tagged (%.0f s)', wG));
RG = load(fullfile(dirs.G, 'lunar_pnt_mission_study.mat'), 'run_meta');
ok = check(ok, strcmp(RG.run_meta.polishSource(3).origin, 'fresh') && ...
           all(strcmp({RG.run_meta.polishSource([1 2 4 5 6]).origin}, 'reused')), '(G) polishSource: sc 3 fresh, the rest reused');

%% (I) G0 FAIL stops the study (seam: G0's own series with slot 6 moved 0.3 ND)
d6 = tulip_phasing();  d6(6) = d6(6) + 0.3;
[txt, ran, wI] = run_study(file, struct('headless', true, 'outDir', dirs.I, 'serviceDays', 3, ...
                                        'reuse', matA, 'g0SlotOverride', d6));
RI = load(fullfile(dirs.I, 'lunar_pnt_mission_study.mat'), 'pole', 'Vpole');
ok = check(ok, ran && has_gate(txt, 'G0', 'FAIL') && contains(txt, 'VERDICT: STOPPED AT VALIDATION (G0 FAIL)') && ...
           isempty(RI.pole) && strcmp(RI.Vpole.status, 'FAIL'), ...
           sprintf('(I) G0 FAIL -> STOPPED AT VALIDATION, saved, no pole series (%.0f s)', wI));

fprintf('test_lunar_pnt_mission_study: %s\n', tern(ok, 'PASS', 'FAIL'));
end

function [txt, ran, w, errId] = run_study(file, opts)
% RUN_STUDY  Run the study script with study_opts = opts inside a throwaway
%   function workspace (exec_script; never the base workspace, which the
%   script's clearvars would wipe); capture its output.
%   INPUTS: file char; opts struct.  OUTPUTS: txt char; ran logical; w [s];
%   errId char (the error identifier when the run threw, else '').
t0 = tic;
ran = true;
errId = '';
try
    txt = capture(@() exec_script(file, opts));
catch err
    txt = getReport(err, 'extended', 'hyperlinks', 'off');
    errId = err.identifier;
    ran = false;
end
w = toc(t0);
fprintf('%s\n', txt);
end

function txt = capture(fn)
% CAPTURE  Console output of fn().  INPUTS: fn function handle.  OUTPUTS: txt char.
assert(isa(fn, 'function_handle'), 'capture needs a function handle');
txt = evalc('fn()');
end

function exec_script(file, study_opts)
% EXEC_SCRIPT  The workspace the script runs in (it clears everything but
%   study_opts, which it reads).  INPUTS: file char; study_opts struct.  OUTPUTS: none.
assert(isstruct(study_opts), 'study_opts must be a struct');
run(file);
end

function tf = has_gate(txt, id, status)
% HAS_GATE  A gate-table row "id status" exists.  INPUTS: txt; id; status regex.  OUTPUTS: tf.
tf = ~isempty(regexp(txt, ['\n\s+' regexptranslate('escape', id) '\s+' status '\s'], 'once'));
end

function rmdir_if(d)
% RMDIR_IF  Remove a folder tree if it exists.  INPUTS: d char.  OUTPUTS: none.
if isfolder(d), rmdir(d, 's'); end
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
