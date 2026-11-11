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
%       no figure left open. HEADLINE PINS: last arrival day 18.799 (1e-3),
%       assignment [2 5 6 1 3 4], total propellant 71.261 kg (1e-3), exactly
%       three rounds, no accepted polish with familyChanged. The saved H is
%       double and passes mission_gdop's slot check again; the fingerprint
%       carries droPhase0 and the orbit-handle MD5s; run header and run_meta
%       carry the git dirty flags; no reused exclusions are named.
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
%   (G) tampered cached entries (sc 3's tArr shifted 1e-6 ND; sc 4's rv0
%       moved 1e-9 ND) are IGNORED and re-polished fresh; the other four
%       are reused.
%   (I) the G0-FAIL stop: G0's own series given a slot moved 0.3 ND
%       (g0SlotOverride seam) -> STOPPED AT VALIDATION (G0 FAIL), saved.
%   (J) gate I2 deliberate failures: i2Km = 1e-6, then i2Mps = 1e-6 (the
%       velocity part alone) -> I2 FAIL, STOPPED AT INPUTS, the knob named
%       in the VERDICT.
%   (E) also checks the REUSED disclosure: VERDICT text (incl. the three
%       reused refusals that decide the exclusions), row tags,
%       run_meta.polishSource.
%   SIDECAR: (E) writes it with exportFullStatus = false (the suite
%       certifies the six once, in test_export_mission_transfers): six
%       fallback entries with the library's alternatives fields, all marked
%       reused with A's results as their origin, six printed rows, the
%       status-run line and the path, run_meta.export written; (B), (C),
%       (I) print the skip line (stopped or unvalidated) and write none;
%       (G) gets an export error (exportFullStatus = 'bad', refused by the
%       export): EXPORT FAILED printed after the VERDICT, recorded in
%       run_meta.export, the run completes. (A) runs with exportTransfers
%       = false.
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
%  M. Casey  headline pins, I2 failures (J), rv0 entry check    (c) 11/10/2026
%  M. Casey  sidecar export: (E) written + reused, (C) skipped  (c) 11/10/2026
%  M. Casey  fix round 1: (E) fallback export, (B)/(I) skipped,
%            (G) export failure survives                      (c) 11/10/2026
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
              'G', fullfile(base, 'G'), 'I', fullfile(base, 'I'), 'J', fullfile(base, 'J'));
for f = fieldnames(dirs).', mkdir(dirs.(f{1})); end
cleanup = onCleanup(@() rmdir_if(base));

msgs = checkcode(file, '-id');
bad  = msgs(arrayfun(@(m) ~isempty(regexp(m.id, '^(SYNER|PARSE)', 'once')) || ...
                     contains(m.message, 'Parse error'), msgs));
ok = check(ok, isempty(bad), sprintf('study parses clean (%d analyzer messages, %d parse errors)', ...
                                     numel(msgs), numel(bad)));

%% (A) default path, 3 service days
figsBefore = numel(findall(groot, 'Type', 'figure'));
[txt, ran, wA] = run_study(file, struct('headless', true, 'outDir', dirs.A, 'serviceDays', 3, ...
                                        'exportTransfers', false));
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
ok = check(ok, numel(R.P) == 6 && all([R.P.ok]) && isa(R.H.r, 'double'), '(A) six accepted polishes, double positions');
% the headline numbers (README, Measured results), pinned
dayNd  = R.fingerprint.tStar/86400;
tLastA = max([R.P.tArr])*dayNd;
propA  = sum([R.P.propKg]);
ok = check(ok, abs(tLastA - 18.799) < 1e-3, sprintf('(A) last arrival day %.6f (pin 18.799)', tLastA));
ok = check(ok, isequal(R.S.assign, [2 5 6 1 3 4]), sprintf('(A) assignment [%s] (pin [2 5 6 1 3 4])', ...
                                                          sprintf('%d ', R.S.assign)));
ok = check(ok, abs(propA - 71.261) < 1e-3, sprintf('(A) total propellant %.6f kg (pin 71.261)', propA));
ok = check(ok, numel(R.rounds) == 3, sprintf('(A) %d schedule rounds (pin 3)', numel(R.rounds)));
ok = check(ok, ~any([R.P.familyChanged]), '(A) no accepted polish changed family');
% the saved history goes back into mission_gdop (its slot check is 1e-9 ND)
Ma = mission_inputs();
try
    Gr = mission_gdop(R.H, Ma, [0 R.H.tEnd], mission_grid('pole'), struct('score', false, 'quiet', true));
    reOk = numel(Gr.t) == numel(R.pole.tDays);
catch err
    reOk = false;
    fprintf('      %s\n', err.message);
end
ok = check(ok, reOk, '(A) the saved H passes mission_gdop''s slot check (usable for a re-run)');
fpNew = {'droPhase0', 'droMd5', 'tulipMd5'};
ok = check(ok, all(isfield(R.fingerprint, fpNew)) && numel(R.fingerprint.droMd5) == 32 && ...
           numel(R.fingerprint.tulipMd5) == 32, '(A) fingerprint carries droPhase0 and the orbit-handle MD5s');
ok = check(ok, islogical(R.run_meta.optimal_controlDirty) && islogical(R.run_meta.proj7Dirty) && ...
           ~isempty(regexp(txt, 'optimal_control \S+ \((clean|DIRTY: \d+ tracked path)', 'once')), ...
           sprintf('(A) run header and run_meta carry the git dirty flags (optimal_control dirty %d)', ...
                   R.run_meta.optimal_controlDirty));
ok = check(ok, ~contains(txt, 'exclusions reused') && isempty(R.run_meta.reusedRefusals), ...
           '(A) a fresh run names no reused exclusions');
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
RBm = load(fullfile(dirs.B, 'lunar_pnt_mission_study.mat'), 'run_meta');
ok = check(ok, contains(txt, 'SIDECAR: not written (the run stopped before six accepted transfers and a PASSED validation)') && ~isfile(fullfile(dirs.B, 'mission_transfers_70mN.mat')) && ...
           strcmp(RBm.run_meta.export.status, 'skipped'), '(B) no sidecar from an unvalidated run (skip line, run_meta.export skipped)');

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
ok = check(ok, contains(txt, 'SIDECAR: not written (the run stopped before six accepted transfers and a PASSED validation)') && ...
           ~isfile(fullfile(dirs.C, 'mission_transfers_70mN.mat')), '(C) sidecar skipped with a printed line, no file');

%% (D) impossible tolerance: C1 FAIL, quick, nothing polished
[txt, ran, wD] = run_study(file, struct('headless', true, 'outDir', dirs.D, 'tol', 1e-7));
ok = check(ok, ran, sprintf('(D) tol = 1e-7 run completed without a throw (%.0f s)', wD));
ok = check(ok, has_gate(txt, 'C1', 'FAIL'), '(D) gate C1 FAIL');
ok = check(ok, ~isempty(regexp(txt, 'VERDICT: [^\n]*C1', 'once')), '(D) the VERDICT names C1');
ok = check(ok, ~contains(txt, 'polish sc') && wD < 60, sprintf('(D) nothing polished, quick (%.0f s)', wD));

%% (E) requirements path, reduced: chunk 5, 3 service days
[txt, ran, wE] = run_study(file, struct('headless', true, 'outDir', dirs.E, 'serviceDays', 3, ...
                                        'scoreRequirements', true, 'chunkSel', 5, 'reuse', matA, ...
                                        'exportFullStatus', false));
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
ok = check(ok, ~isempty(regexp(txt, 'VERDICT: [^\n]*\(exclusions reused: 3 refused polish\(es\) from [^\n]*; not re-polished in this run\)', 'once')), ...
           '(E) the VERDICT names the three reused refusals (exclusions)');
tags = cellfun(@(id) ~isempty(regexp(txt, ['\n\s+' id '\s+\S+\s+\S+\s+REUSED'], 'once')), ...
               {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'PF', 'PR', 'PP', 'PH', 'PC', 'PA'});
ok = check(ok, all(tags), sprintf('(E) P1..P6, PF, PR, PP, PH, PC, PA rows tagged REUSED (%d of 12)', nnz(tags)));
RE = load(fullfile(dirs.E, 'lunar_pnt_mission_study.mat'), 'run_meta');
ok = check(ok, all(strcmp({RE.run_meta.polishSource.origin}, 'reused')) && ...
           all(strcmp({RE.run_meta.polishSource.file}, matA)), '(E) run_meta.polishSource names the source file');
% the sidecar of a REUSED run, fallback status (exportFullStatus = false: the
% suite certifies the six once, in test_export_mission_transfers)
sideE = fullfile(dirs.E, 'mission_transfers_70mN.mat');
okSide = isfile(sideE);
mt = struct('entries', struct('status', {}));
if okSide
    SE = load(sideE);
    mt = SE.mission_transfers;
    okSide = numel(mt.entries) == 6 && all([mt.mission.reused]) && ...
             all(arrayfun(@(q) strcmp(q.polishOrigin.file, matA), mt.mission)) && ...
             all(contains({mt.entries.source}, '(polish reused)')) && ...
             all(startsWith({mt.mission.statusRun}, 'fallback')) && all([mt.entries.status] == 1) && ...
             isequal(fieldnames(mt.entries), fieldnames(Ma.cat.alternatives));
end
rowsE = regexp(txt, '\n\s+\d\s+0\.\d{6}\s+0\.\d{6}\s+\d+\.\d{4}\s+[1-4]\s+lunar_pnt_mission', 'match');
REx = load(fullfile(dirs.E, 'lunar_pnt_mission_study.mat'), 'run_meta');
ok = check(ok, okSide && numel(rowsE) == 6 && contains(txt, ['wrote ' sideE]) && ...
           contains(txt, 'status run: 0 certified by the library''s certifier, 6 fallback; certifier UNFENCED') && ...
           strcmp(REx.run_meta.export.status, 'written'), ...
           sprintf('(E) sidecar of a REUSED run: 6 fallback entries marked reused, rows, status-run line, path (statuses [%s])', ...
                   sprintf('%d ', [mt.entries.status])));

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
% and sc 4's departure state moved 1e-9 ND: the rv0 entry check must catch it
key4 = sprintf('%d_%d_%d_%d_%d', 4, R.S.assign(4), R.S.iD(4), R.S.iA(4), R.S.nRev(4));
ik4 = find(strcmp(polishCache.key, key4));
polishCache.P{ik4}.rv0(1) = polishCache.P{ik4}.rv0(1) + 1e-9;
save(matG, 'polishCache', '-append');
[txt, ran, wG] = run_study(file, struct('headless', true, 'outDir', dirs.G, 'serviceDays', 3, ...
                                        'makeFigures', false, 'reuse', matG, 'exportFullStatus', 'bad'));
ok = check(ok, ran && contains(txt, 'IGNORED cached sc 3') && contains(txt, 'polish sc 3 ->') && ...
           has_gate(txt, 'P3', 'PASS') && isempty(regexp(txt, '\n\s+P3\s+PASS\s+\S+\s+REUSED', 'once')), ...
           sprintf('(G) tampered entry ignored, sc 3 re-polished fresh, P3 not tagged (%.0f s)', wG));
ok = check(ok, ~isempty(regexp(txt, 'IGNORED cached sc 4[^\n]*cached rv0 differs', 'once')) && ...
           contains(txt, 'polish sc 4 ->') && isempty(regexp(txt, '\n\s+P4\s+PASS\s+\S+\s+REUSED', 'once')), ...
           '(G) a cached rv0 off the current DRO state is ignored, sc 4 re-polished fresh');
RG = load(fullfile(dirs.G, 'lunar_pnt_mission_study.mat'), 'run_meta');
ok = check(ok, all(strcmp({RG.run_meta.polishSource([3 4]).origin}, 'fresh')) && ...
           all(strcmp({RG.run_meta.polishSource([1 2 5 6]).origin}, 'reused')), '(G) polishSource: sc 3, 4 fresh, the rest reused');
% an export error (exportFullStatus = 'bad' is refused by the export) is
% printed, recorded, and does not stop the run or lose the VERDICT
RGm = load(fullfile(dirs.G, 'lunar_pnt_mission_study.mat'), 'run_meta');
iV = strfind(txt, 'VERDICT: PASS');  iX = strfind(txt, 'EXPORT FAILED: export_mission_transfers:opts');
ok = check(ok, ran && ~isempty(iV) && ~isempty(iX) && iX(1) > iV(1) && strcmp(RGm.run_meta.export.status, 'failed') && ...
           contains(RGm.run_meta.export.reason, 'export_mission_transfers:opts') && ...
           ~isfile(fullfile(dirs.G, 'mission_transfers_70mN.mat')), ...
           '(G) EXPORT FAILED printed after the VERDICT, recorded in run_meta.export, run completed, no sidecar');

%% (I) G0 FAIL stops the study (seam: G0's own series with slot 6 moved 0.3 ND)
d6 = tulip_phasing();  d6(6) = d6(6) + 0.3;
[txt, ran, wI] = run_study(file, struct('headless', true, 'outDir', dirs.I, 'serviceDays', 3, ...
                                        'reuse', matA, 'g0SlotOverride', d6));
RI = load(fullfile(dirs.I, 'lunar_pnt_mission_study.mat'), 'pole', 'Vpole');
ok = check(ok, ran && has_gate(txt, 'G0', 'FAIL') && contains(txt, 'VERDICT: STOPPED AT VALIDATION (G0 FAIL)') && ...
           isempty(RI.pole) && strcmp(RI.Vpole.status, 'FAIL'), ...
           sprintf('(I) G0 FAIL -> STOPPED AT VALIDATION, saved, no pole series (%.0f s)', wI));
ok = check(ok, contains(txt, 'SIDECAR: not written (the run stopped before six accepted transfers and a PASSED validation)') && ~isfile(fullfile(dirs.I, 'mission_transfers_70mN.mat')), ...
           '(I) no sidecar from a run whose validation FAILED');

%% (J) gate I2 deliberate failures: a position limit, then a velocity limit, too tight
[txt, ran, wJ] = run_study(file, struct('headless', true, 'outDir', dirs.J, 'i2Km', 1e-6));
ok = check(ok, ran && has_gate(txt, 'I2', 'FAIL') && ...
           ~isempty(regexp(txt, 'VERDICT: STOPPED AT INPUTS \(I2[^\n]*\(reduced: i2Km 1e-06\)', 'once')) && ...
           ~contains(txt, '2. SCHEDULE'), sprintf('(J) i2Km 1e-6: I2 FAIL, STOPPED AT INPUTS, knob named (%.0f s)', wJ));
[txt, ran, wJ2] = run_study(file, struct('headless', true, 'outDir', dirs.J, 'i2Mps', 1e-6));
ok = check(ok, ran && has_gate(txt, 'I2', 'FAIL') && ...
           ~isempty(regexp(txt, 'VERDICT: STOPPED AT INPUTS \(I2[^\n]*\(reduced: i2Mps 1e-06\)', 'once')), ...
           sprintf('(J) i2Mps 1e-6: I2 FAIL on velocity alone, STOPPED AT INPUTS (%.0f s)', wJ2));

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
