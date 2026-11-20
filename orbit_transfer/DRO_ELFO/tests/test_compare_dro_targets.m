function ok = test_compare_dro_targets()
%% Purpose:
%
%   Test of the user script compare_dro_targets on three tiny synthetic
%   sheets (tulip, ELFO raan 0, ELFO raan 90; 2 x 2 torus, rungs 3, 2, 1 N)
%   in tempname. Hand-computed expectations:
%     tulip   accepted [4 3 1], t_f rung 2 = 1.0
%     elfo 0  accepted [4 2 0], t_f rung 2 = 1.5   -> d share(2 N) = 0.5 - 0.75 = -0.25, t_f ratio 1.5
%     elfo 90 accepted [4 4 3], t_f 1.0 2.0 4.0   -> d share(1 N) = 0.75 - 0.25 = +0.5, t_f ratio 2
%   Checked: the saved .mat (three scorecards, combined table, differences)
%   and figure; n/a in the refuted column without a _conj file, the
%   refuted share with one (elfo 0 at 2 N: one of two accepted cells has
%   conjPass 0 -> 0.5); the deep-probe table (absent: "not run"; present:
%   only rungs below 0.5 N); the errors missingSheet, sheetMismatch and
%   knob; the override struct lasting one run; the solve = true path with
%   a stub torus function (no real solve).
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true iff every check
%                                                   passed
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here    = fileparts(mfilename('fullpath'));
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(fileparts(here));
setup_dro_elfo_paths();
addpath(fullfile(fileparts(fileparts(here)), 'costate_common'));

ok  = true;
tmp = tempname;  mkdir(tmp);
cleanTmp = onCleanup(@() rmdir_if(tmp));
rungs = [3 2 1];

% --- the three synthetic sheets -------------------------------------------
[okT, tfT] = deal(false(2, 2, 3), nan(2, 2, 3));
okT(:,:,1) = true;  okT(:,:,2) = true;  okT(2,2,2) = false;  okT(1,1,3) = true;
tfT(:,:,1) = 0.5;  tfT(:,:,2) = 1.0;  tfT(:,:,3) = 2.0;  tfT(~okT) = NaN;
okE0 = false(2, 2, 3);  okE0(:,:,1) = true;  okE0(1,:,2) = true;
tfE0 = nan(2, 2, 3);  tfE0(:,:,1) = 0.75;  tfE0(:,:,2) = 1.5;  tfE0(~okE0) = NaN;
okE9 = true(2, 2, 3);  okE9(2,2,3) = false;
tfE9 = nan(2, 2, 3);  tfE9(:,:,1) = 1.0;  tfE9(:,:,2) = 2.0;  tfE9(:,:,3) = 4.0;  tfE9(~okE9) = NaN;
kT  = struct('target', 'tulip', 'depSize', 1, 'raanDeg', NaN);
kE0 = struct('target', 'elfo',  'depSize', 1, 'raanDeg', 0);
kE9 = struct('target', 'elfo',  'depSize', 1, 'raanDeg', 90);
write_sheet(fullfile(tmp, 'dro_tulip_tau1.mat'),          rungs, okT,  tfT,  kT);
write_sheet(fullfile(tmp, 'dro_elfo_tau1_raan000.mat'),   rungs, okE0, tfE0, kE0);
write_sheet(fullfile(tmp, 'dro_elfo_tau1_raan090.mat'),   rungs, okE9, tfE9, kE9);
write_direct(fullfile(tmp, 'dro_tulip_tau1_direct.mat'), 2, 3);          % tulip sidecar: 1.5-turn arcs about the Moon

base = struct('n_phase', 2, 'rungs', rungs, 'out_dir', tmp, 'fig_visible', 'off');

% --- 1. first run: no conj files, no deep probe ---------------------------
[txt, stillThere] = run_script(base);
fprintf('%s', txt);
f_mat = fullfile(tmp, 'compare_dro_targets.mat');
f_png = fullfile(tmp, 'compare_dro_targets.png');
ok = check('saved .mat and figure exist', isfile(f_mat) && isfile(f_png)) && ok;
S = load(f_mat);
ok = check('three scorecards, tulip first', numel(S.cards) == 3 && strcmp(S.cards(1).tag, 'dro_tulip_tau1') && ...
           strcmp(S.cards(3).tag, 'dro_elfo_tau1_raan090')) && ok;
ok = check('scorecard shares tulip [1 .75 .25]', max(abs(S.cards(1).D.perRung.shareAccepted.' - [1 0.75 0.25])) < 1e-15) && ok;
d0 = S.diffs(S.diffs.target == "elfo raan 0" & S.diffs.thrustN == 2, :);
ok = check(sprintf('elfo 0 minus tulip, 2 N: d share %+.4f (-0.25), t_f ratio %.4f (1.5)', d0.dShare, d0.tfRatio), ...
           abs(d0.dShare + 0.25) < 1e-12 && abs(d0.tfRatio - 1.5) < 1e-12) && ok;
d9 = S.diffs(S.diffs.target == "elfo raan 90" & S.diffs.thrustN == 1, :);
ok = check(sprintf('elfo 90 minus tulip, 1 N: d share %+.4f (+0.5), t_f ratio %.4f (2)', d9.dShare, d9.tfRatio), ...
           abs(d9.dShare - 0.5) < 1e-12 && abs(d9.tfRatio - 2) < 1e-12) && ok;
ok = check('no conj file: refuted column is NaN for all targets and prints n/a', ...
           all(isnan(S.combined_table.tulip_refuted)) && all(isnan(S.combined_table.elfo_raan_0_refuted)) && ...
           contains(txt, 'n/a')) && ok;
ok = check('combined table: one row per rung, tulip median t_f at 2 N is 1.0', height(S.combined_table) == 3 && ...
           abs(S.combined_table.tulip_tf(2) - 1) < 1e-12) && ok;
ok = check('printed: group labels, ladder summary, differences block', contains(txt, 'elfo raan 90') && ...
           contains(txt, 'ladder summary') && contains(txt, 'differences (elfo minus tulip)')) && ok;
ok = check('ladder: tulip lowestAny 1, elfo 0 lowestAny 2 (rung 1 N has no accepted cell)', ...
           S.cards(1).D.ladder.lowestRungAnyCell == 1 && S.cards(2).D.ladder.lowestRungAnyCell == 2) && ok;
ok = check('deep probe absent: one line "deep probe: not run"', contains(txt, 'deep probe: not run') && isempty(S.deep)) && ok;
ok = check('override struct cleared by the script', ~stillThere) && ok;
ok = check(sprintf('rev column with real arcs: tulip median revolutions %.4f (1.5) at every rung, others n/a-dash', ...
           S.combined_table.tulip_rev(1)), max(abs(S.combined_table.tulip_rev - 1.5)) < 1e-3 && ...
           all(isnan(S.combined_table.elfo_raan_0_rev))) && ok;
lines = splitlines(txt);
tab = lines(contains(lines, ' |'));
ok = check(sprintf('table lines under 120 characters (longest %d), bars aligned', max(cellfun(@numel, tab))), ...
           max(cellfun(@numel, tab)) < 120 && numel(unique(cellfun(@numel, tab))) == 1) && ok;
ok = check('legend line printed', contains(txt, 'acc share of cells accepted') && contains(txt, 'refut refuted share')) && ok;
ok = check('complete sheets: "4 of 4 cells attempted", no INCOMPLETE, saved complete = [1 1 1]', ...
           contains(txt, 'tulip: 4 of 4 cells attempted') && ~contains(txt, 'INCOMPLETE') && ...
           isequal(S.complete, [true true true]) && isequal(S.attempted, [4 4 4]) && isequal(S.total, [4 4 4])) && ok;
ok = check('ladder columns renamed', contains(txt, 'cells attempted') && contains(txt, 'median attempts per attempted cell')) && ok;

% --- 2. conj file and a deep-probe sheet -----------------------------------
cp = ones(2, 2, 3);  cp(1,2,2) = 0;                                      % elfo 0, 2 N: accepted cells (1,1),(1,2); (1,2) refuted
conjPass = cp;
save(fullfile(tmp, 'dro_elfo_tau1_raan000_conj.mat'), 'conjPass');
mkdir(fullfile(tmp, 'deep_probe'));
okD = true(1, 1, 3);  tfD = nan(1, 1, 3);  tfD(:) = [3 5 8];
write_sheet(fullfile(tmp, 'deep_probe', 'dro_tulip_tau1.mat'), [0.5 0.3 0.2], okD, tfD, kT);
[txt, ~] = run_script(base);
S = load(f_mat);
ok = check(sprintf('with conj: elfo 0 refuted share at 2 N = %.3f (0.5), elfo 0 at 3 N = %.3f (0)', ...
           S.combined_table.elfo_raan_0_refuted(2), S.combined_table.elfo_raan_0_refuted(1)), ...
           abs(S.combined_table.elfo_raan_0_refuted(2) - 0.5) < 1e-12 && S.combined_table.elfo_raan_0_refuted(1) == 0) && ok;
ok = check('with conj: tulip (no file) is still n/a', all(isnan(S.combined_table.tulip_refuted))) && ok;
ok = check('deep probe: table "below 0.5 N" printed, 0.5 N rung excluded, median t_f of 0.3 N is 5', ...
           contains(txt, 'below 0.5 N') && ~contains(txt, 'deep probe: not run') && numel(S.deep) == 1 && ...
           isequal(S.deep(1).thrustN(:).', [0.3 0.2]) && abs(S.deep(1).perRung.tfMedianND(1) - 5) < 1e-12) && ok;

% --- 3. override lasts one run: n_phase and rungs revert to the knobs ------
ok = check_error('second run without the first run''s overrides reads the knobs (6 x 6, 11 rungs)', ...
                 'compare_dro_targets:sheetMismatch', @() run_script(struct('out_dir', tmp, 'fig_visible', 'off'))) && ok;

% --- 4. errors --------------------------------------------------------------
tmp2 = fullfile(tmp, 'partial');  mkdir(tmp2);
copyfile(fullfile(tmp, 'dro_tulip_tau1.mat'), tmp2);
copyfile(fullfile(tmp, 'dro_elfo_tau1_raan000.mat'), tmp2);
o2 = base;  o2.out_dir = tmp2;
ok = check_error('missing sheet, solve = false', 'compare_dro_targets:missingSheet', @() run_script(o2)) && ok;
try
    run_script(o2);
catch err
    ok = check('error names the file and the command', contains(err.message, 'dro_elfo_tau1_raan090.mat') && ...
               contains(err.message, 'run_dro_elfo_torus')) && ok;
end
ok = check('nothing written by the failed run', ~isfile(fullfile(tmp2, 'compare_dro_targets.mat'))) && ok;
o3 = base;  o3.rungs = [3 2];
ok = check_error('sheet with other rungs than the knobs', 'compare_dro_targets:sheetMismatch', @() run_script(o3)) && ok;
o4 = base;  o4.nonsense = 1;
ok = check_error('unknown override field', 'compare_dro_targets:knob', @() run_script(o4)) && ok;

% --- 4b. incomplete sheet -----------------------------------------------------
tmp3 = fullfile(tmp, 'incomplete');  mkdir(tmp3);
copyfile(fullfile(tmp, 'dro_tulip_tau1.mat'), tmp3);
copyfile(fullfile(tmp, 'dro_elfo_tau1_raan000.mat'), tmp3);
write_sheet(fullfile(tmp3, 'dro_elfo_tau1_raan090.mat'), rungs, okE9, tfE9, kE9, [1 1; 0 0]);   % 2 of 4 attempted
o6 = base;  o6.out_dir = tmp3;
lastwarn('');
txt = run_script(o6);
[~, wid] = lastwarn;
S6 = load(fullfile(tmp3, 'compare_dro_targets.mat'));
ok = check('incomplete sheet: line "elfo raan 90: 2 of 4 cells attempted INCOMPLETE" and the warning id', ...
           contains(txt, 'elfo raan 90: 2 of 4 cells attempted INCOMPLETE') && strcmp(wid, 'compare_dro_targets:incomplete')) && ok;
ok = check('incomplete: notice above the combined table and above the differences block', ...
           numel(strfind(txt, 'INCOMPLETE: elfo raan 90 (2 of 4 cells attempted)')) >= 2) && ok;
ok = check('incomplete: saved attempted [4 4 2], total [4 4 4], complete [1 1 0], ladder shows "2 of 4"', ...
           isequal(S6.attempted, [4 4 2]) && isequal(S6.total, [4 4 4]) && isequal(S6.complete, [true true false]) && ...
           contains(txt, '2 of 4')) && ok;

% --- 4c. sheet identity ------------------------------------------------------
tmp4 = fullfile(tmp, 'identity');  mkdir(tmp4);
copyfile(fullfile(tmp, 'dro_tulip_tau1.mat'), tmp4);
copyfile(fullfile(tmp, 'dro_tulip_tau1.mat'), fullfile(tmp4, 'dro_elfo_tau1_raan000.mat'));   % tulip file under an ELFO tag
copyfile(fullfile(tmp, 'dro_elfo_tau1_raan090.mat'), tmp4);
o7 = base;  o7.out_dir = tmp4;
ok = check_error('tulip file copied under an ELFO tag is refused', 'compare_dro_targets:sheetMismatch', @() run_script(o7)) && ok;
try
    run_script(o7);
catch err
    ok = check('message names the tag and the stored key', contains(err.message, 'dro_elfo_tau1_raan000') && ...
               contains(err.message, 'target tulip') && contains(err.message, 'expected target elfo')) && ok;
end
copyfile(fullfile(tmp, 'dro_elfo_tau1_raan090.mat'), fullfile(tmp4, 'dro_elfo_tau1_raan000.mat'));   % raan 90 file under raan 0
ok = check_error('raan 90 file under the raan 0 tag is refused', 'compare_dro_targets:sheetMismatch', @() run_script(o7)) && ok;
copyfile(fullfile(tmp, 'dro_elfo_tau1_raan000.mat'), fullfile(tmp4, 'dro_elfo_tau1_raan000.mat'));
write_sheet(fullfile(tmp4, 'dro_elfo_tau1_raan090.mat'), rungs, okE9, tfE9, []);                 % no sheet_key at all
try
    run_script(o7);
    ok = check('sheet without sheet_key is refused', false) && ok;
catch err
    ok = check('sheet without sheet_key is refused, message says so', strcmp(err.identifier, 'compare_dro_targets:sheetMismatch') ...
               && contains(err.message, 'no sheet_key')) && ok;
end

% --- 4d. a refused field does not poison the next run ---------------------
[msg, left] = run_bad(struct('out_dir', tmp, 'nonsense', 1));
ok = check(sprintf('unknown field: error knob, compare_opts already cleared (%s)', msg), strcmp(msg, 'compare_dro_targets:knob') && ~left) && ok;

% --- 5. solve = true through a stub torus function -------------------------
o5 = base;  o5.out_dir = tmp2;  o5.solve = true;  o5.torus_fcn = @stub_torus;
txt = run_script(o5);
sc = load(fullfile(tmp2, 'stub_called.mat'));
ok = check(sprintf('stub got only the missing sheet (elfo %g), grid %d x %d, rungs %s, outDir', ...
           sc.o.sheets(1).raanDeg, sc.o.nD, sc.o.nA, mat2str(sc.o.rungs)), numel(sc.o.sheets) == 1 && ...
           strcmp(sc.o.sheets(1).target, 'elfo') && sc.o.sheets(1).raanDeg == 90 && sc.o.nD == 2 && ...
           sc.o.nA == 2 && isequal(sc.o.rungs, rungs) && strcmp(sc.o.outDir, tmp2) && isempty(sc.sel)) && ok;
S5 = load(fullfile(tmp2, 'compare_dro_targets.mat'));
ok = check('after the stub run: three scorecards, solving line printed', numel(S5.cards) == 3 && ...
           contains(txt, 'solving 1 missing sheet(s): dro_elfo_tau1_raan090')) && ok;

fprintf('\ntest_compare_dro_targets: %s\n', ternary(ok, 'ALL PASS', 'FAILED'));
if ~ok, error('test_compare_dro_targets:fail', 'a check failed'); end
end

% ------------------------------------------------------------------------
function [txt, stillThere] = run_script(compare_opts)
% RUN_SCRIPT  run the script with the override struct, capture its output.
% INPUTS: compare_opts struct. OUTPUTS: txt char (printed text); stillThere logical (struct still in workspace).
assert(isstruct(compare_opts), 'overrides must be a struct');   % the script reads it from this workspace
txt = evalc('compare_dro_targets');
stillThere = exist('compare_opts', 'var') == 1;
end

function [id, left] = run_bad(compare_opts)
% RUN_BAD  run the script expecting an error; report its id and whether the override struct is still in the workspace.
% INPUTS: compare_opts struct. OUTPUTS: id char; left logical.
id = '';
try
    evalc('compare_dro_targets');
catch err
    id = err.identifier;
end
left = exist('compare_opts', 'var') == 1;
end

function stub_torus(sel, o)
% STUB_TORUS  stands in for run_dro_elfo_torus: records its arguments, writes the requested sheets.
% INPUTS: sel as passed; o struct. OUTPUTS: none (files in o.outDir).
save(fullfile(o.outDir, 'stub_called.mat'), 'sel', 'o');
for ks = 1:numel(o.sheets)
    write_sheet(fullfile(o.outDir, sprintf('dro_elfo_tau%g_raan%03d.mat', o.sheets(ks).droTau, o.sheets(ks).raanDeg)), ...
                o.rungs, true(o.nD, o.nA, numel(o.rungs)), ones(o.nD, o.nA, numel(o.rungs)), ...
                struct('target', 'elfo', 'depSize', o.sheets(ks).droTau, 'raanDeg', o.sheets(ks).raanDeg));
end
end

function write_sheet(f, rungs, OK, TF, key, ATT)
% WRITE_SHEET  Minimal sheet in the engine's format, with the identity variable sheet_key the front door writes.
% INPUTS: f char; rungs [1 x r]; OK logical [nD x nA x r]; TF [nD x nA x r]; key struct sheet_key ([] = none);
%   ATT [nD x nA] attempt counters (optional, default ones).
szD = size(OK, 1);  szA = size(OK, 2);  nR = numel(rungs);
OK = reshape(OK, szD, szA, nR);  TF = reshape(TF, szD, szA, nR);
WALL = 10*ones(szD, szA, nR);
if nargin < 6, ATT = ones(szD, szA); end
if nargin < 5, key = []; end
FLYKM = nan(szD, szA, nR);  ACCDZ = FLYKM;  RES = FLYKM;  Z8 = nan(8, szD, szA, nR);
sD = (0:szD-1)/szD;  sA = (0:szA-1)/szA;
meta = struct('muStar', 0.012150585609624, 'lStar', 389703.264829278, 'tStar', 382981.289129055, ...
              'depFamily', 'dro', 'arrFamily', 'elfo', 'depParams', struct(), 'arrParams', struct(), ...
              'ispS', 1710, 'm0kg', 150, 'N', 400);
save(f, 'rungs', 'sD', 'sA', 'TF', 'FLYKM', 'ACCDZ', 'RES', 'WALL', 'OK', 'Z8', 'ATT', 'meta');
if ~isempty(key)
    sheet_key = key;
    save(f, 'sheet_key', '-append');
end
end

function write_direct(f, n, nR)
% WRITE_DIRECT  Sidecar whose every cell holds a 1.5-turn circular arc about the Moon (inertial frame), written in
%   the rotating frame as the engine stores it, flight time 3 ND.  INPUTS: f char; n grid size; nR rungs.
muStar = 0.012150585609624;  tfD = 3;  nt = 401;  t = linspace(0, tfD, nt);
ang = 2*pi*1.5*t/tfD;
rot = 0.1*[cos(ang); sin(ang); zeros(1, nt)];
for kk = 1:nt
    cc = cos(t(kk));  ss = sin(t(kk));
    rot(:,kk) = [cc ss 0; -ss cc 0; 0 0 1]*rot(:,kk) + [1 - muStar; 0; 0];
end
DX = zeros(7, nt, n, n, nR);
for ka = 1:n
    for kb = 1:n
        for kr = 1:nR
            DX(1:3,:,ka,kb,kr) = rot;
        end
    end
end
DU = zeros(3, nt, n, n, nR);  DTF = tfD*ones(n, n, nR);  DDONE = true(n, n, nR);
save(f, 'DX', 'DU', 'DTF', 'DDONE');
end

function r = check(msg, cond)
% CHECK  print PASS / FAIL for one condition.  INPUTS: msg char, cond logical.  OUTPUTS: r logical.
r = logical(cond);
fprintf('  [%s] %s\n', ternary(r, 'PASS', 'FAIL'), msg);
end

function r = check_error(msg, id, fh)
% CHECK_ERROR  true iff fh() throws the identifier id.  INPUTS: msg char; id char; fh handle.  OUTPUTS: r logical.
r = false;
try
    fh();
    fprintf('  [FAIL] %s (no error)\n', msg);
catch err
    r = strcmp(err.identifier, id);
    fprintf('  [%s] %s (got %s)\n', ternary(r, 'PASS', 'FAIL'), msg, err.identifier);
end
end

function s = ternary(c, a, b)
% TERNARY  c ? a : b.  INPUTS: c logical; a, b any.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

function rmdir_if(d)
% RMDIR_IF  remove a folder when present.  INPUTS: d char.
if isfolder(d), rmdir(d, 's'); end
end
