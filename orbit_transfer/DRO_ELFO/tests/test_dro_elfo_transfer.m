function ok = test_dro_elfo_transfer()
%% Purpose:
%
%   Test of the user script dro_elfo_transfer (run through the override
%   struct dro_elfo_opts, from inside a helper so the script's variables
%   stay out of this workspace). Checks:
%     1  15 N run: T.gates.pass, the .mat (with T and the knobs) and the
%        captured-entries file exist, the captured file holds one entry
%     2  a second run appends: two entries
%     3  the printed verdict table has a PASS row per gate and the
%        t_f / mass / delta-V line
%     4  max_iter = 1: prints the reason (named, no error stack), creates
%        no figure and writes no files
%     5  an unknown override field (also `dep`) and a malformed elfo knob
%        refuse by name; an override lasts ONE run (the second plain run
%        uses the written knobs); entries append onto an older file whose
%        stored entry has an extra field (cell array)
%     6  the ELFO knob semantics: NaN raan_deg / sA with IM satellite 3
%        resolve to its plane RAAN and M0 / 360 (planOnly through the
%        solver seam, no solve)
%   The solves go through a cache in tempdir/dro_elfo_transfer_test_cache
%   (the seam knob `solver`), keyed by the call arguments and the source
%   text of pair_transfer, located_min_altitude, thrust_ladder_library,
%   fly_transfer, ms_tfmin and get_family_orbit: cold about 3.5 min, run it in the background.
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
dre     = fileparts(here);
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(dre);
setup_dro_elfo_paths();

ok = true;
tmp = tempname;
cleanTmp = onCleanup(@() rmdir_if(tmp));
cacheDir = fullfile(tempdir, 'dro_elfo_transfer_test_cache');
if ~isfolder(cacheDir), mkdir(cacheDir); end
solver = @(varargin) cached_pair(cacheDir, varargin{:});

o1 = struct('thrust_N', 15, 'make_plot', false, 'out_dir', tmp, 'solver', solver);
[T, txt] = run_script(o1);
fprintf('%s', txt);
ok = check(sprintf('15 N run: ok %d, gates.pass %d, t_f %.6f ND', T.ok, T.gates.pass, T.tf_nd), T.ok && T.gates.pass) && ok;
ok = check('results .mat exists with T and knobs', isfile(fullfile(tmp, 'dro_elfo_transfer.mat')) && ...
           all(isfield(load(fullfile(tmp, 'dro_elfo_transfer.mat')), {'T', 'knobs'}))) && ok;
capf = fullfile(tmp, 'dro_elfo_captured_entries.mat');
S = load(capf);
ok = check(sprintf('captured file holds %d entry (cell array)', numel(S.entries)), iscell(S.entries) && numel(S.entries) == 1 && ...
           abs(S.entries{1}.tf_nd - T.tf_nd) < 1e-12) && ok;
ok = check('verdict table: PASS rows for flown, accDz, admissible, moon and the summary line', ...
           contains(txt, 'flown') && contains(txt, 'moon') && contains(txt, 'PASS') && ~contains(txt, 'FAIL') && ...
           contains(txt, 'delta-V')) && ok;
[~, ~] = run_script(o1);
S = load(capf);
ok = check(sprintf('second run appends: %d entries', numel(S.entries)), numel(S.entries) == 2) && ok;

%% failure path
tmp2 = tempname;
cleanTmp2 = onCleanup(@() rmdir_if(tmp2));
figsBefore = findall(groot, 'Type', 'figure');
o2 = struct('thrust_N', 15, 'max_iter', 1, 'make_plot', true, 'out_dir', tmp2);
try
    [T2, txt2] = run_script(o2);
    figsAfter = findall(groot, 'Type', 'figure');
    ok = check('max_iter 1: not ok, reason printed, no error', ~T2.ok && contains(txt2, 'NOT SOLVED') && ...
               contains(txt2, T2.reason) && ~isempty(T2.reason)) && ok;
    fprintf('    printed: %s', strtrim(txt2));
    fprintf('\n');
    ok = check('max_iter 1: no figure, no files, out_dir not even created', numel(figsAfter) == numel(figsBefore) && ...
               ~isfile(fullfile(tmp2, 'dro_elfo_transfer.mat')) && ...
               ~isfile(fullfile(tmp2, 'dro_elfo_captured_entries.mat')) && ~isfolder(tmp2)) && ok;
catch err
    ok = check(['failure path threw: ' err.message], false) && ok;
end

%% unknown knob
ok = check_error('unknown override field refuses', 'dro_elfo_transfer:knob', ...
                 @() run_script(struct('thrust_Nn', 5, 'out_dir', tmp2, 'solver', solver))) && ok;

%% overrides last one run; dep is refused; captured entries tolerate dissimilar entries
rec0 = @(dep, tgt, sD, sA, thr, opts) struct('ok', false, 'reason', 'seam', 'thr', thr, 'dep', dep);
[T6, ex1] = run_twice(struct('thrust_N', 15, 'solver', rec0), struct('solver', rec0));
ok = check(sprintf('override read once: run 1 thrust %g N, then variable cleared (%d), plain run 2 thrust %g N', ...
           T6.first.thr, ex1, T6.second.thr), T6.first.thr == 15 && ~ex1 && T6.second.thr == 5) && ok;
ok = check_error('dro_elfo_opts.dep is refused', 'dro_elfo_transfer:knob', ...
                 @() run_script(struct('dep', struct('family', 'dro'), 'out_dir', tmp2))) && ok;
ok = check_error('malformed elfo (no .source) refused', 'dro_elfo_transfer:elfo', ...
                 @() run_script(struct('elfo', struct('sat', 1), 'out_dir', tmp2, 'solver', rec0))) && ok;
ok = check_error('unknown elfo source refused', 'dro_elfo_transfer:elfo', ...
                 @() run_script(struct('elfo', struct('source', 'nope'), 'out_dir', tmp2, 'solver', rec0))) && ok;
ok = check_error('elfo.sat out of range refused', 'dro_elfo_transfer:elfo', ...
                 @() run_script(struct('elfo', struct('source', 'im_elfo_optimum', 'sat', 9), 'out_dir', tmp2, 'solver', rec0))) && ok;
tmp3 = tempname;  mkdir(tmp3);
cleanTmp3 = onCleanup(@() rmdir_if(tmp3));
old = T.entry;  old.extraField = 'older file';                       % a stored entry with a field the new one lacks
S0 = struct('entries', old);
save(fullfile(tmp3, 'dro_elfo_captured_entries.mat'), '-struct', 'S0');
[~, ~] = run_script(struct('thrust_N', 15, 'make_plot', false, 'out_dir', tmp3, 'solver', solver));
S = load(fullfile(tmp3, 'dro_elfo_captured_entries.mat'));
ok = check(sprintf('append onto an older struct-array file: %d cells, first has extraField %d, second %d', numel(S.entries), ...
           isfield(S.entries{1}, 'extraField'), isfield(S.entries{2}, 'extraField')), iscell(S.entries) && ...
           numel(S.entries) == 2 && isfield(S.entries{1}, 'extraField') && ~isfield(S.entries{2}, 'extraField')) && ok;

%% ELFO knob semantics (no solve: the solver seam records its arguments)
[raan_opt, m0_opt, plane_of, oe] = im_elfo_optimum();
rec = @(dep, tgt, sD, sA, thr, opts) struct('ok', false, 'reason', 'seam', 'tgt', tgt, 'sA', sA, 'sD', sD, 'dep', dep);
[T3, ~] = run_script(struct('elfo', struct('source', 'im_elfo_optimum', 'sat', 3), 'out_dir', tmp2, 'solver', rec));
ok = check(sprintf('IM sat 3: raanDeg %.1f (plane %d), sA %.4f, oe %s', T3.tgt.params.raanDeg, plane_of(3), T3.sA, ...
           mat2str(T3.tgt.params.oe)), T3.tgt.params.raanDeg == raan_opt(plane_of(3)) && ...
           abs(T3.sA - m0_opt(3) / 360) < 1e-15 && isequal(T3.tgt.params.oe, oe)) && ok;
[T4, ~] = run_script(struct('elfo', struct('source', 'elements', 'oe', [9000 0.5 60 90]), 'raan_deg', 400, ...
                             'out_dir', tmp2, 'solver', rec));
ok = check(sprintf('raw elements: sA %.2f, raan_deg as given %.0f (wrapped by the locus)', T4.sA, T4.tgt.params.raanDeg), ...
           T4.sA == 0 && T4.tgt.params.raanDeg == 400) && ok;
[T5, ~] = run_script(struct('target', struct('family', 'tulip', 'params', struct('Np', 7, 'pm', -1)), ...
                             'sA', 0.1, 'out_dir', tmp2, 'solver', rec));
ok = check('target override: the tulip struct reaches the solver unchanged', strcmp(T5.tgt.family, 'tulip') && T5.sA == 0.1) && ok;
fprintf('test_dro_elfo_transfer: %s\n', ternary(ok, 'ALL PASS', 'FAILURES'));
end

% ------------------------------------------------------------------------
function [Tt, ex] = run_twice(o1, o2)
% RUN_TWICE  run the script with override o1, then again with override o2 in the same workspace.
%   INPUTS: o1, o2 structs.  OUTPUTS: Tt (.first, .second the solver results); ex true if dro_elfo_opts still
%   exists in the workspace after the first run.
T = [];
dro_elfo_opts = o1;
evalc('dro_elfo_transfer');
Tt.first = T;
ex = exist('dro_elfo_opts', 'var') == 1;
dro_elfo_opts = o2;
evalc('dro_elfo_transfer');
Tt.second = T;
end

function [T, txt] = run_script(dro_elfo_opts)
% RUN_SCRIPT  run the user script with the override struct, capture its
%   printed output.  INPUTS: dro_elfo_opts struct (overrides).  OUTPUTS: T
%   the solver result the script holds; txt char the printed text.
assert(isstruct(dro_elfo_opts), 'overrides must be a struct');   % the script reads it from this workspace
txt = evalc('dro_elfo_transfer');
end

function T = cached_pair(cacheDir, dep, arr, sD, sA, thrustN, opts)
% CACHED_PAIR  pair_transfer through a tempdir cache keyed by every
%   argument and the source text of pair_transfer.  INPUTS: cacheDir
%   char; the pair_transfer arguments.  OUTPUTS: T struct.
src = source_hash({'pair_transfer', 'located_min_altitude', 'thrust_ladder_library', 'fly_transfer', ...
                   'ms_tfmin', 'get_family_orbit'});
key = sha_hex([jsonencode(struct('dep', dep, 'arr', arr, 'sD', sD, 'sA', sA, 'thrustN', thrustN, ...
                                 'opts', opts, 'v', 'v2')) src]);
f = fullfile(cacheDir, ['dt_' key '.mat']);
if isfile(f)
    S = load(f, 'T');  T = S.T;
    fprintf('  (cached: %s)\n', f);
else
    opts.quiet = true;
    T = pair_transfer(dep, arr, sD, sA, thrustN, opts);
    save(f, 'T');
end
end

function h = source_hash(names)
% SOURCE_HASH  SHA-256 over the source text of the named functions.  INPUTS: names cellstr.
%   OUTPUTS: h char.
parts = cell(1, numel(names));
for k = 1:numel(names)
    f = which(names{k});
    if isempty(f), error('test_dro_elfo_transfer:source', 'cannot resolve %s for the cache key', names{k}); end
    parts{k} = [names{k}, ':', sha_hex(fileread(f))];
end
h = sha_hex(strjoin(parts, ';'));
end

function s = sha_hex(str)
% SHA_HEX  SHA-256 of a char array as lower-case hex.  INPUTS: str char.
%   OUTPUTS: s char [1 x 64].
md = java.security.MessageDigest.getInstance('SHA-256');
b  = typecast(md.digest(unicode2native(str, 'UTF-8')), 'uint8');
s  = lower(reshape(dec2hex(b).', 1, []));
end

function rmdir_if(d)
% RMDIR_IF  remove a folder tree when it exists.  INPUTS: d char.
if isfolder(d), rmdir(d, 's'); end
end

function r = check(msg, cond)
% CHECK  print PASS / FAIL for one condition.  INPUTS: msg char, cond
%   logical.  OUTPUTS: r logical.
r = logical(cond);
fprintf('  [%s] %s\n', ternary(r, 'PASS', 'FAIL'), msg);
end

function r = check_error(msg, id, fh)
% CHECK_ERROR  true iff fh() throws the identifier id.  INPUTS: msg char,
%   id char, fh function handle.  OUTPUTS: r logical.
r = false;
try
    fh();
    fprintf('  [FAIL] %s (no error)\n', msg);
catch err
    r = strcmp(err.identifier, id);
    fprintf('  [%s] %s (%s)\n', ternary(r, 'PASS', 'FAIL'), msg, err.identifier);
end
end

function s = ternary(c, a, b)
% TERNARY  a if c else b.  INPUTS: c logical, a, b char.  OUTPUTS: s char.
if c, s = a; else, s = b; end
end
