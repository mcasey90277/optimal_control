function ok = test_plot_pair_transfer()
%% Purpose:
%
%   Test of plot_pair_transfer on a solved DRO -> ELFO transfer (raanDeg 0,
%   phases 0 / 0, 15 N). The solve runs ONCE and is cached in
%   tempdir/plot_pair_transfer_test_cache keyed by a SHA-256 of the solve
%   inputs and of the source text of pair_transfer, located_min_altitude,
%   thrust_ladder_library, fly_transfer, ms_tfmin and get_family_orbit (cold about 3.5 min:
%   run it in the background). Checks: one figure, 2 x 3 tiles, written
%   PNG file when opts.file is set, nothing written when it is not, a
%   figure that stays invisible with visible 'off', the verdict lines of
%   mode 'verdict' (one per blocking gate, PASS on the solved result),
%   T.ok = false errors plot_pair_transfer:notSolved, and a regated
%   failure shows FAIL in the lines.
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
cc      = fileparts(here);
ot      = fileparts(cc);
home    = getenv('HOME');
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(cc, fullfile(fileparts(ot), 'oclib'), fullfile(home, 'casadi-3.7.0'));
if isempty(which('pumpkyn.cr3bp.tfMin'))
    od = cd(fullfile(home, 'Desktop', 'proj7', 'external', 'pumpkynPie'));  startup();  cd(od);
end

ok = true;
dro1  = struct('family', 'dro',  'params', struct('tau', 1));
elfo0 = struct('family', 'elfo', 'params', struct('raanDeg', 0));
T = cached_solve(dro1, elfo0);
fprintf('  solved cell: ok %d, t_f %.6f ND, gates.pass %d\n', T.ok, T.tf_nd, T.gates.pass);
ok = check('cell solved and gated', T.ok && T.gates.pass) && ok;

tmp = tempname;  mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));

%% figure, tiles, file
f1 = fullfile(tmp, 'pair.png');
hf = plot_pair_transfer(T, struct('visible', 'off', 'file', f1, 'title', 'test'));
closeHf = onCleanup(@() close_if_valid(hf));
tl = findobj(hf, 'Type', 'tiledlayout');
nax = numel(findobj(hf, 'Type', 'axes'));
ok = check(sprintf('one figure, tiledlayout 2 x 3 (grid %s), %d axes', mat2str(tl.GridSize), nax), ...
           isgraphics(hf, 'figure') && isequal(tl.GridSize, [2 3]) && nax == 6) && ok;
ok = check('figure stays invisible with visible off', strcmp(hf.Visible, 'off')) && ok;
d = dir(f1);
ok = check(sprintf('PNG written (%d file)', numel(d)), numel(d) == 1 && d(1).bytes > 5e3) && ok;
hf2 = plot_pair_transfer(T, struct('visible', 'off'));
closeHf2 = onCleanup(@() close_if_valid(hf2));
ok = check('no file written when opts.file is empty', numel(dir(fullfile(tmp, '*.png'))) == 1) && ok;

%% tile content: time axes, located minima, verdict text inside its tile
axs = findobj(hf, 'Type', 'axes');
tEnd = T.flight.t(end) * 382981.289129055 / 86400;
tiles = flipud(axs);                              % creation order: 1..6 (6 = text tile)
okx = true;
for kk = [3 4 5]
    okx = okx && isequal(tiles(kk).XLim, [0 tEnd]);
end
ok = check(sprintf('thrust / altitude / mass tiles have xlim [0 %.4f] d', tEnd), okx) && ok;
nv = numel(findobj(tiles(4), 'Type', 'line', 'Marker', 'v'));
ok = check(sprintf('altitude tile marks %d located minima (Moon and Earth finite)', nv), nv == 2) && ok;
yl = tiles(4).YLim;
ok = check(sprintf('altitude y lower limit %.3g <= 0.5 x Moon floor / located minimum', yl(1)), ...
           yl(1) <= 0.5 * min([T.thresholds.floorKm, T.gates.moonMinKm])) && ok;
drawnow;
ht = findobj(hf, 'Tag', 'verdictText');
ex = ht.Extent;                                   % normalized to the text tile's axes
ok = check(sprintf('verdict text fits its tile (right edge %.3f, bottom %.3f of the tile)', ex(1) + ex(3), ex(2)), ...
           ex(1) >= 0 && ex(1) + ex(3) <= 1 && ex(2) >= 0) && ok;

%% verdict lines
L = plot_pair_transfer('verdict', T);
txt = strjoin(L, newline);
fprintf('%s\n', txt);
ok = check(sprintf('%d verdict lines, all gates PASS', numel(L)), ...
           contains(txt, 'flown') && contains(txt, 'moon') && contains(txt, 'PASS') && ~contains(txt, 'FAIL')) && ok;
cases = {struct('floorKm', 1e5), struct('gateKm', 1e-12), struct('accTol', 0), struct('floorEarthKm', 1e7)};
for kk = 1:numel(cases)
    Tr = pair_transfer('regate', T, cases{kk});
    rows = failed_rows(plot_pair_transfer('verdict', Tr));
    ok = check(sprintf('regate %s: FAIL rows {%s} equal gates.failed {%s}', jsonencode(cases{kk}), strjoin(rows, ','), ...
               strjoin(Tr.gates.failed, ',')), ~isempty(rows) && isequal(sort(rows(:)), sort(Tr.gates.failed(:)))) && ok;
end
ok = check_error('verdict mode without T refuses by name', 'plot_pair_transfer:usage', ...
                 @() plot_pair_transfer('verdict')) && ok;
ok = check_error('no arguments refuses by name', 'plot_pair_transfer:usage', @() plot_pair_transfer()) && ok;

%% refuses an unsolved transfer
Tbad = T;  Tbad.ok = false;  Tbad.reason = 'top rung failed';
ok = check_error('T.ok = false refuses by name', 'plot_pair_transfer:notSolved', ...
                 @() plot_pair_transfer(Tbad, struct('visible', 'off'))) && ok;
fprintf('test_plot_pair_transfer: %s\n', ternary(ok, 'ALL PASS', 'FAILURES'));
end

% ------------------------------------------------------------------------
function T = cached_solve(dep, arr)
% CACHED_SOLVE  pair_transfer at 15 N through a tempdir cache keyed by the
%   inputs and the source text of pair_transfer.  INPUTS: dep, arr structs.
%   OUTPUTS: T struct.
cd_ = fullfile(tempdir, 'plot_pair_transfer_test_cache');
if ~isfolder(cd_), mkdir(cd_); end
src = sha_hex([jsonencode(struct('dep', dep, 'arr', arr, 'v', 'v2')), source_hash({'pair_transfer', 'located_min_altitude', ...
               'thrust_ladder_library', 'fly_transfer', 'ms_tfmin', 'get_family_orbit'})]);
f = fullfile(cd_, ['pp_' src '.mat']);
if isfile(f)
    S = load(f, 'T');  T = S.T;
    fprintf('  (cached: %s)\n', f);
else
    T = pair_transfer(dep, arr, 0, 0, 15, struct('quiet', true));
    save(f, 'T');
end
end

function h = source_hash(names)
% SOURCE_HASH  SHA-256 over the source text of the named functions.  INPUTS: names cellstr.
%   OUTPUTS: h char.
parts = cell(1, numel(names));
for k = 1:numel(names)
    f = which(names{k});
    if isempty(f), error('test_plot_pair_transfer:source', 'cannot resolve %s for the cache key', names{k}); end
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

function close_if_valid(h)
% CLOSE_IF_VALID  close a figure handle we opened.  INPUTS: h handle.
if ~isempty(h) && isgraphics(h), close(h); end
end

function rows = failed_rows(L)
% FAILED_ROWS  names of the gate rows of a verdict table that end in FAIL (the BLOCKING summary row excluded).
%   INPUTS: L cellstr.  OUTPUTS: rows cellstr.
rows = {};
for k = 1:numel(L)
    tk = regexp(L{k}, '^(\S+)\s.*FAIL\s*$', 'tokens', 'once');
    if ~isempty(tk) && ~strcmp(tk{1}, 'BLOCKING'), rows = [rows, tk(1)]; end
end
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
