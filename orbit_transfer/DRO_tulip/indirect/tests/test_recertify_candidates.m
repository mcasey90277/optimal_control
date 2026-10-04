function ok = test_recertify_candidates()
% TEST_RECERTIFY_CANDIDATES  recertify_candidates without a solve (fake
% certifier, fake setup and seeder through its seams): the item order
% (legacy candidates, rib stop points, unmatched primaries), the chunk
% share chunk:nChunk:end, resume (a re-run certifies nothing already
% saved, but retries an item that threw), the seeds it builds, and the
% moved rule (a re-polish that leaves the root is recorded moved).
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
tmp = tempname;  mkdir(tmp);
callLog = fullfile(tmp, 'calls.txt');

% ---- a fake harvest: 2 legacy candidates, 1 stop, 2 unmatched primaries ----
z1 = [1; 2; 3; 4; 5; 6; 7; 8];  z2 = 2*z1;  zp5 = 3*z1;  zp7 = 4*z1;
cand = @(z, sD, sA) struct('ok', false, 'reason', 'conjugate test verdict 0', 'z', z, 'Y', z(1)*ones(14, 24), ...
                           'flyKm', 0.1, 'flyVms', 0.01, 'sD', sD, 'sA', sA, 'source', 'fake');
cands = [cand(z1, 0, 0.25), cand(z2, 0, 0.5)];
needRecert = [1 2];
stops = struct('sA', 0.75, 'sDfrom', 0.5, 'sDto', 0.49, 'z', z1, 'Y', ones(14, 24), 'why', 'dense', 'source', 'fake rib');
unmatched = [5 7];
unmatchedInfo = struct('k', {5, 7}, 'sD', {0.125, 0.375}, 'sA', {0.25, 0.5}, 'z8', {zp5, zp7});
harvestMat = fullfile(tmp, 'harvest.mat');
harvestKey = backfill_status_layer('key', struct('cands', {cands}, 'needRecert', needRecert, 'stops', {stops}, ...
                                                 'unmatched', unmatched, 'unmatchedInfo', {unmatchedInfo}));
save(harvestMat, 'cands', 'needRecert', 'stops', 'unmatched', 'unmatchedInfo', 'harvestKey');

B = struct('stateD', @(s) [s; 0; 0; 0; 0; 0], 'stateA', @(s) [0; s; 0; 0; 0; 0], 'Tnd', 1, 'cnd', 1, 'mu', 0.01);
seeder = @(z8, rv0, K) struct('tf', z8(8), 'tGrid', linspace(0, z8(8), K + 1), 'Y', repmat([rv0(:); 1; z8(1:7)], 1, K + 1));
o = struct('B', B, 'seedFromZ8', seeder, 'pool', [], 'allowUnfenced', true, ...
           'certifier', @(seed, rv0, rvf, B_, co) fakeCert(seed, rv0, rvf, B_, co, callLog));

% ---- chunk 1 of 2: items 1, 3, 5 of [cand1 cand2 stop1 prim5 prim7] ---------
out1 = fullfile(tmp, 'recert_1.mat');
recertify_candidates(harvestMat, 1, 2, out1, o);
L = load(out1);  it = L.items;
ok = chk(ok, isequal({it.kind}, {'cand', 'stop', 'prim'}) && isequal([it.index], [1 1 7]), ...
         'chunk 1 of 2 takes items 1, 3, 5 in order (cand, stop, prim)');
ok = chk(ok, nCalls(callLog) == 3, 'three certifier calls');
ok = chk(ok, ~it(1).moved && it(3).moved, 'the unmoved candidate is kept; the primary that left its root is recorded moved');
ok = chk(ok, isequal(it(1).C.seedY1, [0; 0; 0; 0; 0; 0; 1; z1(1:7)]) && it(1).C.sD == 0 && it(1).C.sA == 0.25, ...
         'a candidate is re-seeded from its own junctions, node 1 forced to [rv0; 1; z(1:7)], at its own phases');
ok = chk(ok, it(2).C.sD == 0.49 && it(2).C.sA == 0.75 && isequal(it(2).C.seedY1(8:14), z1(1:7)), ...
         'the stop point is seeded from the rib''s last point and solved at the stop phase');
ok = chk(ok, it(3).C.sD == 0.375 && isequal(it(3).C.seedY1(8:14), zp7(1:7)), 'an unmatched primary is re-polished from its own z8');

% ---- resume: nothing re-certified --------------------------------------------
recertify_candidates(harvestMat, 1, 2, out1, o);
L = load(out1);
ok = chk(ok, nCalls(callLog) == 3 && numel(L.items) == 3, 'a re-run of a finished chunk certifies nothing');

% ---- chunk 2 of 2, one item throws: saved with its error, retried on re-run ----
out2 = fullfile(tmp, 'recert_2.mat');
o2 = o;  o2.certifier = @(seed, rv0, rvf, B_, co) fakeCert(seed, rv0, rvf, B_, co, callLog, 0.5);   % throws at sA 0.5
recertify_candidates(harvestMat, 2, 2, out2, o2);
L = load(out2);  it = L.items;
ok = chk(ok, isequal({it.kind}, {'cand', 'prim'}) && isequal([it.index], [2 5]) && ~isempty(it(1).err) && isempty(it(2).err), ...
         'chunk 2 of 2 takes items 2, 4; a throw is saved with its error');
ok = chk(ok, ~it(2).moved, 'the primary that stayed on its root is not moved');
n0 = nCalls(callLog);
recertify_candidates(harvestMat, 2, 2, out2, o);
L = load(out2);  it = L.items;
ok = chk(ok, nCalls(callLog) == n0 + 1 && numel(it) == 2 && isempty(it(1).err) && it(1).index == 2, ...
         'the re-run retries only the item that threw, replacing it');
lines = splitlines(strtrim(fileread([out2 '.log'])));
ok = chk(ok, numel(lines) == 3, 'one log line per item certified');
% a chunk file of ANOTHER harvest is refused on resume
L = load(out1);  items = L.items;  harvestKey = repmat('0', 1, 32);  save(out1, 'items', 'harvestKey');
threw = '';
try, recertify_candidates(harvestMat, 1, 2, out1, o); catch ME, threw = ME.identifier; end
ok = chk(ok, strcmp(threw, 'recertify_candidates:harvestKey') && nCalls(callLog) == n0 + 1, ...
         'a resume refuses a chunk file made against another harvest');
ok = chk(ok, ~isfile([out1 '.part']) && ~isfile([out2 '.part']), 'no partial file left behind');
rmdir(tmp, 's');
if ok, fprintf('test_recertify_candidates: ALL PASS\n'); else, fprintf('test_recertify_candidates: FAIL\n'); end
end

function C = fakeCert(seed, rv0, rvf, B, co, callLog, throwAtSA)
% FAKECERT  A certifier that returns the seed's own root (moved by 1e-3 for
% the primary at sD 0.375), logs the call, and throws at sA = throwAtSA.
% INPUTS: seed; rv0; rvf; B; co (.sD .sA); callLog; throwAtSA (optional).
% OUTPUTS: C (.z .Y .status .reason .sD .sA .seedY1).
fid = fopen(callLog, 'a');  fprintf(fid, '%g %g\n', co.sD, co.sA);  fclose(fid);
if nargin >= 7 && co.sA == throwAtSA, error('fake:throw', 'fake certifier throws'); end
z = [seed.Y(8:14, 1); seed.tf];
if abs(co.sD - 0.375) < 1e-12, z(1) = z(1)*(1 + 1e-3); end
C = struct('ok', true, 'reason', 'certified', 'z', z, 'Y', seed.Y(:, 1:end-1), 'status', 4, 'stage', 8, ...
           'status_reason', 'full stack passed', 'flyKm', 0, 'flyVms', 0, 'sD', co.sD, 'sA', co.sA, 'seedY1', seed.Y(:, 1));
assert(numel(rv0) == 6 && numel(rvf) == 6 && isstruct(B));
end

function n = nCalls(f)
% NCALLS  Lines in the call log.  INPUTS: f.  OUTPUTS: n.
n = 0;  if isfile(f), n = numel(splitlines(strtrim(fileread(f)))); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
