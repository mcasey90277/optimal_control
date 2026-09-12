function ok = test_second_order_parallel()
%% Purpose:
%
%   The parallel sweep must be the SERIAL sweep, faster -- not a second
%   implementation that happens to agree. This runs a handful of catalog
%   entries both ways, into two throwaway sidecars, and requires every
%   written-back field to match BITWISE.
%
%   The measurement is pure (no shared state, no file access inside
%   measureOne), so any difference would mean the parallel path is feeding
%   it different inputs -- which is exactly the bug worth catching, since
%   the parallel driver gathers each chunk's inputs in the client rather
%   than letting workers index the catalog.
%
%   It also checks the resume contract the chunking changes: the sidecar is
%   written after each CHUNK, so a run stopped by its entry budget leaves a
%   sidecar that a later run continues from without re-measuring.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/12/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
ind = fullfile(fileparts(here), 'DRO_tulip', 'indirect');
addpath(here, ind);
catMat = fullfile(ind, 'results', 'costate_catalog_dro_tulip_70mN.mat');
if ~isfile(catMat)
    fprintf('  SKIP  no 70 mN catalog on disk (%s)\n', catMat);  return
end
nEnt = 4;
sideS = [tempname '_serial.mat'];  sideP = [tempname '_par.mat'];
cleanup = onCleanup(@() cellfun(@(f) delIf(f), {sideS, sideP}));

pool = gcp('nocreate');
ok = chk(ok, ~isempty(pool), sprintf('a pool is available (%d workers)', ...
         tern(isempty(pool), 0, pool.NumWorkers)));

tic;
Ss = second_order_pass(catMat, struct('sideMat', sideS, 'maxEntries', nEnt, 'parallel', false));
tS = toc;
tic;
Sp = second_order_pass(catMat, struct('sideMat', sideP, 'maxEntries', nEnt, 'parallel', true, 'chunk', nEnt));
tP = toc;

ok = chk(ok, Ss.nDone == nEnt && Sp.nDone == nEnt, ...
         sprintf('both paths measured %d entries (serial %.0f s, parallel %.0f s, %.1fx)', ...
                 nEnt, tS, tP, tS/max(tP, eps)));

% BITWISE on every written-back field
f = {'nInterior', 'multiplicity', 'minRelSigma', 'nInteriorCand', 'nNearMiss', ...
     'nZero', 'nUnresolved', 'nEndCand', 'conjClear', 'h6margin', 'h6ok', ...
     'h6clearance', 'liftMargin', 'liftCertified'};
worst = '';  allSame = true;
for k = 1:numel(f)
    a = [Ss.rows(1:nEnt).(f{k})];  b = [Sp.rows(1:nEnt).(f{k})];
    same = isequaln(a, b);
    allSame = allSame && same;
    if ~same, worst = [worst ' ' f{k}]; end %#ok<AGROW>
end
ok = chk(ok, allSame, sprintf('every written-back field is bitwise identical%s', ...
         tern(allSame, sprintf(' (%d fields, %d entries)', numel(f), nEnt), [': DIFFER at' worst])));

% the candidate records too, not just the summary counts
cSame = true;
for q = 1:nEnt
    ca = Ss.rows(q).candidates;  cb = Sp.rows(q).candidates;
    cSame = cSame && numel(ca) == numel(cb);
    if cSame && ~isempty(ca)
        cSame = cSame && isequaln([ca.sigMinRel], [cb.sigMinRel]) && ...
                         isequaln({ca.kind}, {cb.kind});
    end
end
ok = chk(ok, cSame, 'and the per-candidate records match, kind for kind');

% RESUME across the chunk boundary
S2 = second_order_pass(catMat, struct('sideMat', sideP, 'maxEntries', 2, 'parallel', true, 'chunk', 2));
ok = chk(ok, S2.nDone == nEnt + 2, ...
         sprintf('a resumed run continues rather than re-measuring: %d done after +2', S2.nDone));

if ok, fprintf('TEST_SECOND_ORDER_PARALLEL: ALL PASS\n'); else, fprintf('TEST_SECOND_ORDER_PARALLEL: FAIL\n'); end
end

function delIf(f)
% DELIF  Delete a file if it exists.  INPUTS: f.  OUTPUTS: none.
if isfile(f), delete(f); end
end

function v = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
