function out = backfill_status_layer(mode, recordMat, varargin)
%% Purpose:
%
%   The OPTIMALITY-STATUS LAYER for a library of record WITHOUT a rebuild
%   (spec 6). Two stages around the re-certification (recertify_candidates):
%
%     'harvest'  -- read every candidate the builds left on disk (arrival
%                   sheets, rib files, rib checkpoints, hole-filler files),
%                   give each primary of the record the junctions of the
%                   first candidate that is the same root (same_root), list
%                   the legacy candidates whose status is ambiguous
%                   ("conjugate test verdict 0": FAIL or UNDETERMINED) and
%                   the rib stop points worth a re-solve, and the primaries
%                   no source covers. Writes <outDir>/harvest.mat. No solve.
%     'assemble' -- read harvest.mat and every recert_*.mat, and write the v2
%                   catalog: primaries UNCHANGED (content key asserted equal
%                   to the record's), status 4 grid, junctions per entry,
%                   the alternatives table (make_alternative +
%                   dedup_alternatives) and status_key.
%
%  ASSUMPTIONS / NOTES:
%
% • Source layouts (by the variables a file holds): S = arrival sheet
%   (candidates S.cand{j}, phases (S.problem.sD, S.sA(j))); R = one rib or
%   an array/cell of ribs (.pts, .refused when present, .stop), with rec =
%   a hole-filler file (rec(k).others harvested when present); pts without
%   R = a rib checkpoint (.pts, .refused, sA from .identity). Points carry
%   their own .sD/.sA; the rib's .sA fills a missing one.
% • MIXED STRUCT GENERATIONS: candidates from before the status layer have
%   no .stage/.status; they are kept in a cell array and harmonised (missing
%   fields = []) only once, at the end. Candidates below the floor (no z, a
%   rho refusal without flyKm, an inadmissible flight) are skipped early.
% • needRecert: the legacy verdict-0 candidates, one per (phases, root).
%   stops: every rib that STALLED on a refusal the classifier places at 2
%   or 3 (a conjugate or gate finding, not a polish failure): its last
%   certified point is the seed, the stall's target phase the solve
%   (recertify_candidates, kind 'stop'). A rib with no certified point
%   yields no stop. Stalls at the same arrival phase whose targets agree to
%   1e-3 are one stop (two builds walked the same rib): the shortest step
%   from its seed is kept.
% • Junctions are donated to a primary only by a candidate AT its cell
%   (sD, sA within 1e-8, circular) that is same_root: junctions carry the
%   endpoint states, so the same costates at a neighbouring phase are
%   another transfer.
% • THE HARVEST KEY (MD5 of the candidates' z8 and phases, the stops, the
%   unmatched primaries and the work list, in order) binds the indices:
%   harvest.mat and every recert_*.mat carry it; recertify_candidates
%   (resume) and assemble refuse a mismatch, and harvest refuses to run
%   while recert_*.mat files sit in outDir.
% • assemble reads each (kind, index) once across recert files (an overlap:
%   last file wins, logged and counted). A re-certified result replaces
%   EVERY legacy copy of its (phases, root) -- e.g. the 24 x 24 build's
%   copy of a 24 x 48 candidate -- so dedup order is not load-bearing.
% • assemble: a re-certified candidate REPLACES its legacy record unless
%   the re-certification moved to another root -- then the legacy record
%   stays (status inferred) and the new root enters as its own row (spec 7:
%   a moved root is never relabelled). A re-polished primary that moved is
%   counted and NOT written; without .allowUnmatched a primary without
%   junctions makes assemble refuse (backfill_status_layer:unmatched).
% • The primaries' status is 4 by the record (every primary passed the full
%   stack when certified); a re-polish under the new certifier that does
%   not reproduce 4 is listed in status_layer.primNot4, not relabelled.
%
%% Inputs:
%
%  mode                     char                    'harvest' | 'assemble' |
%                                                   'key' (second argument:
%                                                   a harvest struct)
%  recordMat                char                    the library of record (.mat,
%                                                   one variable, one sheet)
%  'harvest':  sources      cell {1 x n} char       source .mat files
%              outDir       char                    where harvest.mat goes
%  'assemble': outDir       char                    harvest.mat + recert_*.mat
%              opts         struct (optional)       .allowUnmatched [false]
%
%% Outputs:
%
%  out                      struct                  'harvest': .cands (struct
%                                                   array, certify_root results
%                                                   + .sD .sA .source) .primJ
%                                                   {1 x nEntries} .needRecert
%                                                   .stops .unmatched
%                                                   .unmatchedInfo .nBySource
%                                                   .sources .recordMat
%                                                   .harvestKey
%                                                   'assemble': the v2 catalog
%                                                   'key': char [1 x 32]
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

switch mode
    case 'harvest'
        out = harvest(recordMat, varargin{1}, varargin{2});
    case 'assemble'
        opts = struct();  if numel(varargin) >= 2, opts = varargin{2}; end
        out = assemble(recordMat, varargin{1}, opts);
    case 'key'                                 % backfill_status_layer('key', H): the harvest key
        out = harvestKey(recordMat);
    otherwise
        error('backfill_status_layer:mode', 'mode must be ''harvest'', ''assemble'' or ''key'', got ''%s''', mode);
end
end

% ============================================================================
function H = harvest(recordMat, sources, outDir)
% HARVEST  Candidates on disk -> harvest.mat (no solve).
% INPUTS:  recordMat [char]; sources {1 x n} char; outDir [char].
% OUTPUTS: H [struct] what harvest.mat holds.
c = loadCatalog(recordMat);  sh = c.sheets(1);
if ischar(sources), sources = {sources}; end
old = dir(fullfile(outDir, 'recert_*.mat'));
if ~isempty(old)
    error('backfill_status_layer:staleRecert', ['%s holds %d recert_*.mat file(s) indexed against an earlier ' ...
          'harvest: move them away before harvesting again'], outDir, numel(old));
end
cl = {};  stops = struct('sA', {}, 'sDfrom', {}, 'sDto', {}, 'z', {}, 'Y', {}, 'why', {}, 'source', {});
nBySource = zeros(1, numel(sources));
for ks = 1:numel(sources)
    [got, st] = readSource(sources{ks});
    cl = [cl, got];
    stops = [stops, st];
    nBySource(ks) = numel(got);
end
cands = harmonise(cl);
stops = dedupStops(stops);

% ---- junctions for the primaries, by root -------------------------------------
nE = size(sh.z8, 2);  primJ = cell(1, nE);
Zc = nan(8, numel(cl));  hasY = false(1, numel(cl));  sDc = nan(1, numel(cl));  sAc = sDc;
for m = 1:numel(cl)
    Zc(:, m) = cl{m}.z(:);  sDc(m) = cl{m}.sD;  sAc(m) = cl{m}.sA;
    hasY(m) = isfield(cl{m}, 'Y') && ~isempty(cl{m}.Y);
end
onCell = @(a, b) abs(mod(a - b + 0.5, 1) - 0.5) < 1e-8;
prim = primaryList(sh);
for q = 1:numel(prim)
    k = prim(q).k;  z8 = sh.z8(:, k);
    % only a candidate AT the primary's cell may donate junctions (junctions
    % hold the endpoint states: the same root at a neighbouring phase is
    % another transfer); prefilter on t_f, then the library's rule decides
    near = find(hasY & onCell(sDc, prim(q).sD) & onCell(sAc, prim(q).sA) & abs(Zc(8, :) - z8(8)) <= 1e-6*abs(z8(8)));
    for m = near
        if same_root(Zc(:, m), z8), primJ{k} = cl{m}.Y;  break, end
    end
end
isPrim = false(1, nE);  isPrim([prim.k]) = true;
unmatched = find(isPrim & cellfun(@isempty, primJ));
unmatchedInfo = prim(ismember([prim.k], unmatched));

% ---- the ambiguous legacy candidates, one per (phases, root) -------------------
needRecert = [];
for m = 1:numel(cl)
    [~, why] = optimality_status(cl{m});
    if ~startsWith(why, 'necessary only (inferred): legacy verdict 0'), continue, end
    dup = false;
    for q = needRecert
        if nearPhase(cl{q}.sD, cl{m}.sD) && nearPhase(cl{q}.sA, cl{m}.sA) && same_root(cl{q}.z, cl{m}.z)
            dup = true;  break
        end
    end
    if ~dup, needRecert(end+1) = m; end
end

H = struct('needRecert', needRecert, 'unmatched', unmatched, 'nBySource', nBySource, ...
           'recordMat', recordMat);
H.cands = cands;  H.primJ = primJ;  H.stops = stops;  H.unmatchedInfo = unmatchedInfo;  H.sources = sources;
H.harvestKey = harvestKey(H);
if ~isfolder(outDir), mkdir(outDir); end
save(fullfile(outDir, 'harvest.mat'), '-struct', 'H');     % v7: v7.3 writes one HDF5 object per field per candidate
end

% ============================================================================
function c2 = assemble(recordMat, outDir, opts)
% ASSEMBLE  harvest.mat + recert_*.mat -> the v2 catalog (saved and returned).
% INPUTS:  recordMat [char]; outDir [char]; opts [struct] .allowUnmatched.
% OUTPUTS: c2 [struct] the catalog.
allowUnmatched = isfield(opts, 'allowUnmatched') && isequal(opts.allowUnmatched, true);
[c, vn] = loadCatalog(recordMat);  sh = c.sheets(1);
H = load(fullfile(outDir, 'harvest.mat'));
key = harvestKey(H);
assert(isfield(H, 'harvestKey') && strcmp(H.harvestKey, key), 'backfill_status_layer:harvestKey', ...
       'harvest.mat does not match its own stamped key: rebuild it');
cl = num2cell(H.cands);
stops = struct([]);  if isfield(H, 'stops'), stops = H.stops; end
primJ = H.primJ;
extra = {};  nRecert = 0;  nMoved = 0;  movedPrim = [];  primNot4 = [];  nErr = 0;

% ---- the re-certifications ---------------------------------------------------
% one item per (kind, index) across files; files are bound to this harvest
files = dir(fullfile(outDir, 'recert_*.mat'));
allItems = {};  ids = {};  nOverlap = 0;
for kf = 1:numel(files)
    R = load(fullfile(outDir, files(kf).name));
    if ~(isfield(R, 'harvestKey') && strcmp(R.harvestKey, key))
        error('backfill_status_layer:harvestKey', ['%s was made against another harvest (key mismatch): its ' ...
              'indices do not address these candidates'], files(kf).name);
    end
    for it = R.items(:)'
        id = sprintf('%s:%d', it.kind, it.index);
        at = find(strcmp(ids, id), 1);
        if isempty(at), allItems{end+1} = it;  ids{end+1} = id;
        else
            nOverlap = nOverlap + 1;  allItems{at} = it;   % last write wins
            fprintf('backfill_status_layer: %s re-certified in more than one file; %s wins\n', id, files(kf).name);
        end
    end
end
recertified = [];
for q = 1:numel(allItems)
    it = allItems{q};
    if isfield(it, 'err') && ~isempty(it.err), nErr = nErr + 1;  continue, end
    switch it.kind
        case 'cand'
            Cr = it.C;  Cr.sD = cl{it.index}.sD;  Cr.sA = cl{it.index}.sA;
            if it.moved
                nMoved = nMoved + 1;
                Cr.source = [cl{it.index}.source ' -> re-certification moved to another root'];
                extra{end+1} = Cr;
            else
                Cr.source = [cl{it.index}.source ' (re-certified)'];
                cl{it.index} = Cr;  nRecert = nRecert + 1;  recertified(end+1) = it.index;
            end
        case 'stop'
            s = stops(it.index);  Cr = it.C;  Cr.sD = s.sDto;  Cr.sA = s.sA;
            Cr.source = sprintf('%s stop at sD %.4f, re-solved from %.4f', s.source, s.sDto, s.sDfrom);
            extra{end+1} = Cr;
        case 'prim'
            if it.moved
                movedPrim(end+1) = it.index;
            else
                primJ{it.index} = it.C.Y;
                if ~(isfield(it.C, 'status') && isequal(it.C.status, 4)), primNot4(end+1) = it.index; end
            end
    end
end
% THE RE-CERTIFIED RESULT REPLACES EVERY COPY of its (phases, root): the
% harvest listed one representative per (phases, root), and another build's
% copy of the same legacy candidate must not survive dedup with the
% inferred status (row order then does not matter)
for q = recertified
    for m = 1:numel(cl)
        if m == q || ~(nearPhase(cl{m}.sD, cl{q}.sD) && nearPhase(cl{m}.sA, cl{q}.sA)), continue, end
        if isfield(cl{m}, 'stage') && ~isempty(cl{m}.stage), continue, end      % already stamped
        if same_root(cl{m}.z, H.cands(q).z)
            Cr = cl{q};  Cr.source = [cl{m}.source ' (re-certified as ' H.cands(q).source ')'];  cl{m} = Cr;
        end
    end
end
cl = [cl(:)', extra];

% ---- every primary has its junctions, or refuse --------------------------------
prim = primaryList(sh);
missing = [prim(cellfun(@isempty, primJ([prim.k]))).k];
if ~isempty(missing) && ~allowUnmatched
    error('backfill_status_layer:unmatched', ['%d primar%s without junctions (first entry %d; %d re-polish(es) ' ...
          'moved off the root): run recertify_candidates over every chunk, or pass allowUnmatched'], ...
          numel(missing), plural(numel(missing)), missing(1), numel(movedPrim));
end

% ---- the status layer of the primaries ----------------------------------------
nE = size(sh.z8, 2);
sh.status = zeros(size(sh.has_solution), 'int8');
sh.status(logical(sh.has_solution)) = 4;
reasons = repmat({''}, 1, nE);  reasons([prim.k]) = {'full stack passed'};
sh.status_reason = reasons;
sh.junctions = primJ;

% ---- the alternatives -----------------------------------------------------------
rows = struct([]);
for m = 1:numel(cl)
    Cm = cl{m};
    iD = idxOf(sh.sD_frac, Cm.sD);  iA = idxOf(sh.sA_frac, Cm.sA);
    a = make_alternative(Cm, Cm.sD, Cm.sA, iD, iA, Cm.source);
    if isempty(a), continue, end
    if isempty(rows), rows = a; else, rows(end+1) = a; end
end
alts = dedup_alternatives(rows, sh);
for k = 1:numel(alts), alts(k).sheet = 1; end

c2 = c;  c2.sheets = sh;
c2.alternatives = alts;
c2.status_key = status_key();
c2.status_layer = struct('built', char(datetime('now', 'Format', 'yyyy-MM-dd')), 'recordMat', recordMat, ...
                         'sources', {H.sources}, 'nCandidates', numel(H.cands), 'nRecert', nRecert, ...
                         'nMoved', nMoved, 'movedPrimaries', movedPrim, 'primNot4', primNot4, ...
                         'nRecertErrors', nErr, 'nRecertOverlap', nOverlap, 'nStops', numel(stops), ...
                         'nUnmatched', numel(missing), 'harvestKey', key, ...
                         'note', ['status layer backfilled from candidates on disk (spec 6); primaries are the ' ...
                                  'record''s, bit-identical (catalog_content_key)']);
kRec = catalog_content_key(c);  kNew = catalog_content_key(c2);
assert(strcmp(kRec, kNew), 'backfill_status_layer:contentKey', ...
       'the primaries changed (content key %s vs the record''s %s)', kNew, kRec);
outMat = fullfile(outDir, 'costate_catalog_dro_tulip_70mN.mat');
assert(~strcmp(char(java.io.File(outMat).getCanonicalPath()), char(java.io.File(recordMat).getCanonicalPath())), ...
       'backfill_status_layer:overwrite', 'refusing to overwrite the record %s', recordMat);
S.(vn) = c2;  save(outMat, '-struct', 'S');
end

% ============================================================================
function [cl, stops] = readSource(f)
% READSOURCE  Every above-floor candidate of one source file.
% INPUTS:  f [char] .mat path.
% OUTPUTS: cl {1 x n} candidate structs (+ .sD .sA .source); stops [struct].
cl = {};  stops = struct('sA', {}, 'sDfrom', {}, 'sDto', {}, 'z', {}, 'Y', {}, 'why', {}, 'source', {});
L = load(f);  base = sourceTag(f);
if isfield(L, 'S')                                            % arrival sheet
    S = L.S;  sD = 0;  if isfield(S, 'problem') && isfield(S.problem, 'sD'), sD = S.problem.sD; end
    for jc = 1:numel(S.cand)
        c = S.cand{jc};
        for m = 1:numel(c)
            cl = addCand(cl, c(m), sD, S.sA(jc), sprintf('%s cand, column %d', base, jc));
        end
    end
elseif isfield(L, 'R')                                        % rib(s) or hole filler
    Rs = L.R;  if ~iscell(Rs), Rs = num2cell(Rs); end
    for kr = 1:numel(Rs)
        R = Rs{kr};  sA = NaN;  if isfield(R, 'sA'), sA = R.sA; end
        tag = sprintf('%s rib sA %.4f', base, sA);
        [cl, last] = addPoints(cl, R, 'pts', sA, [tag ' point']);
        cl = addPoints(cl, R, 'refused', sA, [tag ' refusal']);
        st = stallOf(R, last, sA, tag);
        if ~isempty(st), stops(end+1) = st; end
    end
    if isfield(L, 'rec') && isfield(L.rec, 'others')
        for kc = 1:numel(L.rec)
            o = L.rec(kc).others;
            for m = 1:numel(o)
                cl = addCand(cl, o(m), L.rec(kc).sD, L.rec(kc).sA, sprintf('%s filler cell (%d,%d) other', ...
                             base, L.rec(kc).iD, L.rec(kc).iA));
            end
        end
    end
elseif isfield(L, 'pts')                                      % rib checkpoint
    sA = NaN;  if isfield(L, 'identity') && isfield(L.identity, 'sA'), sA = L.identity.sA; end
    tag = sprintf('%s checkpoint sA %.4f', base, sA);
    cl = addPoints(cl, L, 'pts', sA, [tag ' point']);
    cl = addPoints(cl, L, 'refused', sA, [tag ' refusal']);
else
    error('backfill_status_layer:source', 'unrecognised source layout in %s (variables: %s)', f, ...
          strjoin(fieldnames(L)', ' '));
end
end

function [cl, last] = addPoints(cl, R, f, sA, tag)
% ADDPOINTS  The points of R.(f) as candidates; last = the last one with
% z and Y (the stop seed).  INPUTS: cl; R; f; sA; tag.  OUTPUTS: cl; last.
last = [];
if ~isfield(R, f) || isempty(R.(f)), return, end
P = R.(f);
for m = 1:numel(P)
    p = P(m);  sAp = sA;  if isfield(p, 'sA') && ~isempty(p.sA), sAp = p.sA; end
    sDp = NaN;  if isfield(p, 'sD') && ~isempty(p.sD), sDp = p.sD; end
    [cl, kept] = addCand(cl, p, sDp, sAp, tag);
    if kept && isfield(p, 'Y') && ~isempty(p.Y), last = cl{end}; end
end
end

function [cl, kept] = addCand(cl, C, sD, sA, source)
% ADDCAND  Append one candidate if it is above the floor, with its phases
% and provenance.  INPUTS: cl; C; sD; sA; source.  OUTPUTS: cl; kept.
kept = false;
if ~isstruct(C) || ~isfield(C, 'z'), return, end            % e.g. a rho refusal
if optimality_status(C) < 1, return, end                     % below the floor
if ~(isfinite(sD) && isfinite(sA))
    fprintf('backfill_status_layer: skipped a candidate of %s without phases (sD %g, sA %g)\n', source, sD, sA);
    return
end
C.sD = mod(sD, 1);  C.sA = mod(sA, 1);  C.source = source;
cl{end+1} = C;  kept = true;
end

function st = stallOf(R, last, sA, tag)
% STALLOF  The stop item of a rib that stalled on a conjugate or gate finding.
% INPUTS:  R [struct] rib; last [struct] its last certified point ([] none);
%          sA [scalar]; tag [char].
% OUTPUTS: st [struct] (.sA .sDfrom .sDto .z .Y .why .source) or [].
st = [];
if ~isfield(R, 'stop') || isempty(last), return, end
tok = regexp(R.stop, '^stalled at sD = ([\d.]+) stepping to ([\d.]+): (.*)$', 'tokens', 'once');
if isempty(tok), return, end
% the classifier reads the refusal in the certifier's words; a stand-in
% above the floor lets it place the TEXT (2 conjugate / 3 gates / 1 other)
code = optimality_status(struct('z', ones(8, 1), 'flyKm', 0, 'flyVms', 0, 'reason', tok{3}, 'ok', false));
if code < 2, return, end
st = struct('sA', sA, 'sDfrom', last.sD, 'sDto', str2double(tok{2}), 'z', last.z(:), 'Y', last.Y, ...
            'why', tok{3}, 'source', tag);
end

function key = harvestKey(H)
% HARVESTKEY  MD5 of what the re-certification indexes, in order: every
% candidate's z8 and phases, every stop, every unmatched primary, the work
% list. recert files carry it; resume and assemble refuse a mismatch.
% INPUTS: H [struct] harvest.  OUTPUTS: key [char 1 x 32].
md = java.security.MessageDigest.getInstance('MD5');
put = @(kind, v) md.update([uint8(kind), typecast(double(v(:)'), 'uint8')]);
for m = 1:numel(H.cands), put('c', [H.cands(m).z(:); H.cands(m).sD; H.cands(m).sA]); end
for m = 1:numel(H.stops), s = H.stops(m);  put('s', [s.sA; s.sDfrom; s.sDto; s.z(:)]); end
if isfield(H, 'unmatchedInfo')
    for m = 1:numel(H.unmatchedInfo), p = H.unmatchedInfo(m);  put('p', [p.k; p.sD; p.sA; p.z8(:)]); end
end
put('n', H.needRecert);  put('u', H.unmatched);
key = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end

function st = dedupStops(st)
% DEDUPSTOPS  One stop per stall: stops at the same arrival phase whose
% target phases agree to 1e-3 (the same stall, walked by two builds or by a
% rib and its finer twin) are one; the one with the SHORTEST step from its
% seed is kept.  INPUTS: st [struct] stops.  OUTPUTS: st.
keep = true(1, numel(st));
for a = 1:numel(st)
    for b = 1:numel(st)
        if a == b || ~keep(b) || ~nearPhase(st(a).sA, st(b).sA), continue, end
        if abs(mod(st(a).sDto - st(b).sDto + 0.5, 1) - 0.5) >= 1e-3, continue, end
        ga = abs(mod(st(a).sDto - st(a).sDfrom + 0.5, 1) - 0.5);  gb = abs(mod(st(b).sDto - st(b).sDfrom + 0.5, 1) - 0.5);
        if gb < ga || (gb == ga && b < a), keep(a) = false;  break, end
    end
end
st = st(keep);
end

function t = sourceTag(f)
% SOURCETAG  A source's provenance tag: its path below results/ (two builds
% hold files of the same name), else its file name.  INPUTS: f [char].
% OUTPUTS: t [char].
t = regexprep(f, '^.*[/\\]results[/\\]', '');
if strcmp(t, f), [~, t] = fileparts(f); else, t = regexprep(t, '\.mat$', ''); end
end

function S = harmonise(cl)
% HARMONISE  A struct array from candidates of mixed field sets (missing
% fields = []).  INPUTS: cl {1 x n}.  OUTPUTS: S [1 x n] struct.
if isempty(cl), S = struct([]); return, end
names = {};
for m = 1:numel(cl), names = union(names, fieldnames(cl{m}), 'stable'); end
for m = 1:numel(cl)
    for f = setdiff(names, fieldnames(cl{m}))', cl{m}.(f{1}) = []; end
    cl{m} = orderfields(cl{m}, names);
end
S = [cl{:}];
end

function P = primaryList(sh)
% PRIMARYLIST  The record's primaries.  INPUTS: sh [struct] sheet.
% OUTPUTS: P [struct] (.k .iD .iA .sD .sA .z8).
P = struct('k', {}, 'iD', {}, 'iA', {}, 'sD', {}, 'sA', {}, 'z8', {});
for iA = 1:numel(sh.sA_frac)
    for iD = 1:numel(sh.sD_frac)
        if ~sh.has_solution(iD, iA, 1), continue, end
        k = sh.entry_index(iD, iA, 1);
        P(end+1) = struct('k', k, 'iD', iD, 'iA', iA, 'sD', sh.sD_frac(iD), 'sA', sh.sA_frac(iA), 'z8', sh.z8(:, k));
    end
end
end

function [c, vn] = loadCatalog(f)
% LOADCATALOG  The one catalog variable of a .mat.  INPUTS: f.  OUTPUTS: c; vn.
L = load(f);  vn = char(fieldnames(L));  c = L.(vn);
assert(isscalar(c.sheets), 'backfill_status_layer:sheets', 'one sheet expected');
end

function k = idxOf(grid, v)
% IDXOF  Grid index of phase v (mod 1), 0 when off the grid.  INPUTS: grid; v.
% OUTPUTS: k.
k = find(abs(mod(grid - v + 0.5, 1) - 0.5) < 1e-8, 1);
if isempty(k), k = 0; end
end

function tf = nearPhase(a, b)
% NEARPHASE  Same phase mod 1.  INPUTS: a; b.  OUTPUTS: tf.
tf = abs(mod(a - b + 0.5, 1) - 0.5) < 1e-9;
end

function s = plural(n)
% PLURAL  'y' or 'ies'.  INPUTS: n.  OUTPUTS: s.
if n == 1, s = 'y'; else, s = 'ies'; end
end
