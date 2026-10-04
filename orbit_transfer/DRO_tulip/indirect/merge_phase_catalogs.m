function [M, info] = merge_phase_catalogs(base, donor, opts)
%% Purpose:
%
%   Merge two phase catalogs of the SAME problem cell by cell: the result is
%   the base catalog, except that wherever the donor has a certified entry
%   and the base has none, or the donor's is faster, the donor's entry is
%   taken -- its costates, its note, and every per-cell grid that belongs to
%   it (conjugate verdict, hypothesis gates, second-order margins, family
%   code). Made for adopting the 24 x 48 library over the 24 x 24 record
%   (FINDINGS 93): base = the 24 x 48, donor = the record, so the merged
%   library is no worse than the record on any cell it shares.
%
%  ASSUMPTIONS / NOTES:
%
% • NOT A CERTIFICATE. A taken entry carries the donor's verdicts, which were
%   computed by the donor's build; the merged file must be AUDITED
%   (audit_phase_catalog recomputes every verdict and fails closed) before
%   it is relied on. The merge only moves numbers.
% • Same problem, or refused: engine, rungs, nondimensionalisation, orbits,
%   family numbering (a family code means a label only through
%   .families.labels) and the second-order settings must agree.
% • Cells are matched by PHASE (circular, 1e-8), so the donor may sit on a
%   coarser grid; donor cells off the base grid are ignored, and a donor
%   sharing no cell is refused. Only sheet 1 (one sheet each).
% • "Faster" means by more than .tolDays: two builds of the same root differ
%   at solver tolerance and the base's entry is kept then.
% • A faster donor entry REPLACES the base entry in its z8 column; a filled
%   hole takes a new column at the end. Base cells the donor does not
%   improve are bit-identical.
% • STATUS LAYER (when both catalogs carry .status): a base primary displaced
%   by a faster donor entry is kept as an ALTERNATIVE (status from the base's
%   status grid, source 'displaced primary (merge)'); the donor's alternatives
%   are remapped onto the base grid (off-grid cells -> 0) and merged with the
%   base's; all are deduplicated against the merged primaries. Catalogs
%   without the layer merge exactly as before.
%
%% Inputs:
%
%  base                     struct or char          phase catalog (or its .mat)
%  donor                    struct or char          phase catalog (or its .mat)
%  opts                     struct (optional)
%   .tolDays [1e-3] how much faster the donor must be, .donorTag ['donor']
%   the name written into the notes and the merge record
%
%% Outputs:
%
%  M                        struct                  the merged catalog, with
%                                                   .merge (what was taken)
%  info                     struct                  .nFaster .nFilled .nShared
%                                                   .rows [iD iA jD jA tfBase
%                                                   tfDonor filled] (ND, base
%                                                   then donor indices)
%
%% Revision History:
%  M. Casey                                                   (c) 10/03/2026
%  M. Casey  displaced primaries + alternatives merged                10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3, opts = struct(); end
base = asCatalog(base);  donor = asCatalog(donor);
tolDays = fieldd(opts, 'tolDays', 1e-3);  tag = fieldd(opts, 'donorTag', 'donor');

% ---- the same problem, or refused -----------------------------------------
assert(isscalar(base.sheets) && isscalar(donor.sheets), 'merge_phase_catalogs:sheets', 'one sheet each');
sb = base.sheets(1);  sd = donor.sheets(1);
same = @(a, b, f) isfield(a, f) && isfield(b, f) && isequal(a.(f), b.(f));
for f = {'thruster', 'rungs_N'}
    assert(same(base, donor, f{1}), 'merge_phase_catalogs:problem', '%s differs: not the same problem', f{1});
end
assert(base.constants.tStar_s == donor.constants.tStar_s, 'merge_phase_catalogs:problem', 'tStar differs');
for f = {'tauDRO', 'Np', 'pm'}
    assert(same(sb, sd, f{1}), 'merge_phase_catalogs:problem', 'sheet %s differs: not the same orbit pair', f{1});
end
assert(isequal(base.families.labels, donor.families.labels), 'merge_phase_catalogs:families', ...
       'family labels differ: a family code would mean another family');
for f = {'instruments', 'nSub', 'K', 'relTolPair'}
    assert(same(base.second_order, donor.second_order, f{1}), 'merge_phase_catalogs:secondOrder', ...
           'second_order.%s differs: the verdicts were not made the same way', f{1});
end

% ---- which fields travel with a cell, and which with an entry ------------
nD = numel(sb.sD_frac);  nA = numel(sb.sA_frac);  nCol = size(sb.z8, 2);
names = fieldnames(sb);
isGrid = cellfun(@(f) size(sb.(f), 1) == nD && size(sb.(f), 2) == nA, names);
gridF = setdiff(names(isGrid), {'entry_index'});
isEntry = ~isGrid & cellfun(@(f) size(sb.(f), 2) == nCol, names) & ~ismember(names, {'sD_frac', 'sA_frac'});
entryF = names(isEntry);
for f = [gridF; entryF]'
    assert(isfield(sd, f{1}), 'merge_phase_catalogs:fields', 'the donor lacks %s: a taken cell would carry a stale value', f{1});
end

assert(~(isfield(sd, 'status') && ~isfield(sb, 'status')), 'merge_phase_catalogs:layer', ...
       'the donor carries a status layer the base lacks: rebuild/backfill the base first');

% ---- the shared cells ------------------------------------------------------
onGrid = @(g, s) find(abs(mod(g - s + 0.5, 1) - 0.5) < 1e-8, 1);
mapD = arrayfun(@(s) emptyToZero(onGrid(sb.sD_frac, s)), sd.sD_frac);
mapA = arrayfun(@(s) emptyToZero(onGrid(sb.sA_frac, s)), sd.sA_frac);
nShared = nnz(mapD)*nnz(mapA);
assert(nShared > 0, 'merge_phase_catalogs:grid', 'the donor shares no cell with the base');
tolND = tolDays*86400/base.constants.tStar_s;

hasLayer = isfield(sb, 'status') && isfield(sd, 'status');
alts = struct([]);

% ---- take the donor's entry where it is faster or the base has none -------
rows = zeros(0, 7);
for jA = find(mapA)
    for jD = find(mapD)
        if ~sd.has_solution(jD, jA, 1), continue, end
        iD = mapD(jD);  iA = mapA(jA);
        filled = ~sb.has_solution(iD, iA, 1);
        tfB = sb.tf_nd(iD, iA, 1);  tfDn = sd.tf_nd(jD, jA, 1);
        if ~filled && ~(tfDn < tfB - tolND), continue, end
        kd = sd.entry_index(jD, jA, 1);
        if filled, kb = size(sb.z8, 2) + 1; else, kb = sb.entry_index(iD, iA, 1); end
        if ~filled && hasLayer                   % keep the displaced primary BEFORE it is overwritten
            kb0 = sb.entry_index(iD, iA, 1);
            Cp = struct('ok', true, 'stage', 8, 'status', double(sb.status(iD, iA, 1)), ...
                        'status_reason', sb.status_reason{kb0}, 'z', sb.z8(:, kb0), 'Y', sb.junctions{kb0}, ...
                        'flyKm', 0, 'flyVms', 0, 'overridden', false);
            alts = appendAlt(alts, make_alternative(Cp, sb.sD_frac(iD), sb.sA_frac(iA), iD, iA, 'displaced primary (merge)'));
        end
        for f = entryF', sb.(f{1})(:, kb) = sd.(f{1})(:, kd); end
        if filled, why = 'hole in the base'; else, why = sprintf('faster by %.4f d', (tfB - tfDn)*base.constants.tStar_s/86400); end
        sb.entry_notes{kb} = sprintf('%s | merged from %s: %s', sd.entry_notes{kd}, tag, why);
        for f = gridF', sb.(f{1})(iD, iA, :) = sd.(f{1})(jD, jA, :); end
        sb.entry_index(iD, iA, 1) = kb;
        rows(end+1, :) = [iD iA jD jA tfB tfDn filled];
    end
end

M = base;  M.sheets = sb;
M.n_entries = nnz(sb.has_solution);
if hasLayer
    allA = struct([]);
    if isfield(base, 'alternatives'), allA = appendAlt(allA, base.alternatives); end
    if isfield(donor, 'alternatives')
        dA = donor.alternatives;
        for ka = 1:numel(dA)
            a = dA(ka);
            if a.iD > 0, a.iD = mapD(a.iD); end
            if a.iA > 0, a.iA = mapA(a.iA); end
            allA = appendAlt(allA, a);
        end
    end
    allA = appendAlt(allA, alts);
    M.alternatives = dedup_alternatives(allA, M.sheets);
    if isfield(base, 'status_key'), M.status_key = base.status_key; end
end
info = struct('nFaster', nnz(rows(:, 7) == 0), 'nFilled', nnz(rows(:, 7) == 1), 'nShared', nShared, 'rows', rows);
M.merge = struct('donor', tag, 'date', char(datetime('now', 'Format', 'yyyy-MM-dd')), 'tolDays', tolDays, ...
                 'nFaster', info.nFaster, 'nFilled', info.nFilled, 'rows', rows, ...
                 'columns', 'iD iA (base) jD jA (donor) tfBase tfDonor (ND) filled', ...
                 'note', 'taken entries carry the DONOR''s verdicts; audit the merged file before relying on it');
end

% ---------------------------------------------------------------------------
function A = appendAlt(A, B)
% APPENDALT  Concatenate alternative rows, skipping empties; rows built here
% get .sheet = 1 so the field sets agree.  INPUTS: A; B.  OUTPUTS: A.
if isempty(B), return, end
if ~isfield(B, 'sheet'), for kb = 1:numel(B), B(kb).sheet = 1; end, end
if isempty(A), A = B; return, end
if ~isfield(A, 'sheet'), for ka = 1:numel(A), A(ka).sheet = 1; end, end
A = [A(:); B(:)]';
end

function c = asCatalog(c)
% ASCATALOG  A catalog struct from a struct or a .mat path.
% INPUTS: c (struct or char).  OUTPUTS: c (struct).
if ischar(c) || isstring(c), L = load(c);  c = L.(char(fieldnames(L))); end
end

function v = emptyToZero(v)
% EMPTYTOZERO  [] -> 0.  INPUTS: v.  OUTPUTS: v.
if isempty(v), v = 0; end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
