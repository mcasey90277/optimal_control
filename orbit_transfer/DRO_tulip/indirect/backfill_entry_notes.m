function out = backfill_entry_notes(catMat, opts)
%% Purpose:
%
%   Write provenance notes into a catalog packaged BEFORE entries carried
%   them (the 576-entry 70 mN library of 2026-09-15): for every solved
%   entry, find the record that produced it -- the sheet's winning
%   candidate at the spine, a rib point in one of the rib files the
%   receipt names, or a direct-solved cell -- and compose the note the
%   producers now write live: the family label, how it was found, the
%   certifier's remarks read off the point's own verdict fields, and the
%   sweep's remark from the keyed sidecar. The catalog is backed up, the
%   notes and notes_key are written in, the schema is re-validated.
%
%% Inputs:
%
%  catMat                   char                    the packaged catalog
%  opts                     struct (optional)
%   .sheetMat  [<catalog dir>/arrival_sheet_*.mat]  the sD = 0 sheet
%   .receipt   [<catalog dir>/catalog_receipt.mat]  names the rib files
%   .holesMat  [<catalog dir>/fine_rib_direct_holes.mat] the direct cells
%   .sidecar   [<catalog dir>/second_order_progress_v3.mat] the sweep
%   .dryRun    [false]      compose but do not write
%
%% Outputs:
%
%  out                      struct                  .nEntries .nSpine .nRib
%                                                   .nDirect .nUnknown
%                                                   .notes (cellstr)
%
%% Revision History:
%  M. Casey                                                   (c) 09/15/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
d = @(f, v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
[catDir, ~] = fileparts(catMat);
sheetMat = d('sheetMat', firstMatch(fullfile(catDir, 'arrival_sheet_*.mat')));
receiptF = d('receipt', fullfile(catDir, 'catalog_receipt.mat'));
holesMat = d('holesMat', fullfile(catDir, 'fine_rib_direct_holes.mat'));
sidecarF = d('sidecar', fullfile(catDir, 'second_order_progress_v3.mat'));
dryRun = d('dryRun', false);

L = load(catMat);  fn = fieldnames(L);  cat_ = L.(fn{1});
sh = cat_.sheets(1);
has = sh.has_solution(:, :, 1);  idx = sh.entry_index(:, :, 1);  z8 = sh.z8;
sD = sh.sD_frac(:).';  sA = sh.sA_frac(:).';
tStar = cat_.constants.tStar_s;
nEnt = size(z8, 2);
notes = repmat({''}, 1, nEnt);
src = zeros(1, nEnt);                      % 1 spine, 2 rib, 3 direct, 0 unknown

% ---- the family label ----------------------------------------------------
labelOf = @(code) '';
if isfield(sh, 'family_index') && isfield(cat_, 'families')
    fams = cat_.families.families;
    labelOf = @(code) famLabel(code, fams);
end

% ---- the sheet: spine entries -------------------------------------------
S = load(sheetMat);  S = S.S;
spineRow = find(abs(mod(sD - S.problem.sD + 0.5, 1) - 0.5) < 1e-8, 1);
for j = 1:numel(S.sA)
    iA = find(abs(mod(sA - S.sA(j) + 0.5, 1) - 0.5) < 1e-8, 1);
    if isempty(iA) || ~has(spineRow, iA), continue, end
    c = S.cand{j};
    k = find([c.ok] & abs([c.tfDays] - S.TF(j)) < 1e-9, 1);
    if isempty(k), continue, end
    e = idx(spineRow, iA);
    if norm(z8(:, e) - c(k).z(:)) > 1e-9*max(1, norm(z8(:, e))), continue, end   % the spine was displaced by a rib
    if c(k).arc > 0, how = sprintf('arc crossing: arc %d (%s), level %.4f', c(k).arc, arcName(S, c(k).arc), c(k).level);
    else, how = 'seed: a certified root already in the library (direct-found or an earlier round)';
    end
    notes{e} = compose(labelOf(famCode(sh, spineRow, iA)), how, c(k));  src(e) = 1;
end

% ---- the ribs the receipt names ------------------------------------------
if isfile(receiptF)
    rc = load(receiptF);
    for f = rc.ribFiles(:).'
        if ~isfile(f{1}) || contains(f{1}, 'direct_holes'), continue, end
        Rk = load(f{1});
        [~, roundName] = fileparts(fileparts(f{1}));
        for m = 1:numel(Rk.R)
            r = Rk.R(m);
            for q = 1:numel(r.pts)
                p = r.pts(q);
                iD = find(abs(mod(sD - p.sD + 0.5, 1) - 0.5) < 1e-8, 1);  iA = find(abs(mod(sA - p.sA + 0.5, 1) - 0.5) < 1e-8, 1);
                if isempty(iD) || isempty(iA) || ~has(iD, iA), continue, end
                e = idx(iD, iA);
                if src(e) ~= 0 || norm(z8(:, e) - p.z(:)) > 1e-9*max(1, norm(z8(:, e))), continue, end
                how = sprintf('rib step %d of %d off the sD = %.4f spine at sA %.4f (%s, %s)', q, numel(r.pts), ...
                              S.problem.sD, r.sA, roundName, fileName(f{1}));
                notes{e} = compose(labelOf(famCode(sh, iD, iA)), how, p);  src(e) = 2;
            end
        end
    end
end

% ---- the direct-solved cells --------------------------------------------
if isfile(holesMat)
    H = load(holesMat);
    for m = 1:numel(H.R)
        for q = 1:numel(H.R(m).pts)
            p = H.R(m).pts(q);
            iD = find(abs(mod(sD - p.sD + 0.5, 1) - 0.5) < 1e-8, 1);  iA = find(abs(mod(sA - p.sA + 0.5, 1) - 0.5) < 1e-8, 1);
            if isempty(iD) || isempty(iA) || ~has(iD, iA), continue, end
            e = idx(iD, iA);
            if norm(z8(:, e) - p.z(:)) > 1e-9*max(1, norm(z8(:, e))), continue, end
            how = 'direct cell solve seeded from a certified neighbour';
            if isfield(H, 'rec')
                rr = H.rec([H.rec.iD] == iD & [H.rec.iA] == iA & [H.rec.ok]);
                if ~isempty(rr)
                    rr = rr(end);
                    how = sprintf('direct cell solve seeded from (%.4f, %.4f)', rr.seed);
                    if contains(rr.reason, 'improve: cell holds')
                        t = regexp(rr.reason, 'improve: cell holds ([\d.]+) d', 'tokens', 'once');
                        if ~isempty(t), how = sprintf('%s, replaced %s d', how, t{1}); end
                    end
                    if contains(rr.reason, 'direct solve failed') || contains(rr.reason, 'clearance floor')
                        how = [how ' (another seed diverged)'];
                    end
                end
            end
            notes{e} = compose(labelOf(famCode(sh, iD, iA)), how, p);  src(e) = 3;
        end
    end
end

% ---- the sweep's remark, from the keyed sidecar --------------------------
minRel = nan(size(has));                   % the closest near-miss per entry, for the grid
if isfile(sidecarF)
    P = load(sidecarF);  R = P.R;
    for q = 1:numel(R)
        if ~isfield(R(q), 'key') || isempty(R(q).key) || ~R(q).done, continue, end
        e = idx(R(q).key(1), R(q).key(2));
        if e == 0 || norm(z8(:, e) - R(q).z8(:)) > 1e-9*max(1, norm(z8(:, e))), continue, end
        nu = 0;  if isfield(R(q), 'nUnresolved') && ~isempty(R(q).nUnresolved), nu = R(q).nUnresolved; end
        if nu > 0
            notes{e} = strjoin([notes(e), {sprintf('sweep: %d unresolved conjugate candidate(s)', nu)}], ' | ');
            if startsWith(notes{e}, ' | '), notes{e} = notes{e}(4:end); end
        end
        if isfield(R(q), 'minRelSigma') && ~isempty(R(q).minRelSigma), minRel(R(q).key(1), R(q).key(2)) = R(q).minRelSigma; end
    end
end

for e = find(src == 0)
    [iD, iA] = find(idx == e);
    notes{e} = compose(labelOf(famCode(sh, iD, iA)), 'origin not found among the sheet, the ribs and the direct cells', struct());
end
out = struct('nEntries', nEnt, 'nSpine', nnz(src == 1), 'nRib', nnz(src == 2), 'nDirect', nnz(src == 3), ...
             'nUnknown', nnz(src == 0), 'notes', {notes});
fprintf('backfill_entry_notes: %d entries: %d spine, %d rib, %d direct, %d unknown\n', ...
        nEnt, out.nSpine, out.nRib, out.nDirect, out.nUnknown);
if dryRun, return, end

cat_.sheets(1).entry_notes = notes;
if any(isfinite(minRel(:))), cat_.sheets(1).conj_min_rel = minRel; end   % closest near-miss, x median sigma_6
cat_.notes_key = struct('layout', '<family> | <how found> | <certifier remarks> | <sweep remarks>', ...
    'how_found', {{'arc crossing: arc k (<file>), level L', 'seed: ...', ...
                   'rib step k of n off the sD = .. spine at sA .. (<round>, <file>)', ...
                   'direct cell solve seeded from (sD, sA) [, replaced t_f d]'}}, ...
    'remarks', {{'polish plateaued |R|=.. > tolR ..', 'lift margin ..x (gate 10x)', 'H6 margin ..x', ...
                 'sweep: k unresolved conjugate candidate(s)'}}, ...
    'note', ['Backfilled by backfill_entry_notes on ' char(datetime('now')) ' from the sheet, the rib files and the direct cells; ' ...
             'the numbers that gate decisions live in the numeric grids (lift_margin, h6_margin, conj_*).']);
p = catalog_schema('validate', cat_);
assert(isempty(p), 'backfill_entry_notes: the catalog does not validate after the notes: %s', strjoin(p, '; '));
bak = [catMat '.bak_notes'];
if ~isfile(bak), copyfile(catMat, bak); end
L.(fn{1}) = cat_;
save(catMat, '-struct', 'L');
fprintf('backfill_entry_notes: written to %s (backup %s)\n', catMat, bak);
end

% ==========================================================================
function s = compose(label, how, C)
% COMPOSE  family | how found | certifier remarks (from the point's own
% verdict fields).  INPUTS: label; how; C.  OUTPUTS: s.
parts = {};
if ~isempty(label), parts{end+1} = label; end
parts{end+1} = how;
if isfield(C, 'reason') && contains(C.reason, 'plateaued')
    t = regexp(C.reason, 'plateaued at \|R\| = ([^,]+), above tolR = ([^)]+)', 'tokens', 'once');
    if ~isempty(t), parts{end+1} = sprintf('polish plateaued |R|=%s > tolR %s', t{1}, t{2}); end
end
if isfield(C, 'liftMargin') && isnumeric(C.liftMargin) && isfinite(C.liftMargin) && C.liftMargin < 20
    parts{end+1} = sprintf('lift margin %.1fx (gate 10x)', C.liftMargin);
end
if isfield(C, 'h6Margin') && isnumeric(C.h6Margin) && isfinite(C.h6Margin) && C.h6Margin < 3
    parts{end+1} = sprintf('H6 margin %.2fx', C.h6Margin);
end
s = strjoin(parts, ' | ');
end

function c = famCode(sh, iD, iA)
% FAMCODE  The entry's family_index code, 0 without a map.  INPUTS: sh; iD;
% iA.  OUTPUTS: c.
c = 0;
if isfield(sh, 'family_index'), c = double(sh.family_index(iD, iA, 1)); end
end

function s = famLabel(code, fams)
% FAMLABEL  The family clause for a code.  INPUTS: code; fams.  OUTPUTS: s.
if code >= 1, s = ['family ' fams(code).label];
elseif code == -1, s = 'family: unattached root (no arc passes through it)';
elseif code == -2, s = 'family: spine unidentified';
else, s = '';
end
end

function n = arcName(S, ia)
% ARCNAME  The arc file the sheet's candidate came from.  INPUTS: S; ia.
% OUTPUTS: n.
n = '';
if isfield(S, 'arcs') && ia <= numel(S.arcs), n = S.arcs{ia}; end
end

function f = firstMatch(pat)
% FIRSTMATCH  First file matching a pattern, '' if none.  INPUTS: pat.
% OUTPUTS: f.
q = dir(pat);  f = '';
if ~isempty(q), f = fullfile(q(1).folder, q(1).name); end
end

function n = fileName(p)
% FILENAME  Name + extension.  INPUTS: p.  OUTPUTS: n.
[~, a, b] = fileparts(p);  n = [a b];
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
