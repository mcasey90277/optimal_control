function R = compare_phase_catalogs(newCat, refCat, opts)
%% Purpose:
%
%   Is a rebuilt phase catalog THE SAME LIBRARY as a reference one? Cell by
%   cell over the (departure phase, arrival phase) grid it compares
%
%     coverage   which cells hold a certified transfer
%     t_f        the minimum flight time
%     costates   the solution itself: the seven initial costates (t_f is
%                compared on its own line, so it is left out of this norm)
%     family     which continuation family owns the cell
%
%   and prints a verdict. It is the acceptance test of a reproduction
%   (reproduce_library_70mN, section 4), and it reads nothing but the two
%   catalog files.
%
%  ASSUMPTIONS / NOTES:
%
% • A rebuilt library is not byte-identical to its reference: every entry is
%   re-polished, so t_f and z8 agree to solver tolerance, not to the last
%   bit. The tolerances say how close is "the same root".
% • FAMILIES ARE COMPARED AS A PARTITION, not by index number. The family
%   index of a cell depends on the order the arcs were listed, so the same
%   five families can come back numbered differently. What must agree is
%   which cells belong together. The families of the two catalogs are paired
%   off, largest overlap first; a cell outside its family's pairing is a
%   difference.
% • The two catalogs must be on the same grid; another grid is refused.
%
%% Inputs:
%
%  newCat, refCat           char | struct           catalog .mat file, or the
%                                                   catalog struct (one sheet)
%  opts                     struct (optional)
%   .tolTfDays              double                  t_f agreement, days [1e-6]
%   .tolZRel                double                  |dz8| / |z8| [1e-6]
%   .print                  logical                 print the report [true]
%   .maxList                int                     differing cells listed [20]
%
%% Outputs:
%
%  R                        struct                  .ok (a match) .nRef .nNew
%                                                   .nBoth .nOnlyRef .nOnlyNew
%                                                   .nNonfinite .nTfOver
%                                                   .nZOver .nFamilyDiff
%                                                   .familiesCompared .nAgree
%                                                   .worstTfDays
%                                                   .worstZRel .cells (table:
%                                                   iD iA sD sA what)
%
%% Revision History:
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3, opts = struct(); end
tolTfDays = fieldd(opts, 'tolTfDays', 1e-6);
tolZRel   = fieldd(opts, 'tolZRel', 1e-6);
doPrint   = fieldd(opts, 'print', true);
maxList   = fieldd(opts, 'maxList', 20);

%% 1. The two sheets: ONE problem, on ONE grid:
[sn, tStarN, idN] = oneSheet(newCat);
[sr, tStar,  idR] = oneSheet(refCat);
% identical arrays under another engine, orbit or time scale are not the same
% library: the problems must be the same before a single cell is compared
sameProblem = all(abs(idN - idR) <= 1e-9*max(1, abs(idR))) && abs(tStarN - tStar) <= 1e-9*tStar;
assert(sameProblem, 'compare_phase_catalogs:problem', ...
       ['the two catalogs are libraries of DIFFERENT problems ' ...
        '(thrust N, Isp s, m0 kg, DRO tau, tulip Np, branch: new %s, reference %s)'], mat2str(idN, 6), mat2str(idR, 6));
circ = @(x) abs(mod(x + 0.5, 1) - 0.5);                     % circular phase distance
sameGrid = numel(sn.sD_frac) == numel(sr.sD_frac) && numel(sn.sA_frac) == numel(sr.sA_frac) && ...
           max(circ(sn.sD_frac(:) - sr.sD_frac(:))) < 1e-9 && max(circ(sn.sA_frac(:) - sr.sA_frac(:))) < 1e-9;
assert(sameGrid, 'compare_phase_catalogs:grid', 'the two catalogs are on different phase grids; nothing to compare');

%% 2. Coverage: which cells each catalog holds:
hasN = logical(sn.has_solution(:, :, 1));   hasR = logical(sr.has_solution(:, :, 1));
both = hasN & hasR;
[iDs, iAs] = find(both);  nBoth = numel(iDs);

%% 3. Per common cell: flight time, and the costates on their own:
% t_f is compared by itself, so it is left OUT of the costate norm -- in one
% 8-vector a flight time of order 4 dilutes a change in costates of order 1.
dTf = zeros(nBoth, 1);  dZ = zeros(nBoth, 1);
for q = 1:nBoth
    iD = iDs(q);  iA = iAs(q);
    dTf(q) = abs(sn.tf_nd(iD, iA, 1) - sr.tf_nd(iD, iA, 1))*tStar/86400;
    zn = sn.z8(1:7, sn.entry_index(iD, iA, 1));   zr = sr.z8(1:7, sr.entry_index(iD, iA, 1));
    dZ(q) = sqrt(sum((zn - zr).^2))/max(sqrt(sum(zr.^2)), realmin);
end
% A NON-FINITE ENTRY IS A DIFFERENCE. `NaN > tol` is false in MATLAB, so a
% NaN flight time or costate used to compare EQUAL.
nonfinite = ~isfinite(dTf) | ~isfinite(dZ);
tfOver = dTf > tolTfDays & ~nonfinite;   zOver = dZ > tolZRel & ~nonfinite;

%% 4. Families, as a partition of the common cells:
% Without a family map in BOTH catalogs the comparison is NOT ESTABLISHED --
% which is reported as such and is not a match; it is never "0 differences".
famDiff = false(nBoth, 1);
familiesCompared = isfield(sn, 'family_index') && isfield(sr, 'family_index');
if familiesCompared && nBoth > 0
    fN = double(sn.family_index(:, :, 1));  fR = double(sr.family_index(:, :, 1));
    famDiff = ~inPairedFamily(fR(both), fN(both));
end

%% 5. The differing cells, named:
what = strings(0, 1);  ciD = zeros(0, 1);  ciA = zeros(0, 1);
[a, b] = find(hasR & ~hasN);   [ciD, ciA, what] = addCells(ciD, ciA, what, a, b, "missing from the new catalog");
[a, b] = find(hasN & ~hasR);   [ciD, ciA, what] = addCells(ciD, ciA, what, a, b, "only in the new catalog");
[ciD, ciA, what] = addCells(ciD, ciA, what, iDs(nonfinite), iAs(nonfinite), "a non-finite flight time or costate");
[ciD, ciA, what] = addCells(ciD, ciA, what, iDs(tfOver), iAs(tfOver), compose("t_f differs by %.2e d", dTf(tfOver)));
[ciD, ciA, what] = addCells(ciD, ciA, what, iDs(zOver & ~tfOver), iAs(zOver & ~tfOver), compose("z8 differs by %.2e (relative)", dZ(zOver & ~tfOver)));
[ciD, ciA, what] = addCells(ciD, ciA, what, iDs(famDiff), iAs(famDiff), "owned by a different family");
sDcol = reshape(sr.sD_frac(ciD), [], 1);   sAcol = reshape(sr.sA_frac(ciA), [], 1);     % columns, whatever the grids' orientation
cells = table(ciD, ciA, sDcol, sAcol, what, 'VariableNames', {'iD', 'iA', 'sD', 'sA', 'what'});

%% 6. The verdict:
R = struct('nRef', nnz(hasR), 'nNew', nnz(hasN), 'nBoth', nBoth, ...
           'nOnlyRef', nnz(hasR & ~hasN), 'nOnlyNew', nnz(hasN & ~hasR), ...
           'nNonfinite', nnz(nonfinite), 'nTfOver', nnz(tfOver), 'nZOver', nnz(zOver), ...
           'nFamilyDiff', nnz(famDiff), 'familiesCompared', familiesCompared, ...
           'nAgree', nnz(~(nonfinite | tfOver | zOver | famDiff)), ...        % common CELLS with no difference at all
           'worstTfDays', max([dTf(isfinite(dTf)); 0]), 'worstZRel', max([dZ(isfinite(dZ)); 0]), ...
           'tolTfDays', tolTfDays, 'tolZRel', tolZRel, 'cells', cells);
R.ok = R.nOnlyRef == 0 && R.nOnlyNew == 0 && R.nNonfinite == 0 && R.nTfOver == 0 && R.nZOver == 0 && ...
       R.nFamilyDiff == 0 && R.familiesCompared;
if doPrint, printReport(R, maxList); end
end

% ------------------------------------------------------------------------
function [s, tStar, id] = oneSheet(c)
% ONESHEET  The single sheet of a phase catalog, from a file or a struct, with
% the numbers that say WHICH problem it is a library of.
% INPUTS: c (char | struct).  OUTPUTS: s (sheet struct); tStar (s per ND);
% id [1 x 6] thrust N, Isp s, m0 kg, DRO tau, tulip Np, branch.
if ischar(c) || isstring(c)
    assert(isfile(c), 'compare_phase_catalogs:file', 'no such catalog: %s', c);
    L = load(c);  fn = fieldnames(L);  c = L.(fn{1});
end
assert(isstruct(c) && isfield(c, 'sheets') && isscalar(c.sheets), 'compare_phase_catalogs:catalog', ...
       'a phase catalog with exactly one sheet is expected');
s = c.sheets(1);  tStar = c.constants.tStar_s;
id = [c.rungs_N(1), c.thruster.isp_s, c.thruster.m0_kg, s.tauDRO, s.Np, s.pm];
end

% ------------------------------------------------------------------------
function paired = inPairedFamily(fRef, fNew)
% INPAIREDFAMILY  Pair the reference families with the new ones, largest
% overlap first, each family used once; a cell is `paired` when its two
% indices are a matched pair. Index NUMBERS never matter, only which cells
% sit together.  INPUTS: fRef, fNew [n x 1] family index per common cell.
% OUTPUTS: paired [n x 1 logical].
[uR, ~, kR] = unique(fRef);  [uN, ~, kN] = unique(fNew);
overlap = accumarray([kR, kN], 1, [numel(uR), numel(uN)]);     % cells shared by (ref family, new family)
partner = zeros(numel(uR), 1);                                 % new-family number paired with each ref family
work = overlap;
while any(work(:) > 0)
    [~, big] = max(work(:));  [r, n] = ind2sub(size(work), big);
    partner(r) = n;
    work(r, :) = 0;  work(:, n) = 0;                           % each family is used once
end
paired = partner(kR) == kN;
end

% ------------------------------------------------------------------------
function [ciD, ciA, what] = addCells(ciD, ciA, what, iD, iA, label)
% ADDCELLS  Append cells and their label(s) to the list of differences.
% INPUTS: the three lists; iD, iA (indices); label (one string, or one per
% cell).  OUTPUTS: the three lists.
n = numel(iD);
if n == 0, return, end
if isscalar(label), label = repmat(label, n, 1); end
ciD = [ciD; iD(:)];  ciA = [ciA; iA(:)];  what = [what; label(:)];
end

% ------------------------------------------------------------------------
function printReport(R, maxList)
% PRINTREPORT  The comparison as text: one line per check with PASS/FAIL,
% then the differing cells.  INPUTS: R; maxList.  OUTPUTS: none.
pf = @passFail;
fprintf('CATALOG COMPARISON (new against reference)\n');
fprintf('  coverage  : reference %d cells, new %d, common %d; %d missing, %d extra            %s\n', ...
        R.nRef, R.nNew, R.nBoth, R.nOnlyRef, R.nOnlyNew, pf(R.nOnlyRef == 0 && R.nOnlyNew == 0));
fprintf('  t_f       : worst deviation %.3e d (tolerance %.0e d); %d cell(s) over      %s\n', ...
        R.worstTfDays, R.tolTfDays, R.nTfOver, pf(R.nTfOver == 0));
fprintf('  costates  : worst relative deviation %.3e (tolerance %.0e); %d cell(s) over   %s\n', ...
        R.worstZRel, R.tolZRel, R.nZOver, pf(R.nZOver == 0));
fprintf('  finite    : %d cell(s) hold a non-finite flight time or costate                  %s\n', R.nNonfinite, pf(R.nNonfinite == 0));
if R.familiesCompared
    fprintf('  families  : %d cell(s) owned by a different family (compared as a partition)   %s\n', R.nFamilyDiff, pf(R.nFamilyDiff == 0));
else
    fprintf('  families  : NOT ESTABLISHED -- a catalog carries no family map                  FAIL\n');
end
fprintf('  agreement : %d of %d common cells agree in everything\n', R.nAgree, R.nBoth);
fprintf('  VERDICT   : %s\n', pf(R.ok));
n = height(R.cells);
for q = 1:min(n, maxList)
    fprintf('    cell (%2d,%2d)  sD %.4f  sA %.4f : %s\n', R.cells.iD(q), R.cells.iA(q), R.cells.sD(q), R.cells.sA(q), R.cells.what(q));
end
if n > maxList, fprintf('    ... and %d more (R.cells holds them all)\n', n - maxList); end
end

% ------------------------------------------------------------------------
function t = passFail(c)
% PASSFAIL  'PASS' or 'FAIL'.  INPUTS: c (logical).  OUTPUTS: t (char).
if c, t = 'PASS'; else, t = 'FAIL'; end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
