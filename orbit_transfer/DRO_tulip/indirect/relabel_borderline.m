function [c, info] = relabel_borderline(c, altRows)
%% Purpose:
%
%   Label CONSERVATIVELY the alternatives whose optimality status flipped
%   between runs because a check sits at its threshold (a lift margin of
%   9 against a gate of 10, a witness timing out under load). For every
%   audit row (audit_status_layer .altRows) that is a status flip on the
%   SAME root, the alternative takes status = status_meet(stored, re-audit)
%   -- the claim BOTH runs support on the lattice ordered by what was
%   established (1 < 3 < 2, 1 < 3 < 4; 2 and 4 incomparable), NOT the
%   numeric minimum -- with borderline = true and a status_reason naming
%   both statuses, the result, the audit's why and the old reason. A
%   borderline status is a LOWER BOUND: the audit then accepts any
%   re-certification of the same root whose meet with it is it.
%
%  ASSUMPTIONS / NOTES:
%
% • altRows(m).k is the GLOBAL alternative index (audit_status_layer
%   iterates k over the table, chunked or not). The caller binds the rows to
%   this catalog (content keys); an index outside the table is refused.
% • ONLY a status flip on the SAME root is relabelled: the audit's why must
%   start with 'status not reproduced' or 'borderline status not reached',
%   the root must not have moved (spec 7, Review Focus 5), and both the
%   stored and the re-audit status must be integer codes 1..4. Every other
%   BAD row (moved, below the floor, flight, throw, malformed) is listed in
%   info.notRelabelled with why.
% • Every alternative gets a logical scalar .borderline (false unless
%   relabelled here, an existing true is kept) so the struct array stays
%   uniform. OK rows and rows the audit never saw are otherwise untouched.
% • Pure: no I/O.
%
%% Inputs:
%
%  c                        struct                  catalog with .alternatives
%                                                   and .status_layer
%
%  altRows                  struct array            audit rows (.k .ok .why
%                                                   .statusNow .moved)
%
%% Outputs:
%
%  c                        struct                  .alternatives with
%                                                   .borderline on every row;
%                                                   .status_layer.nBorderline
%                                                   .borderlineNote
%
%  info                     struct                  .rows [n x 4] (k, stored,
%                                                   now, newStatus) .notRelabelled
%                                                   (.k .why)
%
%% Revision History:
%  M. Casey                                                   (c) 10/05/2026
%  M. Casey  fix round 1: same-root flips only (why prefix), codes
%            1..4 for stored and re-audit status                     10/05/2026
%  M. Casey  fix round 2: lattice meet (status_meet), not min    10/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

A = c.alternatives;
for k = 1:numel(A)
    if ~isfield(A, 'borderline') || ~isequal(A(k).borderline, true), A(k).borderline = false; end
end
rows = zeros(0, 4);
notRel = struct('k', {}, 'why', {});
for m = 1:numel(altRows)
    r = altRows(m);
    if isequal(r.ok, true), continue, end
    k = r.k;
    assert(isscalar(k) && k >= 1 && k <= numel(A) && k == round(k), 'relabel_borderline:index', ...
           'audit row %d names alternative %g, outside the table (1..%d)', m, k, numel(A));
    sNow = r.statusNow;
    why = char(r.why);
    stored = A(k).status;
    if isequal(r.moved, true)
        notRel(end+1) = struct('k', k, 'why', ['the root moved: ' why]);
        continue
    end
    if ~(startsWith(why, 'status not reproduced') || startsWith(why, 'borderline status not reached'))
        notRel(end+1) = struct('k', k, 'why', ['not a status flip on the same root: ' why]);
        continue
    end
    if ~isCode(sNow)
        notRel(end+1) = struct('k', k, 'why', sprintf('re-audit status is not an integer code 1..4 (%s): %s', num2str(sNow), why));
        continue
    end
    if ~isCode(stored)
        notRel(end+1) = struct('k', k, 'why', sprintf('stored status is not an integer code 1..4 (%s): %s', num2str(stored), why));
        continue
    end
    stored = double(stored);  sNow = double(sNow);
    newSt = status_meet(stored, sNow);                       % the common ground of the two runs (lattice, not min)
    A(k).status_reason = sprintf('borderline: stored %d, re-audit %d -> %d (the common ground of the two runs) -- %s | earlier: %s', ...
                                 stored, sNow, newSt, why, A(k).status_reason);
    A(k).status = newSt;
    A(k).borderline = true;
    rows(end+1, :) = [k, stored, sNow, newSt];
end
c.alternatives = A;
c.status_layer.nBorderline = nnz([A.borderline]);
c.status_layer.borderlineNote = ['a borderline status (alternatives.borderline) is a lower bound: the status flipped ' ...
                                 'between runs at a threshold, the claim both runs support is kept (status_meet on ' ...
                                 'the lattice 1 < 3 < 2, 1 < 3 < 4), and the audit accepts the same root at any status ' ...
                                 'whose meet with it is it'];
info = struct('rows', rows);
info.notRelabelled = notRel;
end

% ---------------------------------------------------------------------------
function tf = isCode(v)
% ISCODE  v is a real finite integer status code in 1..4.
% INPUTS: v (any).  OUTPUTS: tf [logical].
tf = isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v == round(v) && v >= 1 && v <= 4;
end
