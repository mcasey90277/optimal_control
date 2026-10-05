function [c, info] = relabel_borderline(c, altRows)
%% Purpose:
%
%   Label CONSERVATIVELY the alternatives whose optimality status flipped
%   between runs because a check sits at its threshold (a lift margin of
%   9 against a gate of 10, a witness timing out under load). For every
%   audit row (audit_status_layer .altRows) that is BAD, whose root did NOT
%   move, and whose re-audit status is a finite code >= 1, the alternative
%   takes status = min(stored, re-audit), borderline = true, and a
%   status_reason naming both statuses, the audit's why and the old reason.
%   A borderline status is a LOWER BOUND: the audit then accepts any
%   re-certification of the same root at that status or higher.
%
%  ASSUMPTIONS / NOTES:
%
% • altRows(m).k is the GLOBAL alternative index (audit_status_layer
%   iterates k over the table, chunked or not). The caller binds the rows to
%   this catalog (content keys); an index outside the table is refused.
% • A MOVED row is never relabelled (spec 7, Review Focus 5: another root is
%   not this alternative), nor is a row with no finite re-audit status >= 1
%   (below the floor, a throw): both are listed in info.notRelabelled.
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
    if isequal(r.moved, true)
        notRel(end+1) = struct('k', k, 'why', ['the root moved: ' r.why]);
        continue
    end
    if ~(isnumeric(sNow) && isscalar(sNow) && isreal(sNow) && isfinite(sNow) && sNow >= 1)
        notRel(end+1) = struct('k', k, 'why', sprintf('no finite re-audit status >= 1 (%g): %s', sNow, r.why));
        continue
    end
    stored = double(A(k).status);
    newSt = min(stored, double(sNow));
    A(k).status_reason = sprintf('borderline: stored %d, re-audit %d -- %s | earlier: %s', stored, sNow, r.why, A(k).status_reason);
    A(k).status = newSt;
    A(k).borderline = true;
    rows(end+1, :) = [k, stored, sNow, newSt];
end
c.alternatives = A;
c.status_layer.nBorderline = nnz([A.borderline]);
c.status_layer.borderlineNote = ['a borderline status (alternatives.borderline) is a lower bound: the status flipped ' ...
                                 'between runs at a threshold, the lower one is kept, and the audit accepts the same ' ...
                                 'root at that status or higher'];
info = struct('rows', rows);
info.notRelabelled = notRel;
end
