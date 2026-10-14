function [prim, alt2] = sample_status_pair(cat)
%% Purpose:
%
%   Pick the deliverable's sample pair: the on-grid cell of the library
%   that records the most status-2 alternatives (ties broken by the
%   smallest flight-time gap between the cell's primary and its fastest
%   status-2 alternative), and return that cell's primary row and its
%   fastest status-2 row, both as transfers_at rows.
%
%  ASSUMPTIONS / NOTES:
%
% • Single-sheet catalogs with an alternatives layer (fields iD, iA, status).
% • Off-grid alternatives (iD or iA = 0) are not cells and are not counted.
% • A cell with no primary is skipped. No such cell -> a named error.
%
%% Inputs:
%
%  cat                      struct                  single-sheet costate catalog
%
%% Outputs:
%
%  prim                     struct                  transfers_at row (kind 'primary')
%  alt2                     struct                  transfers_at row (kind 'alternative',
%                                                   status 2, smallest t_f in the cell)
%
%% Revision History:
%  M. Casey                                                   (c) 10/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~isfield(cat, 'alternatives') || isempty(cat.alternatives)
    error('sample_status_pair:noAlternatives', 'the catalog records no alternatives');
end
A  = cat.alternatives;
sh = cat.sheets(1);
sel = find([A.status] == 2 & [A.iD] > 0 & [A.iA] > 0);
if isempty(sel)
    error('sample_status_pair:none', 'no on-grid cell records a status-2 alternative');
end
cells = unique([[A(sel).iD].', [A(sel).iA].'], 'rows');
best = struct('n', -1, 'gap', Inf, 'prim', [], 'alt2', []);
for kc = 1:size(cells, 1)
    T = transfers_at(cat, sh.sD_frac(cells(kc, 1)), sh.sA_frac(cells(kc, 2)));
    ip = find(strcmp({T.kind}, 'primary'), 1);
    ia = find(strcmp({T.kind}, 'alternative') & [T.status] == 2);
    if isempty(ip) || isempty(ia), continue; end
    n   = numel(ia);
    gap = abs(T(ia(1)).tfDays - T(ip).tfDays);          % T is sorted by t_f: ia(1) is the fastest
    if n > best.n || (n == best.n && gap < best.gap)
        best = struct('n', n, 'gap', gap, 'prim', T(ip), 'alt2', T(ia(1)));
    end
end
if isempty(best.prim)
    error('sample_status_pair:none', 'no on-grid cell has both a primary and a status-2 alternative');
end
prim = best.prim;
alt2 = best.alt2;
end
