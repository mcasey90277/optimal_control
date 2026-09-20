function [branch, nBranch] = phase_components(safeD, safeA, has)
%% Purpose:
%
%   CANDIDATE COMPONENTS OF A PHASE SHEET. Cells joined by a chain of accepted
%   edges (for instance time-consistent edges, phase_edge_residuals) form a
%   connected component of the grid graph; this labels them, on the torus
%   (both phase axes wrap). It is graph bookkeeping: a component is only as
%   meaningful as the edge test behind it, it is NOT an established solution
%   branch, and two members of one component need not be interpolable (the
%   component may wind round the torus or follow a fold). Renamed from
%   phase_branches after review, 2026-09-19.
%
%   It replaces the family label for this purpose. A family label says which
%   continuation arc an entry was FOUND on; it attaches by flight time, is
%   not reproducible between builds, and was measured to hide 33 jumps inside
%   one label (FINDINGS 80, 82). A branch is computed from the entries'
%   own flight times and costates, so two builds of the same library give
%   the same branches.
%
%  ASSUMPTIONS / NOTES:
%
% • safeD(i,j) is the edge from cell (i,j) to (i+1,j), safeA(i,j) the edge to
%   (i,j+1); the last row / column wraps to the first.
% • An empty cell (has = false) belongs to no branch (label 0) and its edges
%   do not count, whatever the masks say.
% • A branch here is a statement about THIS grid's resolution: an axis that
%   is under-resolved shows few safe edges, and that is the honest answer --
%   interpolating along it is not safe at this spacing.
%
%% Inputs:
%
%  safeD, safeA             [nD x nA] logical       safe edges along each axis
%  has                      [nD x nA] logical       the cells that hold an entry
%
%% Outputs:
%
%  branch                   [nD x nA]               branch number per cell, 0 =
%                                                   no entry; numbered by size,
%                                                   1 = the largest
%  nBranch                  int                     how many
%
%% Revision History:
%  M. Casey                                                   (c) 09/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

[nD, nA] = size(has);
assert(isequal(size(safeD), [nD nA]) && isequal(size(safeA), [nD nA]), 'phase_components:size', ...
       'the two edge masks and the cell mask must have one size');
id = reshape(1:nD*nA, nD, nA);                       % a node per cell
nextD = circshift(id, -1, 1);   nextA = circshift(id, -1, 2);
okD = safeD & has & circshift(has, -1, 1);           % an edge counts only between two entries
okA = safeA & has & circshift(has, -1, 2);
from = [id(okD); id(okA)];   to = [nextD(okD); nextA(okA)];
keep = from ~= to;                                    % a one-cell axis wraps onto itself
G = graph(from(keep), to(keep), [], nD*nA);
comp = conncomp(G);                                  % one label per node
comp(~has(:).') = 0;
% number the branches by size, largest first
labels = unique(comp(comp > 0));
sizes = arrayfun(@(c) nnz(comp == c), labels);
[~, order] = sort(sizes, 'descend');
branch = zeros(nD, nA);
for k = 1:numel(order), branch(comp == labels(order(k))) = k; end
nBranch = numel(labels);
end
