function r = phase_edge_residuals(T, G, s, dim)
%% Purpose:
%
%   Do neighbouring entries of a phase sheet lie on ONE smooth branch? Along
%   one phase axis, every edge between a cell and the next gets the residual
%   of the trapezoid rule built from the two cells' own sensitivities,
%
%       r = [T(next) - T(here)] - (1/2) [G(here) + G(next)] * ds ,
%
%   where G = dT/ds at the cell (phase_sensitivity). On a smooth branch r is
%   third order in the spacing, |r| <= (ds^3/12) max|T'''|. A residual far
%   above that means the two cells are NOT neighbours on one branch: there is
%   a jump between them (another root, another winding, a fold), whatever
%   their family label says. It costs nothing: T and G are already known.
%
%  ASSUMPTIONS / NOTES:
%
% • The axis is PERIODIC: the last cell's edge joins it to the first, with
%   spacing mod(s(1) - s(end), 1). The list need not be uniform.
% • A NaN in T or G (an empty cell) blanks the two edges that touch it.
% • If the sheet is under-resolved along the axis (structure finer than the
%   spacing) every residual is large and the test says nothing; judge that
%   from the residuals' own distribution before reading single edges.
%
%% Inputs:
%
%  T                        [nD x nA] | vector      the sheet (flight time)
%  G                        same size               dT/ds along `dim`
%  s                        [1 x n]                 the phase list along `dim`,
%                                                   in [0,1), increasing
%  dim                      1 | 2                   the axis the edges run along
%
%% Outputs:
%
%  r                        size of T               r(cell) = residual of the
%                                                   edge from that cell to the
%                                                   NEXT one along `dim`
%
%% Revision History:
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

assert(isequal(size(T), size(G)), 'phase_edge_residuals:size', 'T and G must have one size');
assert(any(dim == [1 2]) && numel(s) == size(T, dim), 'phase_edge_residuals:axis', ...
       'the phase list must have one value per cell along dimension %d', dim);
ds = mod(diff([s(:); s(1) + 1]), 1);                 % circular spacing, cell -> next
if dim == 1, ds = ds(:); else, ds = ds(:).'; end     % broadcast along the chosen axis
Tn = circshift(T, -1, dim);   Gn = circshift(G, -1, dim);
r = (Tn - T) - 0.5*(G + Gn).*ds;
end
