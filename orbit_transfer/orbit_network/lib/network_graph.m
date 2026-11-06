function G = network_graph(E, N)
%% Purpose:
%
%   The directed graph of the catalogs: nodes are orbits (named by key),
%   edges are sheets. MATLAB sorts a digraph's edges by (source, target),
%   so every attribute lives IN the edge table and is read back from
%   G.Edges; sheet_index maps an edge to its row of the sheet table.
%
%% Inputs:
%
%  E                        struct [m x 1]          edge_tables output
%
%  N                        struct [n x 1]          Node table (key,
%                                                   family, period_nd, CJ,
%                                                   sigma, stable)
%
%% Outputs:
%
%  G                        digraph                 Nodes: Name, family,
%                                                   period_nd, CJ, sigma,
%                                                   stable. Edges:
%                                                   EndNodes (names),
%                                                   sheet_index,
%                                                   closed_frac,
%                                                   tf_min_days,
%                                                   tf_med_days,
%                                                   n_accepted, n_solved
%                                                   [m x nR], F_low_N,
%                                                   hypothetical (false)
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the catalogs' graph, degrees per node
     here = fileparts(fileparts(mfilename('fullpath')));
     addpath(here, fullfile(fileparts(here), 'costate_common'));
     T = load_catalog_sheets([]);
     G = network_graph(edge_tables(T), orbit_node_table(T, struct('verbose', false)));
     disp([G.Nodes(:, 'Name'), table(indegree(G), outdegree(G), 'VariableNames', {'in', 'out'})]);
     return
end

keys = {N.key};
dep  = {E.dep_key};  arr = {E.arr_key};
miss = [dep(~ismember(dep, keys)), arr(~ismember(arr, keys))];
if ~isempty(miss)
    error('orbit_network:network_graph:node', 'edge key(s) not in the node table: %s', strjoin(unique(miss), ', '));
end
nodeTab = table(keys(:), {N.family}.', [N.period_nd].', [N.CJ].', [N.sigma].', logical([N.stable]).', ...
                'VariableNames', {'Name', 'family', 'period_nd', 'CJ', 'sigma', 'stable'});
edgeTab = table([dep(:), arr(:)], [E.sheet_index].', vertcat(E.closed_frac), vertcat(E.tf_min_days), vertcat(E.tf_med_days), ...
                vertcat(E.n_accepted), vertcat(E.n_solved), [E.F_low_N].', false(numel(E), 1), ...
                'VariableNames', {'EndNodes', 'sheet_index', 'closed_frac', 'tf_min_days', 'tf_med_days', 'n_accepted', 'n_solved', 'F_low_N', 'hypothetical'});
G = digraph(edgeTab, nodeTab);
end
