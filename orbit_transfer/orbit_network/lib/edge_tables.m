function E = edge_tables(T)
%% Purpose:
%
%   Per-edge, per-rung attributes of the network built from a sheet
%   table, over ACCEPTED cells (certified witnesses): solved and accepted
%   counts, the accepted fraction of the phasing grid, the min and median
%   time of flight over accepted cells (NaN when none), and the lowest
%   accepted rung. An unaccepted rung is unresolved, not infeasible, so
%   the lowest accepted rung is reported as what it is and never as a
%   floor or a bracket.
%
%% Inputs:
%
%  T                        struct [n x 1]          Sheet table
%                                                   (load_catalog_sheets)
%
%% Outputs:
%
%  E                        struct [n x 1]          sheet_index, dep_key,
%                                                   arr_key, dep_family,
%                                                   arr_family, rungs_N
%                                                   [1 x nR], n_solved,
%                                                   n_accepted,
%                                                   closed_frac,
%                                                   tf_min_days,
%                                                   tf_med_days [1 x nR],
%                                                   F_low_N, F_low_str
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the 56 edges at 1 N and 0.5 N
     here = fileparts(fileparts(mfilename('fullpath')));
     addpath(here);
     E = edge_tables(load_catalog_sheets([]));
     r1 = find(E(1).rungs_N == 1);  r05 = find(E(1).rungs_N == 0.5);
     fprintf('%-18s %-14s %7s %8s %7s %8s %20s\n', 'dep', 'arr', 'frac1N', 'tf1N[d]', 'frac.5', 'tf.5[d]', 'lowest accepted');
     for kk = 1:numel(E)
         fprintf('%-18s %-14s %7.3f %8.3f %7.3f %8.3f %20s\n', E(kk).dep_key, E(kk).arr_key, E(kk).closed_frac(r1), E(kk).tf_min_days(r1), E(kk).closed_frac(r05), E(kk).tf_min_days(r05), E(kk).F_low_str);
     end
     return
end

n    = numel(T);
rows = cell(n, 1);
for kk = 1:n
    s     = T(kk);
    nR    = numel(s.rungs_N);
    nCell = size(s.accepted, 1)*size(s.accepted, 2);
    toDay = s.constants.tStar_s/86400;
    e.sheet_index = kk;
    e.dep_key     = s.dep_key;
    e.arr_key     = s.arr_key;
    e.dep_family  = s.dep_family;
    e.arr_family  = s.arr_family;
    e.rungs_N     = s.rungs_N;
    e.n_solved    = zeros(1, nR);
    e.n_accepted  = zeros(1, nR);
    e.closed_frac = zeros(1, nR);
    e.tf_min_days = nan(1, nR);
    e.tf_med_days = nan(1, nR);
    for rr = 1:nR
        a  = s.accepted(:, :, rr);
        tf = s.tf_nd(:, :, rr);
        e.n_solved(rr)    = nnz(s.has_solution(:, :, rr));
        e.n_accepted(rr)  = nnz(a);
        e.closed_frac(rr) = nnz(a)/nCell;
        if any(a(:))
            e.tf_min_days(rr) = min(tf(a))*toDay;
            e.tf_med_days(rr) = median(tf(a))*toDay;
        end
    end
    low = find(e.n_accepted > 0, 1, 'last');              % rungs run shallow -> deep
    if isempty(low)
        e.F_low_N = NaN;  e.F_low_str = 'none';
    elseif low == nR
        e.F_low_N = s.rungs_N(nR);  e.F_low_str = sprintf('%g (lowest tested)', s.rungs_N(nR));
    else
        e.F_low_N = s.rungs_N(low);  e.F_low_str = sprintf('%g', s.rungs_N(low));
    end
    rows{kk} = e;
end
E = vertcat(rows{:});
end
