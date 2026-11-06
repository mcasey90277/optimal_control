function R = relay_closure(T, N, src_prefix, relay_prefix, dst_prefix)
%% Purpose:
%
%   Two-leg routes src -> relay -> dst assembled from catalog sheets
%   without new solves. A route is closed at a rung iff both legs have an
%   accepted cell. Its time is PHASE-COMPATIBLE: leg 1 arrives on the
%   relay at phase sA1(b), leg 2 departs at phase sD2(c), and the coast
%   between them on the relay orbit takes T_relay * mod(sD2(c) - sA1(b),
%   1); the route time is the minimum over accepted cell pairs of
%   t1(a, b) + wait + t2(c, d). The sum of the two legs' independent
%   minima is kept as the optimistic lower bound. Phases of both legs
%   must be measured from the same relay-orbit origin (the same
%   get_family_orbit state), which the grid check enforces. Mass: a leg
%   solved at m0 is flyable at the previous leg's arrival mass by scaling
%   thrust and mass by q = m1/m0 <= 1 (same path, bound respected); the
%   note states this. Relay-versus-direct cost is not a question for
%   this function.
%
%% Inputs:
%
%  T                        struct [n x 1]          Sheet table
%                                                   (load_catalog_sheets)
%
%  N                        struct [m x 1]          Node table (key,
%                                                   period_nd)
%
%  src_prefix               char                    Key prefix of the
%                                                   sources ('halo_L1')
%
%  relay_prefix             char                    Of the relays
%                                                   ('halo_L2')
%
%  dst_prefix               char                    Of the destinations
%                                                   ('tulip')
%
%% Outputs:
%
%  R                        struct [k x 1]          src_key, relay_key,
%                                                   dst_key, s1, s2,
%                                                   closed [1 x nR],
%                                                   t_sum_days,
%                                                   t_route_days,
%                                                   wait_days [1 x nR],
%                                                   cells [4 x nR]
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the 32 L1 halo -> L2 halo -> tulip routes of the catalogs, at 1 N
     here = fileparts(fileparts(mfilename('fullpath')));
     addpath(here, fullfile(fileparts(here), 'costate_common'));
     T = load_catalog_sheets([]);
     R = relay_closure(T, orbit_node_table(T, struct('verbose', false)), 'halo_L1', 'halo_L2', 'tulip');
     r1 = find(T(1).rungs_N == 1);
     for kk = 1:numel(R)
         fprintf('%-18s -> %-16s -> %-10s  closed %d  t_sum %.3f d  t_route %.3f d (wait %.3f d)\n', R(kk).src_key, R(kk).relay_key, R(kk).dst_key, R(kk).closed(r1), R(kk).t_sum_days(r1), R(kk).t_route_days(r1), R(kk).wait_days(r1));
     end
     return
end

dep  = {T.dep_key};  arr = {T.arr_key};
keys = {N.key};
isS  = @(k, p) strncmp(k, p, numel(p));
pairs = zeros(0, 2);
for s1 = 1:numel(T)
    if ~(isS(dep{s1}, src_prefix) && isS(arr{s1}, relay_prefix)), continue; end
    for s2 = 1:numel(T)
        if strcmp(dep{s2}, arr{s1}) && isS(arr{s2}, dst_prefix)
            pairs(end+1, :) = [s1 s2];
        end
    end
end
nRoute = size(pairs, 1);
rows   = cell(nRoute, 1);
for qq = 1:nRoute
    s1 = pairs(qq, 1);  s2 = pairs(qq, 2);
    L1 = T(s1);  L2 = T(s2);
    nR = numel(L1.rungs_N);
    kN = find(strcmp(keys, L1.arr_key), 1);
    if isempty(kN)
        error('orbit_network:relay_closure:node', 'relay %s is not in the node table', L1.arr_key);
    end
    if numel(L1.sA_frac) ~= numel(L2.sD_frac) || max(abs(sort(L1.sA_frac(:)) - sort(L2.sD_frac(:)))) > 1e-12
        error('orbit_network:relay_closure:grid', 'route %s -> %s -> %s: leg 1 arrival phases %s and leg 2 departure phases %s differ', ...
              L1.dep_key, L1.arr_key, L2.arr_key, mat2str(L1.sA_frac, 4), mat2str(L2.sD_frac, 4));
    end
    toDay = L1.constants.tStar_s/86400;
    Trel  = N(kN).period_nd*toDay;
    W     = Trel*mod(L2.sD_frac(:).' - L1.sA_frac(:), 1);  % [nA1 x nD2] wait from arrival phase b to departure phase c
    r.src_key   = L1.dep_key;
    r.relay_key = L1.arr_key;
    r.dst_key   = L2.arr_key;
    r.s1        = s1;
    r.s2        = s2;
    r.closed       = false(1, nR);
    r.t_sum_days   = nan(1, nR);
    r.t_route_days = nan(1, nR);
    r.wait_days    = nan(1, nR);
    r.cells        = zeros(4, nR);
    for rr = 1:nR
        a1 = L1.accepted(:, :, rr);  t1 = L1.tf_nd(:, :, rr)*toDay;
        a2 = L2.accepted(:, :, rr);  t2 = L2.tf_nd(:, :, rr)*toDay;
        if ~any(a1(:)) || ~any(a2(:)), continue; end
        r.closed(rr)     = true;
        r.t_sum_days(rr) = min(t1(a1)) + min(t2(a2));
        best = Inf;
        [aa, bb] = find(a1);  [cc, dd] = find(a2);
        for p1 = 1:numel(aa)
            for p2 = 1:numel(cc)
                tot = t1(aa(p1), bb(p1)) + W(bb(p1), cc(p2)) + t2(cc(p2), dd(p2));
                if tot < best
                    best = tot;
                    r.t_route_days(rr) = tot;
                    r.wait_days(rr)    = W(bb(p1), cc(p2));
                    r.cells(:, rr)     = [aa(p1); bb(p1); cc(p2); dd(p2)];
                end
            end
        end
    end
    rows{qq} = r;
end
if nRoute == 0
    R = struct('src_key', {}, 'relay_key', {}, 'dst_key', {}, 's1', {}, 's2', {}, 'closed', {}, 't_sum_days', {}, 't_route_days', {}, 'wait_days', {}, 'cells', {});
else
    R = vertcat(rows{:});
end
end
