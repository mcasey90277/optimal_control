function ok = test_lib()
%% Purpose:
%
%   The lib/ helpers against hand-computed answers. partial_f_test: a
%   six-point case worked by hand (RSS, F and the p-value from the beta
%   function), an exact-fit case, the :nested refusal when X0 is not in
%   the span of X1 and the :shape refusal.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every check passed
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok   = true;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'lib'));

%% partial_f_test: y = 1 + 2 x + 3 z exactly, X0 drops z -> huge F, p ~ 0
x = [0 1 2 3 4 5].';  z = [1 0 1 0 1 0].';
y = 1 + 2*x + 3*z;
o = partial_f_test(y, [ones(6, 1), x], [ones(6, 1), x, z]);
ok = check(ok, max(abs(o.b1 - [1; 2; 3])) < 1e-10 && o.RSS1 < 1e-18 && o.p < 1e-12, 'exact fit: coefficients recovered, RSS1 ~ 0, p ~ 0');
ok = check(ok, o.df1 == 1 && o.df2 == 3, 'degrees of freedom 1 and 3');

%% hand-worked case: y with noise, compare against the closed forms
y = [1.1; 2.9; 5.2; 6.8; 9.1; 11.0];
X0 = [ones(6, 1), x];  X1 = [ones(6, 1), x, z];
o  = partial_f_test(y, X0, X1);
b0 = X0 \ y;  RSS0 = sum((y - X0*b0).^2);
b1 = X1 \ y;  RSS1 = sum((y - X1*b1).^2);
F  = ((RSS0 - RSS1)/1) / (RSS1/3);
p  = betainc(3/(3 + F), 1.5, 0.5);
ok = check(ok, abs(o.F - F) < 1e-12 && abs(o.p - p) < 1e-12, sprintf('hand-worked F = %.4f, p = %.4f reproduced', F, p));
ok = check(ok, o.p > 0 && o.p < 1 && o.R2_1 >= o.R2_0, 'p in (0, 1) and R2 non-decreasing with the larger model');
ok = check(ok, abs(o.pt1(3) - o.p) < 1e-12 && abs(o.t1(3)^2 - o.F) < 1e-9, 'one added column: its t-test p equals the partial-F p and t^2 = F');
ok = check(ok, isequal(size(o.se1), [3 1]) && all(o.se1 > 0), 'standard errors are a positive [p1 x 1] column');
ok = check(ok, isequal(size(o.ci1), [3 2]) && all(o.ci1(:, 1) < o.b1) && all(o.ci1(:, 2) > o.b1), '95% CI brackets each coefficient');
tq = sqrt(3*(1/betaincinv(0.05, 1.5, 0.5) - 1));
ok = check(ok, abs(tq - 3.1824) < 2e-4 && abs(o.ci1(3, 2) - o.ci1(3, 1) - 2*tq*o.se1(3)) < 1e-12, 'CI half-width is t_{0.975,3} = 3.182 standard errors');
%% exact-fit policies and the rank refusal
o2 = partial_f_test(y, [ones(6, 1), x, z, [1; 0; 0; 1; 0; 1]], [ones(6, 1), x, z, [1; 0; 0; 1; 0; 1], [0; 1; 0; 0; 1; 0]]);
ok = check(ok, o2.F >= 0 && o2.p >= 0 && o2.p <= 1, 'four and five columns on six points: finite F, p in [0, 1]');
yE = 1 + 2*x;
o3 = partial_f_test(yE, [ones(6, 1), x], [ones(6, 1), x, z]);
ok = check(ok, o3.F == 0 && o3.p == 1, 'both models exact (RSS0 = RSS1 = 0): no improvement claimed, p = 1');
try
    partial_f_test(y, [ones(6, 1), x], [ones(6, 1), x, 2*x]);
    ok = check(ok, false, ':rank refusal');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:partial_f_test:rank'), ':rank refusal on a rank-deficient X1');
end

%% the p-value against a tabulated point: F(1, 3) upper tail at F = 10.13 is 0.05
o = struct('F', 10.128, 'df1', 1, 'df2', 3);
ok = check(ok, abs(betainc(o.df2/(o.df2 + o.df1*o.F), o.df2/2, o.df1/2) - 0.05) < 2e-4, 'beta-function tail matches the F(1,3) 5% critical value 10.13');

%% refusals
try
    partial_f_test(y, [ones(6, 1), z], [ones(6, 1), x]);
    ok = check(ok, false, ':nested refusal');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:partial_f_test:nested'), ':nested refusal when X0 is not in the span of X1');
end
try
    partial_f_test(y, X0(1:5, :), X1);
    ok = check(ok, false, ':shape refusal');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:partial_f_test:shape'), ':shape refusal on a row mismatch');
end

%% mirror_state: the hand vector, a matrix, and the refusal
x  = [1; 2; 3; 4; 5; 6];
ok = check(ok, isequal(mirror_state(x), [1; -2; 3; -4; 5; -6]), 'mirror_state flips y, vx, vz');
ok = check(ok, isequal(mirror_state([x, 2*x]), [[1; -2; 3; -4; 5; -6], 2*[1; -2; 3; -4; 5; -6]]), 'mirror_state acts column-wise on [6 x n]');
ok = check(ok, isequal(mirror_state(mirror_state(x)), x), 'mirror_state is an involution');
try
    mirror_state(x(1:5));
    ok = check(ok, false, 'mirror_state :shape refusal');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:mirror_state:shape'), 'mirror_state :shape refusal on a 5-vector');
end

%% symmetry_residual: at an OFF-PLANE point a symmetric orbit maps onto itself; a non-periodic state does not
addpath(fullfile(fileparts(fileparts(here)), 'costate_common'));
[~, ~, rv, info] = evalc('get_family_orbit(''dro'', struct(''tau'', 1))');
x0 = rv(1, :).';  Tp = info.periodND;  mu = 0.012150585609624;
ok = check(ok, max(abs(x0([2 4 6]))) < 1e-12, 'the stored initial state sits on the symmetry plane (so a check there would be vacuous)');
[~, X] = pumpkyn.cr3bp.prop([0 0.37*Tp], x0.', mu);
x1 = X(end, 1:6).';
ok = check(ok, max(abs(x1([2 4 6]))) > 1e-3, 'phi(x0, 0.37 T) is off the plane');
[res, tS] = symmetry_residual(x1, Tp, mu);
ok = check(ok, res < 1e-8, sprintf('DRO tau 1: mirror image of the off-plane point lies on the orbit (residual %.1e at t = %.4f of T = %.4f)', res, tS, Tp));
ok = check(ok, tS >= 0 && tS <= Tp, 'minimizing time inside [0, T]');
xn = x1;  xn(5) = xn(5)*1.02;                              % 2% faster: not on the periodic orbit
resN = symmetry_residual(xn, Tp, mu);
ok = check(ok, resN > 1e-6, sprintf('a non-periodic state is NOT mapped onto an orbit (residual %.1e): the gate can fail', resN));


%% edge_tables: synthetic 2-sheet table; the fastest SOLVED cell is refuted and must not be selected
tStar = 382981.289129055;  d = tStar/86400;
S = struct('dep_key', 'a', 'arr_key', 'b', 'dep_family', 'dro', 'arr_family', 'tulip', 'rungs_N', [2 1 0.5], ...
           'has_solution', false(2, 2, 3), 'accepted', false(2, 2, 3), 'tf_nd', nan(2, 2, 3), ...
           'sD_frac', [0 0.5], 'sA_frac', [0 0.5], ...
           'constants', struct('muStar', 0.012150585609624, 'lStar_km', 389703.264829278, 'tStar_s', tStar));
S.has_solution(:, :, 1) = [true true; true false];  S.accepted(:, :, 1) = [false true; true false];  S.tf_nd(:, :, 1) = [1 2; 3 NaN];
S.has_solution(:, :, 2) = [true false; false false]; S.accepted(:, :, 2) = [true false; false false]; S.tf_nd(:, :, 2) = [4 NaN; NaN NaN];
S2 = S;  S2.arr_key = 'c';  S2.has_solution(:) = false;  S2.accepted(:) = false;  S2.has_solution(1, 1, 3) = true;  S2.accepted(1, 1, 3) = true;  S2.tf_nd(1, 1, 3) = 5;
E = edge_tables([S; S2]);
ok = check(ok, numel(E) == 2 && E(2).sheet_index == 2 && strcmp(E(2).arr_key, 'c'), 'edge_tables: one row per sheet, sheet_index and keys carried');
ok = check(ok, isequal(E(1).n_solved, [3 1 0]) && isequal(E(1).n_accepted, [2 1 0]) && isequal(E(1).closed_frac, [2/4 1/4 0]), 'edge_tables: solved and accepted counts, closed_frac over accepted');
ok = check(ok, abs(E(1).tf_min_days(1) - 2*d) < 1e-12, 'edge_tables: the refuted 1-day cell is NOT the minimum; min over ACCEPTED is 2 days');
ok = check(ok, abs(E(1).tf_med_days(1) - 2.5*d) < 1e-12, 'edge_tables: median over accepted cells (2, 3 days -> 2.5)');
ok = check(ok, isnan(E(1).tf_min_days(3)) && isnan(E(1).tf_med_days(3)), 'edge_tables: NaN on a rung with nothing accepted');
ok = check(ok, E(1).F_low_N == 1 && strcmp(E(1).F_low_str, '1'), 'edge_tables: lowest accepted rung 1 N');
ok = check(ok, E(2).F_low_N == 0.5 && strcmp(E(2).F_low_str, '0.5 (lowest tested)'), 'edge_tables: a sheet accepted at the last rung says so, not a bracket');

%% network_graph: scrambled input order; attributes looked up by (dep, arr) are the right ones
Nt = struct('key', {'b'; 'c'; 'a'}, 'family', {'tulip'; 'tulip'; 'dro'}, 'period_nd', {1; 2; 3}, 'CJ', {3; 3; 3}, 'sigma', {0; 0; 0}, 'stable', {true; true; true});
G  = network_graph([E(2); E(1)], Nt);                       % scrambled: sheet 2 first
ok = check(ok, numnodes(G) == 3 && numedges(G) == 2 && isequal(sort(G.Nodes.Name), {'a'; 'b'; 'c'}), 'network_graph: 3 named nodes, 2 edges');
ia = findedge(G, 'a', 'b');  ic = findedge(G, 'a', 'c');
ok = check(ok, G.Edges.sheet_index(ia) == 1 && G.Edges.sheet_index(ic) == 2, 'network_graph: sheet_index follows the edge, not the input order');
ok = check(ok, abs(G.Edges.tf_min_days(ia, 1) - 2*d) < 1e-12 && isnan(G.Edges.tf_min_days(ic, 1)) && abs(G.Edges.tf_min_days(ic, 3) - 5*d) < 1e-12, 'network_graph: per-rung t_f rows belong to their edges');
ok = check(ok, G.Edges.F_low_N(ia) == 1 && G.Edges.F_low_N(ic) == 0.5 && ~any(G.Edges.hypothetical), 'network_graph: scalar attributes follow too; hypothetical false');
ok = check(ok, strcmp(G.Nodes.family{findnode(G, 'a')}, 'dro') && G.Nodes.period_nd(findnode(G, 'c')) == 2, 'network_graph: node variables by name');
try
    network_graph(E(1), Nt(1:2));
    ok = check(ok, false, 'network_graph :node refusal');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:network_graph:node'), 'network_graph :node refusal when an edge key is missing');
end


%% accessibility_plane: raw coverage vs the known-feasible envelope; monotone in the cap; zeros below the data
P = accessibility_plane([S; S2], d*[0.5 2 2.5 10]);
ok = check(ok, isequal(size(P.pairs), [3 4]) && isequal(size(P.pairs_env), [3 4]) && isequal(size(P.closed), [3 4]), 'accessibility_plane: [nR x nC] arrays');
ok = check(ok, isequal(P.pairs(1, :), [0 1 1 1]) && isequal(P.cells(1, :), [0 1 1 2]), 'accessibility_plane: raw at rung 1 counts ACCEPTED cells (2 and 3 days; the refuted 1-day cell is not there)');
ok = check(ok, isequal(P.pairs(3, :), [0 0 0 1]) && isequal(P.cells(3, :), [0 0 0 1]), 'accessibility_plane: raw at rung 3 only sheet 2 (5 days)');
ok = check(ok, isequal(P.closed{1, 4}, 1) && isequal(P.closed{3, 4}, 2) && isempty(P.closed{1, 1}), 'accessibility_plane: closed sheet lists');
ok = check(ok, isequal(P.pairs_env(1, :), [0 1 1 2]), 'accessibility_plane: the envelope at the 2 N bound includes sheet 2''s 0.5 N witness (feasible under a higher bound)');
ok = check(ok, isequal(P.cells_env(1, :), [0 1 1 4]) && isequal(P.cells_env(3, :), P.cells(3, :)), 'accessibility_plane: envelope cells take the min over deeper rungs; the deepest row equals the raw row');
ok = check(ok, all(all(diff(P.pairs, 1, 2) >= 0)) && all(all(diff(P.cells_env, 1, 2) >= 0)) && all(all(diff(P.pairs_env, 1, 1) <= 0)), 'accessibility_plane: non-decreasing in the cap; envelope non-increasing toward deeper rungs');
P0 = accessibility_plane([S; S2], d*[0.1 0.2]);
ok = check(ok, ~any(P0.pairs(:)) && ~any(P0.cells_env(:)), 'accessibility_plane: caps below the data give all zeros, no error');

%% relay_closure: a toy a -> b -> c where the smallest SUM of leg minima is not the smallest ASSEMBLED time
mk = @(dk, ak, A, tf) struct('dep_key', dk, 'arr_key', ak, 'dep_family', 'x', 'arr_family', 'x', 'rungs_N', [2 1], ...
                           'has_solution', A, 'accepted', A, 'tf_nd', tf, 'sD_frac', [0 0.5], 'sA_frac', [0 0.5], ...
                           'constants', struct('muStar', 0.012150585609624, 'lStar_km', 389703.264829278, 'tStar_s', tStar));
% leg 1 (a -> b): at rung 1 two accepted cells: arrive at phase 0 in 3 days, or at phase 0.5 in 3.2 days; rung 2 open
A1 = false(2, 2, 2);  A1(1, 1, 1) = true;  A1(1, 2, 1) = true;
t1 = nan(2, 2, 2);    t1(1, 1, 1) = 3/d;   t1(1, 2, 1) = 3.2/d;
% leg 2 (b -> c): at rung 1 one accepted cell departing at phase 0.5 in 2 days; at rung 2 one cell departing at phase 0 in 2.5 days
A2 = false(2, 2, 2);  A2(2, 1, 1) = true;  A2(1, 1, 2) = true;
t2 = nan(2, 2, 2);    t2(2, 1, 1) = 2/d;   t2(1, 1, 2) = 2.5/d;
Tt = [mk('L1_a', 'L2_b', A1, t1); mk('L2_b', 'tulip_c', A2, t2); mk('L1_a', 'tulip_c', A1, t1)];   % the third is a direct edge, not a route
Nt = struct('key', {'L1_a'; 'L2_b'; 'tulip_c'}, 'period_nd', {1; 4/d; 3});                        % relay period = 4 days
R  = relay_closure(Tt, Nt, 'L1', 'L2', 'tulip');
ok = check(ok, numel(R) == 1 && strcmp(R.relay_key, 'L2_b') && R.s1 == 1 && R.s2 == 2, 'relay_closure: one route a -> b -> c, sheets indexed into T');
ok = check(ok, isequal(R.closed, [true false]), 'relay_closure: closed iff BOTH legs have an accepted cell (rung 2: leg 1 open)');
ok = check(ok, abs(R.t_sum_days(1) - 5) < 1e-9, 'relay_closure: optimistic sum = 3 + 2 days');
% assembled: arrive phase 0 (3 d) then wait to phase 0.5 = 2 d, +2 d = 7 d; or arrive phase 0.5 (3.2 d), no wait, +2 d = 5.2 d -> choose 5.2
ok = check(ok, abs(R.t_route_days(1) - 5.2) < 1e-9 && abs(R.wait_days(1)) < 1e-9 && isequal(R.cells(:, 1), [1; 2; 2; 1]), 'relay_closure: the assembled minimizer is the phase-compatible 5.2-day pairing, not the 5-day sum');
ok = check(ok, R.t_route_days(1) >= R.t_sum_days(1) && isnan(R.t_route_days(2)) && isnan(R.t_sum_days(2)) && all(R.cells(:, 2) == 0), 'relay_closure: t_route >= t_sum; NaN and zero cells when open');
Tbad = Tt;  Tbad(2).sD_frac = [0.1 0.6];
try
    relay_closure(Tbad, Nt, 'L1', 'L2', 'tulip');
    ok = check(ok, false, 'relay_closure :grid refusal');
catch err
    ok = check(ok, strcmp(err.identifier, 'orbit_network:relay_closure:grid'), 'relay_closure :grid refusal when the two legs'' phase grids on the relay differ');
end
fprintf('test_lib: %s\n', verdict(ok));
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
function s = verdict(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
