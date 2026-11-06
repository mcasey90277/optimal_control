function ok = test_orbit_node_table()
%% Purpose:
%
%   The node table on a two-sheet synthetic table (one DRO-to-tulip sheet,
%   one L2-halo-to-tulip sheet sharing the tulip): 3 nodes in order of
%   first appearance, periods equal to the catalogs' T_dep / T_arr, the
%   symplectic pairing |lambda|max |lambda|min = 1 and det M = 1, the DRO
%   and tulip stable and the tau 2.8 halo unstable, the degrees, and the
%   cache round-trip (hit on equal keys, miss on a different key set,
%   rebuild forced).
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
root = fileparts(here);
addpath(root, fullfile(fileparts(root), 'costate_common'));

T  = load_catalog_sheets([]);
T2 = [T(find(strcmp({T.dep_key}, 'dro_tau1') & strcmp({T.arr_key}, 'tulip_Np7'), 1)); ...
      T(find(strcmp({T.dep_key}, 'halo_L2_tau2.8') & strcmp({T.arr_key}, 'tulip_Np7'), 1))];
ok = check(ok, numel(T2) == 2, 'two synthetic rows found');

cache = fullfile(tempdir, 'test_orbit_node_table', 'nodes.mat');
if exist(cache, 'file') == 2, delete(cache); end
t0 = tic;
N  = orbit_node_table(T2, struct('cache', cache, 'verbose', false));
tBuild = toc(t0);
ok = check(ok, numel(N) == 3 && isequal({N.key}, {'dro_tau1', 'halo_L2_tau2.8', 'tulip_Np7'}), 'three nodes, departures first then the shared arrival');
ok = check(ok, abs(N(1).period_nd - T2(1).T_dep) < 1e-9 && abs(N(2).period_nd - T2(2).T_dep) < 1e-9 && abs(N(3).period_nd - T2(1).T_arr) < 1e-9, 'periods match the catalogs to 1e-9');
ok = check(ok, all(arrayfun(@(r) abs(r.lam_max*r.lam_min - 1) < 1e-6 && abs(r.detM - 1) < 1e-6, N)), 'symplectic pairing and det M = 1 to 1e-6');
ok = check(ok, N(1).stable && N(3).stable && ~N(2).stable, 'DRO and tulip stable, the tau 2.8 L2 halo unstable');
ok = check(ok, N(1).sigma == 0 && N(3).sigma == 0 && N(2).sigma > 1.0 && N(2).sigma < 1.4, 'sigma exactly 0 for the stable pair, about 1.18 for the halo');
ok = check(ok, N(2).CJ > 3.0 && N(2).CJ < 3.05 && N(1).CJ > 3.05 && N(1).CJ < 3.06, 'Jacobi constants in the known ranges');
ok = check(ok, isequal([N.out_degree], [1 1 0]) && isequal([N.in_degree], [0 0 2]), 'degrees over the two edges');
ok = check(ok, isequal(size(N(1).x0), [6 1]), 'x0 is a column');

%% cache
t0 = tic;
Nc = orbit_node_table(T2, struct('cache', cache, 'verbose', false));
tHit = toc(t0);
ok = check(ok, isequal(Nc, N), 'cache hit returns the same table');
ok = check(ok, tHit < 0.25*max(tBuild, 0.1), sprintf('cache hit is fast (%.2f s vs %.2f s build)', tHit, tBuild));
N1 = orbit_node_table(T2(1), struct('cache', cache, 'verbose', false));
ok = check(ok, numel(N1) == 2, 'different key set misses the cache and rebuilds (2 nodes)');
S  = load(cache);
ok = check(ok, numel(S.N) == 2, 'the rebuilt table replaced the cache');
Nr = orbit_node_table(T2, struct('cache', cache, 'verbose', false, 'rebuild', true));
ok = check(ok, numel(Nr) == 3 && max(abs([Nr.lam_max] - [N.lam_max])) < 1e-9, 'forced rebuild reproduces the spectrum');

fprintf('test_orbit_node_table: %s\n', verdict(ok));
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
function s = verdict(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
