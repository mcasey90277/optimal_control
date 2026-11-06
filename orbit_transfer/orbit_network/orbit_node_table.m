function N = orbit_node_table(T, opts)
%% Purpose:
%
%   The distinct orbits of a sheet table, each built once with
%   get_family_orbit and characterized by its period, Jacobi constant and
%   monodromy spectrum. sigma = ln|lambda|_max / T is the paper's
%   instability rate (Braik & Ross 2026), set to exactly 0 on a stable
%   orbit (|lambda|_max < 1 + 1e-3: the continuation tolerance leaves the
%   unit-circle pair a hair outside it). It is computed from
%   pumpkyn.cr3bp.monodromy directly because pumpkynPie's orbitProperties
%   errors on a stable orbit (empty unstable-eigenvalue list).
%
%   A cache (default results/orbit_nodes.mat beside this file) is used when
%   its stored meta (keys, muStar, every parameter struct) equals the
%   request; opts.rebuild forces a rebuild, and publication runs rebuild.
%   The callers' gate G1 compares period_nd with the catalogs' own periods
%   as a further check.
%
%% Inputs:
%
%  T                        struct [n x 1]          Sheet table from
%                                                   load_catalog_sheets
%                                                   (dep_key, arr_key,
%                                                   dep_family, dep_params,
%                                                   arr_family, arr_params,
%                                                   constants)
%
%  opts                     struct (optional)       .cache [results/
%                                                   orbit_nodes.mat];
%                                                   .rebuild [false];
%                                                   .verbose [true] print
%                                                   one line per orbit
%                                                   built
%
%% Outputs:
%
%  N                        struct [m x 1]          One row per distinct
%                                                   key, in order of first
%                                                   appearance (departures
%                                                   of T first, then
%                                                   arrivals): key, family,
%                                                   params, period_nd, x0
%                                                   [6 x 1], CJ, lam_max,
%                                                   lam_min, detM, sigma,
%                                                   stable, in_degree,
%                                                   out_degree
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the 18 orbits of the four catalogs, from the cache if present
     here = fileparts(mfilename('fullpath'));
     addpath(here, fullfile(fileparts(here), 'costate_common'));
     N = orbit_node_table(load_catalog_sheets([]));
     fprintf('%-20s %8s %10s %12s %8s %6s %3s %3s\n', 'key', 'T_nd', 'C_J', 'lam_max', 'sigma', 'stable', 'in', 'out');
     for kk = 1:numel(N)
         fprintf('%-20s %8.4f %10.5f %12.4e %8.4f %6d %3d %3d\n', N(kk).key, N(kk).period_nd, N(kk).CJ, N(kk).lam_max, N(kk).sigma, N(kk).stable, N(kk).in_degree, N(kk).out_degree);
     end
     return
end

if nargin < 2 || isempty(opts), opts = struct(); end
here    = fileparts(mfilename('fullpath'));
cache   = getf(opts, 'cache', fullfile(here, 'results', 'orbit_nodes.mat'));
rebuild = getf(opts, 'rebuild', false);
verbose = getf(opts, 'verbose', true);

%% the distinct keys, departures first, then arrivals, order of first appearance
keys  = [{T.dep_key}, {T.arr_key}];
fams  = [{T.dep_family}, {T.arr_family}];
pars  = [{T.dep_params}, {T.arr_params}];
[ukeys, first] = unique(keys, 'stable');
fams  = fams(first);
pars  = pars(first);

%% cache hit? (keys, muStar and every parameter struct must match; periods are re-checked by the callers' G1)
meta = struct('keys', {ukeys}, 'muStar', T(1).constants.muStar, 'params', {pars}, 'date', char(datetime('now', 'Format', 'yyyy-MM-dd')));
if ~rebuild && exist(cache, 'file') == 2
    S = load(cache);
    if isfield(S, 'N') && isfield(S, 'meta') && isequal(S.meta.keys, ukeys) && isequal(S.meta.muStar, meta.muStar) && isequal(S.meta.params, pars)
        N = S.N;
        N = degrees(N, T);
        return
    end
end

%% build
m    = numel(ukeys);
mu   = T(1).constants.muStar;
rows = cell(m, 1);
for kk = 1:m
    [~, ~, rv, info] = evalc('get_family_orbit(fams{kk}, pars{kk})');   % the loader prints its continuation log
    x0 = rv(1, :).';
    Tp = info.periodND;
    M   = pumpkyn.cr3bp.monodromy(x0.', Tp, mu);
    lam = abs(eig(M));
    CJ  = pumpkyn.cr3bp.jacobi(x0.', mu);
    r.key        = ukeys{kk};
    r.family     = fams{kk};
    r.params     = pars{kk};
    r.period_nd  = Tp;
    r.x0         = x0;
    r.CJ         = CJ;
    r.lam_max    = max(lam);
    r.lam_min    = min(lam);
    r.detM       = det(M);
    r.stable     = r.lam_max < 1 + 1e-3;                    % continuation tolerance puts |lam|max a hair above 1
    if r.stable, r.sigma = 0; else, r.sigma = log(r.lam_max)/Tp; end
    r.in_degree  = 0;
    r.out_degree = 0;
    rows{kk}     = r;
    if verbose
        fprintf('  node %2d/%2d  %-20s T = %.4f  C_J = %.5f  |lam|max = %.4e  sigma = %.4f\n', kk, m, r.key, Tp, CJ, r.lam_max, r.sigma);
    end
end
N = vertcat(rows{:});
N = degrees(N, T);

cdir = fileparts(cache);
if ~isempty(cdir) && ~exist(cdir, 'dir'), mkdir(cdir); end
save(cache, 'N', 'meta');
end

% ------------------------------------------------------------------------
function N = degrees(N, T)
%% Purpose:
%
%   In- and out-degree of each node over the sheet table's directed edges.
%
for kk = 1:numel(N)
    N(kk).out_degree = nnz(strcmp({T.dep_key}, N(kk).key));
    N(kk).in_degree  = nnz(strcmp({T.arr_key}, N(kk).key));
end
end

% ------------------------------------------------------------------------
function v = getf(s, name, default)
%% Purpose:
%
%   Field with a default.
%
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
