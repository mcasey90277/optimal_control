function [C, info] = direct_cell_solve(zSeed, rvSeed, sDcell, sAcell, B, opts)
%% Purpose:
%
%   One cell, one seed: a direct Hermite-Simpson + Sundman minimum-time
%   solve at the cell (sDcell, sAcell) warm-started from another certified
%   root's PMP flight, harvested into a multiple-shooting seed and put
%   through the full gate stack. This is the branch-blind search that found
%   the fourth and fifth families of the 70 mN torus (FINDINGS 61, 68) and
%   the cell-by-cell filler's inner step (fill_holes_direct): the
%   collocation solver does not know which branch its seed is on and lands
%   in the basin of the fastest nearby extremal.
%
%   A solution whose periselene sits on the clearance floor is refused
%   before certification: it is a constrained arc, not an unconstrained
%   extremal (FINDINGS 68).
%
%% Inputs:
%
%  zSeed                    [8 x 1]                 the seed root (lambda(0), t_f)
%  rvSeed                   [6 x 1]                 the seed's departure state
%  sDcell, sAcell           double                  the cell's phases
%  B                        struct                  arclength_arrival('setup')
%                                                   (.stateD .stateA .Tnd .cnd
%                                                   .mu .problem)
%  opts                     struct (optional)
%   .N [800] collocation intervals, .K [24] shooting segments,
%   .clearKm [1900] lunar clearance (radius), .maxCpuSec [300] per direct
%   solve, .wallSec [600] per certification, .pool [] fence pool
%
%% Outputs:
%
%  C                        struct                  certify_root output, or
%                                                   a struct with .ok false
%                                                   and .reason when the
%                                                   direct solve failed
%  info                     struct                  .tfDirectDays .perisKm
%                                                   .success .wall
%
%% Revision History:
%  M. Casey                                                   (c) 09/15/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 6, opts = struct(); end
d = @(f, v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
D = fullfile(fileparts(fileparts(here)), 'DRO_tulip', 'direct');
addpath(D, fullfile(D, 'lib'), fullfile(D, 'certify'), fullfile(fileparts(fileparts(here)), 'costate_common'));
N = d('N', 800);  K = d('K', 24);  clearKm = d('clearKm', 1900);
maxCpu = d('maxCpuSec', 300);  wallSec = d('wallSec', 600);
rMoonKm = 1737.4;  floorKm = clearKm - rMoonKm;
tStar = B.problem.tStar;  lStar = B.problem.lStar;
mu = B.mu;  Tmax = B.Tnd;  c = B.cnd;
t0 = tic;
info = struct('tfDirectDays', NaN, 'perisKm', NaN, 'success', false, 'wall', NaN);
C = struct('ok', false, 'reason', '', 'tfDays', NaN, 'z', nan(8, 1), 'sA', sAcell, 'sD', sDcell);

rv0 = B.stateD(sDcell);  rv0 = rv0(1:6);
rvf = B.stateA(sAcell);
% the seed's PMP flight, re-boundaried to this cell, as the warm start
[tauR, rvR] = pumpkyn.cr3bp.tfMinProp(zSeed(8), [rvSeed(:); 1; zSeed(1:7)], Tmax, c, mu);
sN = linspace(0, tauR(end), N + 1);
X0 = interp1(tauR, rvR(:, 1:7), sN, 'spline').';
LV = interp1(tauR, rvR(:, 11:13), sN, 'spline');
U0 = [(-LV ./ max(vecnorm(LV, 2, 2), eps)).'; ones(1, N + 1)];
o = casadi_mintime_dro(rv0, rvf(1:6), Tmax, c, mu, N, X0, U0, zSeed(8), struct('maxIter', 3000, ...
        'scheme', 'hermite-simpson', 'sundman', true, 'returnModel', true, ...
        'minAltKm', floorKm, 'maxCpuSec', maxCpu));
if isfield(o, 'model'), o = rmfield(o, 'model'); end
info.tfDirectDays = o.tf * tStar/86400;
info.perisKm = min(vecnorm(o.X(1:3, :) - [1 - mu; 0; 0], 2, 1))*lStar - rMoonKm;
info.success = o.success;
if ~o.success
    C.reason = 'direct solve failed';
elseif info.perisKm < floorKm + 5
    C.reason = sprintf('direct solution rides the clearance floor (%.0f km, %.1f d)', info.perisKm, info.tfDirectDays);
else
    seed = harvest_ms_seed(o, K);
    copts = struct('wallSec', wallSec, 'sA', sAcell, 'sD', sDcell);
    if ~isempty(d('pool', [])), copts.pool = opts.pool; end
    C = certify_root(seed, rv0, rvf, B, copts);
end
info.wall = toc(t0);
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
