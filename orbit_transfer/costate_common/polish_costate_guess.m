function [z8, info] = polish_costate_guess(zGuess, rv0, rvf, phys, opts)
%% Purpose:
%
%   Turn a costate GUESS into a root of the minimum-time boundary-value
%   problem, or say that it could not be done. Three steps, the same three
%   every consumer of an interpolated guess needs:
%
%     1. fly the guess from t = 0 (its arrival miss is the first honest
%        measure of how good the guess is);
%     2. cut that flight into multiple-shooting junctions (seed_from_z8);
%     3. polish with ms_tfmin.
%
%   It is ONE function so that the whole of it can run under run_capped's
%   hard wall-clock fence: a poor guess can make step 1 crawl near a primary
%   for hours inside a single integration, where no in-process budget fires.
%
%  ASSUMPTIONS / NOTES:
%
% • An UNCONVERGED polish returns NaN, never ms_tfmin's best iterate: a best
%   iterate fed onward looks like a root and is not one.
% • maxIter defaults to 600, not ms_tfmin's 100. A guess that flies within a
%   few hundred km still needs 90-500 damped Newton iterations, because
%   re-flying it for 26 days amplifies its error into every junction; at 100
%   the measured hit rate of an interpolated guess read 7 of 14 instead of
%   12 of 14 (FINDINGS 88).
% • Convergence says the result is A root. Whether it is the root the guess
%   was aimed at is the caller's question (compare with the guess).
%
%% Inputs:
%
%  zGuess                   [8 x 1]                 [lambda(0) (7); t_f], normal
%                                                   chart, ND
%
%  rv0, rvf                 [6 x 1]                 departure and arrival
%                                                   states, rotating frame ND
%
%  phys                     struct                  .Tnd .cnd .muStar, and
%                                                   optionally .lStar [km,
%                                                   389703.264829278] for the
%                                                   reported miss
%
%  opts                     struct (optional)       .K [24] segments,
%                                                   .maxIter [600],
%                                                   .tolR [3e-11],
%                                                   .wallSec [400],
%                                                   .conjTest [false],
%                                                   .polishMax [5] plain Newton
%                                                   steps ms_bvp takes AFTER
%                                                   its capped search
%
%% Outputs:
%
%  z8                       [8 x 1]                 the root, or NaN(8,1)
%
%  info                     struct                  .converged .normR .iters
%                                                   .wall (s), .guessMissKm
%                                                   (flown miss of the GUESS),
%                                                   .conj (if .conjTest),
%                                                   .Y .tGrid (the polished
%                                                   junctions; empty if not
%                                                   converged)
%
%% Revision History:
%  M. Casey                                                   (c) 09/20/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 5, opts = struct(); end
K = 24;          if isfield(opts, 'K'),        K = opts.K;               end
maxIter = 600;   if isfield(opts, 'maxIter'),  maxIter = opts.maxIter;   end
tolR = 3e-11;    if isfield(opts, 'tolR'),     tolR = opts.tolR;         end
wallSec = 400;   if isfield(opts, 'wallSec'),  wallSec = opts.wallSec;   end
conjTest = false; if isfield(opts, 'conjTest'), conjTest = opts.conjTest; end
polishMax = 5;   if isfield(opts, 'polishMax'), polishMax = opts.polishMax; end
lStar = 389703.264829278;  if isfield(phys, 'lStar'), lStar = phys.lStar; end

zGuess = zGuess(:);  rv0 = rv0(:);  rvf = rvf(:);
assert(numel(zGuess) == 8 && isreal(zGuess) && all(isfinite(zGuess)) && zGuess(8) > 0, ...
       'polish_costate_guess:guess', 'the guess must be eight real finite numbers with t_f > 0');

% 1. the guess's own flight
[~, Yg] = pumpkyn.cr3bp.tfMinProp(zGuess(8), [rv0(1:6); 1; zGuess(1:7)], phys.Tnd, phys.cnd, phys.muStar);
info = struct('converged', false, 'normR', NaN, 'iters', NaN, 'wall', NaN, ...
              'guessMissKm', sqrt(sum((Yg(end, 1:3).' - rvf(1:3)).^2))*lStar, 'conj', [], 'Y', [], 'tGrid', []);

% 2. junctions from the flown guess          3. polish
seed = seed_from_z8(zGuess, rv0(1:6), K, phys.Tnd, phys.cnd, phys.muStar);
[z, it] = ms_tfmin(rv0(1:6), rvf(1:6), seed, phys.Tnd, phys.cnd, phys.muStar, ...
                   struct('tolR', tolR, 'wallSec', wallSec, 'maxIter', maxIter, 'polishMax', polishMax, 'conjTest', conjTest));
info.converged = logical(it.converged);  info.normR = it.normR;  info.iters = it.iters;  info.wall = it.wall;
z8 = NaN(8, 1);
if info.converged && isnumeric(z) && numel(z) == 8 && all(isfinite(z))
    z8 = z(:);  info.Y = it.Y;  info.tGrid = it.tGrid;
    if isfield(it, 'conj'), info.conj = it.conj; end
else
    info.converged = false;
end
end
