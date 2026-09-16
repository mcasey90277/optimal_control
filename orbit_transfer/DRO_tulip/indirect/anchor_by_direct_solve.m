function [C, info] = anchor_by_direct_solve(seedAnchorMat, sD0, sA0, orbits, engine, outMat, opts)
%% Purpose:
%
%   Make the FIRST certified root of a new problem -- another DRO period,
%   another tulip, another engine -- so a phase torus can start: a direct
%   Hermite-Simpson + Sundman solve at (sD0, sA0) warm-started from a known
%   transfer's PMP flight (the seed anchor, typically the 70 mN 7-petal
%   anchor), harvested into a multiple-shooting seed and put through the
%   full gate stack for the NEW orbits and engine (direct_cell_solve). A
%   certified root is saved in the anchor layout arclength_arrival reads
%   (best.z, best.it.Y, best.sA, best.sD, best.tfDays, best.origin, Tnd,
%   cnd). No shipped catalog covers 70 mN, an 8-petal tulip or a DRO
%   period off {0.5, 1, 2, 3}; this is how such a campaign begins.
%
%% Inputs:
%
%  seedAnchorMat            char                    a known certified root
%                                                   (anchor .mat) to warm-
%                                                   start from
%  sD0, sA0                 double                  the new anchor's phases
%  orbits                   struct                  .tauDRO .NpTulip .pmTulip
%  engine                   struct                  .thrustN .ispS .m0kg
%  outMat                   char                    anchor file to write
%  opts                     struct (optional)       direct_cell_solve options
%                                                   (.N .K .clearKm
%                                                   .maxCpuSec .wallSec)
%
%% Outputs:
%
%  C                        struct                  certify_root output (.ok
%                                                   says whether outMat was
%                                                   written)
%  info                     struct                  direct-solve report
%                                                   (.tfDirectDays .perisKm
%                                                   .success .wall)
%
%% Revision History:
%  M. Casey                                                   (c) 09/15/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 7, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
seed = load(seedAnchorMat);
if isfield(seed, 'best'), root = seed.best; else, root = seed; end
assert(isfield(root, 'z') && isfield(root, 'it') && isfield(root.it, 'Y'), 'the seed anchor must carry .z and .it.Y');
so = struct('thrustN', engine.thrustN, 'ispS', engine.ispS, 'm0kg', engine.m0kg, 'tauDRO', orbits.tauDRO, ...
            'NpTulip', orbits.NpTulip, 'pmTulip', orbits.pmTulip, 'sD', sD0, 'physicsOnly', true);
[B, ~] = arclength_arrival('setup', so);            % the NEW problem's physics; no root of it exists yet
% the seed's own departure state (its problem, not ours) for the warm start
if isfield(root, 'rv0') && numel(root.rv0) >= 6, rvSeed = root.rv0(1:6);
else, rvSeed = root.it.Y(1:6, 1);                 % the first junction's state
end
pool = capped_pool();
opts.pool = pool;
[C, info] = direct_cell_solve(root.z(:), rvSeed(:), sD0, sA0, B, opts);
fprintf('anchor_by_direct_solve: DRO tau %.3f -> %d-petal tulip (pm %+d), %.0f mN: direct %.3f d, %s\n', ...
        orbits.tauDRO, orbits.NpTulip, orbits.pmTulip, engine.thrustN*1000, info.tfDirectDays, C.reason);
if ~C.ok, return, end
best = struct('z', C.z(:), 'it', struct('Y', C.Y), 'sA', sA0, 'sD', sD0, 'tfDays', C.tfDays, ...
    'origin', sprintf('anchor_by_direct_solve from %s at (%.4f, %.4f): DRO tau %.3f -> %d-petal tulip pm %+d, %.0f mN / %g s / %g kg, certified %s', ...
                      seedAnchorMat, sD0, sA0, orbits.tauDRO, orbits.NpTulip, orbits.pmTulip, engine.thrustN*1000, engine.ispS, engine.m0kg, char(datetime('now'))));
Tnd = B.Tnd;  cnd = B.cnd;
save(outMat, 'best', 'Tnd', 'cnd');
fprintf('  anchor written: %s (%.3f d)\n', outMat, C.tfDays);
end
