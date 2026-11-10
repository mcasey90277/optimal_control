function seedTop = seed_top_from_entry(sheetMat, iD, iA, kr, thrustTargetN)
% SEED_TOP_FROM_ENTRY  Fly a sheet entry and resample it as an engine seed.
%
%   Flies entry (iD, iA, kr) of a thrust_ladder_library sheet with
%   pumpkyn.cr3bp.tfMinProp, resamples the flight on the engine's N + 1
%   uniform nodes and returns the struct thrust_ladder_library's
%   opts.seedTop expects. The seed is thrust-locked (thr = 1) and its t_f is
%   scaled to the target rung by the engine's own rule (Tprev/TN)^0.6, so a
%   12 N entry seeds a 15 N re-solve (reverse warm start, basin repair) and
%   a neighbour cell's 15 N entry seeds a 15 N solve (phase continuation).
%
% INPUTS:
%   sheetMat      - sheet .mat written by thrust_ladder_library     [char]
%   iD, iA, kr    - departure phase, arrival phase, rung indices    [scalars]
%   thrustTargetN - rung the seed will be used at, N                [scalar]
%
% OUTPUTS:
%   seedTop       - .X [7 x N+1] [r; v; m], .U [4 x N+1] [alpha; thr],
%                   .tf scalar ND, .from struct(iD, iA, kr, thrustN, tf_nd)
%
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(here, 'tests'));
    sheetMat = fixture_demo_sheet();  iD = 1;  iA = 1;  kr = 2;  thrustTargetN = 15;
end
Q = load(sheetMat, 'Z8', 'OK', 'rungs', 'sD', 'meta');
if ~Q.OK(iD, iA, kr)
    error('seed_top_from_entry:notOK', 'entry (%d,%d,rung %g N) of %s is not OK', ...
          iD, iA, Q.rungs(kr), sheetMat);
end
meta   = Q.meta;
z8     = Q.Z8(:, iD, iA, kr);
TN     = Q.rungs(kr);
ndp    = nd_propulsion(TN, meta.ispS, meta.m0kg, meta.lStar, meta.tStar);
[tD, rvD] = get_family_orbit(meta.depFamily, meta.depParams);
dep    = phase_state(tD, rvD);
rv0    = dep(Q.sD(iD));
[tt, yy] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(:); 1; z8(1:7)], ndp.Tnd, ndp.cnd, meta.muStar);
N      = meta.N;
tg     = linspace(0, z8(8), N+1);
[tu, iu] = unique(tt(:), 'stable');          % repeated samples break spline; keep the first of each
if numel(tu) < numel(tt), tt = tu;  yy = yy(iu, :); end   % same index selects the rows of yy
Y      = interp1(tt, yy, tg, 'spline');                       % [N+1 x 14]
if any(~isfinite(Y(:)))
    error('seed_top_from_entry:nonFinite', 'resampled flight of entry (%d,%d,%g N) is not finite', iD, iA, TN);
end
lamV   = Y(:, 11:13).';
alpha  = -lamV ./ max(vecnorm(lamV, 2, 1), eps);
seedTop = struct('X', Y(:, 1:7).', 'U', [alpha; ones(1, N+1)], ...
                 'tf', z8(8)*(TN/thrustTargetN)^0.6, ...
                 'from', struct('iD', iD, 'iA', iA, 'kr', kr, 'thrustN', TN, 'tf_nd', z8(8)));
if nargin == 0
    fprintf('seed_top_from_entry: %g N entry tf %.5f -> seed tf %.5f for %g N, m_f/m0 %.4f\n', ...
            TN, z8(8), seedTop.tf, thrustTargetN, seedTop.X(7, end));
end
end
