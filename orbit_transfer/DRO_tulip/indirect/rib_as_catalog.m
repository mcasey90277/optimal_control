function c = rib_as_catalog(R, S, ref, opts)
%% Purpose:
%
%   Wrap a DEPARTURE RIB -- the certified roots at one arrival phase and many
%   departure phases (build_ribs / rib_from_crossing) -- as a one-COLUMN
%   phase catalog, so that everything that reads a catalog (interp_study,
%   phase_catalog_interp, score_interpolator) can read a rib. The departure
%   twin of arrival_sheet_as_catalog: a rib is cheap to walk at a resolution
%   a full library is not (96 departure phases in hours), which makes it the
%   place to measure interpolation along departure phase.
%
%  ASSUMPTIONS / NOTES:
%
% • The column holds the SPINE CROSSING the rib was hung from (the sheet's
%   winner at its departure phase, column R.j) plus every certified rib
%   point. The rib must be that column's: one arrival phase throughout,
%   equal to the sheet's s_A(R.j).
% • The departure LATTICE: a walker that bisects a hard step accepts points
%   between grid lines. With opts.nD the column keeps only the points on
%   the lattice spine + k/nD (the grid steps then mean what the scorer
%   prints) and counts the rest in .rib.nOffLattice; without it every
%   certified point is kept.
% • One family (code 1, label 'rib'): a rib is one continuation from one
%   crossing. Whether neighbours may be blended is still decided by the
%   interpolator's own rules, not by this label.
% • The problem identity (orbits, engine, constants) is COPIED from a
%   reference catalog, as in arrival_sheet_as_catalog: the caller vouches
%   that the rib was walked for that problem. A VIEW for study, not a
%   deliverable: no verdict fields are carried.
%
%% Inputs:
%
%  R                        struct                  one rib: .j (sheet column)
%                                                   .sA .pts (certify_root
%                                                   outputs: .sD .sA .z [8x1])
%                                                   .stop
%
%  S                        struct                  the arrival sheet it hangs
%                                                   from: .sA [1 x nA] .Z8
%                                                   [8 x nA] .problem.sD
%                                                   (spine phase, [0])
%
%  ref                      struct                  a phase catalog of the same
%                                                   problem (.constants
%                                                   .thruster .rungs_N
%                                                   .sheets(1) orbit keys)
%
%  opts                     struct (optional)       .nD [] departure lattice
%                                                   size .tol [1e-8] phase
%                                                   tolerance .outMat '' save
%                                                   here (variable `catalog`)
%
%% Outputs:
%
%  c                        struct                  .sheets(1) with .sD_frac
%                                                   [1 x n] .sA_frac (scalar)
%                                                   .has_solution [n x 1]
%                                                   .tf_nd .entry_index
%                                                   .family_index [n x 1] .z8
%                                                   [8 x n] + the orbit keys;
%                                                   .constants .thruster
%                                                   .rungs_N .n_entries
%                                                   .families.labels .rib
%                                                   (.j .stop .nPts
%                                                   .nOffLattice .nD) .note
%
%% Revision History:
%  M. Casey                                                   (c) 10/03/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 4, opts = struct(); end
nD = fieldd(opts, 'nD', []);  tol = fieldd(opts, 'tol', 1e-8);
sD0 = 0;  if isfield(S, 'problem') && isfield(S.problem, 'sD'), sD0 = S.problem.sD; end
wrapd = @(x) abs(mod(x + 0.5, 1) - 0.5);                        % circular distance from 0

% ---- the rib must be the sheet column it names ----------------------------
assert(R.j >= 1 && R.j <= numel(S.sA) && wrapd(S.sA(R.j) - R.sA) < tol, 'rib_as_catalog:column', ...
       'the rib (s_A %.6f) is not column %d of the sheet (s_A %.6f)', R.sA, R.j, S.sA(min(max(R.j, 1), numel(S.sA))));
assert(all(isfinite(S.Z8(:, R.j))), 'rib_as_catalog:spine', 'the sheet has no certified crossing in column %d', R.j);
pts = R.pts;
assert(all(arrayfun(@(p) wrapd(p.sA - R.sA) < tol, pts)), 'rib_as_catalog:arrival', ...
       'the rib''s points are not all at its arrival phase %.6f', R.sA);
Zr = [pts.z];
assert(isempty(Zr) || all(isfinite(Zr(:))), 'rib_as_catalog:finite', 'a rib point carries non-finite costates');

% ---- the column: the spine crossing, then the rib's points --------------
sDall = mod([sD0, pts.sD], 1);
Z = [S.Z8(:, R.j), Zr];
keep = true(size(sDall));
if ~isempty(nD)
    k = (sDall - sD0)*nD;
    keep = abs(k - round(k)) < tol*nD;
end
nOff = nnz(~keep(2:end));
sDall = sDall(keep);  Z = Z(:, keep);
[sDall, order] = sort(sDall);  Z = Z(:, order);
dup = [false, diff(sDall) < tol];
assert(~any(dup), 'rib_as_catalog:duplicate', 'two rib points at departure phase %.6f', sDall(find(dup, 1)));
n = numel(sDall);

r = ref.sheets(1);
sheet = struct('sD_frac', sDall, 'sA_frac', R.sA, 'has_solution', true(n, 1), 'tf_nd', Z(8, :).', ...
               'entry_index', (1:n).', 'z8', Z, 'family_index', ones(n, 1), ...
               'tauDRO', r.tauDRO, 'Np', r.Np, 'pm', r.pm, 'period_tulip_nd', r.period_tulip_nd);
c = struct('sheets', sheet, 'constants', ref.constants, 'thruster', ref.thruster, 'rungs_N', ref.rungs_N, ...
           'n_entries', n, 'families', struct('labels', {{'rib'}}), ...
           'rib', struct('j', R.j, 'stop', R.stop, 'nPts', numel(pts), 'nOffLattice', nOff, 'nD', nD), ...
           'note', sprintf(['a departure rib at s_A = %.4f (sheet column %d, spine s_D = %g) wrapped as a one-column ' ...
                            'catalog (rib_as_catalog); %d entries, %d off-lattice point(s) dropped; rib stop: %s; a view for study'], ...
                           R.sA, R.j, sD0, n, nOff, R.stop));
outMat = fieldd(opts, 'outMat', '');
if ~isempty(outMat), catalog = c;  save(outMat, 'catalog'); end
end

% ---------------------------------------------------------------------------
function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
