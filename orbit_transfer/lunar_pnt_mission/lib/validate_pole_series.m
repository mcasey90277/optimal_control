function V = validate_pole_series(M, opts)
%% Purpose:
%
%   VALIDATION GATE of the study's south-pole GDOP series, against an
%   independent construction. With the six satellites on their slots for
%   one tulip period, the pole series this study computes (mission_gdop on
%   the 'pole' grid, slot states from the phase_state handles of
%   get_family_orbit) is compared with a DIRECT gdop_grid_from_states call
%   on the paper's own constellation construction,
%   pumpkyn.cr3bp.tulipConstellation(Np, tau0, 1, pm, dtau6), at the same
%   epochs (the engine's half-open resampling, asserted to 1e-12 ND) and
%   the same point (latitude -90, longitude 0, mask 0).
%
%   The two orbit constructions differ by up to about 3 km (tulipConstellation
%   applies no continuation), so the series are not identical. Rule:
%     compared samples  both finite, both GDOP < tol.gdopMax, and the same
%                       visibility count
%     relative diff     max |ours - ref|/ref over the compared samples must
%                       be < tol.rel
%     edges             a visibility-count or finite/NaN mismatch is allowed
%                       ONLY within one sample of a change of the
%                       reference's visibility count (a satellite crossing
%                       the horizon between two nearly equal geometries)
%     coverage          at least 90 percent of the epochs compared (no
%                       vacuous PASS)
%   PASS iff all four hold; FAIL otherwise. UNRESOLVED when the phasing in
%   M differs from tulip_phasing() (the reference would be a different
%   constellation). tulipConstellation's own muStar, lStar and tStar must
%   equal M's to 1e-12 relative (validate_pole_series:constants otherwise),
%   so the reference states are in the mission's units.
%   Tolerance choice: tol.rel = 1e-3 sits about an order of magnitude
%   above the nominal relative difference (the two orbit constructions
%   differ by a few km) and several orders below what moving one slot by
%   0.3 ND produces; the measured numbers are in tests/
%   test_validate_pole_series (printed) and the Task 6 report, not here.
%
%% Inputs:
%
%  M                        struct                  mission_inputs() output
%  opts                     struct (optional)
%    .dtauOverride          [1 x 6] []              Slot offsets of OUR series
%                                                   only [ND] (deliberate-
%                                                   failure seam; the reference
%                                                   keeps tulip_phasing())
%    .Nt                    scalar [round(tauT*tStar/300)]   Epochs over one period
%    .tol                   struct [rel 1e-3, gdopMax 20]
%
%% Outputs:
%
%  V.status                 char                    'PASS' | 'FAIL' | 'UNRESOLVED'
%  V.reason                 char                    One line with the numbers
%  V.maxRel, V.p99Rel       scalar                  Relative GDOP difference
%  V.nCompared, V.Nt        scalar                  Compared samples, epochs
%  V.nVisMismatch           scalar                  Visibility-count mismatches
%  V.nVisOffEdge            scalar                  ... not next to an edge
%  V.nPatternMismatch       scalar                  finite/NaN mismatches
%  V.nPatternOffEdge        scalar                  ... not next to an edge
%  V.tDays                  [Nt x 1]                Epochs [d]
%  V.ours, V.ref            [Nt x 1]                The two GDOP series
%  V.tol, V.wallS
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  M. Casey  header: no run-specific counts                   (c) 11/10/2026
%  M. Casey  reference constants asserted; override phase wrap (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    V = validate_pole_series(mission_inputs());
    fprintf('%s | %s (%.1f s)\n', V.status, V.reason, V.wallS);
    return;
end
if nargin < 2, opts = struct(); end
wall0 = tic;
g     = mission_grid('pole');
Nt    = getopt(opts, 'Nt', round(M.tauT*M.tStar/g.dtSec));
tol   = getopt(opts, 'tol', struct('rel', 1e-3, 'gdopMax', 20));
V = struct('status', 'UNRESOLVED', 'reason', '', 'maxRel', NaN, 'p99Rel', NaN, ...
           'nCompared', 0, 'Nt', Nt, 'nVisMismatch', NaN, 'nVisOffEdge', NaN, ...
           'nPatternMismatch', NaN, 'nPatternOffEdge', NaN, 'tDays', [], ...
           'ours', [], 'ref', [], 'tol', tol, 'wallS', NaN);

%% Reference identity: the paper's phasing
dtauRef = tulip_phasing();
if numel(M.dtau6) ~= 6 || max(abs(M.dtau6(:).' - dtauRef(:).')) > 1e-12
    V.reason = sprintf('M.dtau6 %s differs from tulip_phasing() %s', ...
                       mat2str(M.dtau6, 15), mat2str(dtauRef, 15));
    V.wallS = toc(wall0);
    return;
end

%% Ours: on-slots history, mission_gdop on the pole grid
Mo = M;
if isfield(opts, 'dtauOverride') && ~isempty(opts.dtauOverride)
    Mo.dtau6      = opts.dtauOverride(:).';
    Mo.slotPhase0 = mod(Mo.dtau6, M.tauT)/M.tauT;   % same wrap as mission_gdop's validate mode
end
H = assemble_histories(Mo, [], 1:6, struct('onSlots', true, 'tEnd', M.tauT));
G = mission_gdop(H, Mo, [0 M.tauT], g, struct('score', false, 'keepSeries', true, ...
                 'quiet', true, 'Nt', Nt));

%% Reference: tulipConstellation states straight into the engine
[tauC, rC, ~, muC, lStarC, tStarC] = pumpkyn.cr3bp.tulipConstellation(M.sh.Np, M.tauT, 1, ...
                                                                      M.sh.pm, dtauRef);
cst = [muC lStarC tStarC];  cstM = [M.muStar M.lStar M.tStar];
if any(abs(cst - cstM) > 1e-12*abs(cstM))
    error('validate_pole_series:constants', ['tulipConstellation constants [%.15g %.15g %.15g] ' ...
          'differ from the mission''s [%.15g %.15g %.15g] (muStar, lStar, tStar)'], cst, cstM);
end
[tauC, iu] = unique(tauC(:));
rC = rC(iu, :, :);
if abs(tauC(1)) > 1e-12 || abs(tauC(end) - M.tauT) > 1e-9
    error('validate_pole_series:reference', 'tulipConstellation grid spans [%g %g], not one period', ...
          tauC(1), tauC(end));
end
p = struct('lat_range', g.chunks{1}, 'lon_range', g.lon_range, 'lat_res', g.lat_res, ...
           'lon_res', g.lon_res, 'elevation_mask', g.mask, 'time_samples', Nt);
R = gdop_grid_from_states(tauC, rC, M.muStar, M.lStar, M.tStar, p);
if max(abs(R.tau_i(:) - G.t(:))) > 1e-12
    error('validate_pole_series:epochs', 'reference epochs differ from ours by %.3g ND', ...
          max(abs(R.tau_i(:) - G.t(:))));
end

%% Compare
a  = G.gdop(:, 1);                 b  = R.gdop_time_series(:, 1);
va = double(G.vis(:, 1));          vb = double(R.visibility_count(:, 1));
chg  = [false; diff(vb) ~= 0];
edge = chg | [chg(2:end); false];                    % sample on either side of a change
visMis = va ~= vb;
patMis = isfinite(a) ~= isfinite(b);
cmp = isfinite(a) & isfinite(b) & a < tol.gdopMax & b < tol.gdopMax & ~visMis;
rel = abs(a(cmp) - b(cmp))./b(cmp);
V.maxRel           = max([rel; 0]);
V.p99Rel           = prctile([rel; 0], 99);
V.nCompared        = nnz(cmp);
V.nVisMismatch     = nnz(visMis);
V.nVisOffEdge      = nnz(visMis & ~edge);
V.nPatternMismatch = nnz(patMis);
V.nPatternOffEdge  = nnz(patMis & ~edge);
V.tDays = G.t*M.tStar/86400;
V.ours  = a;  V.ref = b;
pass = V.maxRel < tol.rel && V.nVisOffEdge == 0 && V.nPatternOffEdge == 0 && ...
       V.nCompared >= 0.9*Nt;
if pass, V.status = 'PASS'; else, V.status = 'FAIL'; end
V.reason = sprintf(['%d of %d epochs compared: max rel %.2e (p99 %.2e, tol %.0e); ' ...
                    'visibility mismatches %d (%d off an edge); finite/NaN mismatches %d ' ...
                    '(%d off an edge)'], V.nCompared, Nt, V.maxRel, V.p99Rel, tol.rel, ...
                   V.nVisMismatch, V.nVisOffEdge, V.nPatternMismatch, V.nPatternOffEdge);
V.wallS = toc(wall0);
end

%% ======================================================================
function v = getopt(s, f, v0)
% GETOPT  s.(f) if present and non-empty else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
