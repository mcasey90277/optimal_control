function ok = test_validate_pole_series()
%% Purpose:
%
%   Verifies the south-pole validation gate and the 'pole' grid:
%   (1) mission_grid('pole') is one point at latitude -90 exactly,
%   longitude 0, 300 s, mask 0, and the engine puts it at the Moon's south
%   pole (finite position 1737 km below the Moon's centre); (2) with the
%   six satellites on their slots for one tulip period, the pole series of
%   mission_gdop (phase_state slots) agrees with a direct
%   gdop_grid_from_states call on pumpkyn tulipConstellation states: status
%   PASS, measured relative difference far below the tolerance, no
%   finite/NaN or visibility mismatch away from a visibility edge;
%   (3) DELIBERATE FAILURE: slot 6 moved 0.3 ND in OUR series only must
%   give FAIL (the reference identity is unchanged).
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true if every assertion
%                                                   held (else it throws)
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));
M = mission_inputs();

% ---- 1. the pole grid ----
g = mission_grid('pole');
assert(isequal(g.chunks, {[-90 -90]}) && isequal(g.lon_range, [0 0]) && ...
       g.dtSec == 300 && g.mask == 0, 'pole grid definition');
pos = pumpkyn.cr3bp.fromLLA([], [-pi/2 0 0], M.muStar, M.lStar, 2);
dKm = (pos(1:3) - [1 - M.muStar 0 0])*M.lStar;
assert(all(isfinite(dKm)) && abs(dKm(3) + 1737.1) < 1 && hypot(dKm(1), dKm(2)) < 1e-6, ...
       'pole position [%g %g %g] km', dKm);
fprintf('  pole grid: one point (-90, 0), Moon-centred position [%.3g %.3g %.1f] km\n', dKm);

% ---- 2. nominal: PASS ----
V = validate_pole_series(M);
fprintf('  nominal: %s | %s\n', V.status, V.reason);
assert(strcmp(V.status, 'PASS'), 'nominal status %s', V.status);
assert(V.maxRel < V.tol.rel/5, 'nominal max rel %.3g not well inside tol %.3g', V.maxRel, V.tol.rel);
assert(V.nPatternOffEdge == 0 && V.nVisOffEdge == 0, 'nominal pattern mismatches');
assert(V.nCompared > 0.9*V.Nt, 'too few compared samples (%d of %d)', V.nCompared, V.Nt);

% ---- 3. deliberate failure: slot 6 moved 0.3 ND in our series ----
d2 = M.dtau6;  d2(6) = d2(6) + 0.3;
V2 = validate_pole_series(M, struct('dtauOverride', d2));
fprintf('  moved slot 6: %s | %s\n', V2.status, V2.reason);
assert(strcmp(V2.status, 'FAIL'), 'moved slot must FAIL (got %s)', V2.status);
ok = true;
fprintf('  PASS\n');
end
