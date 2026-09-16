function ok = test_survey_family_bounds()
%% Purpose:
%
%   Tests survey_family_bounds -- the per-family admissibility survey
%   (periselene >= 500 km, whole orbit within 100 Mm of the Moon) -- on
%   real DRO members from pumpkyn/pumpkynPie, plus one member that cannot
%   be built:
%     1. ROWS: one row per member that could be built; the failing member
%        is reported and EXCLUDED, not recorded.
%     2. METRICS: each row's period, periselene altitude and maximum Moon
%        distance agree with an independent recomputation from
%        get_family_orbit.
%     3. VERDICT: admissible == (periAlt >= 500 km && maxMoon <= 1e5 km),
%        row for row, B.admissible mirrors the rows, and both verdicts
%        occur (tau = 1 is admissible, tau = 4 reaches 133 Mm and is not).
%     4. CRITERIA are recorded in B.
%     5. SAVE: outMat holds the same struct.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);

mu      = 0.012150585609624;
lStar   = 389703.264829278;
tStar   = 382981.289129055;
rMoonKm = 1737.4;

grid_ = struct('tau', {1.0, 4.0, NaN});              % 4: reaches 133 Mm (NOT admissible); NaN: cannot be built
outMat = [tempname '.mat'];
cleaner = onCleanup(@() deleteIfThere(outMat));
B = survey_family_bounds('dro', grid_, outMat);

%% 1. Rows:
ok = chk(ok, numel(B.rows) == 2 && isequal([B.rows.params], grid_(1:2)), ...
         sprintf('two buildable members recorded, the failing one excluded (%d rows)', numel(B.rows)));

%% 2-3. Metrics and verdict against an independent recomputation:
for k = 1:numel(B.rows)
    [tt, rv] = get_family_orbit('dro', grid_(k));
    dM  = vecnorm(rv(:,1:3) - [1-mu 0 0], 2, 2) * lStar;
    pa  = min(dM) - rMoonKm;  mx = max(dM);
    row = B.rows(k);
    ok = chk(ok, abs(row.periodDays - tt(end)*tStar/86400) < 1e-9 && ...
                 abs(row.periAltKm - pa) < 1e-6 && abs(row.maxMoonKm - mx) < 1e-6, ...
             sprintf('tau = %g: period %.3f d, periselene %.0f km, max %.0f km reproduced', ...
                     grid_(k).tau, row.periodDays, row.periAltKm, row.maxMoonKm));
    ok = chk(ok, row.admissible == (pa >= 500 && mx <= 1e5) && B.admissible(k) == row.admissible, ...
             sprintf('tau = %g: admissible = %d follows the criteria', grid_(k).tau, row.admissible));
end

%% 4. Criteria:
ok = chk(ok, isequal(B.admissible(:).', [true false]), ...
         'one admissible member and one rejected (tau = 4 leaves the 100 Mm sphere)');
ok = chk(ok, B.criteria.altFloorKm == 500 && B.criteria.maxDistKm == 1e5 && strcmp(B.family, 'dro'), ...
         'criteria and family recorded');

%% 5. Save:
S = load(outMat);
ok = chk(ok, isequaln(S, B), 'outMat holds the same struct');

if ok, fprintf('TEST_SURVEY_FAMILY_BOUNDS: ALL PASS\n');
else,  fprintf('TEST_SURVEY_FAMILY_BOUNDS: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function deleteIfThere(f)
%% Purpose:
%
%   Remove a temporary file if it exists.
%
if isfile(f), delete(f); end
end

function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
