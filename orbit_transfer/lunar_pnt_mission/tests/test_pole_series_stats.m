function ok = test_pole_series_stats()
%% Purpose:
%
%   Verifies pole_series_stats on hand-built series whose answers are
%   counted by hand: (1) all NaN; (2) one outage of known length at the
%   START, in the MIDDLE and at the END of the service window (and a longer
%   outage BEFORE the window, which must not count); (3) values exactly 6
%   (an outage, counted in 6 <= GDOP < 7) and an Inf (an outage, not
%   finite, excluded from the median); (4) first finite / first < 6 days on
%   the whole series; (5) a length mismatch is refused.
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
dt = 300;
tD = (0:19).'*dt/86400;                       % 20 samples
t0 = tD(6);                                   % service = samples 6..20 (15)

% ---- 1. all NaN ----
s = pole_series_stats(tD, nan(20, 1), t0, dt);
assert(isnan(s.firstFiniteDay) && isnan(s.firstLt6Day), 'all NaN: first days');
assert(s.svcSamples == 15 && abs(s.svcDays - 15*dt/86400) < 1e-15, 'all NaN: service samples');
assert(s.fracLt6 == 0 && s.fracFinite == 0 && isnan(s.medianGdop) && isnan(s.p95Gdop), ...
       'all NaN: fractions / median');
assert(abs(s.longestOutageH - 15*dt/3600) < 1e-12 && s.n6to7 == 0 && s.nLt6 == 0, ...
       'all NaN: outage %.4f h', s.longestOutageH);
fprintf('  all NaN: 15 service samples, outage %.2f h, median NaN\n', s.longestOutageH);

% ---- 2. one outage of known length at start / middle / end ----
% pre-window samples 1..5 are an outage of 5 (NaN); they must not count
where = {6:8, 12:15, 17:20};                  % lengths 3, 4, 4
lens  = [3 4 4];
for kw = 1:3
    g = 5*ones(20, 1);  g(1:5) = NaN;
    g(where{kw}) = 9;
    s = pole_series_stats(tD, g, t0, dt);
    assert(abs(s.longestOutageH - lens(kw)*dt/3600) < 1e-12, 'outage %d: %.4f h', kw, s.longestOutageH);
    assert(abs(s.fracLt6 - (15 - lens(kw))/15) < 1e-15 && s.fracFinite == 1, 'outage %d: fractions', kw);
    assert(s.nLt6 == 15 - lens(kw) && s.n6to7 == 0 && s.medianGdop == 5, 'outage %d: counts', kw);
    assert(s.firstFiniteDay == tD(6) && s.firstLt6Day == tD(6 + (kw == 1)*3), 'outage %d: first days', kw);
end
fprintf('  outages at start / middle / end: %d / %d / %d samples recovered\n', lens);

% ---- 3. exactly 6 and an Inf ----
g = 4*ones(20, 1);
g(3) = 5.5;                                    % pre-window, first < 6
g(10) = 6;  g(11) = 6;                         % exactly 6: outage, 6 <= GDOP < 7
g(14) = Inf;                                   % outage, not finite
g(16) = 6.999;                                 % in [6, 7)
s = pole_series_stats(tD, g, t0, dt);
assert(s.n6to7 == 3 && s.nLt6 == 11, 'six/Inf: n6to7 %d nLt6 %d', s.n6to7, s.nLt6);
assert(abs(s.fracLt6 - 11/15) < 1e-15 && abs(s.fracFinite - 14/15) < 1e-15, 'six/Inf: fractions');
assert(abs(s.longestOutageH - 2*dt/3600) < 1e-12, 'six/Inf: longest outage %.4f h', s.longestOutageH);
gf = [4*ones(11, 1); 6; 6; 6.999];             % the 14 finite service samples
assert(s.medianGdop == median(gf) && s.medianGdop == 4 && s.p95Gdop == prctile(gf, 95), ...
       'six/Inf: median %.4f p95 %.4f', s.medianGdop, s.p95Gdop);
assert(s.firstFiniteDay == 0 && s.firstLt6Day == 0, 'six/Inf: first days');
fprintf('  exactly 6 and Inf: n6to7 %d, nLt6 %d, finite %.4f, median %.2f, p95 %.4f\n', ...
        s.n6to7, s.nLt6, s.fracFinite, s.medianGdop, s.p95Gdop);

% ---- 4. first finite / first < 6 before the window ----
g = nan(20, 1);  g(4) = 8;  g(9:20) = 5;
s = pole_series_stats(tD, g, t0, dt);
assert(s.firstFiniteDay == tD(4) && s.firstLt6Day == tD(9), 'first days: %.6f %.6f', ...
       s.firstFiniteDay, s.firstLt6Day);
assert(abs(s.longestOutageH - 3*dt/3600) < 1e-12, 'first days: outage samples 6..8');

% ---- 5. length mismatch ----
thrown = '';
try
    pole_series_stats(tD, nan(19, 1), t0, dt);
catch err
    thrown = err.identifier;
end
assert(strcmp(thrown, 'pole_series_stats:shape'), 'length mismatch not refused (%s)', thrown);
ok = true;
fprintf('  PASS\n');
end
