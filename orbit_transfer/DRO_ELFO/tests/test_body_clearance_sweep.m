function ok = test_body_clearance_sweep()
%% Purpose:
%
%   Test of body_clearance_sweep (Task 7). Sections:
%     1  SYNTHETIC PERILUNE BETWEEN SAMPLES (deliberate failure of the
%        sample minimum): the flight seam returns, for the top-rung entry, a
%        straight flyby of the Moon sampled at two times 600 km above the
%        surface with the closest approach, 150 km, between them (cubic
%        Hermite is exact on a straight line, propagate = false). The
%        sample minimum (600 km) would pass the 500 km floor; the located
%        minimum 150 km is recorded, its time is the mid-time, and the entry
%        is demoted (OK = false), the second entry (5000 km flyby) kept;
%        the backup <sheet>.bak_clearance is written and sheet_key survives
%        the rewrite
%     2  Z8 IDENTITY: a second call flies nothing; after one z8 is changed
%        on the sheet only that entry is re-flown (the stale record is not
%        reused) and its new record carries the new z8
%     3  floors are applied to the raw record: a 1e6 km Earth floor on an
%        unchanged record flies nothing and counts every entry below
%     4  a flight that throws counts as not located and below; an existing
%        backup is never overwritten; an incomplete census (maxEntries)
%        demotes nothing
%     5  REAL ENTRIES (no seam, the true min-time propagator): a 1 x 1 x 2
%        cut of the matched-experiment sheet; both entries clear 500 km
%        (Moon) and 200 km (Earth), the located minimum is never above the
%        sample minimum, and its time lies inside [0, t_f]
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true iff every check
%                                                   passed (throws on a
%                                                   failure)
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here    = fileparts(mfilename('fullpath'));
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(fileparts(here));
setup_dro_elfo_paths();
addpath(here);
wd = tempname;  mkdir(wd);
cleanWd = onCleanup(@() rmdir(wd, 's'));

realMat = fixture_dro_elfo_sheet(fullfile(wd, 'real', 'dro_elfo_tau1_raan000.mat'), ...
                                 struct('iD', 2, 'iA', 2, 'nR', 2));
Q0 = load(realMat);
mu = Q0.meta.muStar;  L = Q0.meta.lStar;  Rm = 1737.4;
zTop = Q0.Z8(:, 1, 1, 1);

%% 1. synthetic perilune between samples
synMat = fullfile(wd, 'syn', 'dro_elfo_tau1_raan000.mat');
mkdir(fileparts(synMat));  copyfile(realMat, synMat);
nCalls = 0;
flyBy = @(hPeri, hSample) flyby(mu, L, Rm, hPeri, hSample);
seam  = @(z8, rv0, rvf, phys) pick_flight(z8, zTop, flyBy);
o = struct('flyFcn', @counted, 'propagate', false);
S1 = body_clearance_sweep(synMat, o);
E1 = load(S1.sideMat);
assert(S1.done && S1.nFlown == 2 && nCalls == 2, 'both entries flown once');
assert(abs(E1.MOONKM(1, 1, 1) - 150) < 1e-6, 'located perilune 150 km (got %.6f)', E1.MOONKM(1, 1, 1));
assert(abs(E1.TMOON(1, 1, 1) - 0.5) < 1e-9, 'perilune at the mid-time (got %.6f)', E1.TMOON(1, 1, 1));
Fs = flyBy(150, 600);
hs = L*sqrt(sum((Fs.Y(:, 1:3) - [1 - mu 0 0]).^2, 2)) - Rm;
assert(min(hs) > 500 && abs(min(hs) - 600) < 1e-6, 'the SAMPLE minimum (%.3f km) would pass the 500 km floor', min(hs));
assert(abs(E1.MOONKM(1, 1, 2) - 5000) < 1e-6, 'second entry 5000 km');
assert(S1.nBelow == 1 && S1.nBelowMoon == 1 && isequal(S1.demoted, [1 1 1]), 'only the perilune entry demoted');
Q1 = load(synMat);
assert(~Q1.OK(1, 1, 1) && Q1.OK(1, 1, 2), 'OK after demotion');
assert(isfile([synMat(1:end-4) '.bak_clearance']), 'backup beside the sheet');
B1 = load([synMat(1:end-4) '.bak_clearance'], '-mat');
assert(all(B1.OK(:)), 'the backup holds the pre-demotion OK');
assert(isequaln(Q1.sheet_key, Q0.sheet_key), 'sheet_key survives the demotion rewrite');
assert(isfield(Q1.meta, 'clearanceSweep') && Q1.meta.clearanceSweep.nDemoted == 1, 'meta stamped');
fprintf('  synthetic flyby: sample min %.1f km, located %.6f km at t %.6f -> demoted\n', ...
        min(hs), E1.MOONKM(1, 1, 1), E1.TMOON(1, 1, 1));

%% 2. Z8 identity
nCalls = 0;
S2 = body_clearance_sweep(synMat, o);
assert(S2.nFlown == 0 && nCalls == 0 && S2.nBelow == 1, 'resume: nothing re-flown, the demotion record kept');
Qz = load(synMat);  Z8 = Qz.Z8;  Z8(3, 1, 1, 2) = Z8(3, 1, 1, 2)*(1 + 1e-9);
save(synMat, '-append', 'Z8');
S2b = body_clearance_sweep(synMat, o);
E2 = load(S2b.sideMat);
assert(S2b.nFlown == 1 && nCalls == 1, 'exactly the changed entry re-flown (got %d)', S2b.nFlown);
assert(isequaln(E2.Z8SWEPT(:, 1, 1, 2), Z8(:, 1, 1, 2)), 'the new record carries the new z8');
assert(isequaln(E2.Z8SWEPT(:, 1, 1, 1), E1.Z8SWEPT(:, 1, 1, 1)), 'the unchanged record kept');

%% 3. floors on the raw record
nCalls = 0;
S3 = body_clearance_sweep(synMat, struct('flyFcn', @counted, 'propagate', false, ...
                                         'floorEarthKm', 1e6, 'demote', false));
assert(S3.nFlown == 0 && S3.nBelow == 2 && S3.nBelowEarth == 2, 'impossible Earth floor: all below, nothing flown');

%% 4. throwing flight, backup kept, incomplete census
thrMat = fullfile(wd, 'thr', 'dro_elfo_tau1_raan000.mat');
mkdir(fileparts(thrMat));  copyfile(realMat, thrMat);
marker = 'kept';  save([thrMat(1:end-4) '.bak_clearance'], 'marker', '-mat');
S4a = body_clearance_sweep(thrMat, struct('flyFcn', @(varargin) error('t:x', 'boom'), ...
                                          'propagate', false, 'maxEntries', 1));
Q4 = load(thrMat);
assert(~S4a.done && S4a.nTodo == 1 && all(Q4.OK(:)), 'incomplete census demotes nothing');
S4 = body_clearance_sweep(thrMat, struct('flyFcn', @(varargin) error('t:x', 'boom'), 'propagate', false));
Q4 = load(thrMat);
assert(S4.nNonFinite == 2 && S4.nBelow == 2 && ~any(Q4.OK(:)), 'thrown flights: not located, below, demoted');
Bm = load([thrMat(1:end-4) '.bak_clearance'], '-mat');
assert(isfield(Bm, 'marker'), 'an existing backup is never overwritten');

%% 5. real entries, true propagator
S5 = body_clearance_sweep(realMat);
E5 = load(S5.sideMat);
Q5 = load(realMat);
assert(S5.done && S5.nDone == 2 && S5.nBelow == 0 && all(Q5.OK(:)), 'real entries clear both floors');
for kr = 1:2
    ndp = nd_propulsion(Q5.rungs(kr), Q5.meta.ispS, Q5.meta.m0kg, L, Q5.meta.tStar);
    [tD, rvD] = get_family_orbit(Q5.meta.depFamily, Q5.meta.depParams);
    dep = phase_state(tD, rvD);
    rv0 = dep(Q5.sD(1));
    [~, y] = pumpkyn.cr3bp.tfMinProp(Q5.Z8(8, 1, 1, kr), [rv0(:); 1; Q5.Z8(1:7, 1, 1, kr)], ...
                                     ndp.Tnd, ndp.cnd, mu);
    hSm = min(L*sqrt(sum((y(:, 1:3) - [1 - mu 0 0]).^2, 2)) - Rm);
    assert(E5.MOONKM(1, 1, kr) <= hSm + 1e-9, 'located <= sample minimum');
    assert(E5.TMOON(1, 1, kr) >= 0 && E5.TMOON(1, 1, kr) <= Q5.Z8(8, 1, 1, kr) + 1e-12, 'time inside [0, tf]');
    fprintf('  real %g N entry: Moon located %.3f km (sample %.3f) at %.5f ND, Earth %.1f km\n', ...
            Q5.rungs(kr), E5.MOONKM(1, 1, kr), hSm, E5.TMOON(1, 1, kr), E5.EARTHKM(1, 1, kr));
end

ok = true;
fprintf('test_body_clearance_sweep: PASS\n');

    function F = counted(z8, rv0, rvf, phys)
        % COUNTED  The flight seam, counting its calls.  INPUTS: z8 [8 x 1];
        % rv0, rvf [6 x 1]; phys [1 x 1] struct.  OUTPUTS: F [1 x 1] struct.
        nCalls = nCalls + 1;
        F = seam(z8, rv0, rvf, phys);
    end
end

% ------------------------------------------------------------------------
function F = pick_flight(z8, zTop, flyBy)
% PICK_FLIGHT  The perilune flyby for the top-rung entry, a far one else.
% INPUTS: z8 [8 x 1]; zTop [8 x 1]; flyBy [1 x 1] function handle.
% OUTPUTS: F [1 x 1] struct (.t [2 x 1], .Y [2 x 14]).
if isequaln(z8(:), zTop(:)), F = flyBy(150, 600); else, F = flyBy(5000, 6000); end
end

function F = flyby(mu, L, Rm, hPeri, hSample)
% FLYBY  Straight-line Moon flyby: closest approach hPeri km at t = 0.5,
% samples at t = 0 and 1, each hSample km above the surface (ND units).
% INPUTS: mu, L km, Rm km, hPeri km, hSample km, each [1 x 1].  OUTPUTS: F .t [2 x 1],
% .Y [2 x 14].
d0 = (Rm + hPeri)/L;
v  = 2*sqrt((Rm + hSample)^2 - (Rm + hPeri)^2)/L;     % ND distance per ND time
t  = [0; 1];
pos = [1 - mu + d0 + 0*t, v*(t - 0.5), 0*t];
F = struct('t', t, 'Y', [pos, repmat([0 v 0], 2, 1), ones(2, 1), zeros(2, 7)]);
end
