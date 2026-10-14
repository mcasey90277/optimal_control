function ok = test_extremes_movies_short()
% TEST_EXTREMES_MOVIES_SHORT  The extremes movie builder films a PAIR of
% library rows (a cell's primary and its fastest status-2 alternative) with
% pumpkyn's own propagator: both .mp4 and .gif are written, the mp4 frames
% read back at 1280 x 720, the returned stem is the one asked for, both
% panels were flown from their junction states, each flight ARRIVES on the
% arrival orbit at the row's phase (an oracle from the orbit table, not from
% the movie code), a forced z8 flight of the primary lands where the
% junction flight does, each flight's time base runs 0 .. z8(8) and passes
% through every stored junction, and the function's dependency closure
% anywhere in the repo is exactly {costate_catalog_extremes, status_key,
% get_family_orbit} (no flight-layer file). Frame-capped (6) to stay short.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
cc = fullfile(fileparts(fileparts(here)), 'costate_common');
addpath(here, cc);
L = load(fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
[prim, alt2] = sample_status_pair(c);
ok = chk(ok, strcmp(prim.kind, 'primary') && strcmp(alt2.kind, 'alternative') && alt2.status == 2, ...
    sprintf('sample pair: primary t_f %.3f d, status-2 alternative t_f %.3f d', prim.tfDays, alt2.tfDays));
stem = fullfile(tempname, 'pair_movie');  mkdir(fileparts(stem));
opts = struct('pair', {{prim, alt2}}, 'labels', 'PA', 'nFrames', 6, 'fly', 'auto');
[got, info] = costate_catalog_extremes_movies(c, [], 'time', stem, opts);
ok = chk(ok, strcmp(got, stem), 'the returned stem is the one asked for');
ok = chk(ok, isfile([stem '.mp4']) && isfile([stem '.gif']), 'both .mp4 and .gif exist');
vr = VideoReader([stem '.mp4']);
ok = chk(ok, vr.Width == 1280 && vr.Height == 720, sprintf('mp4 frames are %d x %d', vr.Width, vr.Height));
ok = chk(ok, isequal(info.methods, {'junctions', 'junctions'}), sprintf('both panels flown from junctions (%s, %s)', info.methods{:}));
ok = chk(ok, info.nFrames == 6, 'the frame cap was honored');

%% Oracle: each flight ends on the arrival orbit at the row's arrival phase,
%  runs on a time base from 0 to the row's own t_f, and passes through every
%  stored junction (the old engine's defect check, recomputed here)
sh = c.sheets(1);
pA = sh.arr_params;  pA.muStar = c.constants.muStar;
[tT, rvT] = get_family_orbit(sh.arr_family, pA);
lStar = c.constants.lStar_km;  tStar = c.constants.tStar_s;
rows = {prim, alt2};
for k = 1:2
    F = info.F{k};
    rvf = interp1(tT/tT(end), rvT, mod(rows{k}.sA, 1), 'spline');
    missKm = sqrt(sum((F.Y(end, 1:3) - rvf(1:3)).^2))*lStar;
    ok = chk(ok, missKm < 1, sprintf('panel %d arrives on the arrival orbit at sA = %.4f (miss %.3f km)', k, rows{k}.sA, missKm));
    tfRow = rows{k}.z8(8);
    ok = chk(ok, F.t(1) == 0 && issorted(F.t) && abs(F.t(end) - tfRow) < 1e-12 && abs(F.tfNd - tfRow) < 1e-12, ...
        sprintf('panel %d time base runs 0 .. z8(8) = %.6f ND, sorted', k, tfRow));
    ok = chk(ok, abs(F.tfDays - tfRow*tStar/86400) < 1e-9 && abs(F.dvKms - rows{k}.dvKms) < 1e-6, ...
        sprintf('panel %d t_f and Delta-V match the row (%.3f d, %.4f km/s)', k, F.tfDays, F.dvKms));
    J = rows{k}.junctions;  K = size(J, 2);  h = tfRow/K;  worstKm = 0;
    for q = 2:K
        Yq = interp1(F.t, F.Y, (q - 1)*h);
        worstKm = max(worstKm, sqrt(sum((Yq(1:3) - J(1:3, q).').^2))*lStar);
    end
    ok = chk(ok, worstKm < 1e-3, sprintf('panel %d passes through its %d stored junctions (worst %.2e km)', k, K, worstKm));
end

%% A forced z8 flight of the primary lands where its junction flight does
opts.fly = 'z8';  opts.pair = {prim, prim};
[~, infoZ] = costate_catalog_extremes_movies(c, [], 'time', [stem '_z8'], opts);
ok = chk(ok, isequal(infoZ.methods, {'z8', 'z8'}), 'opts.fly = z8 is honored');
dEnd = sqrt(sum((infoZ.F{1}.Y(end, 1:3) - info.F{1}.Y(end, 1:3)).^2))*lStar;
ok = chk(ok, dEnd < 1, sprintf('z8 flight ends within %.3f km of the junction flight', dEnd));

%% Dependency closure inside the repo: no flight-layer file (anywhere in the repo)
repoRoot = fileparts(fileparts(fileparts(here)));
req = matlab.codetools.requiredFilesAndProducts(fullfile(here, 'costate_catalog_extremes_movies.m'));
inRepo = req(startsWith(req, repoRoot));
[~, names] = cellfun(@fileparts, inRepo, 'UniformOutput', false);
names = setdiff(names, 'costate_catalog_extremes_movies');
want = {'costate_catalog_extremes', 'get_family_orbit', 'status_key'};
ok = chk(ok, isequal(sort(names(:)), sort(want(:))), sprintf('repo dependencies are exactly {%s} (got {%s})', ...
    strjoin(want, ', '), strjoin(names, ', ')));
if ok, fprintf('test_extremes_movies_short: ALL PASS\n'); else, fprintf('test_extremes_movies_short: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
