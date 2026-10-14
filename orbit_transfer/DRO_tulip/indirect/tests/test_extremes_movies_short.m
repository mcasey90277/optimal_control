function ok = test_extremes_movies_short()
% TEST_EXTREMES_MOVIES_SHORT  The extremes movie builder films a PAIR of
% library rows (a cell's primary and its fastest status-2 alternative)
% through the fly layer: both .mp4 and .gif are written, the mp4 frames read
% back at 1280 x 720, the returned stem is the one asked for, and both
% panels were flown from their junction states. Frame-capped (6) to stay
% short.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
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
if ok, fprintf('test_extremes_movies_short: ALL PASS\n'); else, fprintf('test_extremes_movies_short: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
