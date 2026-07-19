% TEST_MEE_RES_TO_CART  Task 5: adapter correctness + visualization front door.
%
% Three checks against the certified 10 N min-fuel case
% (results/MEE_M2_10N.mat, m_f=1377.10 kg):
%   (a) mee_res_to_cart_res reconstructs the certified numbers: apogee-start
%       |r|~1.1028, GEO-end |r|~1.000, unit inertial thrust direction
%       (|alpha_ECI|=1 at every node, since beta_RTN is unit by construction
%       and Rrtn2eci is orthonormal), and mass preserved (1377.10 kg).
%   (b) transfer_movie accepts a `res` STRUCT directly (not only a .mat
%       path) -- exercised by actually calling it with the struct in a
%       time-boxed subprocess (kill-guarded so the test doesn't sit through
%       a full 300-frame render) and confirming it does NOT fail on the
%       old `load(matFile)` line (which errors with
%       'MATLAB:string:MustBeStringScalarOrCharacterVector' when matFile is
%       a struct -- confirmed pre-fix in this task's TDD red step).
%   (c) gergaud_plot renders a static trajectory PNG from a Cartesian res.
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/mee_res_to_cart_res.m (adapter under test).
%   [2] earth_elliptic_to_geo/transfer_movie.m (struct-or-path input).
%   [3] earth_elliptic_to_geo/gergaud_plot.m (static trajectory PNG).

here = fileparts(mfilename('fullpath')); cd(here);
scratchDir = fullfile(tempdir, 'gergaud_test_scratch');

%% (a) adapter reconstruction on the certified 10 N case
S = load(fullfile(module_root(),'results','MEE_M2_10N.mat'));
c = mee_res_to_cart_res(S.res.fuel.X, S.res.fuel.U, S.res.fuel.dL, S.res.sigma, 10, 1.5, 1);
r = c.fuel.X(1:3,:); rmag = vecnorm(r); an = vecnorm(c.fuel.U(1:3,:));
assert(abs(rmag(1)-1.1028) < 1e-3, 'apogee start ~1.103');
assert(abs(rmag(end)-1.0) < 1e-3, 'GEO end ~1.000');
assert(all(abs(an-1) < 1e-9), 'unit inertial thrust dir');
assert(abs(1500*c.fuel.X(7,end) - 1377.10) < 0.1, 'mass preserved');
fprintf('test_mee_res_to_cart (a) adapter reconstruction PASSED\n');

%% (b) struct handoff: transfer_movie(structRes, outStem) must not fail
% trying to load() a struct. transfer_movie always renders 300 frames
% internally (hardcoded, not something this task may change), so we run it
% in a detached subprocess and kill it after a few seconds -- "a couple of
% frames is enough" -- rather than sit through the full render here.
if ~exist(scratchDir, 'dir'), mkdir(scratchDir); end
tmpDir  = tempname(scratchDir); mkdir(tmpDir);
probeM  = fullfile(tmpDir, 'tm_probe.m');
logFile = fullfile(tmpDir, 'tm_probe.log');
outStem = fullfile(tmpDir, 'tm_probe_out');
wrapSh  = fullfile(tmpDir, 'run_with_kill.sh');
killSecs = 18;

fidP = fopen(probeM, 'w');
fprintf(fidP, [ ...
    'hereP = ''%s''; cd(hereP);\n' ...
    'Sp = load(fullfile(hereP,''results'',''MEE_M2_10N.mat''));\n' ...
    'cp = mee_res_to_cart_res(Sp.res.fuel.X, Sp.res.fuel.U, Sp.res.fuel.dL, Sp.res.sigma, 10, 1.5, 1);\n' ...
    'fprintf(''CALLING_TRANSFER_MOVIE_WITH_STRUCT\\n'');\n' ...
    'transfer_movie(cp, ''%s'');\n' ...
    'fprintf(''TRANSFER_MOVIE_RETURNED_FULLY\\n'');\n'], module_root(), outStem);
fclose(fidP);

fidS = fopen(wrapSh, 'w');
fprintf(fidS, [ ...
    '#!/bin/bash\n' ...
    '/Applications/MATLAB_R2025b.app/bin/matlab -batch "run(''%s'')" > %s 2>&1 &\n' ...
    'pid=$!\n' ...
    '(sleep %d; kill -TERM $pid 2>/dev/null; sleep 2; kill -9 $pid 2>/dev/null) &\n' ...
    'wd=$!\n' ...
    'wait $pid; rc=$?\n' ...
    'kill $wd 2>/dev/null\n' ...
    'exit $rc\n'], probeM, logFile, killSecs);
fclose(fidS);
system(sprintf('chmod +x %s', wrapSh));
system(wrapSh);   % blocks up to ~killSecs+2s; exit code deliberately unchecked (kill is expected)

logTxt = '';
if isfile(logFile), logTxt = fileread(logFile); end
assert(contains(logTxt, 'CALLING_TRANSFER_MOVIE_WITH_STRUCT'), ...
    'struct-handoff probe did not even reach the transfer_movie call:\n%s', logTxt);
assert(~contains(logTxt, 'MustBeStringScalarOrCharacterVector') && ...
       ~contains(logTxt, 'Argument must be a text scalar') && ...
       ~contains(logTxt, 'Error using load'), ...
    'transfer_movie(struct,...) must not fail trying to load() a struct -- probe log:\n%s', logTxt);

rmdir(tmpDir, 's');  % scratch-only artifacts (probe script, log, partial movie frames)
fprintf('test_mee_res_to_cart (b) transfer_movie struct handoff PASSED\n');

%% (c) gergaud_plot: static trajectory PNG from a Cartesian res
pngOut = fullfile(scratchDir, 'test_gergaud_plot_smoke.png');
if isfile(pngOut), delete(pngOut); end
gergaud_plot(c, pngOut, 'test smoke: MEE M2 10N');
d = dir(pngOut);
assert(~isempty(d) && d.bytes > 0, 'gergaud_plot must write a non-empty PNG');
fprintf('test_mee_res_to_cart (c) gergaud_plot smoke PASSED (%s, %d bytes)\n', pngOut, d.bytes);

fprintf('test_mee_res_to_cart PASSED\n');
