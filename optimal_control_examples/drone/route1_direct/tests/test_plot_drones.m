function ok = test_plot_drones()
%% Purpose:
%
%   plot_drones renders the M = 4 fixture to a PNG headlessly (no visible
%   figure), the file is non-empty, and the figure has four axes;
%   movie_drones writes a short MP4 whose frame size is 1280 x 720 (a
%   multiple of 16: anything else shears into diagonal streaks); the front
%   door script run_drones.m parses (no checkcode error).
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
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
S = load(fullfile(here, 'data', 'drones_direct_M4.mat'));
png = [tempname, '.png'];  mp4 = [tempname, '.mp4'];
cleanup = onCleanup(@() cellfun(@(f) delete(f), {png, mp4}));

fig = plot_drones(S.out, struct('visible', false, 'file', png));
ax = findobj(fig, 'Type', 'axes');
ok = chk(ok, numel(ax) == 4, sprintf('four panels (%d axes)', numel(ax)));
ok = chk(ok, isfile(png) && dir(png).bytes > 1000, 'PNG written');
close(fig);

movie_drones(S.out, mp4, struct('nFrames', 12, 'fps', 6));
ok = chk(ok, isfile(mp4), 'MP4 written');
if isfile(mp4)
    vr = VideoReader(mp4);
    ok = chk(ok, vr.Width == 1280 && vr.Height == 720, sprintf('frame %dx%d (multiple of 16)', vr.Width, vr.Height));
end

msgs = checkcode(fullfile(here, 'run_drones.m'), '-id');
isErr = false(1, numel(msgs));
for km = 1:numel(msgs)
    isErr(km) = strncmp(msgs(km).id, 'E', 1) || contains(msgs(km).message, 'Parse error');
end
ok = chk(ok, ~any(isErr), sprintf('run_drones.m parses (%d checkcode errors)', sum(isErr)));

if ok, fprintf('TEST_PLOT_DRONES: ALL PASS\n'); else, fprintf('TEST_PLOT_DRONES: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
