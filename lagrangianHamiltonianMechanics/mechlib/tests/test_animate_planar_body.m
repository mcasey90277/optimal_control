function ok = test_animate_planar_body()
%% Purpose:
%
%   mech.animate_planar_body writes a 20-frame MP4 with exact 1280 x 720
%   frames from a rectangle sliding on a circle (with a marker, a segment,
%   an arrow and a subtitle), survives a single-point pose, and refuses an
%   unsorted time grid.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every check passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

t = linspace(0, 2, 41);
pose = [cos(pi*t); sin(pi*t); pi*t + pi/2];
body.outline  = [-0.15 0.15 0.15 -0.15; -0.08 -0.08 0.08 0.08];
body.markers  = struct('name', 'nose', 'pos', [0.15; 0], 'color', [1 0.3 0.3]);
body.segments = struct('pts', [-0.1 0.1; 0 0], 'color', [0.5 0.9 1], 'width', 3);
body.arrows   = struct('name', 'v', 'fun', @(tk) [cos(pi*tk); sin(pi*tk); -0.3*sin(pi*tk); 0.3*cos(pi*tk)], 'color', [1 0.75 0.3]);
body.subtitle = @(tk) sprintf('t = %.2f s', tk);
stem = fullfile(tempdir, 'test_animate_planar_body');
file = mech.animate_planar_body(t, pose, body, stem, struct('fps', 10, 'seconds', 2, 'title', 'test'));
ok = check(ok, exist(file, 'file') == 2, sprintf('MP4 written: %s', file));
d = dir(file);
ok = check(ok, d.bytes > 10000, sprintf('MP4 is %d bytes (> 10000)', d.bytes));
vr = VideoReader(file);
fr = readFrame(vr);
ok = check(ok, isequal(size(fr, 1:2), [720 1280]), sprintf('frame is %d x %d (720 x 1280)', size(fr, 1), size(fr, 2)));
ok = check(ok, abs(vr.NumFrames - 20) <= 1, sprintf('%d frames (20 asked)', vr.NumFrames));

file1 = mech.animate_planar_body([0 1], [0 0; 0 0; 0 0], struct('outline', body.outline), fullfile(tempdir, 'test_animate_point'), struct('fps', 5, 'seconds', 0.4));
ok = check(ok, exist(file1, 'file') == 2, 'a single-point pose (zero path extent) still writes a movie with finite limits');

ok = check(ok, throws(@() mech.animate_planar_body([0 2 1], pose(:, 1:3), body, stem), 'mech:animate_planar_body:t'), ...
           'unsorted t refused with :t');
ok = check(ok, throws(@() mech.animate_planar_body(t, pose(1:2, :), body, stem), 'mech:animate_planar_body:shape'), ...
           'pose with two rows refused with :shape');

tNaN = t;  tNaN(5) = NaN;
ok = check(ok, throws(@() mech.animate_planar_body(tNaN, pose, body, stem), 'mech:animate_planar_body:t'), ...
           't containing a NaN refused with :t');
ok = check(ok, throws(@() mech.animate_planar_body(t, pose, body, stem, struct('traces', pose(1:2, 1:end-1))), 'mech:animate_planar_body:shape'), ...
           'opts.traces of size [2 x (N-1)] refused with :shape');

if ok, fprintf('TEST_ANIMATE_PLANAR_BODY: ALL PASS\n'); else, fprintf('TEST_ANIMATE_PLANAR_BODY: FAILURES\n'); end
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
%% Purpose:
%
%   Print one PASS/FAIL line and fold the condition into the verdict.
%
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && logical(cond);
end

% ------------------------------------------------------------------------
function tf = throws(fun, id)
%% Purpose:
%
%   True when fun() throws an error with identifier id.
%
tf = false;
try
    fun();
catch err
    tf = strcmp(err.identifier, id);
end
end
