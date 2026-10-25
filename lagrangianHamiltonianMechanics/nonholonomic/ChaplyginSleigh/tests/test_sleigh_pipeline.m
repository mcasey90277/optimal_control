function ok = test_sleigh_pipeline()
%% Purpose:
%
%   The pieces run_sleigh.m calls, headless and into a temp folder: the
%   body shape, plot_sleigh's two PNGs (also for the degenerate rest
%   case), a 20-frame movie through mech.animate_planar_body with the
%   sleigh's arrows and subtitle, and a parse check of the script itself
%   (scripts take no arguments, so the script is run by hand once).
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
sdir = fileparts(here);
root = fileparts(fileparts(sdir));
addpath(sdir, fullfile(root, 'mechlib'));

p = sleigh_params(1, 0.1, 0.2);
t = 0:0.02:4;
[~, yy] = ode113(@(a, b) sleigh_rhs(a, b, p), t, [0; 0; 0; 1; 2], odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
Y = yy.';
out = sleigh_quantities(t, Y, p);

body = sleigh_body(p);
ok = check(ok, size(body.outline, 1) == 2 && size(body.outline, 2) >= 4, 'sleigh_body: outline is [2 x P], P >= 4');
names = {body.markers.name};
ok = check(ok, any(strcmp(names, 'com')) && any(strcmp(names, 'blade')), 'sleigh_body: markers com and blade');
cm = body.markers(strcmp(names, 'com')).pos;
ok = check(ok, max(abs(cm - [p.a; 0])) < 1e-15, 'sleigh_body: center of mass at (a, 0) in blade coordinates');
ok = check(ok, numel(body.segments) >= 1, 'sleigh_body: the blade segment is present');

outDir = fullfile(tempdir, 'test_sleigh_pipeline');
if ~exist(outDir, 'dir'), mkdir(outDir); end
if ~isempty(dir(fullfile(outDir, 'tsp*'))), delete(fullfile(outDir, 'tsp*')); end
files = plot_sleigh(t, Y, out, p, struct('label', 'tsp', 'visible', false), outDir);
ok = check(ok, numel(files) == 2 && all(cellfun(@(f) exist(f, 'file') == 2, files)), 'plot_sleigh writes its two PNGs');
outR = sleigh_quantities(t, zeros(5, numel(t)), p);
filesR = plot_sleigh(t, zeros(5, numel(t)), outR, p, struct('label', 'tsp_rest', 'visible', false), outDir);
ok = check(ok, numel(filesR) == 2 && all(cellfun(@(f) exist(f, 'file') == 2, filesR)), 'plot_sleigh survives the rest case (C = 0, degenerate ellipse)');

arrows = sleigh_arrows(t, Y, out, p);
body.arrows   = arrows;
body.subtitle = @(tk) sprintf('t = %5.2f s', tk);
file = mech.animate_planar_body(t, [out.rb; Y(3, :)], body, fullfile(outDir, 'tsp'), ...
                                struct('fps', 10, 'seconds', 2, 'traces', cat(3, out.rc, out.rb), 'title', 'test'));
ok = check(ok, exist(file, 'file') == 2, sprintf('movie written through mech.animate_planar_body: %s', file));

msgs = checkcode(fullfile(sdir, 'run_sleigh.m'), '-id');
bad = msgs(arrayfun(@(mm) startsWith(mm.id, 'E') || contains(mm.message, 'Parse error'), msgs));
ok = check(ok, isempty(bad), sprintf('run_sleigh.m parses (%d checkcode errors)', numel(bad)));

if ok, fprintf('TEST_SLEIGH_PIPELINE: ALL PASS\n'); else, fprintf('TEST_SLEIGH_PIPELINE: FAILURES\n'); end
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
