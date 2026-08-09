% TEST_VIZ_SMOKE  Viz functions execute headless and write files.
% Not a beauty contest -- existence + nonzero size (plus, for the
% cinematic movie, the frame-size lock). Both movie smokes are trimmed to
% ~2 seconds of playback to stay fast.
%
% Exercises BOTH halves of the figure-ownership contract added in the
% task-9 fix report (2026-08-09): plot_pdg_solution/plot_footprint are
% called here WITH their output captured, so this test -- not the
% function -- owns the handles and is responsible for closing them
% (nargout>0 means the function itself leaves them open). A bare-statement
% call (nargout==0) is the OTHER half of the contract and is exercised
% implicitly every time either function is called elsewhere in this repo
% without capturing an output.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
solC = solve_pdg_colloc(P, struct('N', 30));
solV = solve_pdg_convex(P, struct('tf', solC.tf, 'Nconv', 60));
ctrl = tvlqr_design(solC, P);
out  = sim_closed_loop(solC, ctrl, P, struct());
mc   = run_monte_carlo(solC, ctrl, P, struct('Nrun', 8));

od = fullfile(tempdir, 'bl_viz_smoke');  if ~exist(od,'dir'), mkdir(od); end
cleanupOd = onCleanup(@() rmdir(od, 's'));   % always remove the tempdir,
                                              % pass or fail
fig1 = plot_pdg_solution(solC, solV, fullfile(od, 'sol.png'));
close(fig1);
fig2 = plot_footprint(mc, P, fullfile(od, 'fp.png'));
close(fig2);
movie_landing(out, solC, P, fullfile(od, 'mov.mp4'), struct('duration', 2));

% CINEMATIC movie (viz/movie_landing_cinematic.m) -- presentation product,
% NOT part of the front door, so nothing else in the suite would touch it.
% Trimmed to ~2 s of playback the only way this movie's timeline can be
% trimmed: it has no 'duration' knob (its length is set by the flight plus
% the title card, fade, slow-motion finale and hold), so the card/hold and
% the slow-motion lead are cut instead. Existence + nonzero + the
% divide-by-16 1920x1080 frame lock, which is the one property of this
% function a silent regression could plausibly break.
movie_landing_cinematic(out, solC, P, fullfile(od, 'cine.mp4'), ...
    struct('titleSec',0.3, 'fadeSec',0.2, 'holdSec',0.3, 'slowLead',0.5, ...
           'slowRamp',0.3, 'fps',20));
for f = {'sol.png','fp.png','mov.mp4','cine.mp4'}
    d = dir(fullfile(od, f{1}));
    assert(~isempty(d) && d.bytes > 1e3, 'missing/empty %s', f{1});
end
vC = VideoReader(fullfile(od, 'cine.mp4'));
assert(vC.Width == 1920 && vC.Height == 1080, ...
    'cinematic frame size %dx%d, expected 1920x1080 (divide-by-16 lock)', ...
    vC.Width, vC.Height);
assert(vC.NumFrames > 10, 'cinematic movie has only %d frames', vC.NumFrames);
fprintf('test_viz_smoke PASS\n');
