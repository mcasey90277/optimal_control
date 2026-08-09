% TEST_VIZ_SMOKE  Viz functions execute headless and write files.
% Not a beauty contest -- existence + nonzero size only. Movie smoke uses
% 2 seconds of frames to stay fast.
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
plot_pdg_solution(solC, solV, fullfile(od, 'sol.png'));
plot_footprint(mc, P, fullfile(od, 'fp.png'));
movie_landing(out, solC, P, fullfile(od, 'mov.mp4'), struct('duration', 2));
for f = {'sol.png','fp.png','mov.mp4'}
    d = dir(fullfile(od, f{1}));
    assert(~isempty(d) && d.bytes > 1e3, 'missing/empty %s', f{1});
end
fprintf('test_viz_smoke PASS\n');
